//! koth_ff: A high-performance LC-MS feature finder for timsTOF and mzML data.
//!
//! ## Pipeline
//!
//! 1. **Hill detection** – groups MS1 peaks across scans into chromatographic hills
//! 2. **Feature finding** – groups hills into isotope envelopes with charge state assignment
//! 3. **Scoring** – scores isotope pattern quality using the averagine model
//!
//! ## Quick start
//!
//! ```no_run
//! use std::path::Path;
//! use koth_ff::{run_hills, run_features, run_scoring, config::KothConfig};
//!
//! let config = KothConfig::default();
//! let input = Path::new("data.mzML");
//!
//! let spectra = koth_ff::read_spectra(input, &config.file).unwrap();
//! let hills = run_hills(&spectra, &config.hills, &config.file);
//! let features = run_features(&hills, &config.features, &config.file).unwrap();
//! let scored = run_scoring(&features, &config.scoring, &config.features);
//! ```

pub mod alignment;
pub mod config;
pub mod input;
pub mod error;
pub mod features;
pub mod hills;
pub mod io;
pub mod lfq;
pub mod mem;
pub mod models;
pub mod output;
pub mod scoring;

use std::path::Path;

use config::{FeaturesConfig, FileConfig, HillsConfig, ScoringConfig, ToleranceType};
use error::KothError;
use models::{Feature, Hill, ScoredFeature, Spectrum};
use rand::seq::SliceRandom;

/// Read MS1 spectra from an mzML file or Bruker .d directory.
pub fn read_spectra(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    io::read_spectra(path, file)
}

/// Stage 1: Detect chromatographic hills from MS1 spectra.
pub fn run_hills(spectra: &[Spectrum], config: &HillsConfig, file: &FileConfig) -> Vec<Hill> {
    if file.decoy_mode {
        log::info!("Decoy mode: shuffling {} spectra before hill detection", spectra.len());
        let mut shuffled = spectra.to_vec();
        shuffled.shuffle(&mut rand::thread_rng());
        hills::detect_hills(&shuffled, config, file)
    } else {
        hills::detect_hills(spectra, config, file)
    }
}

/// Stage 1 (streaming): Detect hills by reading the mzML/Bruker file directly,
/// processing one spectrum at a time without building a Vec<Spectrum>.
/// This is the preferred API for large files — peak memory is O(active_hills)
/// rather than O(total_peaks).
pub fn run_hills_streaming(path: &Path, config: &HillsConfig, file: &FileConfig) -> Result<Vec<Hill>, KothError> {
    hills_streaming_inner(path, config, file)
}

/// Stage 1 (MS2): Detect MS2 hills from an mzML file, partitioned by the
/// precursor isolation window of each MS2 spectrum (DIA channels).
/// Only mzML is supported — Bruker .d MS2 frames are not yet handled.
pub fn run_ms2_hills_streaming(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_lowercase();
    if !(name.ends_with(".mzml") || name.ends_with(".mzml.gz")) {
        log::warn!(
            "MS2 hill detection only supports mzML inputs; skipping '{}'",
            path.display()
        );
        return Ok(Vec::new());
    }
    let iter = io::mzml::stream_mzml_ms2(path)?;
    log::info!("Streaming MS2 hill detection from {}", path.display());
    Ok(hills::detect_ms2_hills_from_iter(iter, config, file))
}

fn hills_streaming_inner(path: &Path, config: &HillsConfig, file: &FileConfig) -> Result<Vec<Hill>, KothError> {
    let effective_file = match maybe_calibrate(path, config, file)? {
        Some(f) => std::borrow::Cow::Owned(f),
        None => std::borrow::Cow::Borrowed(file),
    };
    let file = effective_file.as_ref();

    #[cfg(feature = "tdf")]
    {
        if path.extension().and_then(|e| e.to_str()) == Some("d") || path.is_dir() {
            // Bruker: still requires loading all frames (timsrust doesn't expose a streaming API)
            let mut spectra = io::read_spectra(path, file)?;
            if file.decoy_mode {
                log::info!(
                    "Decoy mode: shuffling {} Bruker spectra before hill detection",
                    spectra.len()
                );
                spectra.shuffle(&mut rand::thread_rng());
            }
            return Ok(hills::detect_hills(&spectra, config, file));
        }
    }

    // mzML path — collect first if decoy mode so we can shuffle
    let iter = io::mzml::stream_mzml(path)?;
    if file.decoy_mode {
        let mut spectra: Vec<_> = iter.collect();
        log::info!(
            "Decoy mode: shuffling {} spectra before hill detection",
            spectra.len()
        );
        spectra.shuffle(&mut rand::thread_rng());
        log::info!("Streaming hill detection from {} (decoy)", path.display());
        Ok(hills::detect_hills_from_iter(spectra.into_iter(), config, file))
    } else {
        log::info!("Streaming hill detection from {}", path.display());
        Ok(hills::detect_hills_from_iter(iter, config, file))
    }
}

/// If `adaptive_mz_tolerance` is enabled and the tolerance is ppm-typed,
/// run a wide-tolerance pass-1 sweep to gather the empirical ppm-delta
/// distribution and return a `FileConfig` whose `mz_tolerance` has been
/// replaced with `min(median + sigma_mult × σ, pass1_ceiling)`.
///
/// Returns `Ok(None)` when calibration is disabled or not applicable
/// (Dalton tolerances, decoy mode where shuffled spectra wouldn't give a
/// real-signal distribution, or when too few samples were recorded to
/// compute a meaningful σ).
fn maybe_calibrate(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Option<FileConfig>, KothError> {
    if !file.adaptive_mz_tolerance {
        return Ok(None);
    }
    if !matches!(file.mz_tolerance_type, ToleranceType::Ppm) {
        log::warn!(
            "adaptive_mz_tolerance is enabled but mz_tolerance_type is not Ppm — skipping calibration"
        );
        return Ok(None);
    }
    if file.decoy_mode {
        log::info!(
            "Decoy mode active — skipping adaptive m/z tolerance calibration (shuffled spectra would skew the distribution)"
        );
        return Ok(None);
    }

    let pass1_multiplier = file.adaptive_mz_tolerance_pass1_multiplier.max(1.0);
    let pass1_ceiling = file.mz_tolerance * pass1_multiplier;
    let mut pass1_file = file.clone();
    pass1_file.mz_tolerance = pass1_ceiling;
    // Disable adaptive in the pass-1 file so we don't recurse.
    pass1_file.adaptive_mz_tolerance = false;

    let calibrated = calibrate_streaming(path, config, &pass1_file)?;
    let Some((calibrated_ppm, n_samples)) = calibrated else {
        log::warn!(
            "Adaptive m/z calibration: insufficient samples to compute σ — keeping user-set tolerance {:.2} ppm",
            file.mz_tolerance
        );
        return Ok(None);
    };

    let final_ppm = calibrated_ppm.min(pass1_ceiling);
    log::info!(
        "Adaptive m/z tolerance: user-set {:.2} ppm, pass-1 ceiling {:.2} ppm, calibrated median+σ → {:.2} ppm (n={} samples)",
        file.mz_tolerance,
        pass1_ceiling,
        final_ppm,
        n_samples
    );

    let mut out = file.clone();
    out.mz_tolerance = final_ppm;
    // Disable adaptive in the returned config so downstream doesn't re-run.
    out.adaptive_mz_tolerance = false;
    Ok(Some(out))
}

/// Pass-1 sweep that returns `(calibrated_ppm, n_samples)` if enough
/// matches landed in the histogram, else `None`.
fn calibrate_streaming(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Option<(f64, u64)>, KothError> {
    use hills::detector::HillDetector;

    let sigma_mult = file.adaptive_mz_tolerance_sigma_mult;
    let noise_sigma = file.noise_filter_sigma;

    // Branch on input format. For Bruker we already buffer the spectra
    // (timsrust has no streaming API), so the pass-1 + pass-2 cost is
    // bounded by `2× detect_hills` over the same Vec. For mzML, pass 1
    // streams the file directly and discards everything except the
    // ppm-delta histogram.
    #[cfg(feature = "tdf")]
    {
        if path.extension().and_then(|e| e.to_str()) == Some("d") || path.is_dir() {
            let spectra = io::read_spectra(path, file)?;
            let mut det = HillDetector::new(config, file).with_calibration_recording();
            for mut spec in spectra {
                if let Some(sigma) = noise_sigma {
                    hills::noise::filter_spectrum(&mut spec, sigma);
                }
                det.process_scan(&spec);
            }
            let n = det.calibration_sample_count();
            return Ok(det.calibrated_tolerance_ppm(sigma_mult).map(|p| (p, n)));
        }
    }

    let iter = io::mzml::stream_mzml(path)?;
    let mut det = HillDetector::new(config, file).with_calibration_recording();
    for mut spec in iter {
        if let Some(sigma) = noise_sigma {
            hills::noise::filter_spectrum(&mut spec, sigma);
        }
        det.process_scan(&spec);
    }
    let n = det.calibration_sample_count();
    Ok(det.calibrated_tolerance_ppm(sigma_mult).map(|p| (p, n)))
}

/// Stage 2: Detect isotope features from hills.
pub fn run_features(hills: &[Hill], config: &FeaturesConfig, file: &FileConfig) -> Result<Vec<Feature>, KothError> {
    Ok(features::detect_features(hills, config, file))
}

/// Stage 3: Score isotope features using the averagine model.
///
/// Applies three AND-ed retention filters from `FeaturesConfig`:
///   - `min_isotope_score` (Bhattacharyya vs averagine)
///   - `min_cosine_score`  (mean chromatographic cosine of isotope hills)
///   - `min_combined_score` (= isotope × cosine)
///
/// All thresholds at 0.0 keep every scored feature.
pub fn run_scoring(
    features: &[Feature],
    scoring_cfg: &ScoringConfig,
    features_cfg: &FeaturesConfig,
) -> Vec<ScoredFeature> {
    let mut scored = scoring::score_features(features, scoring_cfg);
    let any_filter = features_cfg.min_isotope_score > 0.0
        || features_cfg.min_cosine_score > 0.0
        || features_cfg.min_combined_score > 0.0;
    if any_filter {
        let before = scored.len();
        scored.retain(|sf| {
            sf.isotope_score >= features_cfg.min_isotope_score
                && sf.cosine_score >= features_cfg.min_cosine_score
                && sf.combined_score >= features_cfg.min_combined_score
        });
        log::info!(
            "Retained {}/{} features after score filters (iso\u{2265}{:.2}, cos\u{2265}{:.2}, comb\u{2265}{:.2})",
            scored.len(),
            before,
            features_cfg.min_isotope_score,
            features_cfg.min_cosine_score,
            features_cfg.min_combined_score
        );
    }
    scored
}
