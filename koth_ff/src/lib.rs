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
pub mod pipeline;
pub mod scoring;
pub mod stats;

// In-process / streaming API surface. Re-exported at the crate root so callers
// write `koth_ff::run_pipeline` / `koth_ff::PipelineSink` without reaching into
// the module path.
pub use pipeline::{
    group_ms2_hills_by_window, run_pipeline, run_pipeline_from_spectra, run_pipeline_streaming,
    run_pipeline_streaming_from_spectra, run_pipeline_with_ms2, FeatureFindingOutput,
    PipelineOptions, PipelineSink,
};

// Re-export the core in-memory result types at the crate root so a downstream
// crate (koth_tracer, uno) can name them without depending on the module layout.
pub use models::{Feature, Hill, IsolationWindow, ScoredFeature, Spectrum};

use std::path::Path;

use config::{FeaturesConfig, FileConfig, HillsConfig, ScoringConfig};
use error::KothError;
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

/// Stage 1 (MS2): Detect MS2 hills partitioned by the precursor isolation window
/// of each MS2 spectrum (DIA channels).
///
/// Supported inputs: mzML (`.mzML`/`.mzML.gz`), Bruker **diaPASEF** `.d`
/// (`--features tdf`), and Thermo **DIA** `.raw` (`--features thermo`). On a
/// diaPASEF `.d` the fixed isolation windows are read straight from the raw
/// frames' quadrupole settings and each window segment is IM-collapsed into an
/// MS2 spectrum (see [`io::bruker::read_bruker_ms2`]). On a DIA `.raw` each MS2
/// scan becomes one MS2 spectrum stamped with its precursor isolation window
/// (see [`io::thermo::read_thermo_ms2`]). ddaPASEF `.d` and DDA `.raw` are
/// unsupported for MS2 (precursor reconstruction is out of scope): they yield an
/// empty set plus a warning, leaving MS1 unaffected.
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

    #[cfg(feature = "tdf")]
    if name.ends_with(".d")
        || (path.is_dir() && path.extension().and_then(|e| e.to_str()) == Some("d"))
    {
        // Bruker diaPASEF: stream one MS2 spectrum per isolation-window segment
        // out of the raw `.d`, then group by window exactly as the mzML path
        // does. ddaPASEF / unknown acquisitions return an empty Vec (with a
        // warning) from `read_bruker_ms2`, so MS1 is left untouched.
        let spectra = io::bruker::read_bruker_ms2(path, file)?;
        if spectra.is_empty() {
            return Ok(Vec::new());
        }
        log::info!(
            "Streaming MS2 hill detection from {} ({} Bruker diaPASEF MS2 spectra)",
            path.display(),
            spectra.len()
        );
        return Ok(hills::detect_ms2_hills_from_iter(
            spectra.into_iter(),
            config,
            file,
        ));
    }

    #[cfg(feature = "thermo")]
    if name.ends_with(".raw") {
        // Thermo DIA: one MS2 spectrum per MS2 scan, stamped with its precursor
        // isolation window, then grouped by window exactly as the mzML path does.
        // A DDA `.raw` (or one that cannot be confidently classified as DIA)
        // returns an empty Vec (with a warning) from `read_thermo_ms2`, so MS1 is
        // left untouched.
        let spectra = io::thermo::read_thermo_ms2(path)?;
        if spectra.is_empty() {
            return Ok(Vec::new());
        }
        log::info!(
            "Streaming MS2 hill detection from {} ({} Thermo DIA MS2 spectra)",
            path.display(),
            spectra.len()
        );
        return Ok(hills::detect_ms2_hills_from_iter(
            spectra.into_iter(),
            config,
            file,
        ));
    }

    if !(name.ends_with(".mzml") || name.ends_with(".mzml.gz")) {
        log::warn!(
            "MS2 hill detection supports mzML, Bruker diaPASEF .d, and Thermo DIA .raw \
             (.raw needs --features thermo) inputs only; skipping '{}'",
            path.display()
        );
        return Ok(Vec::new());
    }
    let iter = io::mzml::stream_mzml_ms2(path)?;
    log::info!("Streaming MS2 hill detection from {}", path.display());
    Ok(hills::detect_ms2_hills_from_iter(iter, config, file))
}

fn hills_streaming_inner(path: &Path, config: &HillsConfig, file: &FileConfig) -> Result<Vec<Hill>, KothError> {
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

    // Thermo .raw: no streaming API, so batch-load all MS1 spectra (mirrors the
    // Bruker path above). `read_spectra` routes `.raw` to the native reader when
    // built with `--features thermo`, or returns a clear "rebuild with
    // --features thermo" error otherwise — handled here (rather than the mzML
    // fall-through below) so the message is actionable in both builds.
    if path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("raw")) {
        let mut spectra = io::read_spectra(path, file)?;
        if file.decoy_mode {
            log::info!(
                "Decoy mode: shuffling {} Thermo spectra before hill detection",
                spectra.len()
            );
            spectra.shuffle(&mut rand::thread_rng());
        }
        return Ok(hills::detect_hills(&spectra, config, file));
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
        // Decode the mzML (decompress + XML parse + peak extraction) on a
        // background reader thread so it overlaps with `process_scan` on the
        // consumer side. Order is preserved, so hill detection is unchanged.
        let iter = io::prefetch::prefetch(iter, PREFETCH_CAPACITY);
        Ok(hills::detect_hills_from_iter(iter, config, file))
    }
}

/// Bounded look-ahead (in spectra) for the mzML reader thread. At ~1–2k peaks
/// per MS1 scan this caps the prefetch buffer at a few MB while giving the
/// consumer enough slack to stay busy across decode-time variance.
const PREFETCH_CAPACITY: usize = 64;

/// Stage 2: Detect isotope features from hills.
pub fn run_features(hills: &[Hill], config: &FeaturesConfig, file: &FileConfig) -> Result<Vec<Feature>, KothError> {
    if file.mz_recalibration {
        if file.decoy_mode {
            log::info!(
                "Decoy mode active — skipping m/z recalibration (shuffled spectra would skew the surface)"
            );
        } else {
            match features::learn_recal_model(hills, config, file) {
                Some(model) => {
                    log::info!(
                        "m/z recalibration: learned isotope-consistency surface from {} spacing samples (global offset {:.3} ppm, σ {:.3} ppm) — applying in pass 2",
                        model.n_samples,
                        model.global_offset(),
                        model.global_sigma,
                    );
                    return Ok(features::detect_features_with_recal(hills, config, file, Some(&model)));
                }
                None => {
                    log::warn!(
                        "m/z recalibration enabled but too few isotope-spacing samples to build a surface — using uncorrected features"
                    );
                }
            }
        }
    }
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
    let mut scored = scoring::score_features(
        features,
        scoring_cfg,
        &features_cfg.sulfur_offsets,
    );
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
