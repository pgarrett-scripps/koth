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

/// Crate version with the git commit it was built from, e.g. `0.1.0 (a1b2c3d4e5f6)`,
/// or `0.1.0 (a1b2c3d4e5f6-dirty)` when the working tree had uncommitted changes.
///
/// This is what `--version` reports on both binaries, and what the benchmark's
/// pin check reads. `KOTH_GIT_SHA` is set by `build.rs`; it is `unknown` when the
/// crate is built from a source archive with no git directory.
pub const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), " (", env!("KOTH_GIT_SHA"), ")");

pub mod alignment;
pub mod config;
pub mod error;
pub mod features;
pub mod hills;
pub mod input;
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

/// Apply the decoy construction to a full spectrum set: when `file.decoy_mode`
/// is set, shuffle the **entire** set (Fisher–Yates over `rand::thread_rng`)
/// before hill detection; otherwise return it untouched.
///
/// This is the single definition of the "shuffle the whole set, then detect"
/// step. Every entry point — in-memory ([`run_hills`]), streaming
/// ([`hills_streaming_inner`]), and the MS2-buffered pipeline path
/// ([`pipeline::run_pipeline_streaming_from_spectra`]) — routes its decoy shuffle
/// through here so the RNG source and full-set semantics live in exactly one
/// place.
pub(crate) fn shuffle_if_decoy(mut spectra: Vec<Spectrum>, file: &FileConfig) -> Vec<Spectrum> {
    if file.decoy_mode {
        log::info!(
            "Decoy mode: shuffling {} spectra before hill detection",
            spectra.len()
        );
        spectra.shuffle(&mut rand::thread_rng());
    }
    spectra
}

/// Stage 1: Detect chromatographic hills from MS1 spectra.
pub fn run_hills(spectra: &[Spectrum], config: &HillsConfig, file: &FileConfig) -> Vec<Hill> {
    let spectra = shuffle_if_decoy(spectra.to_vec(), file);
    hills::detect_hills(&spectra, config, file)
}

/// Stage 1 (streaming): Detect hills by reading the mzML/Bruker file directly,
/// processing one spectrum at a time without building a `Vec<Spectrum>`.
/// This is the preferred API for large files — peak memory is O(active_hills)
/// rather than O(total_peaks).
pub fn run_hills_streaming(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
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
/// (see `io::thermo::read_thermo_ms2`). ddaPASEF `.d` and DDA `.raw` are
/// unsupported for MS2 (precursor reconstruction is out of scope): they yield an
/// empty set plus a warning, leaving MS1 unaffected.
pub fn run_ms2_hills_streaming(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
    let fmt = io::detect_format(path);

    #[cfg(feature = "tdf")]
    if fmt == io::InputFormat::BrukerD {
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
    if fmt == io::InputFormat::ThermoRaw {
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

    if !matches!(fmt, io::InputFormat::Mzml | io::InputFormat::MzmlGz) {
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

fn hills_streaming_inner(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
    let fmt = io::detect_format(path);

    #[cfg(feature = "tdf")]
    if fmt == io::InputFormat::BrukerD {
        // Bruker: still requires loading all frames (timsrust doesn't expose a streaming API)
        let spectra = shuffle_if_decoy(io::read_spectra(path, file)?, file);
        return Ok(hills::detect_hills(&spectra, config, file));
    }

    // Thermo .raw: no streaming API, so batch-load all MS1 spectra (mirrors the
    // Bruker path above). `read_spectra` routes `.raw` to the native reader when
    // built with `--features thermo`, or returns a clear "rebuild with
    // --features thermo" error otherwise — handled here (rather than the mzML
    // fall-through below) so the message is actionable in both builds.
    if fmt == io::InputFormat::ThermoRaw {
        let spectra = shuffle_if_decoy(io::read_spectra(path, file)?, file);
        return Ok(hills::detect_hills(&spectra, config, file));
    }

    // mzML path — collect first if decoy mode so we can shuffle
    let iter = io::mzml::stream_mzml(path)?;
    if file.decoy_mode {
        let spectra = shuffle_if_decoy(iter.collect(), file);
        log::info!("Streaming hill detection from {} (decoy)", path.display());
        Ok(hills::detect_hills_from_iter(
            spectra.into_iter(),
            config,
            file,
        ))
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
pub fn run_features(
    hills: &[Hill],
    config: &FeaturesConfig,
    file: &FileConfig,
) -> Result<Vec<Feature>, KothError> {
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
                    return Ok(features::detect_features_with_recal(
                        hills,
                        config,
                        file,
                        Some(&model),
                    ));
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
    let mut scored = scoring::score_features(features, scoring_cfg, &features_cfg.sulfur_offsets);
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

#[cfg(test)]
mod shuffle_tests {
    use super::*;
    use config::FileConfig;

    fn spectra(n: usize) -> Vec<Spectrum> {
        (0..n)
            .map(|i| Spectrum {
                scan_index: i,
                retention_time: i as f64 * 0.1,
                peaks: Vec::new(),
                ms_level: 1,
                isolation_window: None,
            })
            .collect()
    }

    /// Decoy off: the set is returned untouched, in the exact input order (the
    /// non-decoy path the byte-parity tests exercise).
    #[test]
    fn no_decoy_is_identity() {
        let file = FileConfig {
            decoy_mode: false,
            ..FileConfig::default()
        };
        let out = shuffle_if_decoy(spectra(64), &file);
        let idx: Vec<usize> = out.iter().map(|s| s.scan_index).collect();
        assert_eq!(idx, (0..64).collect::<Vec<_>>());
    }

    /// Decoy on: the FULL set is permuted with no drops or duplicates (a decoy
    /// run must shuffle every spectrum, not a subset). Order is `thread_rng`, so
    /// we assert the preserved multiset, which is the load-bearing invariant.
    #[test]
    fn decoy_permutes_full_set() {
        let file = FileConfig {
            decoy_mode: true,
            ..FileConfig::default()
        };
        let out = shuffle_if_decoy(spectra(256), &file);
        assert_eq!(out.len(), 256);
        let mut idx: Vec<usize> = out.iter().map(|s| s.scan_index).collect();
        idx.sort_unstable();
        assert_eq!(idx, (0..256).collect::<Vec<_>>());
    }
}
