//! koth_ff: A high-performance LC-MS feature finder for timsTOF and mzML data.
//!
//! The in-memory detection stages live in the [`koth_core`] crate, which has no
//! file I/O or native dependencies. They are re-exported here at their original
//! paths (`koth_ms::hills`, `koth_ms::models`, `koth_ms::config::KothConfig`,
//! `koth_ms::run_features`, ...); this crate adds the mzML, Bruker `.d` and
//! Thermo `.raw` readers, the output writers, cross-run alignment and LFQ, and
//! the `koth_ff` / `koth_align` executables.
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
//! use koth_ms::{run_hills, run_features, run_scoring, config::KothConfig};
//!
//! let config = KothConfig::default();
//! let input = Path::new("data.mzML");
//!
//! let spectra = koth_ms::read_spectra(input, &config.file).unwrap();
//! let hills = run_hills(&spectra, &config.hills, &config.file);
//! let features = run_features(&hills, &config.features, &config.file).unwrap();
//! let scored = run_scoring(&features, &config.scoring, &config.features, config.file.polarity);
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
pub mod input;
pub mod io;
pub mod lfq;
pub mod output;
pub mod pipeline;

// The in-memory stages, re-exported from koth-core at their original paths.
pub use koth_core::{features, hills, mem, models, scoring, stats};

// Stage-by-stage entry points, re-exported from koth-core.
pub use koth_core::{run_features, run_hills, run_scoring};

// In-process / streaming API surface. Re-exported at the crate root so callers
// write `koth_ms::run_pipeline` / `koth_ms::PipelineSink` without reaching into
// the module path.
pub use pipeline::{
    group_ms2_hills_by_window, run_pipeline, run_pipeline_from_spectra, run_pipeline_streaming,
    run_pipeline_streaming_from_spectra, run_pipeline_with_ms2, FeatureFindingOutput,
    PipelineOptions, PipelineSink,
};

// Re-export the core in-memory result types at the crate root so a downstream
// crate (koth_tracer, uno) can name them without depending on the module layout.
pub use models::{Feature, Hill, IsolationWindow, Polarity, ScoredFeature, Spectrum};

use std::path::Path;

use config::{FileConfig, HillsConfig};
use error::KothError;
use koth_core::shuffle_if_decoy;

/// Read MS1 spectra from an mzML file or Bruker .d directory.
pub fn read_spectra(path: &Path, file: &FileConfig) -> Result<Vec<Spectrum>, KothError> {
    io::read_spectra(path, file)
}

/// Stage 1 (streaming): Detect hills from mzML, Bruker `.d`, or Thermo `.raw`
/// using bounded spectrum buffers. Completed hills are retained for assembly;
/// native readers also retain scan metadata. Decoy shuffling and optional TIC
/// normalization still collect spectra. Gzipped mzML is decompressed as it is
/// parsed, never held whole in memory.
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
    let (iter, status) = io::mzml::stream_mzml_ms2_checked(path)?;
    log::info!("Streaming MS2 hill detection from {}", path.display());
    let hills = hills::detect_ms2_hills_from_iter(iter, config, file);
    status.check()?;
    Ok(hills)
}

fn hills_streaming_inner(
    path: &Path,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
    let iter = io::stream_spectra(path, file)?;
    if file.decoy_mode {
        // Full-run shuffling deliberately retains its collecting semantics.
        let spectra = shuffle_if_decoy(iter.collect::<Result<Vec<_>, _>>()?, file);
        return Ok(hills::detect_hills_from_iter(
            spectra.into_iter(),
            config,
            file,
        ));
    }

    log::info!("Streaming MS1 hill detection from {}", path.display());
    // Bruker already owns a bounded producer. Thermo and mzML are decoded on
    // this prefetch thread, overlapping I/O with hill detection.
    let iter: io::SpectrumStream = if io::detect_format(path) == io::InputFormat::BrukerD {
        iter
    } else {
        Box::new(io::prefetch::prefetch(iter, PREFETCH_CAPACITY))
    };
    detect_hills_from_results(iter, config, file)
}

/// Stop on the first reader error and discard partial hills. The detector's
/// infallible iterator API must never turn an I/O failure into a successful run.
fn detect_hills_from_results(
    iter: io::SpectrumStream,
    config: &HillsConfig,
    file: &FileConfig,
) -> Result<Vec<Hill>, KothError> {
    let mut error = None;
    let spectra = iter.map_while(|item| match item {
        Ok(spectrum) => Some(spectrum),
        Err(e) => {
            error = Some(e);
            None
        }
    });
    let hills = hills::detect_hills_from_iter(spectra, config, file);
    match error {
        Some(error) => Err(error),
        None => Ok(hills),
    }
}

/// Bounded look-ahead in spectra for the mzML and Thermo reader thread.
const PREFETCH_CAPACITY: usize = 64;

#[cfg(test)]
mod ms1_streaming_tests {
    use super::*;

    #[test]
    fn reader_error_discards_partial_detection() {
        let items = vec![
            Ok(Spectrum {
                scan_index: 0,
                retention_time: 1.0,
                peaks: vec![],
                ms_level: 1,
                isolation_window: None,
                faims_cv: None,
            }),
            Err(KothError::ThermoError("decode failure".into())),
        ];
        let result = detect_hills_from_results(
            Box::new(items.into_iter()),
            &HillsConfig::default(),
            &FileConfig::default(),
        );
        assert!(matches!(result, Err(KothError::ThermoError(_))));
    }
}
