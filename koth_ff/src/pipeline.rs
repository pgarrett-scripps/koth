//! In-process / streaming feature-finding API.
//!
//! The in-memory half ([`PipelineSink`], [`PipelineOptions`],
//! [`FeatureFindingOutput`], [`run_pipeline_from_spectra`],
//! [`run_pipeline_streaming_from_spectra`], [`group_ms2_hills_by_window`]) lives in
//! `koth_core::pipeline` and is re-exported here; this module adds the
//! file-path entry points.
//!
//! This module exposes the `koth_ff` MS1 pipeline as a **library seam** so that
//! downstream tools can obtain the hills and isotope-chain features directly in
//! memory, without `koth_ff` writing intermediate `hills.parquet` /
//! `features.parquet` files that they then re-read.
//!
//! It is an *additional* entry point layered on top of the same staged
//! functions the `koth_ff` binary uses ([`crate::run_hills_streaming`],
//! [`crate::run_features`], [`crate::run_scoring`]). The binary's file-writing
//! path is unchanged, and the output produced here is bit-for-bit what that
//! path would serialize (see the parity tests in `pipeline_tests.rs`).
//!
//! # Streaming granularity — and why
//!
//! The pipeline has a hard data dependency: **isotope-chain assembly needs the
//! complete hill set** (it builds an m/z-sorted index over every hill, then
//! grows chains across it), and hills themselves are only finalized after the
//! whole run is read (co-elution *splitting*, baseline *filtering* and stable
//! *`hill_id` assignment* are all global post-passes). So there is no honest
//! "emit a finalized hill mid-scan" granularity: a hill's identity isn't known
//! until end-of-read.
//!
//! The API therefore reflects the algorithm rather than pretending otherwise:
//!
//! 1. With default settings, spectra arrive through bounded buffers during hill
//!    detection. Native readers retain scan metadata; completed hills are kept.
//!    See [`crate::run_hills_streaming`] for options that still collect spectra.
//! 2. When detection completes, the sink receives the **whole finalized hill
//!    set once**, by value ([`PipelineSink::on_hills`]).
//! 3. Features are then assembled and handed to the sink **one at a time, by
//!    value** ([`PipelineSink::on_feature`]) and freed as they are consumed — so
//!    a streaming consumer folds each feature into its own structures without
//!    ever holding the full feature vector.
//!
//! # Consumers
//!
//! * **`koth_tracer`** (DIA pseudo-spectra) implements [`PipelineSink`]: it
//!   takes ownership of the hill set in `on_hills` to build its trace index,
//!   then consumes features incrementally in `on_feature`. With
//!   [`PipelineOptions::emit_ms2`] it *also* receives the DIA fragment hills in
//!   `on_ms2_hills`, partitioned by precursor isolation window, so it can map a
//!   precursor feature to the window(s) that isolated it.
//! * **`uno`** (MS1 search) just wants the final feature list: it calls
//!   [`run_pipeline`] / [`run_pipeline_from_spectra`] and reads
//!   [`FeatureFindingOutput::features`].
//!
//! # MS2 (DIA fragment hills) — opt-in
//!
//! MS2 fragment hill detection ([`crate::run_ms2_hills_streaming`]) is exposed
//! through the *same* seam but is **off by default**. When
//! [`PipelineOptions::emit_ms2`] is set, the MS2 hills are delivered once via
//! [`PipelineSink::on_ms2_hills`] (or [`FeatureFindingOutput::ms2_hills`]) after
//! the whole MS1 side. Each MS2 hill carries its precursor isolation window in
//! [`Hill::isolation_window`](crate::models::Hill::isolation_window), and [`group_ms2_hills_by_window`] groups them per
//! DIA channel. MS2 is supported on **mzML**, **Bruker diaPASEF `.d`**
//! (`--features tdf`), and **Thermo DIA `.raw`** (`--features thermo`); DDA
//! acquisitions (ddaPASEF `.d`, DDA `.raw`) degrade to an empty MS2 set + a
//! warning, leaving MS1 untouched. With `emit_ms2` false the MS1 path is
//! byte-for-byte identical to before.

use std::path::Path;

use crate::config::KothConfig;
use crate::error::KothError;

pub use koth_core::pipeline::*;

/// Run the full pipeline from a file path, streaming results to `sink`.
/// No output files are written.
///
/// Hill detection routes through [`crate::run_hills_streaming`], using bounded
/// MS1 spectrum buffers for mzML, Bruker `.d`, and Thermo `.raw`. Decoy shuffling
/// and optional TIC normalization collect spectra before detection.
pub fn run_pipeline_streaming<S: PipelineSink>(
    path: &Path,
    config: &KothConfig,
    opts: &PipelineOptions,
    sink: &mut S,
) -> Result<(), KothError> {
    let hills = crate::run_hills_streaming(path, &config.hills, &config.file)?;
    run_pipeline_from_hills(hills, config, opts, sink)?;
    // MS2 is opt-in and independent of MS1: only after the whole MS1 side has
    // been delivered do we (re)read the input for DIA fragment hills. When
    // `emit_ms2` is false this branch is skipped entirely, so the MS1-only path
    // is byte-for-byte unchanged. mzML-only: `.d`/`.raw` return an empty set +
    // warning via `run_ms2_hills_streaming`.
    if opts.emit_ms2 {
        // MS2 hill detection uses the MS2-resolved config: `[hills]` with any
        // `[hills_ms2]` overrides applied. Absent `[hills_ms2]` ⇒ same as MS1.
        let ms2 = crate::run_ms2_hills_streaming(path, &config.ms2_hills(), &config.file)?;
        sink.on_ms2_hills(ms2);
    }
    Ok(())
}

/// Convenience: run the pipeline from a file path and collect all hills and
/// features into owned vectors. The all-in-memory shape `uno` uses.
pub fn run_pipeline(
    path: &Path,
    config: &KothConfig,
    opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError> {
    let mut out = FeatureFindingOutput::default();
    run_pipeline_streaming(path, config, opts, &mut out)?;
    Ok(out)
}

/// Convenience: run the full pipeline **plus** MS2 (DIA fragment) hill detection
/// from a file path, collecting MS1 hills, MS1 features, and MS2 hills into owned
/// vectors. The combined shape a `koth_tracer` consumer uses when it wants both
/// sides at once rather than a streaming sink.
///
/// This forces MS2 emission on (equivalent to setting
/// [`PipelineOptions::emit_ms2`] `true`); the rest of `opts` is honored as-is.
/// mzML-only for the MS2 side: Bruker `.d` / Thermo `.raw` yield an empty
/// `ms2_hills` plus a warning, leaving the MS1 result unaffected.
///
/// MS2 fragment hills are detected with [`KothConfig::ms2_hills`] — `[hills]`
/// plus any `[hills_ms2]` overrides — so an in-process caller drives MS1 and MS2
/// hill settings independently from one `KothConfig`. With no `[hills_ms2]` the
/// MS2 side uses `[hills]` verbatim.
pub fn run_pipeline_with_ms2(
    path: &Path,
    config: &KothConfig,
    opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError> {
    let opts = PipelineOptions {
        emit_ms2: true,
        ..opts.clone()
    };
    let mut out = FeatureFindingOutput::default();
    run_pipeline_streaming(path, config, &opts, &mut out)?;
    Ok(out)
}

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
