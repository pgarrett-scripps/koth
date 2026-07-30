//! In-process / streaming feature-finding API.
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
//! 1. Spectra stream in one at a time during hill detection (peak memory is
//!    `O(active hills)`, not `O(total peaks)` — the reader is already streaming).
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
//! [`Hill::isolation_window`], and [`group_ms2_hills_by_window`] groups them per
//! DIA channel. MS2 is supported on **mzML**, **Bruker diaPASEF `.d`**
//! (`--features tdf`), and **Thermo DIA `.raw`** (`--features thermo`); DDA
//! acquisitions (ddaPASEF `.d`, DDA `.raw`) degrade to an empty MS2 set + a
//! warning, leaving MS1 untouched. With `emit_ms2` false the MS1 path is
//! byte-for-byte identical to before.

use std::path::Path;

use rand::seq::SliceRandom;

use crate::config::KothConfig;
use crate::error::KothError;
use crate::models::{Feature, Hill, IsolationWindow, ScoredFeature, Spectrum};

/// Knobs that mirror the `koth_ff` binary's optional stages. Everything else is
/// taken from the [`KothConfig`] passed in explicitly.
#[derive(Debug, Clone)]
pub struct PipelineOptions {
    /// Run the averagine isotope-pattern scoring stage. `true` mirrors the
    /// binary's default; `false` mirrors `--no-scoring` (features are still
    /// emitted, wrapped as [`ScoredFeature`] with zeroed scores).
    pub scoring: bool,
    /// Also run MS2 hill detection (DIA, partitioned by precursor isolation
    /// window) over the same input, and deliver the finalized MS2 hill set via
    /// [`PipelineSink::on_ms2_hills`] (and [`FeatureFindingOutput::ms2_hills`]).
    ///
    /// **Default `false`.** When `false`, the MS1 feature path is byte-for-byte
    /// unchanged — no MS2 reader is opened and no MS2 hook fires; this exactly
    /// reproduces the historic MS1-only [`run_pipeline`] behavior.
    ///
    /// MS2 detection mirrors [`crate::run_ms2_hills_streaming`]: mzML, Bruker
    /// diaPASEF `.d` (`--features tdf`), and Thermo DIA `.raw` (`--features
    /// thermo`) are supported; DDA acquisitions (ddaPASEF `.d`, DDA `.raw`) yield
    /// an empty MS2 set plus a warning (the MS1 result is unaffected). MS2
    /// detection is independent of MS1 feature finding — they share only the input.
    pub emit_ms2: bool,
}

impl Default for PipelineOptions {
    fn default() -> Self {
        Self { scoring: true, emit_ms2: false }
    }
}

/// A consumer notified as pipeline stages finalize.
///
/// Both methods have defaults, so a consumer implements only what it needs.
/// The call contract is: `on_hills` **exactly once**, then `on_feature`
/// **once per feature**, in detection order (the same set and order the file
/// writer receives before it applies its own charge>0 filter and
/// intensity-descending sort).
pub trait PipelineSink {
    /// Called once, when hill detection has fully completed. The hills carry
    /// their final `hill_id`s (post-split, post-filter) and are identical to
    /// the rows of `hills.{tsv,parquet}`. Ownership is transferred: feature
    /// structs already hold their own (cheap, `Arc`-backed) hill clones, so the
    /// pipeline no longer needs this vector.
    fn on_hills(&mut self, hills: Vec<Hill>) {
        let _ = hills;
    }

    /// Called once per finalized feature, by value, in detection order.
    fn on_feature(&mut self, feature: ScoredFeature) {
        let _ = feature;
    }

    /// Called **at most once**, with the whole finalized MS2 hill set, by value,
    /// **only when [`PipelineOptions::emit_ms2`] is `true`**. Default no-op.
    ///
    /// MS2 hills are the DIA fragment traces produced by
    /// [`crate::run_ms2_hills_streaming`] — detected independently of MS1, and
    /// partitioned by precursor isolation window. Each hill carries its window
    /// in [`Hill::isolation_window`] (`Some(..)` for every MS2 hill), so a
    /// consumer can map an MS1 precursor feature to the window(s) that isolated
    /// it by testing `window.lower <= precursor_mz <= window.upper`. Group them
    /// per window with [`group_ms2_hills_by_window`].
    ///
    /// Called after `on_hills` and all `on_feature` calls, so a consumer already
    /// holds the MS1 side before the MS2 hills arrive. Not called at all when
    /// `emit_ms2` is `false`, nor (with a warning) when the input is not mzML.
    fn on_ms2_hills(&mut self, hills: Vec<Hill>) {
        let _ = hills;
    }
}

/// Group a finalized MS2 hill set by its precursor isolation window, in a
/// deterministic order (sorted by [`IsolationWindow::key`]).
///
/// MS1 hills (those with `isolation_window == None`) are skipped, so this is
/// safe to pass any hill slice. This is the convenience shape a `koth_tracer`
/// consumer uses to walk one DIA channel at a time; the same information is
/// available per hill via [`Hill::isolation_window`] without grouping.
pub fn group_ms2_hills_by_window(hills: Vec<Hill>) -> Vec<(IsolationWindow, Vec<Hill>)> {
    use std::collections::BTreeMap;
    // BTreeMap keyed on the stable integer window key → deterministic order.
    let mut groups: BTreeMap<(i64, i64, i64), (IsolationWindow, Vec<Hill>)> = BTreeMap::new();
    for h in hills {
        let Some(iw) = h.isolation_window else { continue };
        groups
            .entry(iw.key())
            .or_insert_with(|| (iw, Vec::new()))
            .1
            .push(h);
    }
    groups.into_values().collect()
}

/// Owned, in-memory result of a full run — the convenience shape for consumers
/// (like `uno`) that just want the final lists rather than a streaming callback.
///
/// `features` is exactly the vector the binary hands to its feature writer
/// **before** the writer's own presentation step (which drops `charge == 0`
/// features and sorts by descending intensity). Consumers that want the
/// on-disk row set can replicate that with
/// `features.iter().filter(|f| f.feature.charge > 0)`.
#[derive(Debug, Default)]
pub struct FeatureFindingOutput {
    pub hills: Vec<Hill>,
    pub features: Vec<ScoredFeature>,
    /// Finalized MS2 (DIA fragment) hills, partitioned by precursor isolation
    /// window — populated **only** when [`PipelineOptions::emit_ms2`] is `true`
    /// (empty otherwise). Each carries its window in [`Hill::isolation_window`];
    /// group with [`group_ms2_hills_by_window`].
    pub ms2_hills: Vec<Hill>,
}

/// Wrap an unscored [`Feature`] as a [`ScoredFeature`] with zeroed scores —
/// byte-for-byte what the binary's `--no-scoring` path emits.
pub(crate) fn wrap_unscored(f: Feature) -> ScoredFeature {
    ScoredFeature {
        cosine_score: f.cosine_score,
        feature: f,
        neutron_offset: 0,
        isotope_score: 0.0,
        combined_score: 0.0,
        theoretical_pattern: Vec::new(),
    }
}

/// Core: given a finalized hill set, detect + score features and drive `sink`.
/// Hills are handed to the sink *after* feature detection (which borrows them);
/// features already own their hill clones, so this transfers ownership rather
/// than cloning.
fn drive_from_hills<S: PipelineSink>(
    hills: Vec<Hill>,
    config: &KothConfig,
    opts: &PipelineOptions,
    sink: &mut S,
) -> Result<(), KothError> {
    // Stage 2: isotope-chain features (needs the whole hill set). Reuses the
    // exact crate-level entry the binary calls, including the m/z-recalibration
    // branch keyed on `file.mz_recalibration`.
    let features = crate::run_features(&hills, &config.features, &config.file)?;

    // Stage 3: scoring (optional) — produced before we release the hills so the
    // ordering the writer would see is preserved.
    let scored: Vec<ScoredFeature> = if opts.scoring {
        crate::run_scoring(&features, &config.scoring, &config.features)
    } else {
        features.into_iter().map(wrap_unscored).collect()
    };

    // Hand the finalized hill set to the sink (by value) now that features have
    // been built and own their own hill clones.
    sink.on_hills(hills);

    // Stream features out one at a time, freeing each as it is consumed.
    for sf in scored {
        sink.on_feature(sf);
    }
    Ok(())
}

/// Run the full pipeline from a file path, streaming results to `sink`.
/// No output files are written.
///
/// Hill detection routes exactly as the binary does (streaming mzML with a
/// prefetch reader thread; batch-load for Bruker `.d` / Thermo `.raw`; decoy
/// shuffle when `file.decoy_mode`), via [`crate::run_hills_streaming`].
pub fn run_pipeline_streaming<S: PipelineSink>(
    path: &Path,
    config: &KothConfig,
    opts: &PipelineOptions,
    sink: &mut S,
) -> Result<(), KothError> {
    let hills = crate::run_hills_streaming(path, &config.hills, &config.file)?;
    drive_from_hills(hills, config, opts, sink)?;
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

/// Run the full pipeline from an in-memory iterator of spectra, streaming
/// results to `sink`. No files are written and no `Vec<Spectrum>` is retained
/// (unless `hills.tic_norm_window > 0` or `file.decoy_mode`, which inherently
/// need the whole set — matching the binary).
///
/// This is the file-free entry point: callers that already hold spectra (or
/// synthesize them, as the parity tests do) use this.
pub fn run_pipeline_streaming_from_spectra<I, S>(
    spectra: I,
    config: &KothConfig,
    opts: &PipelineOptions,
    sink: &mut S,
) -> Result<(), KothError>
where
    I: Iterator<Item = Spectrum>,
    S: PipelineSink,
{
    if opts.emit_ms2 {
        // MS2 hills are detected over the same spectra as MS1, but a streaming
        // iterator can only be consumed once — so when MS2 is requested we
        // buffer the input once and run both detectors over it (MS2 is
        // independent of MS1). This mirrors what the file path does with two
        // reads, without a second pass over the source. Detect MS2 from the
        // *unshuffled* buffer (isolation-window hills need RT order); the MS1
        // decoy shuffle, if any, applies only to the MS1 path below.
        let buf: Vec<Spectrum> = spectra.collect();
        // MS2 detector filters to `ms_level == 2` internally, so hand it the
        // whole buffer. `detect_ms2_hills_from_iter` groups by isolation window.
        let ms2 = crate::hills::detect_ms2_hills_from_iter(
            buf.iter().cloned(),
            &config.ms2_hills(),
            &config.file,
        );
        // The MS1 detector does *not* filter by level (it processes every
        // spectrum it is given), so we must exclude MS2 scans here — otherwise
        // fragment peaks would pollute the MS1 hills. On the file path the MS1
        // and MS2 readers are already separate streams; this reproduces that.
        let ms1: Vec<Spectrum> = buf.into_iter().filter(|s| s.ms_level != 2).collect();
        let hills = if config.file.decoy_mode {
            let mut ms1 = ms1;
            ms1.shuffle(&mut rand::thread_rng());
            crate::hills::detect_hills_from_iter(ms1.into_iter(), &config.hills, &config.file)
        } else {
            crate::hills::detect_hills_from_iter(ms1.into_iter(), &config.hills, &config.file)
        };
        drive_from_hills(hills, config, opts, sink)?;
        sink.on_ms2_hills(ms2);
        return Ok(());
    }

    let hills = if config.file.decoy_mode {
        // Decoy mode needs the full set to shuffle before detection, matching
        // `crate::run_hills` / the binary's mzML decoy branch.
        let mut buf: Vec<Spectrum> = spectra.collect();
        buf.shuffle(&mut rand::thread_rng());
        crate::hills::detect_hills_from_iter(buf.into_iter(), &config.hills, &config.file)
    } else {
        crate::hills::detect_hills_from_iter(spectra, &config.hills, &config.file)
    };
    drive_from_hills(hills, config, opts, sink)
}

/// A [`PipelineSink`] that simply collects everything into owned vectors.
#[derive(Debug, Default)]
struct CollectingSink {
    out: FeatureFindingOutput,
}

impl PipelineSink for CollectingSink {
    fn on_hills(&mut self, hills: Vec<Hill>) {
        self.out.hills = hills;
    }
    fn on_feature(&mut self, feature: ScoredFeature) {
        self.out.features.push(feature);
    }
    fn on_ms2_hills(&mut self, hills: Vec<Hill>) {
        self.out.ms2_hills = hills;
    }
}

/// Convenience: run the pipeline from a file path and collect all hills and
/// features into owned vectors. The all-in-memory shape `uno` uses.
pub fn run_pipeline(
    path: &Path,
    config: &KothConfig,
    opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError> {
    let mut sink = CollectingSink::default();
    run_pipeline_streaming(path, config, opts, &mut sink)?;
    Ok(sink.out)
}

/// Convenience: run the pipeline from an in-memory spectrum iterator and
/// collect all hills and features into owned vectors.
pub fn run_pipeline_from_spectra<I>(
    spectra: I,
    config: &KothConfig,
    opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError>
where
    I: Iterator<Item = Spectrum>,
{
    let mut sink = CollectingSink::default();
    run_pipeline_streaming_from_spectra(spectra, config, opts, &mut sink)?;
    Ok(sink.out)
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
    let opts = PipelineOptions { emit_ms2: true, ..opts.clone() };
    let mut sink = CollectingSink::default();
    run_pipeline_streaming(path, config, &opts, &mut sink)?;
    Ok(sink.out)
}

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
