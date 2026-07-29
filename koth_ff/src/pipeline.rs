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
//!   then consumes features incrementally in `on_feature`.
//! * **`uno`** (MS1 search) just wants the final feature list: it calls
//!   [`run_pipeline`] / [`run_pipeline_from_spectra`] and reads
//!   [`FeatureFindingOutput::features`].

use std::path::Path;

use rand::seq::SliceRandom;

use crate::config::KothConfig;
use crate::error::KothError;
use crate::models::{Feature, Hill, ScoredFeature, Spectrum};

/// Knobs that mirror the `koth_ff` binary's optional stages. Everything else is
/// taken from the [`KothConfig`] passed in explicitly.
#[derive(Debug, Clone)]
pub struct PipelineOptions {
    /// Run the averagine isotope-pattern scoring stage. `true` mirrors the
    /// binary's default; `false` mirrors `--no-scoring` (features are still
    /// emitted, wrapped as [`ScoredFeature`] with zeroed scores).
    pub scoring: bool,
}

impl Default for PipelineOptions {
    fn default() -> Self {
        Self { scoring: true }
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
    drive_from_hills(hills, config, opts, sink)
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

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
