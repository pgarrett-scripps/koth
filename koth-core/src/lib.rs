//! koth-core: the in-memory LC-MS feature detection stages of
//! [koth](https://github.com/pgarrett-scripps/koth).
//!
//! This crate reads and writes no files and has no native dependencies. A
//! caller brings its own spectra (for example from its own mzML, Thermo or
//! Bruker reader) as [`Spectrum`] values and runs:
//!
//! 1. **Hill detection** – groups MS1 peaks across scans into chromatographic hills
//! 2. **Feature finding** – groups hills into isotope envelopes with charge state assignment
//! 3. **Scoring** – scores isotope pattern quality using the averagine model
//!
//! either stage by stage ([`run_hills`], [`run_features`], [`run_scoring`]) or
//! in one call ([`run_pipeline_from_spectra`], or the streaming
//! [`run_pipeline_streaming_from_spectra`] with a [`PipelineSink`]).
//!
//! The `koth-ms` crate (the `koth_ff` and `koth_align` executables) adds the
//! file readers and writers and cross-run LFQ, and runs exactly these stages;
//! it re-exports this crate's modules at their original `koth_ms::` paths.
//!
//! # Example
//!
//! ```
//! use koth_core::{run_pipeline_from_spectra, KothConfig, Peak, PipelineOptions, Spectrum};
//!
//! let spectra: Vec<Spectrum> = (0..20)
//!     .map(|i| Spectrum {
//!         scan_index: i,
//!         retention_time: i as f64 * 0.05,
//!         peaks: vec![Peak { mz: 500.0, intensity: 1.0e5, ion_mobility: 0.0 }],
//!         ms_level: 1,
//!         isolation_window: None,
//!         faims_cv: None,
//!     })
//!     .collect();
//!
//! let out = run_pipeline_from_spectra(
//!     spectra.into_iter(),
//!     &KothConfig::default(),
//!     &PipelineOptions::default(),
//! )?;
//! println!("{} hills, {} features", out.hills.len(), out.features.len());
//! # Ok::<(), koth_core::Error>(())
//! ```

pub mod config;
pub mod error;
pub mod features;
pub mod hills;
pub mod mem;
pub mod models;
pub mod pipeline;
pub mod scoring;
pub mod stats;

pub use config::KothConfig;
pub use error::{Error, Result};
pub use models::{Feature, Hill, IsolationWindow, Peak, Polarity, ScoredFeature, Spectrum};
pub use pipeline::{
    group_ms2_hills_by_window, run_pipeline_from_hills, run_pipeline_from_spectra,
    run_pipeline_streaming_from_spectra, FeatureFindingOutput, PipelineOptions, PipelineSink,
};

use std::collections::HashMap;

use config::{FeaturesConfig, FileConfig, HillsConfig, ScoringConfig};
use rand::seq::SliceRandom;

/// Apply the decoy construction to a full spectrum set: when `file.decoy_mode`
/// is set, shuffle the **entire** set (Fisher–Yates over `rand::thread_rng`)
/// before hill detection; otherwise return it untouched.
///
/// This is the single definition of the "shuffle the whole set, then detect"
/// step. Every entry point — in-memory ([`run_hills`]), `koth-ms`'s streaming
/// file reader, and the MS2-buffered pipeline path
/// ([`pipeline::run_pipeline_streaming_from_spectra`]) — routes its decoy shuffle
/// through here so the RNG source and full-set semantics live in exactly one
/// place.
pub fn shuffle_if_decoy(mut spectra: Vec<Spectrum>, file: &FileConfig) -> Vec<Spectrum> {
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

/// Stage 2: Detect isotope features from hills.
pub fn run_features(
    hills: &[Hill],
    config: &FeaturesConfig,
    file: &FileConfig,
) -> Result<Vec<Feature>> {
    if has_multiple_faims_channels(hills.iter().map(|h| h.faims_cv)) {
        let mut groups: HashMap<Option<u32>, (Option<f32>, Vec<Hill>)> = HashMap::new();
        for hill in hills {
            let cv = hills::canonical_faims_cv(hill.faims_cv);
            groups
                .entry(cv.map(f32::to_bits))
                .or_insert_with(|| (cv, Vec::new()))
                .1
                .push(hill.clone());
        }
        let mut groups: Vec<_> = groups.into_values().collect();
        groups.sort_by(|a, b| hills::cmp_faims_cv(a.0, b.0));

        let mut out = Vec::new();
        for (cv, group) in groups {
            log::info!(
                "Detecting features for FAIMS CV {} ({} hills)",
                cv.map_or_else(|| "missing".to_string(), |v| v.to_string()),
                group.len()
            );
            out.extend(run_features_one_channel(&group, config, file)?);
        }
        return Ok(out);
    }
    run_features_one_channel(hills, config, file)
}

fn run_features_one_channel(
    hills: &[Hill],
    config: &FeaturesConfig,
    file: &FileConfig,
) -> Result<Vec<Feature>> {
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
    polarity: Polarity,
) -> Vec<ScoredFeature> {
    if has_multiple_faims_channels(features.iter().map(Feature::faims_cv)) {
        let mut groups: HashMap<Option<u32>, (Option<f32>, Vec<Feature>)> = HashMap::new();
        for feature in features {
            let cv = hills::canonical_faims_cv(feature.faims_cv());
            debug_assert!(
                feature
                    .hills
                    .iter()
                    .all(|h| hills::canonical_faims_cv(h.faims_cv) == cv),
                "feature contains hills from different FAIMS CV channels"
            );
            groups
                .entry(cv.map(f32::to_bits))
                .or_insert_with(|| (cv, Vec::new()))
                .1
                .push(feature.clone());
        }
        let mut groups: Vec<_> = groups.into_values().collect();
        groups.sort_by(|a, b| hills::cmp_faims_cv(a.0, b.0));

        let mut out = Vec::new();
        for (_cv, group) in groups {
            out.extend(run_scoring_one_channel(
                &group,
                scoring_cfg,
                features_cfg,
                polarity,
            ));
        }
        return out;
    }
    run_scoring_one_channel(features, scoring_cfg, features_cfg, polarity)
}

fn run_scoring_one_channel(
    features: &[Feature],
    scoring_cfg: &ScoringConfig,
    features_cfg: &FeaturesConfig,
    polarity: Polarity,
) -> Vec<ScoredFeature> {
    let mut scored = scoring::score_features(
        features,
        scoring_cfg,
        &features_cfg.sulfur_offsets,
        &features_cfg.isotope_model.model(),
        polarity,
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

fn has_multiple_faims_channels(values: impl Iterator<Item = Option<f32>>) -> bool {
    let mut first: Option<Option<u32>> = None;
    for value in values {
        let key = hills::canonical_faims_cv(value).map(f32::to_bits);
        match first {
            None => first = Some(key),
            Some(first_key) if first_key != key => return true,
            Some(_) => {}
        }
    }
    false
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
                faims_cv: None,
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

/// Compiles and runs the README's example as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
