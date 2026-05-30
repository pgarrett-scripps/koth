pub mod anchors;
pub mod drift;
pub mod warp;

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::models::{Hill, ScoredFeature};

use anchors::{find_anchors, AnchorPair};
use drift::DriftFit;
use warp::RtWarp;

/// Configuration for multi-run alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentConfig {
    /// PPM tolerance for finding anchor feature pairs
    pub anchor_mass_ppm: f64,
    /// Normalised RT window [0, 1] for anchor matching (e.g. 0.05 = ±5% of gradient)
    pub rt_anchor_window: f64,
    /// Ion mobility tolerance for anchor matching (1/K0 units)
    pub im_tolerance: f64,
    /// Minimum `combined_score` (isotope × chromato cosine) for a feature to
    /// be eligible as an alignment anchor.
    pub min_anchor_combined_score: f64,
    /// Minimum anchor count required to fit a warp; falls back to identity if below
    pub min_anchor_count: usize,
    /// Bandwidth for the sliding-window median (fraction of normalised RT range)
    pub rt_warp_bandwidth: f64,
    /// Sigma threshold for anchor sigma-clipping
    pub rt_warp_sigma_clip: f64,
    /// Number of sigma-clip iterations
    pub rt_warp_clip_iters: usize,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            anchor_mass_ppm: 10.0,
            rt_anchor_window: 0.05,
            im_tolerance: 0.05,
            min_anchor_combined_score: 0.5,
            min_anchor_count: 10,
            rt_warp_bandwidth: 0.1,
            rt_warp_sigma_clip: 3.0,
            rt_warp_clip_iters: 5,
        }
    }
}

/// All input data for a single LC-MS run.
pub struct RunInput {
    pub name: String,
    /// Scored features — used for anchor finding and as consensus source for the reference run
    pub features: Vec<ScoredFeature>,
    /// All hills for this run — used for LFQ grid extraction
    pub hills: Vec<Hill>,
    /// Retention time (minutes) indexed by absolute scan index.
    /// If empty, hill RT positions are interpolated from rt_start / rt_end.
    pub scan_times: Vec<f64>,
}

impl RunInput {
    /// Overall RT span of this run, derived from scan_times if available.
    pub fn rt_range(&self) -> (f64, f64) {
        if !self.scan_times.is_empty() {
            let min = self.scan_times.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = self.scan_times.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            return (min, max);
        }
        // Derive from feature RT extents as fallback
        let min = self
            .features
            .iter()
            .map(|f| f.feature.rt_start())
            .fold(f64::INFINITY, f64::min);
        let max = self
            .features
            .iter()
            .map(|f| f.feature.rt_end())
            .fold(f64::NEG_INFINITY, f64::max);
        if min.is_finite() && max.is_finite() && max > min {
            (min, max)
        } else {
            (0.0, 1.0)
        }
    }
}

/// Alignment parameters for a single non-reference run.
pub struct RunAlignment {
    pub n_anchors: usize,
    pub rt_warp: RtWarp,
    pub mass_drift: DriftFit,
    pub im_drift: DriftFit,
    pub run_rt_range: (f64, f64),
    pub ref_rt_range: (f64, f64),
    /// All anchor pairs used for this run's alignment.
    /// Retained for diagnostic reporting; aligned 1:1 with the masks below.
    pub anchors: Vec<AnchorPair>,
    /// Per-anchor "kept" mask from the RT-warp sigma-clip iterations.
    pub rt_active: Vec<bool>,
    /// Per-anchor "kept" mask from the mass-drift sigma-clip iterations.
    pub mass_active: Vec<bool>,
    /// Per-anchor "kept" mask from the IM-drift fit.
    /// `false` for anchors that lacked IM on either side OR were sigma-clipped.
    pub im_active: Vec<bool>,
}

impl RunAlignment {
    /// Map an absolute run RT → absolute reference RT.
    pub fn warp_rt(&self, run_rt: f64) -> f64 {
        self.rt_warp.apply(run_rt, self.run_rt_range, self.ref_rt_range)
    }

    /// Apply mass-PPM drift correction to a run mz value.
    /// `run_rt` should be the uncorrected run retention time.
    pub fn correct_mz(&self, mz: f64, run_rt: f64) -> f64 {
        let rt_norm = anchors::normalize_rt(run_rt, self.run_rt_range);
        let ppm = self.mass_drift.predict(rt_norm);
        mz / (1.0 + ppm / 1e6)
    }

    /// Apply ion-mobility drift correction.
    pub fn correct_im(&self, im: f64, run_rt: f64) -> f64 {
        let rt_norm = anchors::normalize_rt(run_rt, self.run_rt_range);
        im - self.im_drift.predict(rt_norm)
    }

    /// Predict the expected RT of a reference feature in this run's native RT space.
    /// Uses an approximate inverse warp (valid for the small-to-moderate shifts typical in LC-MS).
    pub fn predict_run_rt(&self, ref_rt: f64) -> f64 {
        let ref_norm = anchors::normalize_rt(ref_rt, self.ref_rt_range);
        let delta = self.rt_warp.delta_at_ref_norm(ref_norm);
        let run_norm = (ref_norm - delta).clamp(0.0, 1.0);
        run_norm * (self.run_rt_range.1 - self.run_rt_range.0) + self.run_rt_range.0
    }

    /// Predict the expected m/z of a reference feature in this run's native m/z space.
    pub fn predict_run_mz(&self, ref_mz: f64, ref_rt: f64) -> f64 {
        let rt_norm = anchors::normalize_rt(ref_rt, self.ref_rt_range);
        let ppm = self.mass_drift.predict(rt_norm);
        ref_mz * (1.0 + ppm / 1e6)
    }

    /// Predict the expected ion mobility of a reference feature in this run's native IM space.
    pub fn predict_run_im(&self, ref_im: f64, ref_rt: f64) -> f64 {
        let rt_norm = anchors::normalize_rt(ref_rt, self.ref_rt_range);
        ref_im + self.im_drift.predict(rt_norm)
    }
}

/// Full alignment result for a set of runs.
pub struct AlignmentResult {
    /// Index of the chosen reference run in the input slice
    pub reference_idx: usize,
    pub reference_name: String,
    /// Alignments for every non-reference run, keyed by run name
    pub alignments: HashMap<String, RunAlignment>,
}

/// Compute alignment parameters for all runs relative to an auto-selected reference.
///
/// The reference is the run with the most features whose `combined_score`
/// >= `config.min_anchor_combined_score`.
pub fn align_runs(runs: &[RunInput], config: &AlignmentConfig) -> AlignmentResult {
    let ref_idx = runs
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let n = r
                .features
                .iter()
                .filter(|f| f.combined_score >= config.min_anchor_combined_score)
                .count();
            (i, n)
        })
        .max_by_key(|(_, n)| *n)
        .map(|(i, _)| i)
        .unwrap_or(0);

    let reference = &runs[ref_idx];
    let ref_rt_range = reference.rt_range();

    let mut alignments = HashMap::new();

    for (i, run) in runs.iter().enumerate() {
        if i == ref_idx {
            continue;
        }

        let run_rt_range = run.rt_range();
        let anchors = find_anchors(run, reference, config);

        log::info!(
            "Alignment: run '{}' → '{}': {} anchors",
            run.name,
            reference.name,
            anchors.len()
        );

        let (rt_warp, rt_active) = if anchors.len() >= config.min_anchor_count {
            warp::fit_rt_warp(&anchors, config)
        } else {
            log::warn!(
                "Run '{}': only {} anchors (< {}), using identity RT warp",
                run.name,
                anchors.len(),
                config.min_anchor_count
            );
            (warp::identity_warp(), vec![false; anchors.len()])
        };

        let (mass_drift, mass_active) = drift::fit_mass_drift(&anchors, config);
        let (im_drift, im_active) = drift::fit_im_drift(&anchors, config);

        let n_anchors = anchors.len();
        alignments.insert(
            run.name.clone(),
            RunAlignment {
                n_anchors,
                rt_warp,
                mass_drift,
                im_drift,
                run_rt_range,
                ref_rt_range,
                anchors,
                rt_active,
                mass_active,
                im_active,
            },
        );
    }

    AlignmentResult {
        reference_idx: ref_idx,
        reference_name: reference.name.clone(),
        alignments,
    }
}
