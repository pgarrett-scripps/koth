pub mod anchors;
pub mod drift;
pub mod warp;

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::models::{Hill, ScoredFeature};

use anchors::{find_anchors, AnchorPair};
use drift::DriftFit;
use warp::RtWarp;

/// RT-warp model. `Ransac` is the only supported model: it tolerates the ~25%
/// wrong-match contamination in coordinate-only anchors and degrades far less
/// than a sliding-median or affine fit as contamination rises, with fewer
/// tuning knobs.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WarpKind {
    /// RANSAC inlier selection → median piecewise-linear fit on the inliers →
    /// isotonic monotonicity. Robust to a contaminated anchor set (wrong
    /// coordinate-matches): the diagonal is found from a global line consensus
    /// rather than refitting on a contaminated seed.
    Ransac,
}

fn default_warp_kind() -> WarpKind {
    WarpKind::Ransac
}

fn default_ransac_thresh() -> f64 {
    0.01
}

/// Configuration for multi-run alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AlignmentConfig {
    /// Explicit reference run name; None selects the highest-quality run.
    pub reference_run: Option<String>,
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
    /// used by the RANSAC warp's inlier median piecewise-linear fit.
    pub rt_warp_bandwidth: f64,
    /// Which warp model to fit. Only `ransac` (RANSAC-inlier piecewise warp) is
    /// supported.
    #[serde(default = "default_warp_kind")]
    pub rt_warp_kind: WarpKind,
    /// Inlier half-band for the RANSAC warp, in normalised RT units
    /// (0.01 ≈ 1.4 min on a 142-min gradient). An anchor is an inlier if its
    /// reference RT is within this distance of the consensus diagonal.
    #[serde(default = "default_ransac_thresh")]
    pub rt_warp_ransac_thresh: f64,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            reference_run: None,
            anchor_mass_ppm: 10.0,
            rt_anchor_window: 0.05,
            im_tolerance: 0.05,
            min_anchor_combined_score: 0.5,
            min_anchor_count: 10,
            rt_warp_bandwidth: 0.1,
            rt_warp_kind: WarpKind::Ransac,
            rt_warp_ransac_thresh: default_ransac_thresh(),
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
    /// Optional measured RT bounds when hills are streamed and features are absent.
    /// This is metadata, not an indexed scan-time vector.
    pub rt_bounds: Option<(f64, f64)>,
}

impl RunInput {
    /// Overall RT span of this run, derived from scan_times if available.
    pub fn rt_range(&self) -> (f64, f64) {
        if !self.scan_times.is_empty() {
            let min = self
                .scan_times
                .iter()
                .cloned()
                .fold(f64::INFINITY, f64::min);
            let max = self
                .scan_times
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            return (min, max);
        }
        if let Some((min, max)) = self.rt_bounds {
            if min.is_finite() && max.is_finite() && max > min {
                return (min, max);
            }
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
            let min = self
                .hills
                .iter()
                .map(|h| h.rt_start)
                .fold(f64::INFINITY, f64::min);
            let max = self
                .hills
                .iter()
                .map(|h| h.rt_end)
                .fold(f64::NEG_INFINITY, f64::max);
            if min.is_finite() && max.is_finite() && max > min {
                (min, max)
            } else {
                (0.0, 1.0)
            }
        }
    }
}

/// Per-region estimate of the post-warp RT-residual spread σ(run_rt_norm),
/// in normalised RT units. Built from the RANSAC inlier anchors: for each
/// inlier, the residual `ref_rt_norm − (run_rt_norm + warp_delta)` is the RT
/// scatter the warp could not remove. Binned over run-normalised RT with
/// empirical-Bayes shrinkage toward the global σ (so sparse bins don't produce
/// a wild estimate) and a floor. `sigma_norm_at` returns `None` when the model
/// is degenerate (no inliers / zero spread), signalling callers to fall back to
/// the raw RT term.
#[derive(Debug, Clone, Serialize)]
pub struct RtSigmaModel {
    /// Shrunk σ per RT bin, in normalised RT units. Empty ⇒ degenerate.
    bins: Vec<f64>,
    /// Global (all-inlier) σ in normalised RT units.
    global: f64,
}

impl RtSigmaModel {
    /// σ at a run-normalised RT position, or `None` if the model is degenerate.
    pub fn sigma_norm_at(&self, run_norm: f64) -> Option<f64> {
        if self.bins.is_empty() || self.global <= 0.0 {
            return None;
        }
        let nb = self.bins.len();
        let idx = ((run_norm.clamp(0.0, 1.0) * nb as f64) as usize).min(nb - 1);
        Some(self.bins[idx])
    }
}

/// Number of RT bins for the residual-spread model.
const RT_SIGMA_BINS: usize = 12;
/// Empirical-Bayes pseudo-count: a bin's σ is shrunk toward the global σ as if
/// it had this many extra global-variance observations. Larger ⇒ more shrinkage.
const RT_SIGMA_SHRINK_N0: f64 = 20.0;
/// Floor on any bin's σ, as a fraction of the global σ, to avoid an
/// over-peaked RT likelihood in a bin that happens to look near-zero.
const RT_SIGMA_FLOOR_FRAC: f64 = 0.25;

fn build_rt_sigma_model(anchors: &[AnchorPair], active: &[bool], warp: &RtWarp) -> RtSigmaModel {
    // Post-warp residuals of inlier anchors, keyed by run-normalised position.
    let mut xs: Vec<f64> = Vec::new();
    let mut res: Vec<f64> = Vec::new();
    for (a, &ok) in anchors.iter().zip(active.iter()) {
        if !ok {
            continue;
        }
        let pred = a.run_rt_norm + warp.delta_at_run_norm(a.run_rt_norm);
        xs.push(a.run_rt_norm);
        res.push(a.ref_rt_norm - pred);
    }
    if res.len() < 2 {
        return RtSigmaModel {
            bins: Vec::new(),
            global: 0.0,
        };
    }

    let var = |v: &[f64]| -> f64 {
        let n = v.len();
        if n < 2 {
            return 0.0;
        }
        let mean = v.iter().sum::<f64>() / n as f64;
        v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64
    };
    let global_var = var(&res);
    let global = global_var.sqrt();
    if global <= 0.0 {
        return RtSigmaModel {
            bins: Vec::new(),
            global: 0.0,
        };
    }
    let floor = RT_SIGMA_FLOOR_FRAC * global;

    let nb = RT_SIGMA_BINS;
    let mut bin_res: Vec<Vec<f64>> = vec![Vec::new(); nb];
    for (&x, &r) in xs.iter().zip(res.iter()) {
        let idx = ((x.clamp(0.0, 1.0) * nb as f64) as usize).min(nb - 1);
        bin_res[idx].push(r);
    }

    let bins: Vec<f64> = bin_res
        .iter()
        .map(|b| {
            let n = b.len() as f64;
            // Empirical-Bayes shrink the bin variance toward the global variance.
            let shrunk_var =
                (n * var(b) + RT_SIGMA_SHRINK_N0 * global_var) / (n + RT_SIGMA_SHRINK_N0);
            shrunk_var.sqrt().max(floor)
        })
        .collect();

    RtSigmaModel { bins, global }
}

/// Alignment parameters for a single non-reference run.
pub struct RunAlignment {
    pub n_anchors: usize,
    pub rt_warp: RtWarp,
    pub mass_drift: DriftFit,
    pub im_drift: DriftFit,
    pub run_rt_range: (f64, f64),
    pub ref_rt_range: (f64, f64),
    /// Per-region post-warp RT-residual spread, in normalised RT units.
    pub rt_sigma: RtSigmaModel,
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
        self.rt_warp
            .apply(run_rt, self.run_rt_range, self.ref_rt_range)
    }

    /// Apply mass-PPM drift correction to a run mz value.
    /// `run_rt` should be the uncorrected run retention time.
    pub fn correct_mz(&self, mz: f64, run_rt: f64) -> f64 {
        // `mass_drift` is fit with x = ref_rt_norm (see `drift::fit_mass_drift`),
        // so it must be evaluated on the reference-frame normalised RT — the same
        // axis `predict_run_mz` uses. Warp the uncorrected run RT into reference
        // space before normalising, rather than normalising against the run's own
        // range (which would evaluate the fit off its own axis).
        let ref_rt = self.warp_rt(run_rt);
        let rt_norm = anchors::normalize_rt(ref_rt, self.ref_rt_range);
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
/// \>= `config.min_anchor_combined_score`.
pub fn align_runs(runs: &[RunInput], config: &AlignmentConfig) -> AlignmentResult {
    let auto_ref_idx = runs
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

    let ref_idx = config.reference_run.as_ref().map_or(auto_ref_idx, |name| {
        runs.iter()
            .position(|r| &r.name == name)
            .expect("configured reference run must exist")
    });
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
            match config.rt_warp_kind {
                WarpKind::Ransac => warp::fit_ransac_warp(&anchors, config),
            }
        } else {
            log::warn!(
                "Run '{}': only {} anchors (< {}), using identity RT warp",
                run.name,
                anchors.len(),
                config.min_anchor_count
            );
            (warp::identity_warp(), vec![false; anchors.len()])
        };

        let (mass_drift, mass_active) = drift::fit_mass_drift(&anchors);
        let (im_drift, im_active) = drift::fit_im_drift(&anchors);

        let rt_sigma = build_rt_sigma_model(&anchors, &rt_active, &rt_warp);

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
                rt_sigma,
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
