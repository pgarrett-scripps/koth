use std::collections::HashSet;

use serde::Serialize;

use crate::models::ScoredFeature;

use super::{AlignmentConfig, RunInput};

/// A matched pair of features between a run and the reference.
#[derive(Debug, Clone, Serialize)]
pub struct AnchorPair {
    /// Reference feature RT in normalized [0, 1] space
    pub ref_rt_norm: f64,
    /// Run feature RT in normalized [0, 1] space
    pub run_rt_norm: f64,
    /// Reference monoisotopic mz
    pub ref_mz: f64,
    /// Run monoisotopic mz (uncorrected)
    pub run_mz: f64,
    /// Reference ion mobility (0.0 = not available)
    pub ref_im: f64,
    /// Run ion mobility (0.0 = not available)
    pub run_im: f64,
}

impl AnchorPair {
    pub fn ppm_error(&self) -> f64 {
        (self.run_mz - self.ref_mz) / self.ref_mz * 1e6
    }

    pub fn im_delta(&self) -> f64 {
        self.run_im - self.ref_im
    }
}

pub fn normalize_rt(rt: f64, range: (f64, f64)) -> f64 {
    let (min, max) = range;
    let span = max - min;
    if span.abs() < 1e-9 {
        return 0.5;
    }
    ((rt - min) / span).clamp(0.0, 1.0)
}

/// Find anchor pairs between a query run and the reference run.
///
/// Only high-confidence features (combined_score >= min_anchor_combined_score)
/// from both runs are considered. Matching is 1:1 (greedy by lowest PPM error).
pub fn find_anchors(
    run: &RunInput,
    reference: &RunInput,
    config: &AlignmentConfig,
) -> Vec<AnchorPair> {
    let run_rt_range = run.rt_range();
    let ref_rt_range = reference.rt_range();

    // Build sorted index of high-confidence reference features by (charge, mz)
    let mut ref_sorted: Vec<(usize, &ScoredFeature)> = reference
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            f.combined_score >= config.min_anchor_combined_score && f.feature.charge > 0
        })
        .collect();
    ref_sorted.sort_by(|a, b| {
        a.1.feature.charge.cmp(&b.1.feature.charge).then(
            a.1.monoisotopic_mz()
                .partial_cmp(&b.1.monoisotopic_mz())
                .unwrap(),
        )
    });

    // Same for run features
    let mut run_sorted: Vec<(usize, &ScoredFeature)> = run
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            f.combined_score >= config.min_anchor_combined_score && f.feature.charge > 0
        })
        .collect();
    run_sorted.sort_by(|a, b| {
        a.1.feature.charge.cmp(&b.1.feature.charge).then(
            a.1.monoisotopic_mz()
                .partial_cmp(&b.1.monoisotopic_mz())
                .unwrap(),
        )
    });

    let mut anchors = Vec::new();
    let mut claimed_run: HashSet<usize> = HashSet::new();

    for (_, ref_feat) in &ref_sorted {
        let ref_mz = ref_feat.monoisotopic_mz();
        let ref_rt_norm = normalize_rt(ref_feat.feature.rt_apex(), ref_rt_range);
        let charge = ref_feat.feature.charge;
        let mz_tol = ref_mz * config.anchor_mass_ppm / 1e6;

        // Find start of matching charge+mz band in run_sorted
        let lo = run_sorted.partition_point(|(_, f)| {
            f.feature.charge < charge
                || (f.feature.charge == charge && f.monoisotopic_mz() < ref_mz - mz_tol)
        });

        let mut best: Option<(f64, usize, &ScoredFeature)> = None;

        for (slot, (orig_idx, run_feat)) in run_sorted[lo..].iter().enumerate() {
            if run_feat.feature.charge != charge {
                break;
            }
            let run_mz = run_feat.monoisotopic_mz();
            if run_mz > ref_mz + mz_tol {
                break;
            }
            if claimed_run.contains(orig_idx) {
                continue;
            }
            let run_rt_norm = normalize_rt(run_feat.feature.rt_apex(), run_rt_range);
            if (run_rt_norm - ref_rt_norm).abs() > config.rt_anchor_window {
                continue;
            }
            // Optional IM filter
            let ref_im = ref_feat.feature.im_apex();
            let run_im = run_feat.feature.im_apex();
            if ref_im != 0.0 && run_im != 0.0 && (run_im - ref_im).abs() > config.im_tolerance {
                continue;
            }
            let ppm = (run_mz - ref_mz).abs() / ref_mz * 1e6;
            if best.is_none() || ppm < best.as_ref().unwrap().0 {
                best = Some((ppm, lo + slot, run_feat));
            }
        }

        if let Some((_, run_slot, run_feat)) = best {
            let (orig_idx, _) = run_sorted[run_slot];
            claimed_run.insert(orig_idx);
            anchors.push(AnchorPair {
                ref_rt_norm,
                run_rt_norm: normalize_rt(run_feat.feature.rt_apex(), run_rt_range),
                ref_mz,
                run_mz: run_feat.monoisotopic_mz(),
                ref_im: ref_feat.feature.im_apex(),
                run_im: run_feat.feature.im_apex(),
            });
        }
    }

    anchors
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(ref_mz: f64, run_mz: f64, ref_im: f64, run_im: f64) -> AnchorPair {
        AnchorPair {
            ref_rt_norm: 0.0,
            run_rt_norm: 0.0,
            ref_mz,
            run_mz,
            ref_im,
            run_im,
        }
    }

    #[test]
    fn normalize_rt_maps_range_to_unit_interval() {
        assert_eq!(normalize_rt(10.0, (10.0, 20.0)), 0.0);
        assert_eq!(normalize_rt(15.0, (10.0, 20.0)), 0.5);
        assert_eq!(normalize_rt(20.0, (10.0, 20.0)), 1.0);
    }

    #[test]
    fn normalize_rt_clamps_out_of_range() {
        assert_eq!(normalize_rt(5.0, (10.0, 20.0)), 0.0);
        assert_eq!(normalize_rt(25.0, (10.0, 20.0)), 1.0);
    }

    #[test]
    fn normalize_rt_degenerate_range_is_midpoint() {
        assert_eq!(normalize_rt(42.0, (10.0, 10.0)), 0.5);
    }

    #[test]
    fn ppm_error_is_signed_relative_mass_shift() {
        // run 10 ppm heavy than a 1000 Da reference.
        let a = pair(1000.0, 1000.0 * (1.0 + 10.0 / 1e6), 0.0, 0.0);
        assert!((a.ppm_error() - 10.0).abs() < 1e-9);
        // and signed the other way.
        let b = pair(1000.0, 1000.0 * (1.0 - 4.0 / 1e6), 0.0, 0.0);
        assert!((b.ppm_error() + 4.0).abs() < 1e-9);
    }

    #[test]
    fn im_delta_is_run_minus_ref() {
        let a = pair(1000.0, 1000.0, 1.0, 1.25);
        assert!((a.im_delta() - 0.25).abs() < 1e-12);
    }
}
