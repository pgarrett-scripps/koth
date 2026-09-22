use super::{grid::XicGrid, score::ColumnScores, LfqConfig};

/// Result of integrating a single peak from an XIC grid.
#[derive(Debug, Clone)]
pub struct PeakResult {
    /// Summed intensity over all isotopologue rows within [start_bin, end_bin]
    pub intensity: f64,
    /// Peak height: total isotopologue intensity at the apex column only
    /// (the MBR/grid analog of the feature finder's `intensityApex`). Used by
    /// the `apex` quant estimator to avoid integration-window variance.
    pub apex_intensity: f64,
    pub apex_bin: usize,
    pub start_bin: usize,
    pub end_bin: usize,
    pub hybrid_score: f32,
    /// Spectral Bhattacharyya at the apex column: observed-vs-theoretical isotope
    /// pattern match (penalises missing peaks). NOT a cosine — the inter-isotope
    /// cosine is `coelution` below.
    pub bhattacharyya_at_apex: f32,
    /// Individual hybrid components at the apex column, exposed so a downstream
    /// rescorer can weight them independently instead of using the composite.
    pub rt_score_at_apex: f32,
    pub int_score_at_apex: f32,
    pub coelution: f32,
}

impl PeakResult {
    pub fn empty() -> Self {
        Self {
            intensity: 0.0,
            apex_intensity: 0.0,
            apex_bin: 0,
            start_bin: 0,
            end_bin: 0,
            hybrid_score: 0.0,
            bhattacharyya_at_apex: 0.0,
            rt_score_at_apex: 0.0,
            int_score_at_apex: 0.0,
            coelution: 1.0,
        }
    }
}

/// Find the best peak in the grid and integrate it.
///
/// `col_totals` must be the pre-computed per-column intensity sums from `score_grid`;
/// they are reused here for apex finding to avoid recomputing them.
///
/// Algorithm:
/// 1. Locate the apex column (highest raw intensity sum).
/// 2. Expand left and right, stopping when any of:
///    - 20 bins have been added on that side,
///    - the spectral Bhattacharyya drops below `min_spectral_bhattacharyya`,
///    - the hybrid score drops below half the apex score.
/// 3. Sum all intensities in [start_bin, end_bin] across all isotopologue rows.
pub fn integrate(
    grid: &XicGrid,
    scores: &ColumnScores,
    col_totals: &[f32],
    config: &LfqConfig,
) -> PeakResult {
    if grid.is_empty() {
        return PeakResult::empty();
    }
    let n_cols = grid.n_cols();
    if n_cols == 0 {
        return PeakResult::empty();
    }

    // Find apex column by raw intensity (unbiased by RT-score penalty).
    // col_totals was computed by score_grid — no recomputation needed.
    let apex = col_totals[..n_cols]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(n_cols / 2);

    // Gate on raw signal, not quality score — a peak at the grid edge can have
    // rt_score=0 and hybrid=0 despite valid signal, which would falsely zero it out.
    if col_totals[apex] <= 0.0 {
        return PeakResult::empty();
    }

    let max_score = scores.hybrid[apex];
    // half_max for expansion: if hybrid is zero (edge peak), fall back to
    // a nominal threshold so expansion still uses min_spectral_bhattacharyya as gate.
    let half_max = if max_score > 0.0 {
        max_score * 0.5
    } else {
        0.0
    };
    let spec_min = config.min_spectral_bhattacharyya as f32;

    // Expand outward from the apex. The quality gates below stop expansion when
    // the isotope pattern or hybrid score decays, but they do not notice a
    // second peak whose envelope resembles this one: the profile falls into a
    // valley, climbs back out, and expansion follows it. `cut` watches the
    // running minimum and stops at the valley once the profile rises out of it
    // by `peak_cut_discrimination` of the current height.
    // Bound the expansion in window terms, not bins, so the integration span
    // does not change when the grid resolution does.
    let max_steps = ((config.peak_max_halfwidth_frac * n_cols as f64 / 2.0).ceil() as usize).max(1);

    // Expand left
    let mut start = apex;
    let mut left_steps = 0usize;
    while start > 0 && left_steps < max_steps {
        let candidate = start - 1;
        if scores.bhattacharyya[candidate] < spec_min {
            break;
        }
        if scores.hybrid[candidate] < half_max {
            break;
        }
        start = candidate;
        left_steps += 1;
    }

    // Expand right
    let mut end = apex;
    let mut right_steps = 0usize;
    while end + 1 < n_cols && right_steps < max_steps {
        let candidate = end + 1;
        if scores.bhattacharyya[candidate] < spec_min {
            break;
        }
        if scores.hybrid[candidate] < half_max {
            break;
        }
        end = candidate;
        right_steps += 1;
    }

    // Reported cell intensity. Default: raw box-sum over [start, end] across all
    // isotopologue rows. When `config.averagine_projection` is set: the summed
    // averagine-projected (on-pattern) intensity instead, which `score_grid`
    // precomputed per column. Apex selection and window expansion above stay on
    // the raw signal / hybrid score either way — only the reported value changes.
    let (intensity, apex_intensity): (f64, f64) = if config.averagine_projection {
        let sum: f64 = scores.projected[start..=end]
            .iter()
            .map(|&v| v as f64)
            .sum();
        (sum, scores.projected[apex] as f64)
    } else {
        let sum: f64 = grid
            .intensities
            .iter()
            .flat_map(|row| row[start..=end].iter())
            .map(|&v| v as f64)
            .sum();
        (sum, col_totals[apex] as f64)
    };

    PeakResult {
        intensity,
        apex_intensity,
        apex_bin: apex,
        start_bin: start,
        end_bin: end,
        hybrid_score: max_score,
        bhattacharyya_at_apex: scores.bhattacharyya[apex],
        rt_score_at_apex: scores.rt[apex],
        int_score_at_apex: scores.intensity[apex],
        coelution: scores.coelution,
    }
}

#[cfg(test)]
mod scan_density_tests {
    use crate::models::Hill;
    use std::sync::Arc;

    fn hill(rt_start: f64, rt_end: f64, scan_start: usize, scan_end: usize) -> Hill {
        Hill {
            hill_id: 0,
            mz: 500.0,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: (rt_start + rt_end) / 2.0,
            rt_start,
            rt_end,
            rt_width: rt_end - rt_start,
            im: 0.0,
            im_std: 0.0,
            scan_start,
            scan_apex: scan_start,
            scan_end,
            n_scans: scan_end.saturating_sub(scan_start) + 1,
            skipped_scans: 0,
            intensity_sum: 1.0,
            intensity_max: 1.0,
            hill_score: 1.0,
            intensity_profile: Arc::from(&[1.0f32][..]),
            isolation_window: None,
            faims_cv: None,
        }
    }

    #[test]
    fn scan_spacing_is_the_median_over_hills_that_span_several_scans() {
        // Ten scans across ten seconds is one second per scan.
        let hills: Vec<Hill> = (0..9).map(|_| hill(0.0, 10.0 / 60.0, 0, 10)).collect();
        let spacing = crate::lfq::median_scan_spacing(&hills).expect("usable hills");
        assert!(
            (spacing * 60.0 - 1.0).abs() < 1e-6,
            "got {} s",
            spacing * 60.0
        );
    }

    #[test]
    fn no_usable_timing_leaves_the_configured_count_in_place() {
        // A hill confined to one scan cannot report a spacing.
        assert!(crate::lfq::median_scan_spacing(&[hill(0.0, 5.0, 3, 3)]).is_none());
        assert!(crate::lfq::median_scan_spacing(&[]).is_none());
    }
}
