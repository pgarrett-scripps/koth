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
    let cut = config.peak_cut_discrimination as f32;
    // Too few populated columns and a single noisy bin looks like a valley, so
    // short peaks are never cut. FlashLFQ applies the same minimum.
    let populated = col_totals[..n_cols].iter().filter(|&&v| v > 0.0).count();
    let may_cut = cut > 0.0 && populated >= 5;
    let expand = |from: usize, step: isize| -> usize {
        let mut edge = from;
        let mut valley = col_totals[from];
        let mut valley_bin = from;
        for _ in 0..20usize {
            let next = edge as isize + step;
            if next < 0 || next as usize >= n_cols {
                break;
            }
            let candidate = next as usize;
            if scores.bhattacharyya[candidate] < spec_min || scores.hybrid[candidate] < half_max {
                break;
            }
            let height = col_totals[candidate];
            if height < valley {
                valley = height;
                valley_bin = candidate;
            } else if may_cut && height > 0.0 && (height - valley) / height > cut {
                // A rise above the running minimum can be a neighbouring peak or
                // one noisy bin during the descent. Require the rise to clear the
                // column just past the valley as well before cutting, so a single
                // dip cannot truncate a real peak.
                let beyond = valley_bin as isize + step;
                let confirmed = if beyond >= 0 && (beyond as usize) < n_cols {
                    let second = col_totals[beyond as usize];
                    (height - second) / height > cut
                } else {
                    true
                };
                if confirmed {
                    return valley_bin;
                }
            }
            edge = candidate;
        }
        edge
    };
    let start = expand(apex, -1);
    let end = expand(apex, 1);

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
mod tests {
    use super::*;
    use crate::lfq::{grid::XicGrid, score::ColumnScores, LfqConfig};

    /// Two peaks of similar shape separated by a shallow valley. The quality
    /// gates cannot separate them, because a neighbour with a similar envelope
    /// keeps the isotope and hybrid scores high across the valley.
    fn two_peaks() -> (XicGrid, ColumnScores, Vec<f32>) {
        let n = 40;
        let mut grid = XicGrid::empty(1, n, 0.0, 1.0);
        grid.n_slots_filled = 1;
        let mut scores = ColumnScores::new(n);
        let mut totals = vec![0.0f32; n];
        for i in 0..n {
            let x = i as f32;
            let a = 100.0 * (-((x - 10.0) / 3.0).powi(2) / 2.0).exp();
            let b = 70.0 * (-((x - 26.0) / 3.0).powi(2) / 2.0).exp();
            let height = a + b;
            grid.intensities[0][i] = height;
            totals[i] = height;
            // Both peaks look equally good to the quality gates.
            scores.bhattacharyya[i] = 0.9;
            scores.hybrid[i] = 0.9;
            scores.rt[i] = 0.9;
            scores.intensity[i] = 0.9;
        }
        (grid, scores, totals)
    }

    #[test]
    fn expansion_without_a_cut_runs_into_the_neighbouring_peak() {
        let (grid, scores, totals) = two_peaks();
        let config = LfqConfig::default();
        assert_eq!(config.peak_cut_discrimination, 0.0);
        let peak = integrate(&grid, &scores, &totals, &config);
        assert_eq!(peak.apex_bin, 10);
        // Reaches past the valley at bin 18 and into the second peak.
        assert!(peak.end_bin > 20, "end_bin was {}", peak.end_bin);
    }

    #[test]
    fn a_valley_cut_stops_at_the_boundary_between_the_two_peaks() {
        let (grid, scores, totals) = two_peaks();
        let config = LfqConfig {
            peak_cut_discrimination: 0.6,
            ..LfqConfig::default()
        };
        let peak = integrate(&grid, &scores, &totals, &config);
        assert_eq!(peak.apex_bin, 10);
        assert!(
            (14..=22).contains(&peak.end_bin),
            "expected the cut near the valley, got {}",
            peak.end_bin
        );
        let uncut = integrate(&grid, &scores, &totals, &LfqConfig::default());
        assert!(peak.intensity < uncut.intensity);
    }
}
