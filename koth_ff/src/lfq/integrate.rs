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
    pub spectral_cosine_at_apex: f32,
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
            spectral_cosine_at_apex: 0.0,
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
///    - the spectral cosine drops below `spectral_cosine_min`,
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
    // a nominal threshold so expansion still uses spectral_cosine_min as gate.
    let half_max = if max_score > 0.0 { max_score * 0.5 } else { 0.0 };
    let spec_min = config.spectral_cosine_min as f32;

    // Expand left
    let mut start = apex;
    let mut left_steps = 0usize;
    while start > 0 && left_steps < 20 {
        let candidate = start - 1;
        if scores.spectral[candidate] < spec_min {
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
    while end + 1 < n_cols && right_steps < 20 {
        let candidate = end + 1;
        if scores.spectral[candidate] < spec_min {
            break;
        }
        if scores.hybrid[candidate] < half_max {
            break;
        }
        end = candidate;
        right_steps += 1;
    }

    // Integrate across all isotopologue rows in [start, end]
    let intensity: f64 = grid
        .intensities
        .iter()
        .flat_map(|row| row[start..=end].iter())
        .map(|&v| v as f64)
        .sum();

    PeakResult {
        intensity,
        apex_intensity: col_totals[apex] as f64,
        apex_bin: apex,
        start_bin: start,
        end_bin: end,
        hybrid_score: max_score,
        spectral_cosine_at_apex: scores.spectral[apex],
    }
}
