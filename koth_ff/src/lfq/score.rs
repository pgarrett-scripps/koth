use super::{grid::XicGrid, LfqConfig, ScoreMode};

/// Per-column scores for an XIC grid.
#[derive(Debug, Clone)]
pub struct ColumnScores {
    pub rt: Vec<f32>,
    pub intensity: Vec<f32>,
    pub spectral: Vec<f32>,
    pub hybrid: Vec<f32>,
}

impl ColumnScores {
    /// Allocate a zeroed scratch instance for reuse across (feature, run) pairs.
    pub fn new(n_cols: usize) -> Self {
        Self {
            rt: vec![0.0f32; n_cols],
            intensity: vec![0.0f32; n_cols],
            spectral: vec![0.0f32; n_cols],
            hybrid: vec![0.0f32; n_cols],
        }
    }
}

/// Score every column of `grid` against the theoretical isotope pattern.
///
/// Results are written into the pre-allocated `scores` buffer.
/// `col_totals` is filled with per-column intensity sums (reused by `integrate`).
/// `obs` is a scratch slice of length ≥ n_rows used for per-column spectral scoring.
///
/// RT score penalises columns far from the centre (bin 50 of 100).
/// Intensity score is sqrt(col_total / max_total).
/// Spectral score is the cosine similarity between the observed isotopologue
/// intensities at that column and the theoretical averagine pattern.
/// Hybrid = cbrt(rt × intensity × spectral).
/// Note: apex detection uses raw col_total (not hybrid), so the RT score
/// affects only peak-expansion gating and TDC ranking, not apex selection.
pub fn score_grid(
    grid: &XicGrid,
    theoretical_pattern: &[f64],
    config: &LfqConfig,
    scores: &mut ColumnScores,
    col_totals: &mut [f32],
    obs: &mut [f64],
) {
    let n_cols = grid.n_cols();
    let n_rows = grid.n_rows();
    let center = (n_cols as f64 - 1.0) / 2.0;

    // Fill col_totals; reused by integrate() for apex finding.
    for c in 0..n_cols {
        col_totals[c] = grid.intensities.iter().map(|row| row[c]).sum();
    }
    let max_total = col_totals[..n_cols].iter().cloned().fold(0.0f32, f32::max);

    // Theoretical pattern trimmed to n_rows and L1-normalised.
    // norm_theory is pre-computed once to avoid recomputing it per column.
    let mut theory = [0.0f64; 8]; // 8 isotopes is well above any practical limit
    let theory = &mut theory[..n_rows];
    {
        let mut sum = 0.0f64;
        for (i, &v) in theoretical_pattern.iter().take(n_rows).enumerate() {
            theory[i] = v;
            sum += v;
        }
        if sum > 0.0 {
            for v in theory.iter_mut() {
                *v /= sum;
            }
        } else {
            let uniform = 1.0 / n_rows as f64;
            for v in theory.iter_mut() {
                *v = uniform;
            }
        }
    }
    let norm_theory: f64 = theory.iter().map(|x| x * x).sum::<f64>().sqrt();

    let single_isotope = n_rows == 1 || grid.n_slots_filled <= 1;

    for col in 0..n_cols {
        // RT score: 1 - cbrt(|col - center| / center)
        let rt_dev = if center > 0.0 {
            (col as f64 - center).abs() / center
        } else {
            0.0
        };
        let rt_score = (1.0 - rt_dev.cbrt()).max(0.0) as f32;
        scores.rt[col] = rt_score;

        // Intensity score: sqrt(col_total / max_total)
        let int_score = if max_total > 0.0 {
            (col_totals[col] / max_total).sqrt()
        } else {
            0.0
        };
        scores.intensity[col] = int_score;

        // Spectral score: cosine similarity between observed and theoretical pattern.
        // Fill the obs scratch buffer in-place to avoid per-column allocation.
        let spec_score = if single_isotope {
            1.0f32
        } else {
            for (r, v) in obs[..n_rows].iter_mut().enumerate() {
                *v = grid.intensities[r][col] as f64;
            }
            let dot: f64 = obs[..n_rows].iter().zip(theory.iter()).map(|(x, y)| x * y).sum();
            let norm_a: f64 = obs[..n_rows].iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm_a < 1e-12 || norm_theory < 1e-12 {
                0.0f32
            } else {
                (dot / (norm_a * norm_theory)).clamp(0.0, 1.0) as f32
            }
        };
        scores.spectral[col] = spec_score;

        // Hybrid / mode-selected score.
        let hybrid = match config.score_mode {
            ScoreMode::Hybrid => {
                let p = rt_score as f64 * int_score as f64 * spec_score as f64;
                p.cbrt() as f32
            }
            ScoreMode::Rt => rt_score,
            ScoreMode::Intensity => int_score,
            ScoreMode::Spectral => spec_score,
        };
        scores.hybrid[col] = hybrid;
    }
}
