use super::{grid::XicGrid, LfqConfig, ScoreMode};
use crate::scoring::averagine::{bhattacharyya_score, lookup_template};
use crate::scoring::elements::K_PATTERN;

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
    neutral_mass: f64,
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

    // Full-length averagine template (sums ~1 over all K_PATTERN positions),
    // used only by the Bhattacharyya path. Unlike `theory` above — which is
    // trimmed to n_rows and re-normalised, discarding the missed-mass tail —
    // this keeps the full distribution so `bhattacharyya_score` can penalise a
    // column for expected isotope peaks beyond n_rows that are absent.
    let bc_template: [f64; K_PATTERN] = if config.spectral_bhattacharyya {
        lookup_template(neutral_mass)
    } else {
        [0.0; K_PATTERN]
    };

    // Grid-level isotope co-elution scalar (constant across columns), folded
    // into the hybrid below when enabled. Distinct from the per-column spectral
    // term: this asks whether the isotope rows rise and fall *together* across
    // RT. Computed as the theory-weighted mean cosine of the monoisotope row
    // against each higher isotope row that carries signal, over all columns.
    // Rows with no signal are skipped (absence is the spectral term's job, not
    // co-elution's); when fewer than two rows carry signal it defaults to 1.0.
    let coelution: f64 = if config.spectral_coelution {
        isotope_coelution(&grid.intensities, theory, n_rows, n_cols)
    } else {
        1.0
    };

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

        // Spectral score: agreement between the observed isotopologue
        // intensities at this column and the theoretical averagine pattern.
        // Fill the obs scratch buffer in-place to avoid per-column allocation.
        let spec_score = if config.spectral_bhattacharyya {
            // Bhattacharyya path: penalises absent-but-expected isotope peaks
            // via the missed-mass term, so a lone monoisotope no longer scores
            // a free 1.0. Note: no single_isotope short-circuit — the penalty
            // is exactly the behaviour we want for single-peak columns.
            for (r, v) in obs[..n_rows].iter_mut().enumerate() {
                *v = grid.intensities[r][col] as f64;
            }
            bhattacharyya_score(&obs[..n_rows], &bc_template) as f32
        } else if single_isotope {
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
                let base = rt_score as f64 * int_score as f64 * spec_score as f64;
                if config.spectral_coelution {
                    (base * coelution).powf(0.25) as f32
                } else {
                    base.cbrt() as f32
                }
            }
            ScoreMode::Rt => rt_score,
            ScoreMode::Intensity => int_score,
            ScoreMode::Spectral => spec_score,
        };
        scores.hybrid[col] = hybrid;
    }
}

/// Theory-weighted mean cosine of the monoisotope row against each higher
/// isotope row that carries signal, over all RT columns — the grid's
/// isotope-to-isotope co-elution quality in [0, 1]. Rows with no signal are
/// skipped (absence is the spectral term's concern, not co-elution's); returns
/// 1.0 when fewer than two rows carry signal (nothing to co-elute).
fn isotope_coelution(intensities: &[Vec<f32>], theory: &[f64], n_rows: usize, n_cols: usize) -> f64 {
    if n_rows < 2 {
        return 1.0;
    }
    let row0 = &intensities[0];
    let (mut wsum, mut csum) = (0.0f64, 0.0f64);
    for k in 1..n_rows {
        let w = theory.get(k).copied().unwrap_or(0.0);
        if w <= 0.0 {
            continue;
        }
        let rowk = &intensities[k];
        let (mut dot, mut n0, mut nk) = (0.0f64, 0.0f64, 0.0f64);
        for c in 0..n_cols {
            let a = row0[c] as f64;
            let b = rowk[c] as f64;
            dot += a * b;
            n0 += a * a;
            nk += b * b;
        }
        if nk <= 0.0 {
            continue; // M+k absent → not a co-elution failure
        }
        let cos = if n0 > 1e-12 {
            (dot / (n0.sqrt() * nk.sqrt())).clamp(0.0, 1.0)
        } else {
            0.0
        };
        csum += w * cos;
        wsum += w;
    }
    if wsum > 0.0 {
        csum / wsum
    } else {
        1.0
    }
}

#[cfg(test)]
mod spectral_tests {
    use super::*;
    use crate::lfq::grid::XicGrid;

    /// 3 isotope rows × 3 RT bins; only the monoisotope row carries signal,
    /// at the centre column. This is the canonical lone-monoisotope decoy: an
    /// m/z-shifted decoy that happened to grab a single hill.
    fn lone_mono_grid() -> XicGrid {
        let mut g = XicGrid::empty(3, 3, 0.0, 1.0);
        g.intensities[0][1] = 1.0e7;
        g.n_slots_filled = 1;
        g
    }

    /// The cosine path hands a lone monoisotope a free 1.0 spectral score
    /// (`single_isotope` short-circuit); the Bhattacharyya path must instead
    /// penalise it for the absent M+1/M+2 that averagine predicts at this mass.
    /// This pins the loophole closed so the freebie can't silently return.
    #[test]
    fn bhattacharyya_penalises_lone_monoisotope_decoy() {
        let grid = lone_mono_grid();
        let neutral_mass = 1500.0; // averagine predicts substantial M+1/M+2 here
        let theory = vec![0.5, 0.3, 0.2]; // consumed only by the cosine path
        let mut scores = ColumnScores::new(3);
        let mut col_totals = vec![0.0f32; 3];
        let mut obs = vec![0.0f64; 3];

        let mut cfg = LfqConfig::default();
        cfg.score_mode = ScoreMode::Spectral;

        cfg.spectral_bhattacharyya = false;
        score_grid(&grid, &theory, neutral_mass, &cfg, &mut scores, &mut col_totals, &mut obs);
        let cosine_center = scores.spectral[1];

        cfg.spectral_bhattacharyya = true;
        score_grid(&grid, &theory, neutral_mass, &cfg, &mut scores, &mut col_totals, &mut obs);
        let bc_center = scores.spectral[1];

        assert!(
            (cosine_center - 1.0).abs() < 1e-6,
            "cosine freebie should score a lone monoisotope 1.0, got {cosine_center}"
        );
        assert!(
            bc_center > 0.0 && bc_center < 0.85,
            "bhattacharyya should penalise a lone monoisotope (0 < s < 0.85), got {bc_center}"
        );
    }

    /// Co-elution: isotope rows that rise/fall together score ~1; rows peaking
    /// at different columns (interference) score low; a lone row defaults to 1.
    #[test]
    fn coelution_separates_coeluting_from_interfering() {
        let theory = [0.5, 0.3, 0.2];
        // M and M+1 both peak at the centre column → co-elute.
        let coel = vec![
            vec![1.0f32, 9.0, 1.0],
            vec![0.5f32, 4.5, 0.5],
            vec![0.2f32, 1.8, 0.2],
        ];
        let good = isotope_coelution(&coel, &theory, 3, 3);
        // M+1 peaks at a different column than M → interference.
        let interfere = vec![
            vec![9.0f32, 1.0, 0.0],
            vec![0.0f32, 1.0, 9.0],
            vec![0.0f32, 0.0, 0.0],
        ];
        let bad = isotope_coelution(&interfere, &theory, 3, 3);
        assert!(good > 0.95, "co-eluting rows should score ~1, got {good}");
        assert!(bad < 0.5, "interfering rows should score low, got {bad}");
        assert!(good > bad);

        // A lone monoisotope (higher rows empty) cannot be judged → 1.0.
        let lone = vec![vec![0.0f32, 5.0, 0.0], vec![0.0f32; 3], vec![0.0f32; 3]];
        assert_eq!(isotope_coelution(&lone, &theory, 3, 3), 1.0);
    }
}
