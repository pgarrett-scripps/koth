use super::{grid::XicGrid, LfqConfig, ScoreMode};
use crate::scoring::averagine::bhattacharyya_score;
use crate::scoring::elements::K_PATTERN;

/// Per-column scores for an XIC grid.
#[derive(Debug, Clone)]
pub struct ColumnScores {
    pub rt: Vec<f32>,
    pub intensity: Vec<f32>,
    /// **Bhattacharyya** coefficient per column: observed isotope abundances vs
    /// the theoretical averagine pattern, with a penalty for expected-but-missing
    /// peaks. This is the "does the isotope pattern match theory?" signal.
    pub bhattacharyya: Vec<f32>,
    pub hybrid: Vec<f32>,
    /// **Averagine-projected intensity** per column: the observed isotopologue
    /// vector passed through a matched filter for the L2-normalised theoretical
    /// pattern (`<observed, pattern_hat>`) — the on-pattern signal component,
    /// with the part orthogonal to the fingerprint rejected. Only filled when
    /// `LfqConfig.averagine_projection` is set; left zeroed otherwise (the raw
    /// box-sum path never reads it). Summed over the integration window by
    /// `integrate` to form the reported cell intensity.
    pub projected: Vec<f32>,
    /// **Co-elution cosine** (grid-level scalar): cosine similarity between the
    /// monoisotope's XIC trace and each higher isotope's XIC trace across RT —
    /// i.e. "do the matched isotopes rise and fall together?". Orthogonal to
    /// the Bhattacharyya pattern match. 1.0 when <2 isotope rows carry signal.
    pub coelution: f32,
}

impl ColumnScores {
    /// Allocate a zeroed scratch instance for reuse across (feature, run) pairs.
    pub fn new(n_cols: usize) -> Self {
        Self {
            rt: vec![0.0f32; n_cols],
            intensity: vec![0.0f32; n_cols],
            bhattacharyya: vec![0.0f32; n_cols],
            hybrid: vec![0.0f32; n_cols],
            projected: vec![0.0f32; n_cols],
            coelution: 1.0,
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
    // Full-length averagine template (sums ~1 over all K_PATTERN positions) for
    // the Bhattacharyya score. Depends only on the consensus feature's neutral
    // mass, so it is precomputed once per feature by the caller rather than
    // rebuilt on every (run, target/decoy) invocation of this hot function.
    bc_template: &[f64; K_PATTERN],
    config: &LfqConfig,
    scores: &mut ColumnScores,
    col_totals: &mut [f32],
    obs: &mut [f64],
    // When `Some(σ)`, the RT term becomes a Gaussian likelihood `exp(−½ z²)`
    // with `z = (col − center) / σ` (σ in grid-column units); when `None`,
    // the raw `1 − cbrt(rt_dev)` term is used. See `LfqConfig.rt_spread_scoring`.
    rt_sigma_cols: Option<f32>,
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
    // `bc_template` (passed in) keeps the full averagine distribution, unlike
    // `theory` above which is trimmed to n_rows and re-normalised, so
    // `bhattacharyya_score` can penalise a column for expected isotope peaks
    // beyond n_rows that are absent.

    // Co-elution cosine (always computed, grid-level scalar): the theory-weighted
    // mean cosine of the monoisotope XIC trace against each higher isotope's XIC
    // trace across RT — "do the matched isotopes rise and fall together?". This
    // is the true inter-isotope cosine; it is orthogonal to the Bhattacharyya
    // pattern match. Rows with no signal are skipped; <2 signal rows → 1.0.
    let coelution: f64 =
        isotope_coelution(&grid.intensities, theory, n_rows, n_cols, config.lone_coelution);
    scores.coelution = coelution as f32;

    // Averagine projection (opt-in, LfqConfig.averagine_projection): L2 norm of
    // the (L1-normalised) theory pattern, precomputed once so each column's
    // observed vector can be projected onto the pattern direction below. A zero
    // norm (degenerate/empty pattern) leaves the projection disabled and the
    // raw box-sum path untouched.
    let theory_l2: f64 = if config.averagine_projection {
        theory.iter().map(|&t| t * t).sum::<f64>().sqrt()
    } else {
        0.0
    };

    for col in 0..n_cols {
        // RT score: 1 - cbrt(|col - center| / center)
        let rt_dev = if center > 0.0 {
            (col as f64 - center).abs() / center
        } else {
            0.0
        };
        let rt_score = match rt_sigma_cols {
            Some(sigma) if sigma > 0.0 => {
                // σ-normalised Gaussian RT likelihood, region-aware.
                let z = (col as f64 - center) / sigma as f64;
                (-0.5 * z * z).exp() as f32
            }
            // Raw closeness term (original behaviour).
            _ => (1.0 - rt_dev.cbrt()).max(0.0) as f32,
        };
        scores.rt[col] = rt_score;

        // Intensity score: sqrt(col_total / max_total)
        let int_score = if max_total > 0.0 {
            (col_totals[col] / max_total).sqrt()
        } else {
            0.0
        };
        scores.intensity[col] = int_score;

        // Spectral pattern match: Bhattacharyya coefficient of the observed
        // isotope abundances vs the theoretical averagine pattern (penalises
        // expected-but-missing peaks). This is the ONLY observed-vs-theory
        // metric — a redundant observed-vs-theory cosine was removed; the
        // cosine that matters is the inter-isotope co-elution term above.
        for (r, v) in obs[..n_rows].iter_mut().enumerate() {
            *v = grid.intensities[r][col] as f64;
        }
        let bhattacharyya = bhattacharyya_score(&obs[..n_rows], bc_template) as f32;
        scores.bhattacharyya[col] = bhattacharyya;

        // Averagine-projected intensity (opt-in): keep the component of this
        // column's observed isotope vector that lies along the pattern
        // direction (`<observed, pattern_hat>`), discarding the orthogonal
        // off-pattern part (chemical noise, co-isobars, lone monoisotopes).
        // `obs[..n_rows]` was just filled above; `theory` is L1-normalised.
        if theory_l2 > 0.0 {
            let dot: f64 = obs[..n_rows]
                .iter()
                .zip(theory.iter())
                .map(|(&o, &t)| o * t)
                .sum();
            scores.projected[col] = (dot / theory_l2) as f32;
        }

        // Hybrid = geometric mean of the four quality signals:
        //   rt × intensity × bhattacharyya(pattern) × coelution(inter-isotope cosine).
        let hybrid = match config.score_mode {
            ScoreMode::Hybrid => {
                let base = rt_score as f64 * int_score as f64 * bhattacharyya as f64 * coelution;
                base.powf(0.25) as f32
            }
            ScoreMode::Rt => rt_score,
            ScoreMode::Intensity => int_score,
            ScoreMode::Spectral => bhattacharyya,
        };
        scores.hybrid[col] = hybrid;
    }
}

/// Theory-weighted mean cosine of the monoisotope row against each higher
/// isotope row that carries signal, over all RT columns — the grid's
/// isotope-to-isotope co-elution quality in [0, 1]. Rows with no signal are
/// skipped (absence is the spectral term's concern, not co-elution's); returns
/// 1.0 when fewer than two rows carry signal (nothing to co-elute).
fn isotope_coelution(
    intensities: &[Vec<f32>],
    theory: &[f64],
    n_rows: usize,
    n_cols: usize,
    lone_value: f64,
) -> f64 {
    if n_rows < 2 {
        return lone_value;
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
        lone_value
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
        let bc_template = crate::scoring::averagine::lookup_template(neutral_mass);
        let theory = vec![0.5, 0.3, 0.2];
        let mut scores = ColumnScores::new(3);
        let mut col_totals = vec![0.0f32; 3];
        let mut obs = vec![0.0f64; 3];

        let mut cfg = LfqConfig::default();
        cfg.score_mode = ScoreMode::Spectral;

        score_grid(&grid, &theory, &bc_template, &cfg, &mut scores, &mut col_totals, &mut obs, None);
        let bc_center = scores.bhattacharyya[1];

        // The pattern match must penalise a lone monoisotope for the absent
        // M+1/M+2 that averagine predicts — no free perfect score.
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
        let good = isotope_coelution(&coel, &theory, 3, 3, 1.0);
        // M+1 peaks at a different column than M → interference.
        let interfere = vec![
            vec![9.0f32, 1.0, 0.0],
            vec![0.0f32, 1.0, 9.0],
            vec![0.0f32, 0.0, 0.0],
        ];
        let bad = isotope_coelution(&interfere, &theory, 3, 3, 1.0);
        assert!(good > 0.95, "co-eluting rows should score ~1, got {good}");
        assert!(bad < 0.5, "interfering rows should score low, got {bad}");
        assert!(good > bad);

        // A lone monoisotope (higher rows empty) cannot be judged → 1.0.
        let lone = vec![vec![0.0f32, 5.0, 0.0], vec![0.0f32; 3], vec![0.0f32; 3]];
        assert_eq!(isotope_coelution(&lone, &theory, 3, 3, 1.0), 1.0);
        // A4: the lone-cell freebie is configurable — a penalised value flows through.
        assert_eq!(isotope_coelution(&lone, &theory, 3, 3, 0.5), 0.5);
    }

    use crate::lfq::integrate::integrate;

    /// Build a 3-isotope × 3-column grid carrying signal only in the centre
    /// column (so the integration window is exactly the apex — no expansion).
    fn centre_grid(mono: f32, m1: f32, m2: f32) -> XicGrid {
        let mut g = XicGrid::empty(3, 3, 0.0, 1.0);
        g.intensities[0][1] = mono;
        g.intensities[1][1] = m1;
        g.intensities[2][1] = m2;
        g.n_slots_filled = 3;
        g
    }

    fn integrate_grid(grid: &XicGrid, theory: &[f64], projection: bool) -> f64 {
        let bc_template = crate::scoring::averagine::lookup_template(1500.0);
        let mut cfg = LfqConfig::default();
        cfg.averagine_projection = projection;
        let mut scores = ColumnScores::new(3);
        let mut col_totals = vec![0.0f32; 3];
        let mut obs = vec![0.0f64; 3];
        score_grid(grid, theory, &bc_template, &cfg, &mut scores, &mut col_totals, &mut obs, None);
        integrate(grid, &scores, &col_totals, &cfg).intensity
    }

    /// With the flag off, the reported intensity is exactly the raw box-sum;
    /// with it on it is the (positive) averagine matched-filter projection.
    #[test]
    fn projection_off_is_raw_sum_on_is_projection() {
        let theory = [0.5, 0.3, 0.2];
        let grid = centre_grid(50.0, 30.0, 20.0);

        let raw = integrate_grid(&grid, &theory, false);
        assert_eq!(raw, 100.0, "flag off must report the raw box-sum (50+30+20)");

        let proj = integrate_grid(&grid, &theory, true);
        // <[50,30,20], [0.5,0.3,0.2]> / ||[0.5,0.3,0.2]||₂ = 38 / 0.6164 ≈ 61.65
        assert!((proj - 61.65).abs() < 0.5, "projection ≈ 61.65, got {proj}");
    }

    /// The matched filter is less inflated by off-pattern contamination than the
    /// raw box-sum: adding intensity in the low-weight M+2 channel raises the
    /// projected value by proportionally less than it raises the raw sum.
    #[test]
    fn projection_downweights_off_pattern_contamination() {
        let theory = [0.5, 0.3, 0.2];
        let clean = centre_grid(50.0, 30.0, 20.0); // raw 100, on-pattern
        let contam = centre_grid(50.0, 30.0, 60.0); // raw 140, +40 off-pattern in M+2

        let raw_ratio = integrate_grid(&contam, &theory, false)
            / integrate_grid(&clean, &theory, false);
        let proj_ratio = integrate_grid(&contam, &theory, true)
            / integrate_grid(&clean, &theory, true);

        assert!(raw_ratio > 1.35, "raw sum tracks the contamination (140/100), got {raw_ratio}");
        assert!(
            proj_ratio < raw_ratio,
            "projection must be less inflated by off-pattern signal (proj {proj_ratio} vs raw {raw_ratio})"
        );
    }
}
