//! Averagine-based theoretical isotope distribution + Bhattacharyya scoring.
//!
//! The averagine model (Senko, 1995) approximates the average amino-acid
//! composition: C₄.₉₃₈₄ H₇.₇₅₈₃ N₁.₃₅₇₇ O₁.₄₇₇₃ S₀.₀₄₁₇ per 111.1254 Da.
//!
//! For a molecule of mass M, element counts are scaled, rounded to integers,
//! and the isotopologue distribution is computed by convolving the five
//! per-element distributions from [`super::elements`]. This module is a
//! thin wrapper that hands off all the arithmetic to the element cache —
//! the cache makes sulfur-aware scoring (vary `n_S` alone) cheap and exact.

use super::elements::{cache, K_PATTERN};

/// Averagine atomic ratios per 111.1254 Da of neutral mass. The single
/// "averagine residue" mass used by Senko's model.
const AVG_RESIDUE_MASS: f64 = 111.1254;
const AVG_C_PER: f64 = 4.9384;
const AVG_H_PER: f64 = 7.7583;
const AVG_N_PER: f64 = 1.3577;
const AVG_O_PER: f64 = 1.4773;
const AVG_S_PER: f64 = 0.0417;

/// Estimate integer element counts (C, H, N, O, S) from a neutral mass.
/// Half-up rounding from the averagine-scaled floats.
pub fn averagine_counts(neutral_mass: f64) -> (u32, u32, u32, u32, u32) {
    let scale = neutral_mass.max(0.0) / AVG_RESIDUE_MASS;
    let round = |x: f64| -> u32 { (x + 0.5).floor().max(0.0) as u32 };
    (
        round(AVG_C_PER * scale),
        round(AVG_H_PER * scale),
        round(AVG_N_PER * scale),
        round(AVG_O_PER * scale),
        round(AVG_S_PER * scale),
    )
}

/// Theoretical isotope distribution for the averagine-estimated composition
/// of a neutral mass. Returns a normalized `[p0, p1, …, p9]`.
pub fn averagine_distribution(neutral_mass: f64) -> [f64; K_PATTERN] {
    let (c, h, n, o, s) = averagine_counts(neutral_mass);
    cache().distribution(c, h, n, o, s)
}

/// Backward-compatible alias used by existing call sites. The return type is
/// now owned (cheap — 10 doubles) rather than `&'static`: pass-through is
/// the same since callers immediately use the array.
pub fn lookup_template(neutral_mass: f64) -> [f64; K_PATTERN] {
    averagine_distribution(neutral_mass)
}

/// Theoretical isotope distribution with an overridden sulfur count.
/// Useful for sulfur-aware scoring: hold C/H/N/O fixed, vary S.
pub fn averagine_distribution_with_sulfur(neutral_mass: f64, n_s: u32) -> [f64; K_PATTERN] {
    let (c, h, n, o, _) = averagine_counts(neutral_mass);
    cache().distribution(c, h, n, o, n_s)
}

/// Sulfur-aware Bhattacharyya scoring.
///
/// Builds isotope-pattern templates for a small set of sulfur counts spanning
/// the realistic biological range, scores observed against each, returns the
/// best BC and the winning `n_S`. This corrects the systematic bias on
/// peptides with 2+ Cys/Met where ³⁴S (4.25 %, +2 Da) elevates M+2 well
/// above the averagine prediction.
///
/// Template set per call: `[0, max(1, averagine_S), averagine_S + 2, averagine_S + 4]`,
/// deduplicated and clamped to `MAX_S` by the element cache.
pub fn bhattacharyya_score_best_sulfur(obs: &[f64], neutral_mass: f64) -> (f64, u32) {
    let (c, h, n, o, s_avg) = averagine_counts(neutral_mass);
    let variants = sulfur_variants(s_avg);

    let mut best_bc = 0.0f64;
    let mut best_s = s_avg;
    for &n_s in &variants {
        let template = cache().distribution(c, h, n, o, n_s);
        let bc = bhattacharyya_score(obs, &template);
        if bc > best_bc {
            best_bc = bc;
            best_s = n_s;
        }
    }
    (best_bc, best_s)
}

/// Deduplicated sulfur-count variants we score against:
/// `[0, max(1, s_avg), s_avg + 2, s_avg + 4]`. For small peptides (s_avg = 0)
/// this collapses to `{0, 1, 2, 4}`. For larger ones (s_avg = 2) it's
/// `{0, 2, 4, 6}`. Always includes both the no-S and a high-S extreme.
fn sulfur_variants(s_avg: u32) -> [u32; 4] {
    [0, s_avg.max(1), s_avg + 2, s_avg + 4]
}

/// Compute Bhattacharyya coefficient between observed and theoretical distributions.
///
/// BC = Σ sqrt(p_i * q_i) over all `k = min(obs.len(), 10)` positions, then
/// multiplied by `(1 - missed_penalty)` where `missed_penalty` is the
/// theoretical mass beyond position k.
///
/// No active-position mask. Positions where the template predicts ~nothing
/// still appear in both `obs_sum` and `theo_sum`. Effect: an over-extended
/// chain that grafted a noise hill at M+8 (where theory is ~0) dilutes
/// `obs_sum`, shrinking p_i at the real positions; the noise position itself
/// contributes ~0 to BC because q_i ≈ 0. Net: BC drops, as it should.
///
/// Returns BC in [0, 1].
pub fn bhattacharyya_score(obs: &[f64], template: &[f64; K_PATTERN]) -> f64 {
    let k = obs.len().min(K_PATTERN);
    if k == 0 {
        return 0.0;
    }

    let obs_sum: f64 = obs[..k].iter().sum();
    if obs_sum <= 0.0 {
        return 0.0;
    }

    let template_k_sum: f64 = template[..k].iter().sum();
    if template_k_sum <= 0.0 {
        return 0.0;
    }
    let missed_penalty = 1.0 - template_k_sum;

    let bc: f64 = (0..k)
        .map(|i| {
            let p = obs[i] / obs_sum;
            let q = template[i] / template_k_sum;
            (p * q).sqrt()
        })
        .sum();

    (bc * (1.0 - missed_penalty)).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averagine_counts_scales_linearly() {
        // 1500 Da peptide → averagine integer composition
        let (c, h, n, o, s) = averagine_counts(1500.0);
        let scale = 1500.0 / AVG_RESIDUE_MASS;
        assert!((c as f64 - 4.9384 * scale).abs() <= 1.0);
        assert!((h as f64 - 7.7583 * scale).abs() <= 1.0);
        assert!((n as f64 - 1.3577 * scale).abs() <= 1.0);
        assert!((o as f64 - 1.4773 * scale).abs() <= 1.0);
        assert!((s as f64 - 0.0417 * scale).abs() <= 1.0);
    }

    #[test]
    fn averagine_distribution_sums_to_one() {
        for mass in [500.0, 1500.0, 3000.0, 5000.0] {
            let d = averagine_distribution(mass);
            let sum: f64 = d.iter().sum();
            assert!((sum - 1.0).abs() < 1e-9, "{mass}: sum={sum}");
        }
    }

    /// Well-aligned 5-hill chain — BC should be high.
    #[test]
    fn well_aligned_chain_scores_high() {
        let template: [f64; K_PATTERN] = [0.50, 0.30, 0.15, 0.04, 0.008, 0.001, 0.0005, 0.0, 0.0, 0.0];
        let obs = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
        let bc = bhattacharyya_score(&obs, &template);
        assert!(bc > 0.95, "BC for well-aligned chain should be > 0.95, got {bc}");
    }

    /// Over-extended chain — same well-aligned first 5 positions plus 3 noise hills
    /// where the template predicts essentially nothing. New BC must drop noticeably.
    #[test]
    fn over_extended_chain_drops_score() {
        let template: [f64; K_PATTERN] = [0.50, 0.30, 0.15, 0.04, 0.008, 0.001, 0.0005, 0.0001, 0.0, 0.0];
        let obs_real = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
        let bc_real = bhattacharyya_score(&obs_real, &template);
        let obs_overextended = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5, 5.0e6, 5.0e6, 5.0e6];
        let bc_over = bhattacharyya_score(&obs_overextended, &template);
        assert!(
            bc_over < bc_real - 0.05,
            "over-extended BC ({bc_over:.3}) should drop at least 0.05 below clean BC ({bc_real:.3})"
        );
    }

    /// Empty / zero input safe.
    #[test]
    fn empty_input_returns_zero() {
        let template: [f64; K_PATTERN] = [0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(bhattacharyya_score(&[], &template), 0.0);
        assert_eq!(bhattacharyya_score(&[0.0, 0.0], &template), 0.0);
    }

    /// Sanity: sulfur-override path returns a different template than the default.
    #[test]
    fn sulfur_override_changes_pattern() {
        let m = 1500.0;
        let avg = averagine_distribution(m);
        let high_s = averagine_distribution_with_sulfur(m, 3);
        assert!(
            (high_s[2] - avg[2]).abs() > 0.005,
            "3-S override should noticeably shift M+2 (avg={}, hi-S={})", avg[2], high_s[2]
        );
    }
}
