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

/// Averagine sulfur count rounded **up** (`ceil`) instead of half-up.
/// Only the sulfur term changes; C/H/N/O keep `averagine_counts`' rounding.
///
/// `s_float = 0.0417·M / 111.1254 = 3.7524e-4·M`, so the ceiling steps at
/// M ≈ 2665 / 5330 / 7995 Da.
pub fn averagine_sulfur_ceil(neutral_mass: f64) -> u32 {
    let scale = neutral_mass.max(0.0) / AVG_RESIDUE_MASS;
    (AVG_S_PER * scale).ceil().max(0.0) as u32
}

/// Resolve user-declared sulfur `offsets` against a peptide's expected count.
///
/// Each offset is applied to `ceil(expected n_S)`, negatives saturate at 0, and
/// the result is deduplicated (order-preserving) so a config like
/// `[-2, -1, 0, 1]` on a small peptide collapses to `{0, 1, 2}` rather than
/// scoring `n_S = 0` three times. Counts above `MAX_S` are clamped by the
/// element cache.
///
/// Empty `offsets` means "no sulfur awareness" and yields a single template at
/// the plain averagine (half-up rounded) count — see
/// [`bhattacharyya_score_best_sulfur`].
pub fn resolve_sulfur_counts(neutral_mass: f64, offsets: &[i8]) -> Vec<u32> {
    let s_ceil = averagine_sulfur_ceil(neutral_mass) as i32;
    let mut out: Vec<u32> = Vec::with_capacity(offsets.len());
    for &off in offsets {
        let n_s = (s_ceil + off as i32).max(0) as u32;
        if !out.contains(&n_s) {
            out.push(n_s);
        }
    }
    out
}

/// Sulfur-aware Bhattacharyya scoring.
///
/// Scores `obs` against one averagine template per user-declared sulfur offset
/// (see [`resolve_sulfur_counts`]) and returns the best BC with its winning
/// `n_S`. This corrects the systematic bias on Cys/Met-rich peptides, where
/// ³⁴S (4.25 %, +2 Da) elevates M+2 well above the averagine prediction.
///
/// Offsets are relative to the **ceiling** of the expected count, so the
/// default `[-1, 0, 1]` gives `{0, 1, 2}` up to 2665 Da, `{1, 2, 3}` to 5330,
/// `{2, 3, 4}` to 7995. Declare `[-1, 0, 1, 2]` (or include a `0`-reaching
/// offset) if you want the no-sulfur template retained on large peptides —
/// ~32 % of 3 kDa tryptic peptides have no sulfur at all.
///
/// An empty offset list disables sulfur awareness: a single plain averagine
/// template at the half-up rounded count, which is the pre-2026-07 `false`
/// behaviour of `sulfur_aware_scoring`.
///
/// Note this is a **max over templates**, so widening the offset list can only
/// raise scores — including for decoys and mis-assembled chains. Judge a change
/// to the list on discrimination (recall / target-decoy separation), never on
/// the isotope-score distribution alone.
pub fn bhattacharyya_score_best_sulfur(
    obs: &[f64],
    neutral_mass: f64,
    offsets: &[i8],
) -> (f64, u32) {
    let (c, h, n, o, s_avg) = averagine_counts(neutral_mass);

    if offsets.is_empty() {
        let template = cache().distribution(c, h, n, o, s_avg);
        return (bhattacharyya_score(obs, &template), s_avg);
    }

    let mut best_bc = 0.0f64;
    let mut best_s = averagine_sulfur_ceil(neutral_mass);
    for n_s in resolve_sulfur_counts(neutral_mass, offsets) {
        let template = cache().distribution(c, h, n, o, n_s);
        let bc = bhattacharyya_score(obs, &template);
        if bc > best_bc {
            best_bc = bc;
            best_s = n_s;
        }
    }
    (best_bc, best_s)
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
#[path = "averagine_tests.rs"]
mod tests;
