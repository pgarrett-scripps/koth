//! Isotope-pattern scoring: Bhattacharyya coefficient against a theoretical
//! distribution from an [`IsotopeModel`].
//!
//! The composition itself lives in [`super::model`]; this module is the scoring
//! layer over it, including the sulfur-aware maximum used for peptides. All
//! arithmetic is handed off to the element cache in [`super::elements`], which
//! makes varying `n_S` alone cheap and exact.

use super::elements::{cache, K_PATTERN};
use super::model::IsotopeModel;

/// Resolve user-declared sulfur `offsets` against a molecule's expected count.
///
/// Each offset is applied to `ceil(expected n_S)`, negatives saturate at 0, and
/// the result is deduplicated (order-preserving) so a config like
/// `[-2, -1, 0, 1]` on a small peptide collapses to `{0, 1, 2}` rather than
/// scoring `n_S = 0` three times. Counts above `MAX_S` are clamped by the
/// element cache.
///
/// Empty `offsets` means "no sulfur awareness" and yields a single template at
/// the plain half-up rounded count — see [`bhattacharyya_score_best_sulfur`].
pub fn resolve_sulfur_counts(model: &IsotopeModel, neutral_mass: f64, offsets: &[i8]) -> Vec<u32> {
    let s_ceil = model.sulfur_ceil(neutral_mass) as i32;
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
/// Scores `obs` against one template per declared sulfur offset (see
/// [`resolve_sulfur_counts`]) and returns the best BC with its winning `n_S`.
/// For peptides this corrects the systematic bias on Cys/Met-rich sequences,
/// where ³⁴S (4.25 %, +2 Da) elevates M+2 well above the prediction.
///
/// Offsets are relative to the **ceiling** of the expected count, so on the
/// peptide model the default `[-1, 0, 1]` gives `{0, 1, 2}` up to 2665 Da,
/// `{1, 2, 3}` to 5330, `{2, 3, 4}` to 7995. Declare `[-1, 0, 1, 2]` (or include
/// a `0`-reaching offset) to keep the no-sulfur template on large peptides —
/// ~32 % of 3 kDa tryptic peptides have no sulfur at all.
///
/// An empty offset list disables sulfur awareness: a single plain template at
/// the half-up rounded count. **A model with no sulfur** (RNA, DNA) does the
/// same whatever the offsets say — see [`IsotopeModel::has_sulfur`].
///
/// Note this is a **max over templates**, so widening the offset list can only
/// raise scores — including for decoys and mis-assembled chains. Judge a change
/// to the list on discrimination (recall / target-decoy separation), never on
/// the isotope-score distribution alone.
pub fn bhattacharyya_score_best_sulfur(
    obs: &[f64],
    neutral_mass: f64,
    offsets: &[i8],
    model: &IsotopeModel,
) -> (f64, u32) {
    let (c, h, n, o, s_avg) = model.counts(neutral_mass);

    if offsets.is_empty() || !model.has_sulfur() {
        let template = cache().distribution(c, h, n, o, s_avg);
        return (bhattacharyya_score(obs, &template), s_avg);
    }

    let mut best_bc = 0.0f64;
    let mut best_s = model.sulfur_ceil(neutral_mass);
    for n_s in resolve_sulfur_counts(model, neutral_mass, offsets) {
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
