//! Exhaustive / non-destructive isotope assembler.
//!
//! An over-complete `(seed, charge)` hypothesis pool is generated in parallel
//! (`build_charge_candidate`) and resolved non-destructively (`resolve_exhaustive`)
//! by claiming each candidate's best valid prefix in descending evidence order.
//! After conflicts, the best remaining free prefix is rescored and requeued.
//!
//! **Chains extend upward only.** The seed IS the monoisotope hypothesis: a
//! chain runs seed → M+1 → M+2 → …, never downward. Every hill is tried as a
//! seed, so a hill a downward walk could have reached is itself a seed that
//! builds the same envelope upward — with the averagine template indexed from
//! its own position rather than from a seed that later turns out not to be the
//! monoisotope. This also keeps `max_isotopes` a bound on the envelope rather
//! than on each direction separately.
//!
//! Consequence to be aware of: nothing downstream reassigns the monoisotope.
//! If a true monoisotope's own upward chain fails to form, its M+1 may be
//! emitted as a separate feature whose reported mass is one neutron heavy;
//! `scoring::isotope_offset_enabled` only re-labels peaks already in a chain
//! and cannot recover a monoisotope that was never chained.

use std::cmp::Ordering;

use rayon::prelude::*;

use crate::config::{CosineAnchor, FeaturesConfig, FileConfig, ImToleranceType, ToleranceType};
use crate::models::Hill;
use crate::models::Polarity;
use crate::scoring::averagine;
use crate::scoring::model::IsotopeModel;

use super::cosine::cosine_similarity;
use super::recalibration::MzRecalModel;

/// Immutable per-detection context threaded through the chain-assembly
/// functions (`resolve_exhaustive`, `build_charge_candidate`, `find_neighbors`).
/// Built once in `detect_features_with_recal` and passed by shared reference, so
/// the assembler's parallel per-seed generation reads only shared immutable data.
pub(super) struct ChainCtx<'a> {
    pub(super) sorted_hills: &'a [&'a Hill],
    pub(super) mz_array: &'a [f64],
    pub(super) im_array: &'a [f64],
    pub(super) scan_starts: &'a [usize],
    pub(super) scan_ends: &'a [usize],
    pub(super) config: &'a FeaturesConfig,
    pub(super) file: &'a FileConfig,
    pub(super) use_im: bool,
    pub(super) min_intensity: f64,
    pub(super) recal: Option<&'a MzRecalModel>,
}

pub(super) struct Candidate {
    pub(super) hill_indices: Vec<usize>, // indices into sorted_hills, lowest mz first
    pub(super) charge: u8,
}

/// Build the best isotope-chain candidate for one `(seed, charge)` pair by
/// extending upward (M+1, M+2, …) from the seed under chromatographic-cosine
/// and isotope-pattern gates; the resolver (`resolve_exhaustive`) scores and
/// claims. Returns `None` when no isotope partner is found (a lone seed —
/// charge 0).
pub(super) fn build_charge_candidate(
    ctx: &ChainCtx,
    seed_idx: usize,
    charge: u8,
) -> Option<Candidate> {
    let ChainCtx {
        sorted_hills,
        config,
        file,
        min_intensity,
        recal,
        ..
    } = *ctx;
    let seed_hill = sorted_hills[seed_idx];
    let seed_mz = seed_hill.mz;
    let seed_rt = seed_hill.rt;
    // Isotope-consistency recalibration: shift an expected isotope m/z by the
    // learned per-region proportional offset. All isotopes of a feature
    // co-elute, so the seed RT keys the whole chain. No-op when `recal` is
    // `None`.
    let recal_mz = |expected: f64| -> f64 {
        match recal {
            Some(m) => expected * (1.0 + m.predict(expected, seed_rt) / 1e6),
            None => expected,
        }
    };
    let use_ppm_tol = matches!(file.mz_tolerance_type, ToleranceType::Ppm);
    // Region-adaptive isotope-match tolerance: when recalibration is on, derive
    // the ppm window from the surface's per-region residual spread σ instead of
    // the fixed `mz_tolerance`, clamped to [tol_floor_ppm, mz_tolerance].
    // Isotopes of one feature span only a few Da, so a single σ read at the seed
    // applies to the whole chain. No-op for Dalton tolerances or when `recal` is
    // `None` (recalibration disabled).
    let effective_tol_ppm = match (use_ppm_tol, recal) {
        (true, Some(m)) => {
            // `f64::clamp` panics if min > max; a user-configured floor above
            // mz_tolerance would otherwise crash every feature-detection task.
            // Capping the floor at the ceiling makes mz_tolerance the true
            // upper bound regardless of how the floor is set.
            let floor = file.mz_recalibration_tol_floor_ppm.min(file.mz_tolerance);
            (file.mz_recalibration_tol_sigma_mult * m.predict_sigma(seed_mz, seed_rt))
                .clamp(floor, file.mz_tolerance)
        }
        _ => file.mz_tolerance,
    };
    let mz_tol = if use_ppm_tol {
        seed_mz * effective_tol_ppm / 1e6
    } else {
        file.mz_tolerance
    };
    // Which hill the chromatographic-cosine gate is anchored to. `Seed` is the
    // default; `Adjacent` reproduces the legacy predecessor anchor. Either way
    // only the cosine *reference* hill changes — the m/z step target,
    // `find_neighbors` predecessor, and the intensity-ratio predecessor all stay
    // the immediate chain predecessor.
    let cosine_anchor = config.cosine_anchor_mode();

    let step = config.neutron_mass / charge as f64;
    {
        // Per-charge averagine template, used to early-stop the chain when the
        // next theoretical isotope falls below the noise floor. The seed is
        // taken as monoisotopic — that is the assembler's contract, not a
        // provisional guess, so template[iso] indexes the chain directly.
        let neutral_mass_seed_mono = file.polarity.neutral_mass(seed_mz, charge);
        let template = config
            .isotope_model
            .model()
            .distribution(neutral_mass_seed_mono);
        let template_mono = template[0].max(1e-12);
        let seed_intensity = seed_hill.intensity_max;

        // Upward chain: M+1, M+2, ... (the only direction)
        let mut right_chain: Vec<usize> = Vec::new();
        for iso in 1..=config.max_isotopes {
            // Early-stop: bail when theory predicts a peak below noise. Stops
            // the chain from absorbing column-bleed / random hills that happen
            // to fall on the M+k ladder but have no isotope-pattern basis.
            // (Cosine ≥ min_chain_cosine on adjacent pairs isn't enough — two
            // neighbouring contaminants can trivially co-elute with each other.)
            if iso >= template.len() {
                break;
            }
            // Predicted-intensity early-stop (legacy, gated). Bails before the
            // hill is searched for when averagine says the M+k peak would land
            // below the noise floor. Cheap tail-trim, but a dim seed's predicted
            // M+1 is itself sub-floor, so this also blocks legitimate
            // low-abundance monos from ever chaining -> they collapse to
            // charge 0. When disabled, extension is terminated purely by
            // evidence (find_neighbors empty / cosine / intensity-ratio gate /
            // min_isotope_step_ratio / template caps) — see `find_neighbors` and the
            // intensity-ratio gate below.
            if config.chain_predicted_intensity_gate {
                let expected = seed_intensity * (template[iso] / template_mono);
                if expected < min_intensity {
                    break;
                }
            }

            let target_mz = recal_mz(seed_mz + step * iso as f64);
            let ref_hill = right_chain
                .last()
                .map(|&i| sorted_hills[i])
                .unwrap_or(seed_hill);

            let mut exclude = Vec::with_capacity(right_chain.len() + 1);
            exclude.push(seed_idx);
            exclude.extend_from_slice(&right_chain);

            let cands = find_neighbors(
                ctx,
                target_mz,
                mz_tol,
                ref_hill,
                config.min_isotope_step_ratio,
                &exclude,
            );
            if cands.is_empty() {
                break;
            }

            // Cosine reference hill. `Seed` (the default) anchors every
            // isotope's chromatographic cosine to the mono seed, matching the
            // biosaur2 / AlphaPept / Dinosaur convention. It was long assumed
            // that seed-anchoring would over-reject far isotopes (chains drift
            // in S/N from the seed as you go out), which is why koth originally
            // anchored to the immediate predecessor (`Adjacent`). The full
            // 20-run PXD003881 cohort refuted that: seed-anchoring recovered
            // +1565 covered PSMs (+0.31 pp recall) with no quant regression —
            // the recall gain is entirely in the M+2/M+3 isotopes (M+1's
            // predecessor IS the seed, so it is unaffected). `Adjacent` anchors
            // to the immediate predecessor (former default, byte-identical to
            // the pre-2026-07 paper output). Only the cosine reference changes —
            // the m/z step and intensity-ratio predecessor stay `ref_hill`.
            let cos_ref_hill = match cosine_anchor {
                CosineAnchor::Adjacent => ref_hill,
                CosineAnchor::Seed => seed_hill,
            };
            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| {
                    (
                        c,
                        cosine_similarity(cos_ref_hill, sorted_hills[c], config.min_scan_overlap),
                    )
                })
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();

            if best_cos < config.min_chain_cosine {
                break;
            }

            // Intensity-ratio gate: candidate's apex intensity vs predecessor
            // should match the averagine ratio template[iso]/template[iso-1]
            // within ± `max_isotope_log2_ratio` log2 units. Catches contaminant
            // hills that co-elute (so pass cosine) but have wrong intensity
            // for a real isotope peak.
            if config.max_isotope_log2_ratio.is_finite() {
                let pred_int = ref_hill.intensity_max;
                let cand_int = sorted_hills[best_c].intensity_max;
                let theo_pred = template[iso - 1].max(1e-12);
                let theo_cand = template[iso].max(1e-12);
                if pred_int > 0.0 && cand_int > 0.0 {
                    let observed_ratio = cand_int / pred_int;
                    let expected_ratio = theo_cand / theo_pred;
                    if (observed_ratio / expected_ratio).log2().abs()
                        > config.max_isotope_log2_ratio
                    {
                        break;
                    }
                }
            }

            right_chain.push(best_c);
        }

        if right_chain.is_empty() {
            return None; // no partners found — single-hill fallback stays charge=0
        }

        // Assemble chain: seed + right (already low→high mz). The seed IS the
        // monoisotope hypothesis — there is no downward extension. Every hill is
        // tried as a seed, so a hill that a downward walk could have reached is
        // itself a seed that builds the same chain upward, with the averagine
        // template correctly indexed from its own position. See the module docs.
        let chain: Vec<usize> = std::iter::once(seed_idx)
            .chain(right_chain.iter().copied())
            .collect();

        Some(Candidate {
            hill_indices: chain,
            charge,
        })
    }
}

/// Summed log evidence for every monoisotope-anchored prefix (index = length - 1).
/// Each M+k contributes a seed-relative apex-ratio and chromatographic-cosine
/// log-likelihood ratio. Ratio residuals are N(0, signal_sigma²) against a broad
/// N(0, 2²) noise null; cosine is Beta(cosine_shape, 1) against Uniform(0, 1).
/// These working models are not calibrated feature probabilities or FDRs.
///
/// Accumulate each sulfur template separately, then maximize each whole-prefix
/// sum over templates. Retained observations never depend on tail length or a
/// combined apex, so removing a tail cannot alter the retained evidence terms.
fn prefix_log_evidence(
    hills: &[&Hill],
    charge: u8,
    config: &FeaturesConfig,
    polarity: Polarity,
) -> Vec<f64> {
    let mut scores = vec![f64::NEG_INFINITY; hills.len()];
    if hills.len() < 2 || charge == 0 {
        return scores;
    }
    let seed_apex = hills[0].intensity_max;
    if !seed_apex.is_finite() || seed_apex <= 0.0 {
        return scores;
    }
    let model = config.isotope_model.model();
    let neutral_mass = polarity.neutral_mass(hills[0].mz, charge);
    let (c, h, n, o, s) = model.counts(neutral_mass);
    let sulfur_counts = if config.sulfur_offsets.is_empty() || !model.has_sulfur() {
        vec![s]
    } else {
        averagine::resolve_sulfur_counts(&model, neutral_mass, &config.sulfur_offsets)
    };
    let observations: Vec<(f64, f64)> = hills[1..]
        .iter()
        .map(|hill| {
            let log_ratio = hill.intensity_max.log2() - seed_apex.log2();
            let cosine = cosine_similarity(hills[0], hill, config.min_scan_overlap);
            (log_ratio, cosine)
        })
        .collect();
    for s in sulfur_counts {
        let template = crate::scoring::elements::cache().distribution(c, h, n, o, s);
        let mut sum = 0.0;
        for (k, &(observed, cosine)) in observations.iter().enumerate() {
            let expected = template.get(k + 1).copied().unwrap_or(0.0);
            let contribution = if expected <= 0.0 || template[0] <= 0.0 {
                f64::NEG_INFINITY
            } else {
                let residual = observed - (expected.log2() - template[0].log2());
                isotope_log_evidence(residual, cosine, config)
            };
            sum += contribution;
            scores[k + 1] = scores[k + 1].max(sum);
        }
    }
    scores
}

/// Select the best valid prefix, both before heap insertion and after conflicts.
/// The heap key is therefore an upper bound on every remaining valid prefix:
/// claiming hills can only remove choices, never expose a higher-scoring one.
/// Tied prefixes choose the smaller set of hills; length never adds evidence.
fn rescore_chain(
    chain: &[usize],
    charge: u8,
    sorted_hills: &[&Hill],
    config: &FeaturesConfig,
    polarity: Polarity,
) -> Option<HeapItem> {
    let hills: Vec<&Hill> = chain.iter().map(|&i| sorted_hills[i]).collect();
    let scores = prefix_log_evidence(&hills, charge, config, polarity);
    let mut ends: Vec<usize> = (1..chain.len()).collect();
    ends.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
    for end in ends {
        if !scores[end].is_finite() {
            continue;
        }
        // Evaluate the unchanged Bhattacharyya claim gate only when enabled.
        // A failing best prefix must not hide a lower-scoring valid prefix.
        if config.exhaustive_min_isotope_score > 0.0
            && score_chain(
                &hills[..=end],
                charge,
                &config.sulfur_offsets,
                &config.isotope_model.model(),
                polarity,
            ) < config.exhaustive_min_isotope_score
        {
            continue;
        }
        return Some(HeapItem {
            chain: chain[..=end].to_vec(),
            charge,
            log_evidence: scores[end],
        });
    }
    None
}

/// Natural-log likelihood ratio for one seed-anchored isotope observation.
fn isotope_log_evidence(log2_ratio_residual: f64, cosine: f64, config: &FeaturesConfig) -> f64 {
    if !log2_ratio_residual.is_finite() || !cosine.is_finite() || cosine <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let signal_sigma = config.isotope_evidence_ratio_sigma;
    const NOISE_SIGMA: f64 = 2.0;
    let ratio_llr = (NOISE_SIGMA / signal_sigma).ln()
        - 0.5 * log2_ratio_residual.powi(2) * (signal_sigma.powi(-2) - NOISE_SIGMA.powi(-2));
    let shape = config.isotope_evidence_cosine_shape;
    let cosine_llr = shape.ln() + (shape - 1.0) * cosine.min(1.0).ln();
    ratio_llr + cosine_llr
}

/// Priority-queue item for the exhaustive resolver.
///
/// Ordering (max-heap): higher summed log evidence first, with smaller
/// monoisotope index then smaller charge as deterministic ties. Envelope
/// length and the diagnostic isotope score never override the evidence.
struct HeapItem {
    chain: Vec<usize>, // m/z-sorted, monoisotope first
    charge: u8,
    log_evidence: f64,
}

impl HeapItem {
    fn mono(&self) -> usize {
        self.chain[0]
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for HeapItem {}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.log_evidence
            .total_cmp(&other.log_evidence)
            // smaller mono index should pop first => compare as GREATER (reverse).
            .then_with(|| other.mono().cmp(&self.mono()))
            // smaller charge pops first => reverse.
            .then_with(|| other.charge.cmp(&self.charge))
    }
}

/// Non-destructive, confidence-ordered conflict resolution.
///
/// Generates one candidate per `(seed, charge)` with ≥1 isotope partner, then
/// resolves contested hills by the maximum summed log evidence over valid
/// prefixes. Conflicts restrict the feasible set to still-free prefixes, whose
/// maximum cannot exceed the previously queued score. Drop only when the mono
/// is taken or no finite-scoring prefix with an isotope partner passes the gate.
pub(super) fn resolve_exhaustive(ctx: &ChainCtx) -> Vec<Candidate> {
    use std::collections::BinaryHeap;

    let sorted_hills = ctx.sorted_hills;
    let config = ctx.config;
    let polarity = ctx.file.polarity;

    // Phase 1+2: over-complete hypothesis pool — one candidate per (seed, charge),
    // built by the SAME `build_charge_candidate` used throughout. Every
    // (seed, charge) is independent and reads only shared immutable arrays, so
    // generation runs in PARALLEL over seeds (like the main detection loop's
    // `into_par_iter`), then the results are heapified once via an O(n)
    // `BinaryHeap::from`. The heap's `Ord` orders by log evidence, then
    // monoisotope index and charge, so the pop sequence is independent of
    // insertion order. `flat_map_iter` keeps
    // the per-seed charge loop sequential inside each parallel task.
    let items: Vec<HeapItem> = (0..sorted_hills.len())
        .into_par_iter()
        .flat_map_iter(|seed_idx| {
            (config.min_charge..=config.max_charge).filter_map(move |charge| {
                build_charge_candidate(ctx, seed_idx, charge).and_then(|c| {
                    rescore_chain(&c.hill_indices, c.charge, sorted_hills, config, polarity)
                })
            })
        })
        .collect();
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::from(items);

    // Phase 3: non-destructive resolution.
    let mut claimed = vec![false; sorted_hills.len()];
    let mut accepted: Vec<Candidate> = Vec::new();

    while let Some(item) = heap.pop() {
        // Monoisotope must be free, else the whole candidate is dead.
        if claimed[item.mono()] {
            continue;
        }
        // Longest monoisotope-anchored prefix whose hills are all still free.
        let mut k = 0usize;
        while k < item.chain.len() && !claimed[item.chain[k]] {
            k += 1;
        }
        if k == item.chain.len() {
            for &i in &item.chain {
                claimed[i] = true;
            }
            accepted.push(Candidate {
                hill_indices: item.chain,
                charge: item.charge,
            });
        } else if k >= 2 {
            // Re-select among all still-free prefixes. The queued score was
            // already their maximum, so rescoring cannot increase priority.
            if let Some(prefix) = rescore_chain(
                &item.chain[..k],
                item.charge,
                sorted_hills,
                config,
                polarity,
            ) {
                debug_assert!(prefix.log_evidence <= item.log_evidence);
                heap.push(prefix);
            }
        }
        // k == 1: only the monoisotope survives — drop (a lone hill is not a
        // charged feature; a lone monoisotope is emitted as charge-0, which output filters).
    }

    accepted
}

fn score_chain(
    chain_hills: &[&Hill],
    charge: u8,
    sulfur_offsets: &[i8],
    model: &IsotopeModel,
    polarity: Polarity,
) -> f64 {
    if chain_hills.len() <= 1 || charge == 0 {
        return 0.0;
    }

    let mut sorted: Vec<&Hill> = chain_hills.to_vec();
    sorted.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

    let neutral_mass = polarity.neutral_mass(sorted[0].mz, charge);

    let min_scan = sorted.iter().map(|h| h.scan_start).min().unwrap_or(0);
    let max_scan = sorted.iter().map(|h| h.scan_end).max().unwrap_or(0);
    if max_scan < min_scan {
        return 0.0;
    }
    let n = max_scan - min_scan + 1;
    let mut combined = vec![0.0f64; n];
    for hill in &sorted {
        for (i, &v) in hill.intensity_profile.iter().enumerate() {
            let idx = hill.scan_start + i - min_scan;
            if idx < n {
                combined[idx] += v as f64;
            }
        }
    }
    let apex_scan = min_scan
        + combined
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);

    let obs: Vec<f64> = sorted
        .iter()
        .map(|h| h.intensity_at_scan(apex_scan))
        .collect();
    let k = obs.len().min(10);

    if sulfur_offsets.is_empty() {
        let template = model.distribution(neutral_mass);
        averagine::bhattacharyya_score(&obs[..k], &template)
    } else {
        averagine::bhattacharyya_score_best_sulfur(&obs[..k], neutral_mass, sulfur_offsets, model).0
    }
}

/// Find candidate isotope partners near `target_mz`.
///
/// Filters: mz tolerance (binary search), not in `exclude_indices` (current
/// chain members), scan overlap with `ref_hill`, IM tolerance, intensity floor.
fn find_neighbors(
    ctx: &ChainCtx,
    target_mz: f64,
    mz_tol: f64,
    ref_hill: &Hill,
    max_decrease: f64,
    exclude_indices: &[usize],
) -> Vec<usize> {
    let ChainCtx {
        mz_array,
        im_array,
        scan_starts,
        scan_ends,
        sorted_hills,
        file,
        use_im,
        ..
    } = *ctx;
    let lo = mz_array.partition_point(|&x| x < target_mz - mz_tol);
    let hi = mz_array.partition_point(|&x| x <= target_mz + mz_tol);

    if lo >= hi {
        return Vec::new();
    }

    let min_intensity = ref_hill.intensity_max * max_decrease;

    (lo..hi)
        .filter(|&i| {
            if exclude_indices.contains(&i) {
                return false;
            }
            if scan_ends[i] < ref_hill.scan_start || scan_starts[i] > ref_hill.scan_end {
                return false;
            }
            if use_im {
                let im_tol = match file.im_tolerance_type {
                    ImToleranceType::Absolute => file.im_tolerance,
                    ImToleranceType::Relative => ref_hill.im.max(im_array[i]) * file.im_tolerance,
                };
                if (im_array[i] - ref_hill.im).abs() > im_tol {
                    return false;
                }
            }
            sorted_hills[i].intensity_max >= min_intensity
        })
        .collect()
}

#[cfg(test)]
#[path = "assemble_tests.rs"]
mod tests;
