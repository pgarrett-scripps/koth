pub mod cosine;
pub mod recalibration;

use std::cmp::Ordering;

use rayon::prelude::*;

use crate::config::{
    CosineAnchor, FeaturesConfig, FileConfig, ImToleranceType, MzUncertaintyMode, ToleranceType,
};
use crate::models::{Feature, Hill};
use crate::scoring::averagine;
use cosine::cosine_similarity;
use recalibration::{MzRecalBuilder, MzRecalModel};

/// When `mz_uncertainty_mode = Kish`, the binary-search window is widened
/// by this factor so candidates with high `mz_se` that pass the exact
/// combined-tolerance check below aren't pre-filtered out. The exact
/// check inside `find_neighbors` then rejects anything outside the
/// per-candidate Kish-combined window.
const KISH_SEARCH_EXPANSION: f64 = 4.0;

struct Candidate {
    hill_indices: Vec<usize>, // indices into sorted_hills, lowest mz first
    charge: u8,
    mean_cosine: f64,
    mean_ppm: f64,
    composite_score: f64,
}

/// Detect isotope features from a list of chromatographic hills.
///
/// 1. Build mz-sorted index.
/// 2. Every hill is tried as a seed; for each charge state the best isotope
///    chain is built via chromatographic-cosine filtering against the seed
///    (gated by `FeaturesConfig.min_chain_cosine`). Each (seed, charge) pair
///    is scored as `isotope_score × cosine_score` (Bhattacharyya isotope-pattern
///    times mean chromatographic cosine) and the charge with the highest
///    combined score is kept as that seed's candidate.
/// 3. All candidates are sorted by composite_score (desc); greedy conflict
///    resolution accepts the highest-scoring candidate whose hills are
///    unclaimed.
/// 4. Feature structs are built from accepted candidates.
pub fn detect_features(hills: &[Hill], config: &FeaturesConfig, file: &FileConfig) -> Vec<Feature> {
    detect_features_with_recal(hills, config, file, None)
}

/// As [`detect_features`], but shifts the expected isotope position during
/// chain extension by the per-region offset from an isotope-consistency
/// recalibration surface (see [`recalibration`]). Pass `None` for the
/// uncorrected behaviour.
pub fn detect_features_with_recal(
    hills: &[Hill],
    config: &FeaturesConfig,
    file: &FileConfig,
    recal: Option<&MzRecalModel>,
) -> Vec<Feature> {
    if hills.is_empty() {
        return Vec::new();
    }

    log::info!("Detecting isotope features from {} hills", hills.len());

    let mut order: Vec<usize> = (0..hills.len()).collect();
    order.sort_by(|&a, &b| hills[a].mz.partial_cmp(&hills[b].mz).unwrap_or(Ordering::Equal));
    let sorted_hills: Vec<&Hill> = order.iter().map(|&i| &hills[i]).collect();
    let mz_array: Vec<f64> = sorted_hills.iter().map(|h| h.mz).collect();
    let im_array: Vec<f64> = sorted_hills.iter().map(|h| h.im).collect();
    let scan_starts: Vec<usize> = sorted_hills.iter().map(|h| h.scan_start).collect();
    let scan_ends: Vec<usize> = sorted_hills.iter().map(|h| h.scan_end).collect();

    let use_im = im_array.iter().any(|&x| x != 0.0);

    // Noise floor for Bhattacharyya scoring (mirrors scoring::score_features).
    let min_intensity = {
        let mut sums: Vec<f64> = hills.iter().map(|h| h.intensity_sum).collect();
        sums.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if sums.is_empty() {
            0.0
        } else {
            let idx = (sums.len() as f64 * 0.05) as usize;
            sums[idx] * 0.8
        }
    };

    let accepted: Vec<Candidate> = if config.exhaustive_assembly {
        // Experimental non-destructive assembler (biosaur2 / AlphaPept style):
        // over-complete (seed, charge) pool, confidence-ordered claiming, and
        // truncation instead of all-or-nothing drop. Chains are built by the
        // SAME `build_charge_candidate` the greedy path uses, so ON vs OFF
        // differ only in conflict resolution, not chain construction.
        resolve_exhaustive(
            &sorted_hills,
            &mz_array,
            &im_array,
            &scan_starts,
            &scan_ends,
            config,
            file,
            use_im,
            min_intensity,
            recal,
        )
    } else {
        // Legacy greedy path (default). Phase 1+2: generate and score one
        // candidate per seed. Each seed's candidate is computed independently
        // from read-only shared arrays, so this is embarrassingly parallel.
        // `into_par_iter().collect()` is order-preserving (rayon's indexed
        // collect), so the downstream sort + greedy claim stay deterministic.
        let candidates: Vec<Candidate> = (0..sorted_hills.len())
            .into_par_iter()
            .map(|seed_idx| {
                generate_best_candidate(
                    seed_idx,
                    &sorted_hills,
                    &mz_array,
                    &im_array,
                    &scan_starts,
                    &scan_ends,
                    config,
                    file,
                    use_im,
                    min_intensity,
                    recal,
                )
            })
            .collect();

        // Phase 3: sort by composite_score desc, greedy conflict resolution.
        let mut candidates = candidates;
        candidates.sort_by(|a, b| {
            b.composite_score
                .partial_cmp(&a.composite_score)
                .unwrap_or(Ordering::Equal)
        });

        let mut claimed = vec![false; sorted_hills.len()];
        candidates
            .into_iter()
            .filter(|c| {
                if c.hill_indices.iter().all(|&i| !claimed[i]) {
                    for &i in &c.hill_indices {
                        claimed[i] = true;
                    }
                    true
                } else {
                    false
                }
            })
            .collect()
    };

    // Phase 4: build Feature structs.
    let features: Vec<Feature> = accepted
        .into_iter()
        .map(|c| {
            let mut chain_hills: Vec<Hill> =
                c.hill_indices.iter().map(|&i| (*sorted_hills[i]).clone()).collect();
            chain_hills.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

            let cosine_sims: Vec<f64> = (0..chain_hills.len().saturating_sub(1))
                .map(|k| {
                    cosine_similarity(
                        &chain_hills[k],
                        &chain_hills[k + 1],
                        config.cosine_intersection,
                        config.min_scan_overlap,
                    )
                })
                .collect();

            let ppm_errors: Vec<f64> = if c.charge > 0 && chain_hills.len() > 1 {
                let step = config.neutron_mass / c.charge as f64;
                (0..chain_hills.len() - 1)
                    .map(|k| {
                        let theo = chain_hills[k].mz + step;
                        let expt = chain_hills[k + 1].mz;
                        (expt - theo) / theo * 1e6
                    })
                    .collect()
            } else {
                Vec::new()
            };

            let mean_cosine = if cosine_sims.is_empty() {
                0.0
            } else {
                cosine_sims.iter().sum::<f64>() / cosine_sims.len() as f64
            };
            let mean_ppm = if ppm_errors.is_empty() {
                0.0
            } else {
                ppm_errors.iter().sum::<f64>() / ppm_errors.len() as f64
            };

            let charge = if c.hill_indices.len() > 1 { c.charge } else { 0 };

            Feature {
                hills: chain_hills,
                charge,
                cosine_score: mean_cosine,
                ppm_error: mean_ppm,
            }
        })
        .collect();

    log::info!("Detected {} isotope features", features.len());
    features
}

/// Learn an isotope-consistency recalibration surface from a pass-1 feature
/// detection over `hills`. Collects, for every accepted charged feature with
/// ≥2 isotope hills, the signed ppm deviation of each adjacent isotope-hill
/// spacing from the theoretical `neutron_mass / z` step, keyed at the higher
/// hill's (m/z, RT). Returns `None` when too few samples are collected to be
/// meaningful (see `mz_recalibration_min_samples`).
pub fn learn_recal_model(
    hills: &[Hill],
    config: &FeaturesConfig,
    file: &FileConfig,
) -> Option<MzRecalModel> {
    let features = detect_features_with_recal(hills, config, file, None);

    // Generous residual cap so a spacing that skipped an isotope can't poison
    // the marginal fallback: 4× the user tolerance (ppm), or a fixed 40 ppm
    // for Dalton tolerances where a ppm cap isn't defined.
    let max_abs_ppm = if matches!(file.mz_tolerance_type, ToleranceType::Ppm) {
        file.mz_tolerance * 4.0
    } else {
        40.0
    };
    let mut builder = MzRecalBuilder::new(max_abs_ppm);

    for f in &features {
        if f.charge < 1 || f.hills.len() < 2 {
            continue;
        }
        let step = config.neutron_mass / f.charge as f64;
        // `f.hills` is sorted by m/z ascending and chain extension stops at the
        // first missing isotope, so adjacent pairs are consecutive isotopes.
        // The 0.5·step guard is a cheap belt-and-braces against any gap.
        for w in f.hills.windows(2) {
            let observed = w[1].mz - w[0].mz;
            if (observed - step).abs() > 0.5 * step {
                continue;
            }
            let resid_ppm = (observed - step) / w[1].mz * 1e6;
            builder.add(w[1].mz, w[1].rt, resid_ppm);
        }
    }

    builder.finalize(
        file.mz_recalibration_mz_bins,
        file.mz_recalibration_rt_bins,
        file.mz_recalibration_min_samples,
    )
}

/// Build the best candidate for `seed_idx` across all charge states.
///
/// For each charge, extends isotope chains left and right using chromatographic
/// cosine filtering against the seed (gated by `min_chain_cosine`). Scores each
/// complete chain as `isotope_score × cosine_score` (Bhattacharyya isotope-pattern
/// × mean chromatographic cosine) and keeps the charge that maximises this.
fn generate_best_candidate(
    seed_idx: usize,
    sorted_hills: &[&Hill],
    mz_array: &[f64],
    im_array: &[f64],
    scan_starts: &[usize],
    scan_ends: &[usize],
    config: &FeaturesConfig,
    file: &FileConfig,
    use_im: bool,
    min_intensity: f64,
    recal: Option<&MzRecalModel>,
) -> Candidate {
    // Legacy behaviour: keep the single highest-composite candidate across all
    // charges (iterating high→low so ties favour the higher charge, matching the
    // original `for charge in (..).rev()` loop). The per-charge chain build is
    // factored into `build_charge_candidate`, shared with the exhaustive resolver
    // so both paths build chains identically — greedy vs exhaustive differ only
    // in the conflict-resolution step.
    let mut best = Candidate {
        hill_indices: vec![seed_idx],
        charge: 0,
        mean_cosine: 0.0,
        mean_ppm: 0.0,
        composite_score: 0.0,
    };
    for charge in (config.min_charge..=config.max_charge).rev() {
        if let Some(c) = build_charge_candidate(
            seed_idx,
            charge,
            sorted_hills,
            mz_array,
            im_array,
            scan_starts,
            scan_ends,
            config,
            file,
            use_im,
            min_intensity,
            recal,
        ) {
            if c.composite_score > best.composite_score {
                best = c;
            }
        }
    }
    best
}

/// Build the best isotope-chain candidate for one `(seed, charge)` pair, using
/// the left/right chromatographic-cosine chain extension shared by both the
/// greedy (`generate_best_candidate`) and exhaustive (`resolve_exhaustive`)
/// assemblers. Returns `None` when no isotope partner is found in either
/// direction (a lone seed — charge 0).
#[allow(clippy::too_many_arguments)]
fn build_charge_candidate(
    seed_idx: usize,
    charge: u8,
    sorted_hills: &[&Hill],
    mz_array: &[f64],
    im_array: &[f64],
    scan_starts: &[usize],
    scan_ends: &[usize],
    config: &FeaturesConfig,
    file: &FileConfig,
    use_im: bool,
    min_intensity: f64,
    recal: Option<&MzRecalModel>,
) -> Option<Candidate> {
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
    // Region-adaptive isotope-match tolerance (#2): when recalibration is on
    // and `mz_recalibration_adaptive_tol` is set, derive the ppm window from
    // the surface's per-region residual spread σ instead of the fixed
    // `mz_tolerance`, clamped to [tol_floor_ppm, mz_tolerance]. Isotopes of one
    // feature span only a few Da, so a single σ read at the seed applies to
    // the whole chain. No-op for Dalton tolerances or when `recal` is `None`.
    let effective_tol_ppm = match (use_ppm_tol, file.mz_recalibration_adaptive_tol, recal) {
        (true, true, Some(m)) => {
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
    let kish_on = matches!(file.mz_uncertainty_mode, MzUncertaintyMode::Kish);
    let sigma_mult = file.mz_uncertainty_sigma_mult;
    // Which hill the chromatographic-cosine gate is anchored to. `Adjacent`
    // (default) reproduces legacy koth exactly; `Seed` / `Hybrid` only change
    // the cosine *reference* hill — the m/z step target, `find_neighbors`
    // predecessor, and the intensity-ratio predecessor all stay the immediate
    // chain predecessor in every mode.
    let cosine_anchor = config.cosine_anchor_mode();

    const PROTON_MASS: f64 = 1.007_276_466_621;

    let step = config.neutron_mass / charge as f64;
    {

        // Per-charge averagine template, used to early-stop the right chain
        // when the next theoretical isotope falls below the noise floor.
        // We assume the seed is monoisotopic for this check (best guess at
        // chain-build time; the offset search in scoring/mod.rs corrects later).
        let neutral_mass_seed_mono = seed_mz * charge as f64 - charge as f64 * PROTON_MASS;
        let template = averagine::lookup_template(neutral_mass_seed_mono);
        let template_mono = template[0].max(1e-12);
        let seed_intensity = seed_hill.intensity_max as f64;

        // Right chain: M+1, M+2, ...
        let mut right_chain: Vec<usize> = Vec::new();
        let mut right_cosines: Vec<f64> = Vec::new();
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
            // right_max_decrease / template caps) — see `find_neighbors` and the
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
                target_mz,
                mz_tol,
                ref_hill,
                config.right_max_decrease,
                mz_array,
                im_array,
                scan_starts,
                scan_ends,
                sorted_hills,
                &exclude,
                file,
                use_im,
                kish_on,
                sigma_mult,
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
            // the pre-2026-07 paper output); `Hybrid` uses the seed for the
            // first `cosine_hybrid_depth` isotopes then falls back to adjacent.
            // Only the cosine reference changes — the m/z step and
            // intensity-ratio predecessor stay `ref_hill`.
            let cos_ref_hill = match cosine_anchor {
                CosineAnchor::Adjacent => ref_hill,
                CosineAnchor::Seed => seed_hill,
                CosineAnchor::Hybrid => {
                    if iso <= config.cosine_hybrid_depth {
                        seed_hill
                    } else {
                        ref_hill
                    }
                }
            };
            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| {
                    (
                        c,
                        cosine_similarity(
                            cos_ref_hill,
                            sorted_hills[c],
                            config.cosine_intersection,
                            config.min_scan_overlap,
                        ),
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
                let pred_int = ref_hill.intensity_max as f64;
                let cand_int = sorted_hills[best_c].intensity_max as f64;
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
            right_cosines.push(best_cos);
        }

        // Left chain: M-1, M-2, ...
        let mut left_chain: Vec<usize> = Vec::new();
        let mut left_cosines: Vec<f64> = Vec::new();
        for iso in 1..=config.max_isotopes {
            let target_mz = recal_mz(seed_mz - step * iso as f64);
            let ref_hill = left_chain
                .last()
                .map(|&i| sorted_hills[i])
                .unwrap_or(seed_hill);

            let mut exclude =
                Vec::with_capacity(right_chain.len() + left_chain.len() + 1);
            exclude.push(seed_idx);
            exclude.extend_from_slice(&right_chain);
            exclude.extend_from_slice(&left_chain);

            let cands = find_neighbors(
                target_mz,
                mz_tol,
                ref_hill,
                config.left_max_decrease,
                mz_array,
                im_array,
                scan_starts,
                scan_ends,
                sorted_hills,
                &exclude,
                file,
                use_im,
                kish_on,
                sigma_mult,
            );
            if cands.is_empty() {
                break;
            }

            // Cosine reference hill (left direction). Same anchor semantics as
            // the right chain: `Adjacent` = the immediate predecessor
            // (`ref_hill`, legacy, byte-identical), `Seed` = the mono seed for
            // every isotope, `Hybrid` = seed for the first `cosine_hybrid_depth`
            // then adjacent. Only the cosine reference changes.
            let cos_ref_hill = match cosine_anchor {
                CosineAnchor::Adjacent => ref_hill,
                CosineAnchor::Seed => seed_hill,
                CosineAnchor::Hybrid => {
                    if iso <= config.cosine_hybrid_depth {
                        seed_hill
                    } else {
                        ref_hill
                    }
                }
            };
            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| {
                    (
                        c,
                        cosine_similarity(
                            cos_ref_hill,
                            sorted_hills[c],
                            config.cosine_intersection,
                            config.min_scan_overlap,
                        ),
                    )
                })
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();

            if best_cos < config.min_chain_cosine {
                break;
            }

            // Intensity-ratio gate (left direction). Each left step shifts the
            // hypothesis: we now assume the previous leftmost (ref_hill) sat at
            // template[1] and the new candidate sits at template[0]. The
            // expected ratio is constant across left steps: template[0]/template[1]
            // (new mono / old "near-mono"). Catches cases where the candidate
            // is dimmer than the predecessor — the wrong direction for a mono.
            if config.max_isotope_log2_ratio.is_finite() {
                let pred_int = ref_hill.intensity_max as f64;
                let cand_int = sorted_hills[best_c].intensity_max as f64;
                let theo_pred = template[1].max(1e-12);
                let theo_cand = template[0].max(1e-12);
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

            left_chain.push(best_c);
            left_cosines.push(best_cos);
        }

        if left_chain.is_empty() && right_chain.is_empty() {
            return None; // no partners found — single-hill fallback stays charge=0
        }

        // Assemble chain: left reversed (low→high mz) + seed + right.
        let chain: Vec<usize> = left_chain
            .iter()
            .rev()
            .copied()
            .chain(std::iter::once(seed_idx))
            .chain(right_chain.iter().copied())
            .collect();

        let all_cosines: Vec<f64> = left_cosines
            .iter()
            .rev()
            .copied()
            .chain(right_cosines.iter().copied())
            .collect();
        let mean_cosine = all_cosines.iter().sum::<f64>() / all_cosines.len() as f64;

        let chain_hills: Vec<&Hill> = chain.iter().map(|&i| sorted_hills[i]).collect();
        let isotope_score = score_chain(&chain_hills, charge, config.sulfur_aware_scoring);
        let composite = isotope_score * mean_cosine.max(0.0);

        let step_da = config.neutron_mass / charge as f64;
        let mean_ppm = if chain_hills.len() > 1 {
            let s: f64 = (0..chain_hills.len() - 1)
                .map(|k| {
                    let theo = chain_hills[k].mz + step_da;
                    let expt = chain_hills[k + 1].mz;
                    (expt - theo) / theo * 1e6
                })
                .sum();
            s / (chain_hills.len() - 1) as f64
        } else {
            0.0
        };

        Some(Candidate {
            hill_indices: chain,
            charge,
            mean_cosine,
            mean_ppm,
            composite_score: composite,
        })
    }
}

/// Bhattacharyya score for an isotope chain.  Returns 0.0 for single-peak
/// or charge=0 chains.
// ---------------------------------------------------------------------------
// Experimental exhaustive / non-destructive assembler (behind
// `features.exhaustive_assembly`). See `FeaturesConfig::exhaustive_assembly`.
// ---------------------------------------------------------------------------

/// Recompute `(composite, mean_cosine, isotope_score)` of a (possibly truncated)
/// m/z-sorted chain. Mirrors the generation formula:
/// `isotope_score × mean_adjacent_cosine`.
fn rescore_chain(
    chain: &[usize],
    charge: u8,
    sorted_hills: &[&Hill],
    config: &FeaturesConfig,
) -> (f64, f64, f64) {
    if chain.len() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let hills: Vec<&Hill> = chain.iter().map(|&i| sorted_hills[i]).collect();
    let mut cos_sum = 0.0;
    for k in 0..hills.len() - 1 {
        cos_sum += cosine_similarity(
            hills[k],
            hills[k + 1],
            config.cosine_intersection,
            config.min_scan_overlap,
        );
    }
    let mean_cosine = cos_sum / (hills.len() - 1) as f64;
    let isotope_score = score_chain(&hills, charge, config.sulfur_aware_scoring);
    let composite = isotope_score * mean_cosine.max(0.0);
    (composite, mean_cosine, isotope_score)
}

/// Priority-queue item for the exhaustive resolver.
///
/// Ordering (max-heap: the "best" candidate pops first):
///   1. more hills (longer envelope) first,
///   2. if `iso_priority`: higher isotope-pattern score, then higher composite;
///      otherwise: higher composite score,
///   3. deterministic tie-break: smaller monoisotope index, then smaller charge.
///
/// The deterministic tie-break is essential: iteration order of the seed loop or
/// any hash structure must never influence which candidate wins a contested hill.
/// `iso_priority` is copied from `FeaturesConfig::exhaustive_isotope_priority` and
/// is identical for every item in a heap, so the order stays total.
struct HeapItem {
    chain: Vec<usize>, // m/z-sorted, monoisotope first
    charge: u8,
    composite: f64,
    isotope_score: f64,
    iso_priority: bool,
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
        let f = |x: f64, y: f64| x.partial_cmp(&y).unwrap_or(Ordering::Equal);
        // length ascending in Ord => longer is "greater" => pops first.
        self.chain
            .len()
            .cmp(&other.chain.len())
            // #2: when iso_priority, best averagine fit wins the contested hill
            // within a length class, before falling back to composite.
            .then_with(|| {
                if self.iso_priority {
                    f(self.isotope_score, other.isotope_score)
                } else {
                    Ordering::Equal
                }
            })
            .then_with(|| f(self.composite, other.composite))
            // smaller mono index should pop first => compare as GREATER (reverse).
            .then_with(|| other.mono().cmp(&self.mono()))
            // smaller charge pops first => reverse.
            .then_with(|| other.charge.cmp(&self.charge))
    }
}

/// Non-destructive, confidence-ordered conflict resolution.
///
/// Generates one candidate per `(seed, charge)` with ≥1 isotope partner, then
/// resolves contested hills by envelope length (longest first). A candidate whose
/// isotope hills are partly claimed is truncated to its free monoisotope-anchored
/// prefix, re-scored, and re-queued — only dropped when its monoisotope hill is
/// taken, nothing but the mono survives, or it fails the isotope-score claim gate.
#[allow(clippy::too_many_arguments)]
fn resolve_exhaustive(
    sorted_hills: &[&Hill],
    mz_array: &[f64],
    im_array: &[f64],
    scan_starts: &[usize],
    scan_ends: &[usize],
    config: &FeaturesConfig,
    file: &FileConfig,
    use_im: bool,
    min_intensity: f64,
    recal: Option<&MzRecalModel>,
) -> Vec<Candidate> {
    use std::collections::BinaryHeap;

    // Phase 1+2: over-complete hypothesis pool — one candidate per (seed, charge),
    // built by the SAME `build_charge_candidate` the greedy path uses. Every
    // (seed, charge) is independent and reads only shared immutable arrays, so
    // generation runs in PARALLEL over seeds (mirroring the greedy path's
    // `into_par_iter`), then the results are heapified once via an O(n)
    // `BinaryHeap::from`. The heap's `Ord` is a total order (envelope length,
    // isotope/composite score, then monoisotope-index and charge tie-breaks), so
    // the pop sequence is independent of insertion order — this stays
    // byte-identical to the former sequential push loop. `flat_map_iter` keeps
    // the per-seed charge loop sequential inside each parallel task.
    let items: Vec<HeapItem> = (0..sorted_hills.len())
        .into_par_iter()
        .flat_map_iter(|seed_idx| {
            (config.min_charge..=config.max_charge).filter_map(move |charge| {
                build_charge_candidate(
                    seed_idx, charge, sorted_hills, mz_array, im_array, scan_starts,
                    scan_ends, config, file, use_im, min_intensity, recal,
                )
                .map(|c| {
                    let (composite, _mc, isotope_score) =
                        rescore_chain(&c.hill_indices, c.charge, sorted_hills, config);
                    HeapItem {
                        chain: c.hill_indices,
                        charge: c.charge,
                        composite,
                        isotope_score,
                        iso_priority: config.exhaustive_isotope_priority,
                    }
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
            // #1 quality gate: drop (without claiming) a candidate that can't
            // clear the isotope-pattern bar, so its hills stay free for a
            // better-fitting feature. No-op when the bar is 0.
            if item.isotope_score < config.exhaustive_min_isotope_score {
                continue;
            }
            for &i in &item.chain {
                claimed[i] = true;
            }
            let (composite, mean_cosine, _iso) =
                rescore_chain(&item.chain, item.charge, sorted_hills, config);
            accepted.push(Candidate {
                hill_indices: item.chain,
                charge: item.charge,
                mean_cosine,
                mean_ppm: 0.0,
                composite_score: composite,
            });
        } else if k >= 2 {
            // Truncate to the free prefix, re-score, re-queue at reduced
            // confidence. Drop the remnant if it fails the #1 gate (freeing its
            // hills) rather than re-queuing junk.
            let prefix: Vec<usize> = item.chain[..k].to_vec();
            let (composite, _mc, isotope_score) =
                rescore_chain(&prefix, item.charge, sorted_hills, config);
            if isotope_score >= config.exhaustive_min_isotope_score {
                heap.push(HeapItem {
                    chain: prefix,
                    charge: item.charge,
                    composite,
                    isotope_score,
                    iso_priority: config.exhaustive_isotope_priority,
                });
            }
        }
        // k == 1: only the monoisotope survives — drop (a lone hill is not a
        // charged feature; the greedy path emits charge-0, which output filters).
    }

    accepted
}

fn score_chain(chain_hills: &[&Hill], charge: u8, sulfur_aware: bool) -> f64 {
    if chain_hills.len() <= 1 || charge == 0 {
        return 0.0;
    }

    let mut sorted: Vec<&Hill> = chain_hills.to_vec();
    sorted.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));

    const PROTON_MASS: f64 = 1.007_276_466_621;
    let neutral_mass = sorted[0].mz * charge as f64 - charge as f64 * PROTON_MASS;

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

    let obs: Vec<f64> = sorted.iter().map(|h| h.intensity_at_scan(apex_scan)).collect();
    let k = obs.len().min(10);

    if sulfur_aware {
        averagine::bhattacharyya_score_best_sulfur(&obs[..k], neutral_mass).0
    } else {
        let template = averagine::lookup_template(neutral_mass);
        averagine::bhattacharyya_score(&obs[..k], &template)
    }
}

/// Find candidate isotope partners near `target_mz`.
///
/// Filters: mz tolerance (binary search), not in `exclude_indices` (current
/// chain members), scan overlap with `ref_hill`, IM tolerance, intensity floor.
#[allow(clippy::too_many_arguments)]
fn find_neighbors(
    target_mz: f64,
    mz_tol: f64,
    ref_hill: &Hill,
    max_decrease: f64,
    mz_array: &[f64],
    im_array: &[f64],
    scan_starts: &[usize],
    scan_ends: &[usize],
    sorted_hills: &[&Hill],
    exclude_indices: &[usize],
    file: &FileConfig,
    use_im: bool,
    kish_on: bool,
    sigma_mult: f64,
) -> Vec<usize> {
    // Binary-search window: the flat tolerance by default; widened by
    // KISH_SEARCH_EXPANSION when Kish is on so candidates whose own SE
    // pushes the combined window beyond `mz_tol` aren't pre-filtered out.
    // The exact Kish check inside the loop still rejects anything outside
    // the per-candidate combined tolerance.
    let search_window = if kish_on {
        mz_tol * KISH_SEARCH_EXPANSION
    } else {
        mz_tol
    };

    let lo = mz_array.partition_point(|&x| x < target_mz - search_window);
    let hi = mz_array.partition_point(|&x| x <= target_mz + search_window);

    if lo >= hi {
        return Vec::new();
    }

    let min_intensity = ref_hill.intensity_max * max_decrease;
    let ref_se_term = if kish_on { sigma_mult * ref_hill.mz_se } else { 0.0 };
    let base_tol_sq = mz_tol * mz_tol;
    let ref_se_sq = ref_se_term * ref_se_term;

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
            if kish_on {
                // Combined-quadrature tolerance:
                //   tol² = mz_tol² + (σ·se_ref)² + (σ·se_cand)²
                let cand_se_term = sigma_mult * sorted_hills[i].mz_se;
                let tol_sq = base_tol_sq + ref_se_sq + cand_se_term * cand_se_term;
                let delta = mz_array[i] - target_mz;
                if delta * delta > tol_sq {
                    return false;
                }
            }
            sorted_hills[i].intensity_max >= min_intensity
        })
        .collect()
}

#[cfg(test)]
mod exhaustive_tests {
    use super::*;
    use crate::models::Hill;
    use std::sync::Arc;

    /// Build a co-eluting hill at a given m/z with a bell-shaped profile.
    fn hill(mz: f64, scan_start: usize, profile: Vec<f32>) -> Hill {
        let n = profile.len();
        Hill {
            hill_id: 0,
            mz,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: 0.5,
            rt_start: 0.0,
            rt_end: 1.0,
            rt_width: 1.0,
            im: 0.0,
            im_std: 0.0,
            scan_start,
            scan_apex: scan_start + n / 2,
            scan_end: scan_start + n.saturating_sub(1),
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: profile.iter().map(|&x| x as f64).sum(),
            intensity_max: profile.iter().copied().fold(0.0f32, f32::max) as f64,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
        }
    }

    fn z2_config() -> FeaturesConfig {
        FeaturesConfig {
            min_charge: 2,
            max_charge: 2,
            exhaustive_assembly: true,
            ..Default::default()
        }
    }

    /// A clean, un-contested 3-hill charge-2 envelope must be assembled
    /// identically (same charge, same hill set) by the greedy and the
    /// exhaustive resolvers.
    #[test]
    fn uncontested_envelope_matches_greedy() {
        let step = FeaturesConfig::default().neutron_mass / 2.0;
        let base = 600.0;
        let prof = vec![2.0f32, 6.0, 10.0, 6.0, 2.0];
        let hills = vec![
            hill(base, 20, prof.clone()),
            hill(base + step, 20, vec![1.5, 4.5, 7.5, 4.5, 1.5]),
            hill(base + 2.0 * step, 20, vec![1.0, 3.0, 5.0, 3.0, 1.0]),
        ];
        let file = FileConfig::default();

        let mut cfg_off = z2_config();
        cfg_off.exhaustive_assembly = false;
        let greedy = detect_features(&hills, &cfg_off, &file);

        let cfg_on = z2_config();
        let exhaustive = detect_features(&hills, &cfg_on, &file);

        assert_eq!(greedy.len(), 1, "greedy should find one feature");
        assert_eq!(exhaustive.len(), 1, "exhaustive should find one feature");
        assert_eq!(greedy[0].charge, 2);
        assert_eq!(exhaustive[0].charge, 2);
        assert_eq!(exhaustive[0].hills.len(), 3);
        assert_eq!(greedy[0].hills.len(), exhaustive[0].hills.len());
    }

    /// The exhaustive resolver must be deterministic: repeated runs (and a
    /// reversed input order) produce the identical feature set — guards against
    /// the past HashMap-seed nondeterminism bug.
    #[test]
    fn exhaustive_is_deterministic() {
        let step = FeaturesConfig::default().neutron_mass / 2.0;
        let prof = vec![2.0f32, 6.0, 10.0, 6.0, 2.0];
        let mut hills = Vec::new();
        for base in [500.0, 500.0 + step, 700.0, 900.25] {
            for iso in 0..3 {
                hills.push(hill(base + iso as f64 * step, 20, prof.clone()));
            }
        }
        let file = FileConfig::default();
        let cfg = z2_config();

        let key = |fs: &[Feature]| -> Vec<(u8, Vec<u64>)> {
            let mut v: Vec<(u8, Vec<u64>)> = fs
                .iter()
                .map(|f| {
                    let mut mzs: Vec<u64> = f.hills.iter().map(|h| h.mz.to_bits()).collect();
                    mzs.sort_unstable();
                    (f.charge, mzs)
                })
                .collect();
            v.sort();
            v
        };

        let a = detect_features(&hills, &cfg, &file);
        let b = detect_features(&hills, &cfg, &file);
        assert_eq!(key(&a), key(&b), "exhaustive output must be deterministic");
        hills.reverse();
        let c = detect_features(&hills, &cfg, &file);
        assert_eq!(key(&a), key(&c), "output must not depend on input hill order");
    }

    /// Non-destructive property: no hill is ever claimed by more than one
    /// accepted feature (every hill owned at most once).
    #[test]
    fn no_hill_claimed_twice() {
        let step = FeaturesConfig::default().neutron_mass / 2.0;
        let prof = vec![2.0f32, 6.0, 10.0, 6.0, 2.0];
        let mut hills = Vec::new();
        for base in [400.0, 400.0 + step, 650.0] {
            for iso in 0..4 {
                hills.push(hill(base + iso as f64 * step, 15, prof.clone()));
            }
        }
        let file = FileConfig::default();
        let cfg = z2_config();
        let feats = detect_features(&hills, &cfg, &file);

        let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for f in &feats {
            for h in &f.hills {
                assert!(
                    seen.insert(h.mz.to_bits()),
                    "hill at mz {} claimed by two features",
                    h.mz
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // cosine_anchor (features.cosine_anchor): "adjacent" (legacy) vs "seed".
    // -----------------------------------------------------------------------

    fn anchor_config(anchor: &str) -> FeaturesConfig {
        let mut c = FeaturesConfig {
            min_charge: 2,
            max_charge: 2,
            ..Default::default()
        };
        c.cosine_anchor = anchor.to_string();
        c
    }

    /// A clean, perfectly co-eluting 3-isotope charge-2 envelope is assembled
    /// identically (same charge, same hill count) under `adjacent` and `seed`:
    /// when all isotopes share the seed's elution profile the anchor choice
    /// makes no difference.
    #[test]
    fn clean_envelope_identical_under_adjacent_and_seed() {
        let step = FeaturesConfig::default().neutron_mass / 2.0;
        let base = 600.0;
        // All three hills share scans 20..=24 with the same bell shape, so
        // every pairwise cosine (seed-anchored or adjacent) is ~1.
        let hills = vec![
            hill(base, 20, vec![2.0, 6.0, 10.0, 6.0, 2.0]),
            hill(base + step, 20, vec![1.5, 4.5, 7.5, 4.5, 1.5]),
            hill(base + 2.0 * step, 20, vec![1.0, 3.0, 5.0, 3.0, 1.0]),
        ];
        let file = FileConfig::default();

        let adj = detect_features(&hills, &anchor_config("adjacent"), &file);
        let seed = detect_features(&hills, &anchor_config("seed"), &file);

        assert_eq!(adj.len(), 1, "adjacent should find one feature");
        assert_eq!(seed.len(), 1, "seed should find one feature");
        assert_eq!(adj[0].charge, 2);
        assert_eq!(seed[0].charge, 2);
        assert_eq!(adj[0].hills.len(), 3, "adjacent keeps all 3 isotopes");
        assert_eq!(seed[0].hills.len(), 3, "seed keeps all 3 isotopes");
    }

    /// A far isotope (M+2) whose apex has drifted so it co-elutes with its
    /// neighbour (M+1) but NOT with the mono seed (M0): accepted under
    /// `adjacent` (cosine vs M+1), rejected under `seed` (cosine vs M0).
    /// Demonstrates that the two anchors genuinely differ.
    #[test]
    fn drifting_far_isotope_kept_by_adjacent_dropped_by_seed() {
        let step = FeaturesConfig::default().neutron_mass / 2.0;
        let base = 600.0;
        // Same bell SHAPE, apex drifting +2 scans per isotope, amplitude tapered
        // ~0.6× each step (so the intensity-ratio gate is comfortably satisfied;
        // cosine is scale-invariant so the shape correlations are unchanged).
        // Adjacent neighbours overlap at lag 2 (autocorr cos ≈ 0.75 ≥ the 0.5
        // min_chain_cosine); the mono↔M+2 pair overlaps at lag 4 (cos ≈ 0.31,
        // below the gate).
        let shape = [1.0f32, 3.0, 6.0, 9.0, 10.0, 9.0, 6.0, 3.0, 1.0];
        let scaled = |k: f32| -> Vec<f32> { shape.iter().map(|&v| v * k).collect() };
        let hills = vec![
            hill(base, 20, scaled(1.0)),               // M0:  scans 20..=28, apex 24
            hill(base + step, 22, scaled(0.6)),        // M+1: scans 22..=30, apex 26
            hill(base + 2.0 * step, 24, scaled(0.36)), // M+2: scans 24..=32, apex 28
        ];
        let file = FileConfig::default();

        // Hills are already m/z-ascending, so the mono seed is index 0. Drive
        // build_charge_candidate directly from that seed so the greedy
        // re-seeding (which could pick M+1 as a fresh mono) can't mask the
        // anchor difference.
        let sorted: Vec<&Hill> = hills.iter().collect();
        let mz: Vec<f64> = sorted.iter().map(|h| h.mz).collect();
        let im: Vec<f64> = sorted.iter().map(|h| h.im).collect();
        let ss: Vec<usize> = sorted.iter().map(|h| h.scan_start).collect();
        let se: Vec<usize> = sorted.iter().map(|h| h.scan_end).collect();

        // Disable the predicted-intensity noise-floor early-stop so the cosine
        // anchor is the ONLY thing that can differ between the two runs.
        let mut cfg_adj = anchor_config("adjacent");
        cfg_adj.chain_predicted_intensity_gate = false;
        let mut cfg_seed = anchor_config("seed");
        cfg_seed.chain_predicted_intensity_gate = false;

        let cand_adj = build_charge_candidate(
            0, 2, &sorted, &mz, &im, &ss, &se, &cfg_adj, &file, false, 0.0, None,
        )
        .expect("adjacent should build a chain from the mono seed");
        let cand_seed = build_charge_candidate(
            0, 2, &sorted, &mz, &im, &ss, &se, &cfg_seed, &file, false, 0.0, None,
        )
        .expect("seed should build a chain from the mono seed");

        assert_eq!(
            cand_adj.hill_indices.len(),
            3,
            "adjacent: M+2 co-elutes with its neighbour M+1 → full envelope"
        );
        assert_eq!(
            cand_seed.hill_indices.len(),
            2,
            "seed: M+2 scored against the mono (cos ≈ 0.31 < 0.5) → dropped"
        );
    }

    /// An unrecognised `cosine_anchor` value is rejected by config validation.
    #[test]
    fn unknown_cosine_anchor_is_a_config_error() {
        assert!(
            anchor_config("banana").validate().is_err(),
            "unknown cosine_anchor must be a config error"
        );
        for v in ["adjacent", "seed", "hybrid", "SEED", "Hybrid"] {
            assert!(
                anchor_config(v).validate().is_ok(),
                "`{v}` should be an accepted cosine_anchor"
            );
        }
    }
}
