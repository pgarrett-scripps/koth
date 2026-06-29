pub mod cosine;

use std::cmp::Ordering;

use rayon::prelude::*;

use crate::config::{FeaturesConfig, FileConfig, ImToleranceType, MzUncertaintyMode, ToleranceType};
use crate::models::{Feature, Hill};
use crate::scoring::averagine;
use cosine::cosine_similarity;

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
    if hills.is_empty() {
        return Vec::new();
    }

    log::info!("Detecting isotope features from {} hills", hills.len());

    let mut order: Vec<usize> = (0..hills.len()).collect();
    order.sort_by(|&a, &b| hills[a].mz.partial_cmp(&hills[b].mz).unwrap());
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

    // Phase 1+2: generate and score one candidate per seed (all hills).
    // Each seed's candidate is computed independently from read-only shared
    // arrays, so this is embarrassingly parallel. `into_par_iter().collect()`
    // is order-preserving (rayon's indexed collect), so `candidates` ends up
    // in the exact same order as the sequential version — the downstream sort
    // + greedy claim are unchanged and the output stays deterministic.
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
    let accepted: Vec<Candidate> = candidates
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
        .collect();

    // Phase 4: build Feature structs.
    let features: Vec<Feature> = accepted
        .into_iter()
        .map(|c| {
            let mut chain_hills: Vec<Hill> =
                c.hill_indices.iter().map(|&i| (*sorted_hills[i]).clone()).collect();
            chain_hills.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

            let cosine_sims: Vec<f64> = (0..chain_hills.len().saturating_sub(1))
                .map(|k| cosine_similarity(&chain_hills[k], &chain_hills[k + 1]))
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
) -> Candidate {
    let seed_hill = sorted_hills[seed_idx];
    let seed_mz = seed_hill.mz;
    let mz_tol = if matches!(file.mz_tolerance_type, ToleranceType::Ppm) {
        seed_mz * file.mz_tolerance / 1e6
    } else {
        file.mz_tolerance
    };
    let kish_on = matches!(file.mz_uncertainty_mode, MzUncertaintyMode::Kish);
    let sigma_mult = file.mz_uncertainty_sigma_mult;

    let mut best = Candidate {
        hill_indices: vec![seed_idx],
        charge: 0,
        mean_cosine: 0.0,
        mean_ppm: 0.0,
        composite_score: 0.0,
    };

    const PROTON_MASS: f64 = 1.007_276_466_621;

    for charge in (config.min_charge..=config.max_charge).rev() {
        let step = config.neutron_mass / charge as f64;

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
            let expected = seed_intensity * (template[iso] / template_mono);
            if expected < min_intensity {
                break;
            }

            let target_mz = seed_mz + step * iso as f64;
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

            // Cosine vs the immediate predecessor in the chain (ref_hill),
            // not vs the seed. For M+1 the predecessor IS the seed; for
            // M+k≥2 it's the previously-claimed isotope hill. Verifies
            // adjacent isotopes co-elute, which is the actual physical
            // constraint — chains drift in S/N from the seed as you go
            // out, so a seed-anchored cosine over-rejects far isotopes.
            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| (c, cosine_similarity(ref_hill, sorted_hills[c])))
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
            let target_mz = seed_mz - step * iso as f64;
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

            // Cosine vs the immediate predecessor in the chain (ref_hill),
            // not vs the seed. For M+1 the predecessor IS the seed; for
            // M+k≥2 it's the previously-claimed isotope hill. Verifies
            // adjacent isotopes co-elute, which is the actual physical
            // constraint — chains drift in S/N from the seed as you go
            // out, so a seed-anchored cosine over-rejects far isotopes.
            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| (c, cosine_similarity(ref_hill, sorted_hills[c])))
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
            continue; // no partners found — single-hill fallback stays charge=0
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

        if composite > best.composite_score {
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

            best = Candidate {
                hill_indices: chain,
                charge,
                mean_cosine,
                mean_ppm,
                composite_score: composite,
            };
        }
    }

    best
}

/// Bhattacharyya score for an isotope chain.  Returns 0.0 for single-peak
/// or charge=0 chains.
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
