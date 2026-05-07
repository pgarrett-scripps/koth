pub mod cosine;

use std::cmp::Ordering;

use crate::config::{FeaturesConfig, FileConfig, ImToleranceType, ToleranceType};
use crate::models::{Feature, Hill};
use crate::scoring::averagine;
use cosine::cosine_similarity;

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
///    chain is built via cosine-similarity filtering.  Each (seed, charge) pair
///    is scored (Bhattacharyya × mean_cosine) and the charge with the highest
///    composite score is kept as that seed's candidate.
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
    let candidates: Vec<Candidate> = (0..sorted_hills.len())
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
                cosine_similarity: mean_cosine,
                ppm_error: mean_ppm,
            }
        })
        .collect();

    log::info!("Detected {} isotope features", features.len());
    features
}

/// Build the best candidate for `seed_idx` across all charge states.
///
/// For each charge, extends isotope chains left and right using cosine-similarity
/// filtering (same logic as before).  Scores each complete chain with
/// Bhattacharyya × mean_cosine and keeps the charge that maximises this.
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

    let mut best = Candidate {
        hill_indices: vec![seed_idx],
        charge: 0,
        mean_cosine: 0.0,
        mean_ppm: 0.0,
        composite_score: 0.0,
    };

    for charge in (config.min_charge..=config.max_charge).rev() {
        let step = config.neutron_mass / charge as f64;

        // Right chain: M+1, M+2, ...
        let mut right_chain: Vec<usize> = Vec::new();
        let mut right_cosines: Vec<f64> = Vec::new();
        for iso in 1..=config.max_isotopes {
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
            );
            if cands.is_empty() {
                break;
            }

            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| (c, cosine_similarity(sorted_hills[seed_idx], sorted_hills[c])))
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();

            if best_cos < config.min_cosine_similarity {
                break;
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
            );
            if cands.is_empty() {
                break;
            }

            let (best_c, best_cos) = cands
                .iter()
                .map(|&c| (c, cosine_similarity(sorted_hills[seed_idx], sorted_hills[c])))
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();

            if best_cos < config.min_cosine_similarity {
                break;
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
        let bhat = score_chain(&chain_hills, charge, min_intensity);
        let composite = bhat * mean_cosine.max(0.0);

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
fn score_chain(chain_hills: &[&Hill], charge: u8, min_intensity: f64) -> f64 {
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
    let template = averagine::lookup_template(neutral_mass);
    averagine::bhattacharyya_score(&obs[..k], template, min_intensity)
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
) -> Vec<usize> {
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
