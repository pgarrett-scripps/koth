pub mod cosine;

use crate::config::{FeaturesConfig, FileConfig, ImToleranceType, ToleranceType};
use crate::models::{Feature, Hill};
use cosine::cosine_similarity;

/// Detect isotope features from a list of chromatographic hills.
///
/// Algorithm (mirrors zenith_feature_finder/features.py detect_features()):
/// 1. Sort hills by mz; build arrays for binary search
/// 2. Process seeds in descending intensity_max order
/// 3. For each unassigned seed, try all charge states (high → low)
/// 4. Search right (M+1, M+2...) and left (M-1, M-2...) for isotope partners
/// 5. Keep the longest chain; assign all chain hills
/// 6. Build Feature from chain
pub fn detect_features(hills: &[Hill], config: &FeaturesConfig, file: &FileConfig) -> Vec<Feature> {
    if hills.is_empty() {
        return Vec::new();
    }

    log::info!("Detecting isotope features from {} hills", hills.len());

    // Build mz-sorted index
    let mut order: Vec<usize> = (0..hills.len()).collect();
    order.sort_by(|&a, &b| hills[a].mz.partial_cmp(&hills[b].mz).unwrap());
    let sorted_hills: Vec<&Hill> = order.iter().map(|&i| &hills[i]).collect();
    let mz_array: Vec<f64> = sorted_hills.iter().map(|h| h.mz).collect();
    let im_array: Vec<f64> = sorted_hills.iter().map(|h| h.im).collect();
    let scan_starts: Vec<usize> = sorted_hills.iter().map(|h| h.scan_start).collect();
    let scan_ends: Vec<usize> = sorted_hills.iter().map(|h| h.scan_end).collect();

    let use_im = im_array.iter().any(|&x| x != 0.0);
    let use_ppm = matches!(file.mz_tolerance_type, ToleranceType::Ppm);

    let mut assigned = vec![false; sorted_hills.len()];

    // Process by descending intensity_max
    let mut intensity_order: Vec<usize> = (0..sorted_hills.len()).collect();
    intensity_order.sort_by(|&a, &b| {
        sorted_hills[b]
            .intensity_max
            .partial_cmp(&sorted_hills[a].intensity_max)
            .unwrap()
    });

    let mut features: Vec<Feature> = Vec::new();

    for &seed_idx in &intensity_order {
        if assigned[seed_idx] {
            continue;
        }

        let seed_hill = sorted_hills[seed_idx];
        let seed_mz = seed_hill.mz;
        let mz_tol = if use_ppm {
            seed_mz * file.mz_tolerance / 1e6
        } else {
            file.mz_tolerance
        };

        let mut best_chain: Vec<usize> = vec![seed_idx];
        let mut best_charge: u8 = 0;

        for charge in (config.min_charge..=config.max_charge).rev() {
            let step = config.neutron_mass / charge as f64;

            // Search right (M+1, M+2, ...)
            let mut right_chain: Vec<usize> = Vec::new();
            for iso in 1..=config.max_isotopes {
                let target_mz = seed_mz + step * iso as f64;
                let ref_hill = right_chain
                    .last()
                    .map(|&i| sorted_hills[i])
                    .unwrap_or(seed_hill);

                let candidates = find_neighbors(
                    target_mz,
                    mz_tol,
                    seed_idx,
                    ref_hill,
                    config.right_max_decrease,
                    &mz_array,
                    &im_array,
                    &scan_starts,
                    &scan_ends,
                    &sorted_hills,
                    &assigned,
                    file,
                    use_im,
                );

                if candidates.is_empty() {
                    break;
                }

                // Compute cosine for each candidate locally — no global cache
                let mut best_c = candidates[0];
                let mut best_cos = f64::NEG_INFINITY;
                for &c in &candidates {
                    let cos = cosine_similarity(sorted_hills[seed_idx], sorted_hills[c]);
                    if cos > best_cos {
                        best_cos = cos;
                        best_c = c;
                    }
                }
                if best_cos < config.min_cosine_similarity {
                    break;
                }
                right_chain.push(best_c);
            }

            // Search left (M-1, M-2, ...)
            let mut left_chain: Vec<usize> = Vec::new();
            for iso in 1..=config.max_isotopes {
                let target_mz = seed_mz - step * iso as f64;
                let ref_hill = left_chain
                    .last()
                    .map(|&i| sorted_hills[i])
                    .unwrap_or(seed_hill);

                let candidates = find_neighbors(
                    target_mz,
                    mz_tol,
                    seed_idx,
                    ref_hill,
                    config.left_max_decrease,
                    &mz_array,
                    &im_array,
                    &scan_starts,
                    &scan_ends,
                    &sorted_hills,
                    &assigned,
                    file,
                    use_im,
                );

                if candidates.is_empty() {
                    break;
                }

                // Compute cosine for each candidate locally — no global cache
                let mut best_c = candidates[0];
                let mut best_cos = f64::NEG_INFINITY;
                for &c in &candidates {
                    let cos = cosine_similarity(sorted_hills[seed_idx], sorted_hills[c]);
                    if cos > best_cos {
                        best_cos = cos;
                        best_c = c;
                    }
                }
                if best_cos < config.min_cosine_similarity {
                    break;
                }
                left_chain.push(best_c);
            }

            // Full chain: left (reversed, low→high mz) + seed + right
            let chain: Vec<usize> = left_chain
                .iter()
                .rev()
                .copied()
                .chain(std::iter::once(seed_idx))
                .chain(right_chain.iter().copied())
                .collect();

            if chain.len() > best_chain.len() {
                best_chain = chain;
                best_charge = charge;
            }
        }

        // Mark all hills in chain as assigned
        for &i in &best_chain {
            assigned[i] = true;
        }

        // Sort chain hills by mz
        let mut chain_hills: Vec<&Hill> = best_chain.iter().map(|&i| sorted_hills[i]).collect();
        chain_hills.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

        // Compute per-adjacent-pair cosine similarities and ppm errors
        let cosine_sims: Vec<f64> = (0..chain_hills.len().saturating_sub(1))
            .map(|k| cosine_similarity(chain_hills[k], chain_hills[k + 1]))
            .collect();

        let ppm_errors: Vec<f64> = if best_charge > 0 && chain_hills.len() > 1 {
            let step = config.neutron_mass / best_charge as f64;
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

        let charge = if best_chain.len() > 1 { best_charge } else { 0 };

        features.push(Feature {
            hills: chain_hills.into_iter().cloned().collect(),
            charge,
            cosine_similarity: mean_cosine,
            ppm_error: mean_ppm,
        });
    }

    log::info!("Detected {} isotope features", features.len());
    features
}

/// Find candidate isotope partners near target_mz.
///
/// Filters by:
/// - within mz_tol of target_mz (binary search)
/// - not already assigned
/// - not the seed itself
/// - scan overlap with ref_hill
/// - IM within tolerance (if IM is in use)
/// - intensity >= ref_hill.intensity_max * max_decrease
#[allow(clippy::too_many_arguments)]
fn find_neighbors<'a>(
    target_mz: f64,
    mz_tol: f64,
    seed_idx: usize,
    ref_hill: &Hill,
    max_decrease: f64,
    mz_array: &[f64],
    im_array: &[f64],
    scan_starts: &[usize],
    scan_ends: &[usize],
    sorted_hills: &[&'a Hill],
    assigned: &[bool],
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
            if assigned[i] || i == seed_idx {
                return false;
            }
            // Scan overlap
            if scan_ends[i] < ref_hill.scan_start || scan_starts[i] > ref_hill.scan_end {
                return false;
            }
            // IM tolerance
            if use_im {
                let im_tol = match file.im_tolerance_type {
                    ImToleranceType::Absolute => file.im_tolerance,
                    ImToleranceType::Relative => {
                        ref_hill.im.max(im_array[i]) * file.im_tolerance
                    }
                };
                if (im_array[i] - ref_hill.im).abs() > im_tol {
                    return false;
                }
            }
            // Intensity decrease
            sorted_hills[i].intensity_max >= min_intensity
        })
        .collect()
}
