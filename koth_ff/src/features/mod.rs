pub mod cosine;
pub mod recalibration;

mod assemble;

use std::cmp::Ordering;

use crate::config::{FeaturesConfig, FileConfig, ToleranceType};
use crate::models::{Feature, Hill};
use assemble::{resolve_exhaustive, ChainCtx};
use cosine::cosine_similarity;
use recalibration::{MzRecalBuilder, MzRecalModel};

/// Detect isotope features from a list of chromatographic hills.
///
/// 1. Build an m/z-sorted index of the hills.
/// 2. Every `(seed hill, charge)` pair builds an isotope chain by extending
///    upward only (M+1, M+2, …) with chromatographic-cosine filtering (gated by
///    `FeaturesConfig.min_chain_cosine`); the seed is the monoisotope
///    hypothesis — see `build_charge_candidate`.
/// 3. The over-complete candidate pool is resolved non-destructively by
///    `resolve_exhaustive`: contested hills are claimed longest-envelope-first
///    and a partly-claimed candidate is truncated to its free prefix and
///    re-queued rather than dropped.
/// 4. Feature structs are built from the accepted candidates.
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

    // Non-destructive assembler (biosaur2 / AlphaPept style): an over-complete
    // (seed, charge) hypothesis pool, generated in parallel and resolved by
    // claiming contested hills longest-envelope-first, truncating a
    // partly-claimed candidate to its free monoisotope-anchored prefix rather
    // than dropping it. See `resolve_exhaustive`.
    let ctx = ChainCtx {
        sorted_hills: &sorted_hills,
        mz_array: &mz_array,
        im_array: &im_array,
        scan_starts: &scan_starts,
        scan_ends: &scan_ends,
        config,
        file,
        use_im,
        min_intensity,
        recal,
    };
    let accepted = resolve_exhaustive(&ctx);

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
