pub mod averagine;

use crate::config::ScoringConfig;
use crate::models::{Feature, ScoredFeature};
use averagine::{bhattacharyya_score, lookup_template};

/// Score features using averagine theoretical isotope patterns.
///
/// For each feature:
/// 1. Look up theoretical isotope pattern by neutral mass
/// 2. Try neutron offsets in [offset_min, offset_max]
/// 3. Score each offset with Bhattacharyya coefficient + zero-offset bonus
/// 4. Keep best offset; if below min_score_threshold, fall back to offset=0
pub fn score_features(features: &[Feature], config: &ScoringConfig) -> Vec<ScoredFeature> {
    log::info!("Scoring {} features", features.len());

    // Compute min_intensity from 5th percentile of hill intensity_sum
    let mut all_intensities: Vec<f64> = features
        .iter()
        .flat_map(|f| f.hills.iter().map(|h| h.intensity_sum))
        .collect();
    all_intensities.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min_intensity = if all_intensities.is_empty() {
        0.0
    } else {
        let p5_idx = (all_intensities.len() as f64 * 0.05) as usize;
        all_intensities[p5_idx] * 0.8
    };
    eprintln!("[koth_ff] score min_intensity: {:.3e} (from {} hill intensities, p5_idx={})",
        min_intensity, all_intensities.len(),
        if all_intensities.is_empty() { 0 } else { (all_intensities.len() as f64 * 0.05) as usize });

    let all_scored: Vec<ScoredFeature> = features
        .iter()
        .map(|feature| score_one(feature, config, min_intensity))
        .collect();

    let charged: Vec<&ScoredFeature> = all_scored.iter().filter(|sf| sf.feature.charge > 0).collect();
    if !charged.is_empty() {
        let max_score = charged.iter().map(|sf| sf.score).fold(f64::NEG_INFINITY, f64::max);
        let mean_score = charged.iter().map(|sf| sf.score).sum::<f64>() / charged.len() as f64;
        let above = charged.iter().filter(|sf| sf.score >= config.min_score_threshold).count();
        eprintln!(
            "[koth_ff] score diag: {} charged features, max={:.3}, mean={:.3}, above_thresh({}): {}",
            charged.len(), max_score, mean_score, config.min_score_threshold, above
        );
        // Print first 5 charged features for inspection
        for sf in charged.iter().take(5) {
            let obs = sf.feature.isotope_profile_apex();
            eprintln!("  charge={} n_isotopes={} obs={:?} score={:.4}",
                sf.feature.charge, sf.feature.hills.len(),
                obs.iter().map(|x| format!("{:.0}", x)).collect::<Vec<_>>(),
                sf.score);
        }
    }

    // Return all charged features; output filtering by score is applied by
    // the caller via FeaturesConfig::min_score.
    let scored: Vec<ScoredFeature> = all_scored
        .into_iter()
        .filter(|sf| sf.feature.charge > 0)
        .collect();

    log::info!("Scored {} features", scored.len());
    scored
}

fn score_one(feature: &Feature, config: &ScoringConfig, min_intensity: f64) -> ScoredFeature {
    let neutral_mass = match feature.monoisotopic_neutral_mass() {
        Some(m) => m,
        None => {
            // Unknown charge — return unscored
            return ScoredFeature {
                feature: feature.clone(),
                neutron_offset: 0,
                score: 0.0,
                theoretical_pattern: Vec::new(),
            };
        }
    };

    let template = lookup_template(neutral_mass);
    let obs = feature.isotope_profile_apex();
    let k = obs.len().min(10);

    let mut best_combined = f64::NEG_INFINITY;
    let mut best_bc = 0.0f64;
    let mut best_offset: i8 = 0;

    for o in config.isotope_offset_min..=config.isotope_offset_max {
        // Shift obs: obs_aligned[j] = obs[j + o] (shift left by o)
        let mut obs_aligned = vec![0.0f64; k];
        for j in 0..k {
            let src = j as i32 + o as i32;
            if src >= 0 && (src as usize) < k {
                obs_aligned[j] = obs[src as usize];
            }
        }

        let sc = bhattacharyya_score(&obs_aligned, template, min_intensity);
        // Bonus is used only to prefer offset=0 when scores are close; never stored.
        let combined = sc + if o == 0 { config.offset_zero_bonus } else { 0.0 };
        if combined > best_combined {
            best_combined = combined;
            best_bc = sc;
            best_offset = o;
        }
    }

    // If below threshold, reset to offset=0
    if best_combined < config.min_score_threshold {
        best_offset = 0;
        best_bc = bhattacharyya_score(&obs, template, min_intensity);
    }

    // Build normalized theoretical pattern for output
    let theo_slice = &template[..k];
    let theo_sum: f64 = theo_slice.iter().sum();
    let theoretical_pattern: Vec<f64> = if theo_sum > 0.0 {
        theo_slice.iter().map(|&x| x / theo_sum).collect()
    } else {
        vec![0.0; k]
    };

    ScoredFeature {
        feature: feature.clone(),
        neutron_offset: best_offset,
        score: best_bc.clamp(0.0, 1.0),
        theoretical_pattern,
    }
}
