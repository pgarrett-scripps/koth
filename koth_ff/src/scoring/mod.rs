pub mod averagine;
pub mod elements;

use crate::config::ScoringConfig;
use crate::models::{Feature, ScoredFeature};
use averagine::{bhattacharyya_score, lookup_template};

/// Score features using averagine theoretical isotope patterns.
///
/// For each feature:
/// 1. Look up theoretical isotope pattern by neutral mass
/// 2. Try neutron offsets in [offset_min, offset_max]
/// 3. Score each offset with Bhattacharyya coefficient + zero-offset bonus
/// 4. Keep best offset; if below `min_isotope_score_for_offset`, fall back to offset=0
pub fn score_features(
    features: &[Feature],
    config: &ScoringConfig,
    sulfur_aware: bool,
) -> Vec<ScoredFeature> {
    log::info!("Scoring {} features", features.len());

    let all_scored: Vec<ScoredFeature> = features
        .iter()
        .map(|feature| score_one(feature, config, sulfur_aware))
        .collect();

    let charged: Vec<&ScoredFeature> = all_scored.iter().filter(|sf| sf.feature.charge > 0).collect();
    if !charged.is_empty() {
        let max_score = charged.iter().map(|sf| sf.combined_score).fold(f64::NEG_INFINITY, f64::max);
        let mean_score = charged.iter().map(|sf| sf.combined_score).sum::<f64>() / charged.len() as f64;
        let above = charged
            .iter()
            .filter(|sf| sf.isotope_score >= config.min_isotope_score_for_offset)
            .count();
        eprintln!(
            "[koth_ff] score diag: {} charged features, max_combined={:.3}, mean_combined={:.3}, isotope>={}: {}",
            charged.len(), max_score, mean_score, config.min_isotope_score_for_offset, above
        );
        // Print first 5 charged features for inspection
        for sf in charged.iter().take(5) {
            let obs = sf.feature.isotope_profile_apex();
            eprintln!(
                "  charge={} n_isotopes={} obs={:?} iso={:.4} cos={:.4} comb={:.4}",
                sf.feature.charge,
                sf.feature.hills.len(),
                obs.iter().map(|x| format!("{:.0}", x)).collect::<Vec<_>>(),
                sf.isotope_score,
                sf.cosine_score,
                sf.combined_score,
            );
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

fn score_one(feature: &Feature, config: &ScoringConfig, sulfur_aware: bool) -> ScoredFeature {
    let neutral_mass = match feature.monoisotopic_neutral_mass() {
        Some(m) => m,
        None => {
            // Unknown charge — return unscored
            return ScoredFeature {
                feature: feature.clone(),
                neutron_offset: 0,
                isotope_score: 0.0,
                cosine_score: 0.0,
                combined_score: 0.0,
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

    let score_obs = |o: &[f64]| -> f64 {
        if sulfur_aware {
            averagine::bhattacharyya_score_best_sulfur(o, neutral_mass).0
        } else {
            bhattacharyya_score(o, &template)
        }
    };

    for o in config.isotope_offset_min..=config.isotope_offset_max {
        // Shift obs: obs_aligned[j] = obs[j + o] (shift left by o)
        let mut obs_aligned = vec![0.0f64; k];
        for j in 0..k {
            let src = j as i32 + o as i32;
            if src >= 0 && (src as usize) < k {
                obs_aligned[j] = obs[src as usize];
            }
        }

        let sc = score_obs(&obs_aligned);
        // Bonus is used only to prefer offset=0 when scores are close; never stored.
        let combined = sc + if o == 0 { config.offset_zero_bonus } else { 0.0 };
        if combined > best_combined {
            best_combined = combined;
            best_bc = sc;
            best_offset = o;
        }
    }

    // If below the isotope-score floor, reset to offset=0 (no neutron reassignment).
    if best_combined < config.min_isotope_score_for_offset {
        best_offset = 0;
        best_bc = score_obs(&obs);
    }

    // Build normalized theoretical pattern for output
    let theo_slice = &template[..k];
    let theo_sum: f64 = theo_slice.iter().sum();
    let theoretical_pattern: Vec<f64> = if theo_sum > 0.0 {
        theo_slice.iter().map(|&x| x / theo_sum).collect()
    } else {
        vec![0.0; k]
    };

    // Three quality scores carried independently so downstream filters can
    // gate on any one of them:
    //   isotope_score  — Bhattacharyya isotope-pattern match (this stage)
    //   cosine_score   — chromatographic co-elution of isotope hills (feature stage)
    //   combined_score — product, the default "quality" knob
    //
    // The chromatographic cosine is now computed over the mutual scan
    // OVERLAP only (see features::cosine::cosine_similarity), which
    // removes the systematic penalty against short low-abundance hills
    // that the previous union-padded form imposed. With that fix the
    // multiplicative composition cleanly separates real features from
    // noise on the PXD003881 ground-truth set.
    let isotope_score = best_bc.clamp(0.0, 1.0);
    let cosine_score = feature.cosine_score.clamp(0.0, 1.0);
    let combined_score = (isotope_score * cosine_score).clamp(0.0, 1.0);

    ScoredFeature {
        feature: feature.clone(),
        neutron_offset: best_offset,
        isotope_score,
        cosine_score,
        combined_score,
        theoretical_pattern,
    }
}
