use super::*;
use crate::config::FileConfig;
use crate::features::detect_features;
use crate::models::{Feature, Hill};
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
        faims_cv: None,
    }
}

fn z2_config() -> FeaturesConfig {
    FeaturesConfig {
        min_charge: 2,
        max_charge: 2,
        ..Default::default()
    }
}

#[test]
fn public_pipeline_never_assembles_isotopes_across_faims_channels() {
    let step = FeaturesConfig::default().neutron_mass / 2.0;
    let profile = vec![100.0, 500.0, 1_000.0, 500.0, 100.0];
    let mut mono = hill(500.0, 0, profile.clone());
    mono.faims_cv = Some(-50.0);
    let mut isotope = hill(500.0 + step, 0, profile);
    isotope.faims_cv = Some(-65.0);

    let features = crate::run_features(
        &[mono.clone(), isotope.clone()],
        &z2_config(),
        &FileConfig::default(),
    )
    .expect("feature detection");
    assert!(
        features.iter().all(|f| f.hills.len() == 1),
        "isotope-spaced hills from different CVs formed a feature"
    );

    isotope.faims_cv = mono.faims_cv;
    let same_cv = crate::run_features(&[mono, isotope], &z2_config(), &FileConfig::default())
        .expect("same-CV feature detection");
    assert!(
        same_cv.iter().any(|f| f.charge == 2 && f.hills.len() == 2),
        "same-CV isotope pair should still assemble"
    );
}

/// A clean, un-contested 3-hill charge-2 envelope is assembled into a single
/// charge-2 feature spanning all three hills.
#[test]
fn uncontested_envelope_assembled() {
    let step = FeaturesConfig::default().neutron_mass / 2.0;
    let base = 600.0;
    let prof = vec![2.0f32, 6.0, 10.0, 6.0, 2.0];
    let hills = vec![
        hill(base, 20, prof.clone()),
        hill(base + step, 20, vec![1.5, 4.5, 7.5, 4.5, 1.5]),
        hill(base + 2.0 * step, 20, vec![1.0, 3.0, 5.0, 3.0, 1.0]),
    ];
    let file = FileConfig::default();

    let features = detect_features(&hills, &z2_config(), &file);

    assert_eq!(features.len(), 1, "should find one feature");
    assert_eq!(features[0].charge, 2);
    assert_eq!(features[0].hills.len(), 3);
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
    assert_eq!(
        key(&a),
        key(&c),
        "output must not depend on input hill order"
    );
}

/// Non-destructive property: no hill is ever claimed by more than one
/// accepted feature (every hill owned at most once).
#[test]
fn no_hill_claimed_twice() {
    let step = FeaturesConfig::default().neutron_mass / 2.0;
    let prof = vec![2.0f32, 6.0, 10.0, 6.0, 2.0];
    // The 400.0 and 400.0 + step groups overlap, so several distinct hills
    // share an m/z. Identity must therefore be `hill_id`, not m/z: two
    // features legitimately holding two *different* hills that happen to sit
    // at the same m/z is not a double-claim.
    let mut hills = Vec::new();
    for base in [400.0, 400.0 + step, 650.0] {
        for iso in 0..4 {
            let mut h = hill(base + iso as f64 * step, 15, prof.clone());
            h.hill_id = hills.len() as u64;
            hills.push(h);
        }
    }
    let file = FileConfig::default();
    let cfg = z2_config();
    let feats = detect_features(&hills, &cfg, &file);

    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for f in &feats {
        for h in &f.hills {
            assert!(
                seen.insert(h.hill_id),
                "hill {} (mz {}) claimed by two features",
                h.hill_id,
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
    // build_charge_candidate directly from that seed so the resolver's
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

    let ctx_adj = ChainCtx {
        sorted_hills: &sorted,
        mz_array: &mz,
        im_array: &im,
        scan_starts: &ss,
        scan_ends: &se,
        config: &cfg_adj,
        file: &file,
        use_im: false,
        min_intensity: 0.0,
        recal: None,
    };
    let ctx_seed = ChainCtx {
        sorted_hills: &sorted,
        mz_array: &mz,
        im_array: &im,
        scan_starts: &ss,
        scan_ends: &se,
        config: &cfg_seed,
        file: &file,
        use_im: false,
        min_intensity: 0.0,
        recal: None,
    };

    let cand_adj = build_charge_candidate(&ctx_adj, 0, 2)
        .expect("adjacent should build a chain from the mono seed");
    let cand_seed = build_charge_candidate(&ctx_seed, 0, 2)
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
    assert!(
        anchor_config("hybrid").validate().is_err(),
        "removed `hybrid` cosine_anchor must now be a config error"
    );
    for v in ["adjacent", "seed", "SEED", "Seed"] {
        assert!(
            anchor_config(v).validate().is_ok(),
            "`{v}` should be an accepted cosine_anchor"
        );
    }
}

// Model-specific regression fixtures use the original fixed likelihood model;
// production defaults may change after the predefined tuning experiment.
fn evidence_config() -> FeaturesConfig {
    FeaturesConfig {
        isotope_evidence_ratio_sigma: 0.5,
        isotope_evidence_cosine_shape: 4.0,
        ..z2_config()
    }
}

/// Synthetic envelopes follow one fixed template so that each extension has
/// an independently known contribution, including contaminants at the tail.
fn evidence_envelope(residuals: &[f64]) -> Vec<Hill> {
    let config = evidence_config();
    let mz = 1_000.0;
    let template = config
        .isotope_model
        .model()
        .distribution(Polarity::Positive.neutral_mass(mz, 2));
    let shape = [0.1, 0.5, 1.0, 0.5, 0.1];
    residuals
        .iter()
        .enumerate()
        .map(|(k, &residual)| {
            let apex = 10_000.0 * template[k] / template[0] * residual.exp2();
            hill(
                mz + k as f64 * config.neutron_mass / 2.0,
                20,
                shape.iter().map(|v| (v * apex) as f32).collect(),
            )
        })
        .collect()
}

fn evidence_item(hills: &[Hill], len: usize) -> HeapItem {
    let sorted: Vec<_> = hills.iter().collect();
    let chain: Vec<_> = (0..len).collect();
    let mut config = evidence_config();
    config.sulfur_offsets.clear();
    let log_evidence = prefix_log_evidence(&sorted[..len], 2, &config, Polarity::Positive)[len - 1];
    HeapItem {
        chain,
        charge: 2,
        log_evidence,
    }
}

#[test]
fn evidence_rewards_good_isotopes_and_penalizes_bad_ones() {
    assert!(isotope_log_evidence(0.0, 1.0, &evidence_config()) > 0.0);
    assert!(isotope_log_evidence(1.5, 1.0, &evidence_config()) < 0.0);
    assert!(isotope_log_evidence(0.0, 0.2, &evidence_config()) < 0.0);
    assert!(isotope_log_evidence(1.0, 0.5, &evidence_config()) < 0.0);
    assert_eq!(
        isotope_log_evidence(0.0, 0.0, &evidence_config()),
        f64::NEG_INFINITY
    );
}

#[test]
fn clean_extensions_add_fixed_evidence_and_outscore_shorter_chain() {
    let hills = evidence_envelope(&[0.0; 6]);
    let short = evidence_item(&hills, 3);
    let long = evidence_item(&hills, 6);
    assert!(long > short);
    let one_isotope = isotope_log_evidence(0.0, 1.0, &evidence_config());
    assert!((short.log_evidence - 2.0 * one_isotope).abs() < 1e-6);
    assert!((long.log_evidence - 5.0 * one_isotope).abs() < 1e-6);
}

#[test]
fn six_hills_with_two_bad_isotopes_lose_to_clean_three() {
    let hills = evidence_envelope(&[0.0, 0.0, 0.0, 0.0, 1.5, 1.5]);
    let short = evidence_item(&hills, 3);
    let long = evidence_item(&hills, 6);
    assert!(short > long, "clean three must beat contaminated six");
    // Removing negative evidence raises the score: truncation must recompute,
    // not scale by length or assume a prefix always scores lower.
    assert!(evidence_item(&hills, 4) > long);
}

#[test]
fn rescored_prefix_reenters_heap_below_remaining_evidence() {
    let hills = evidence_envelope(&[0.0; 6]);
    let long = evidence_item(&hills, 6);
    let mut competitor = evidence_item(&hills, 4);
    competitor.chain = vec![10, 11, 12, 13];
    let mut heap = std::collections::BinaryHeap::from(vec![long, competitor]);
    assert_eq!(heap.pop().unwrap().chain.len(), 6);
    // Simulate the first contested isotope being M+3: the free prefix has
    // only two pieces of evidence, so the four-hill competitor must pop next.
    heap.push(evidence_item(&hills, 3));
    assert_eq!(heap.pop().unwrap().mono(), 10);
    assert_eq!(heap.pop().unwrap().chain.len(), 3);
}

#[test]
fn heap_ties_ignore_length_but_remain_deterministic() {
    let item = |chain, charge| HeapItem {
        chain,
        charge,
        log_evidence: 2.0,
    };
    let mut heap = std::collections::BinaryHeap::from(vec![
        item(vec![3, 4, 5, 6, 7, 8], 2),
        item(vec![0, 1, 2], 3),
        item(vec![0, 2], 2),
    ]);
    let first = heap.pop().unwrap();
    assert_eq!((first.mono(), first.charge), (0, 2));
    assert_eq!(heap.pop().unwrap().charge, 3);
    assert_eq!(heap.pop().unwrap().mono(), 3);
}

#[test]
fn evidence_handles_zero_and_nonfinite_intensities_without_nan() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut hills = evidence_envelope(&[0.0; 3]);
        hills[2].intensity_max = bad;
        assert_eq!(evidence_item(&hills, 3).log_evidence, f64::NEG_INFINITY);
        hills[0].intensity_max = bad;
        assert_eq!(evidence_item(&hills, 3).log_evidence, f64::NEG_INFINITY);
    }
}

#[test]
fn sulfur_evidence_uses_one_template_for_the_entire_chain() {
    let config = evidence_config();
    let model = config.isotope_model.model();
    let mass = Polarity::Positive.neutral_mass(1_000.0, 2);
    let (c, h, n, o, _) = model.counts(mass);
    let mut hills = evidence_envelope(&[0.0; 5]);
    // Deliberately give alternating isotopes different sulfur preferences.
    let t0 = crate::scoring::elements::cache().distribution(c, h, n, o, 0);
    let t2 = crate::scoring::elements::cache().distribution(c, h, n, o, 2);
    for (k, hill) in hills.iter_mut().enumerate().skip(1) {
        let t = if k % 2 == 0 { t2 } else { t0 };
        hill.intensity_max = 10_000.0 * t[k] / t[0];
    }
    let sorted: Vec<_> = hills.iter().collect();
    let actual = prefix_log_evidence(&sorted, 2, &config, Polarity::Positive)[4];
    let expected = averagine::resolve_sulfur_counts(&model, mass, &config.sulfur_offsets)
        .into_iter()
        .map(|s| {
            let t = crate::scoring::elements::cache().distribution(c, h, n, o, s);
            (1..5)
                .map(|k| {
                    let residual = (hills[k].intensity_max / hills[0].intensity_max).log2()
                        - (t[k] / t[0]).log2();
                    isotope_log_evidence(residual, 1.0, &evidence_config())
                })
                .sum::<f64>()
        })
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((actual - expected).abs() < 1e-6);
}

#[test]
fn best_prefix_discards_negative_tail_before_competing() {
    let hills = evidence_envelope(&[0.0, 0.0, 0.0, 0.0, 1.5, 1.5]);
    let sorted: Vec<_> = hills.iter().collect();
    let mut config = evidence_config();
    config.sulfur_offsets.clear();
    let selected =
        rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).unwrap();
    assert_eq!(selected.chain, vec![0, 1, 2, 3]);
    // Before selection, the six-hill chain loses to a clean three. Selecting
    // its good four-hill prefix gives it the correct priority immediately.
    assert!(selected > evidence_item(&hills, 3));
    assert!(evidence_item(&hills, 3) > evidence_item(&hills, 6));
}

#[test]
fn every_conflict_rescore_is_bounded_by_the_queued_best_prefix() {
    for residuals in [
        [0.0, 0.0, 1.5, 0.0, 0.0, 1.5],
        [0.0, 0.0, 0.0, 0.0, 1.5, 1.5],
        [0.0, 1.5, 0.0, 1.5, 0.0, 0.0],
    ] {
        let hills = evidence_envelope(&residuals);
        let sorted: Vec<_> = hills.iter().collect();
        // Includes max-over-sulfur switching: the invariant holds for all
        // templates because conflicts only shrink the feasible prefix set.
        let config = evidence_config();
        let queued =
            rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).unwrap();
        for n in 2..queued.chain.len() {
            let prefix =
                rescore_chain(&queued.chain[..n], 2, &sorted, &config, Polarity::Positive).unwrap();
            assert!(prefix.log_evidence <= queued.log_evidence);
        }
    }
}

#[test]
fn prefix_search_does_not_stop_at_the_first_negative_contribution() {
    let hills = evidence_envelope(&[0.0, 0.0, 1.5, 0.0, 0.0]);
    let sorted: Vec<_> = hills.iter().collect();
    let mut config = evidence_config();
    config.sulfur_offsets.clear();
    let selected =
        rescore_chain(&[0, 1, 2, 3, 4], 2, &sorted, &config, Polarity::Positive).unwrap();
    assert_eq!(
        selected.chain.len(),
        5,
        "later good evidence can outweigh an interior penalty"
    );
}

#[test]
fn best_prefix_respects_claim_gate_and_checks_alternatives() {
    let hills = evidence_envelope(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let sorted: Vec<_> = hills.iter().collect();
    let mut config = evidence_config();
    config.sulfur_offsets.clear();
    config.exhaustive_min_isotope_score = 1.01;
    assert!(rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).is_none());
    // The top-evidence two-hill prefix can fail the missing-isotope penalty
    // while a longer prefix still passes. Search all valid prefixes.
    let hills = evidence_envelope(&[0.0, 0.0, 1.5, 1.5, 1.5, 1.5]);
    let sorted: Vec<_> = hills.iter().collect();
    config.exhaustive_min_isotope_score = 0.0;
    let top = rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).unwrap();
    let top_bc = score_chain(
        &sorted[..top.chain.len()],
        2,
        &[],
        &config.isotope_model.model(),
        Polarity::Positive,
    );
    let full_bc = score_chain(
        &sorted,
        2,
        &[],
        &config.isotope_model.model(),
        Polarity::Positive,
    );
    assert!(
        full_bc > top_bc,
        "fixture must offer a gate-passing alternative"
    );
    config.exhaustive_min_isotope_score = (top_bc + full_bc) / 2.0;
    let alternative =
        rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).unwrap();
    assert!(alternative.chain.len() > top.chain.len());
    assert!(alternative.log_evidence < top.log_evidence);
}

#[test]
fn evidence_config_rejects_invalid_distribution_parameters() {
    for sigma in [0.0, -1.0, 2.0, f64::INFINITY, f64::NAN] {
        let mut config = evidence_config();
        config.isotope_evidence_ratio_sigma = sigma;
        assert!(config.validate().is_err());
    }
    for shape in [0.0, 1.0, f64::INFINITY, f64::NAN] {
        let mut config = evidence_config();
        config.isotope_evidence_cosine_shape = shape;
        assert!(config.validate().is_err());
    }
}

#[test]
fn resolver_claims_good_prefix_without_absorbing_its_bad_tail() {
    let mut hills = evidence_envelope(&[0.0, 0.0, 0.0, 0.0, 1.5, 1.5]);
    for (i, h) in hills.iter_mut().enumerate() {
        h.hill_id = i as u64;
    }
    let sorted: Vec<_> = hills.iter().collect();
    let mz: Vec<_> = hills.iter().map(|h| h.mz).collect();
    let im: Vec<_> = hills.iter().map(|h| h.im).collect();
    let starts: Vec<_> = hills.iter().map(|h| h.scan_start).collect();
    let ends: Vec<_> = hills.iter().map(|h| h.scan_end).collect();
    let mut config = evidence_config();
    config.sulfur_offsets.clear();
    config.max_isotope_log2_ratio = 2.0;
    let file = FileConfig::default();
    let ctx = ChainCtx {
        sorted_hills: &sorted,
        mz_array: &mz,
        im_array: &im,
        scan_starts: &starts,
        scan_ends: &ends,
        config: &config,
        file: &file,
        use_im: false,
        min_intensity: 0.0,
        recal: None,
    };
    assert_eq!(
        build_charge_candidate(&ctx, 0, 2)
            .unwrap()
            .hill_indices
            .len(),
        6
    );
    let resolved = resolve_exhaustive(&ctx);
    let mono = resolved.iter().find(|c| c.hill_indices[0] == 0).unwrap();
    assert_eq!(mono.hill_indices, vec![0, 1, 2, 3]);
    assert!(resolved
        .iter()
        .filter(|c| c.hill_indices[0] != 0)
        .all(|c| c.hill_indices.iter().all(|i| *i >= 4)));
}

#[test]
fn tuned_likelihood_penalizes_clear_mismatches_and_selects_clean_prefix() {
    // Frozen training winner; this regression also documents the less severe
    // penalty relative to the original 0.5/Beta(4,1) pilot model.
    let config = FeaturesConfig {
        sulfur_offsets: vec![],
        ..z2_config()
    };
    assert!(isotope_log_evidence(0.0, 1.0, &config) > 0.0);
    assert!(isotope_log_evidence(2.0, 1.0, &config) < 0.0);
    assert!(isotope_log_evidence(0.0, 0.1, &config) < 0.0);
    let hills = evidence_envelope(&[0.0, 0.0, 0.0, 0.0, 2.0, 2.0]);
    let sorted: Vec<_> = hills.iter().collect();
    let scores = prefix_log_evidence(&sorted, 2, &config, Polarity::Positive);
    assert!(
        scores[2] > scores[5],
        "a clean three beats a six with two clear contaminants"
    );
    let selected =
        rescore_chain(&[0, 1, 2, 3, 4, 5], 2, &sorted, &config, Polarity::Positive).unwrap();
    assert_eq!(selected.chain, vec![0, 1, 2, 3]);
}
