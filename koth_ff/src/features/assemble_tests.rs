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
        }
    }

    fn z2_config() -> FeaturesConfig {
        FeaturesConfig {
            min_charge: 2,
            max_charge: 2,
            ..Default::default()
        }
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
        assert_eq!(key(&a), key(&c), "output must not depend on input hill order");
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

