    use super::*;
    use crate::models::Hill;
    use std::sync::Arc;

    fn make_hill(profile: Vec<f32>) -> Hill {
        let n = profile.len();
        let int_sum: f64 = profile.iter().map(|&x| x as f64).sum();
        let int_max: f64 = profile.iter().copied().fold(0.0f32, f32::max) as f64;
        Hill {
            hill_id: 0,
            mz: 500.0,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: 0.5,
            rt_start: 0.0,
            rt_end: 1.0,
            rt_width: 1.0,
            im: 0.0,
            im_std: 0.0,
            scan_start: 0,
            scan_apex: n / 2,
            scan_end: n.saturating_sub(1),
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: int_sum,
            intensity_max: int_max,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
        }
    }

    /// Small hills (below threshold) are always kept regardless of shape.
    #[test]
    fn small_hill_always_kept() {
        let flat = make_hill(vec![1000.0; 20]); // 20 < 40
        assert!(keep_hill(&flat, 40, 2.0));
    }

    /// A large hill with a clear chromatographic apex (max in the middle,
    /// low intensity at endpoints) is kept.
    #[test]
    fn large_peaked_hill_kept() {
        // 50-scan Gaussian-ish profile peaking at center 25.
        let n = 50;
        let profile: Vec<f32> = (0..n)
            .map(|i| {
                let d = (i as f32 - 25.0) / 8.0;
                10_000.0 * (-d * d).exp() + 50.0 // small baseline at edges
            })
            .collect();
        let h = make_hill(profile);
        assert!(keep_hill(&h, 40, 2.0));
    }

    /// A large flat hill (constant intensity across all scans) — max ≈
    /// endpoints — is dropped. This is the "column-bleed contaminant"
    /// case the filter is designed to catch.
    #[test]
    fn large_flat_hill_dropped() {
        let flat = make_hill(vec![1000.0; 50]);
        assert!(!keep_hill(&flat, 40, 2.0));
    }

    /// A large hill with a monotonic ramp shape (max at one end, low at
    /// the other) is dropped — the apex-ratio test requires both
    /// endpoints to be low.
    #[test]
    fn large_ramp_hill_dropped() {
        let n = 50;
        let profile: Vec<f32> = (0..n).map(|i| 100.0 + i as f32 * 200.0).collect();
        let h = make_hill(profile);
        assert!(!keep_hill(&h, 40, 2.0));
    }

    /// Boundary check: a hill with exactly `min_scans` scans is evaluated
    /// (not skipped as "small").
    #[test]
    fn boundary_at_min_scans_is_evaluated() {
        let flat = make_hill(vec![1000.0; 40]);
        assert!(!keep_hill(&flat, 40, 2.0));
    }

    /// Degenerate case: smoothed endpoint is zero (e.g., heavy gap fill).
    /// The filter conservatively keeps such hills rather than dividing by
    /// zero.
    #[test]
    fn zero_endpoint_hill_kept_conservatively() {
        // Heavy zeros at both ends → smoothed endpoints will be zero.
        let mut profile = vec![0.0f32; 50];
        for i in 20..30 {
            profile[i] = 1000.0;
        }
        let h = make_hill(profile);
        assert!(keep_hill(&h, 40, 2.0));
    }

    /// End-to-end: a mixed input of 3 peaked + 2 flat large hills + 1
    /// small hill is filtered down to 4 (3 peaked + 1 small).
    #[test]
    fn filter_drops_flat_keeps_peaked_and_small() {
        let peaked = |n: usize| {
            (0..n)
                .map(|i| {
                    let d = (i as f32 - n as f32 / 2.0) / (n as f32 / 6.0);
                    10_000.0 * (-d * d).exp() + 50.0
                })
                .collect::<Vec<f32>>()
        };
        let hills = vec![
            make_hill(peaked(45)),
            make_hill(peaked(50)),
            make_hill(peaked(60)),
            make_hill(vec![1000.0; 45]),
            make_hill(vec![1000.0; 50]),
            make_hill(vec![1000.0; 20]), // small — pass-through
        ];
        let filtered = filter_large_baseline_hills(hills, 40, 2.0);
        assert_eq!(filtered.len(), 4);
    }
