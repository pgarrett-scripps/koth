    use super::*;
    use crate::lfq::LfqEntry;

    /// Minimal `LfqEntry` for TDC tests — only the four fields `compute_qvalues`
    /// reads (feature_idx, run_idx, hybrid_score, is_decoy) are meaningful.
    fn entry(feature_idx: usize, run_idx: usize, score: f32, is_decoy: bool) -> LfqEntry {
        LfqEntry {
            feature_idx,
            run_idx,
            intensity: 1.0,
            hybrid_score: score,
            spectral_bhattacharyya: 0.0,
            n_isotopes_found: 0,
            rt_score: 0.0,
            int_score: 0.0,
            coelution: 0.0,
            is_decoy,
            is_mbr: false,
            expected_rt: 0.0,
            apex_rt: 0.0,
            peak_width_rt: 0.0,
            expected_mz: 0.0,
            observed_mz: 0.0,
            expected_im: 0.0,
            observed_im: 0.0,
        }
    }

    /// With no decoys the running FDR is 0 at every target, so every q-value is 0.
    #[test]
    fn all_targets_get_zero_q() {
        let entries = vec![
            entry(0, 0, 0.9, false),
            entry(1, 0, 0.8, false),
            entry(2, 0, 0.7, false),
        ];
        let q = compute_qvalues(&entries);
        assert_eq!(q.len(), 3);
        for v in q.values() {
            assert_eq!(*v, 0.0);
        }
    }

    /// Hand-computed interleaving: scores desc T(1.0) D(0.9) T(0.8) T(0.7).
    /// Forward FDR: T@1.0 → 0/1=0, T@0.8 → 1/2=0.5, T@0.7 → 1/3. Backward
    /// monotonisation (min at-or-below rank): 0 / 1/3 / 1/3.
    #[test]
    fn running_fdr_counts_decoys_over_targets_and_monotonises() {
        let entries = vec![
            entry(0, 0, 1.0, false),
            entry(9, 0, 0.9, true),
            entry(1, 0, 0.8, false),
            entry(2, 0, 0.7, false),
        ];
        let q = compute_qvalues(&entries);
        assert!((q[&(0, 0)] - 0.0).abs() < 1e-9);
        assert!((q[&(1, 0)] - 1.0 / 3.0).abs() < 1e-9, "q was {}", q[&(1, 0)]);
        assert!((q[&(2, 0)] - 1.0 / 3.0).abs() < 1e-9, "q was {}", q[&(2, 0)]);
        // decoy entries are not reported
        assert_eq!(q.len(), 3);
    }

    /// Q-values must be monotone non-decreasing as score falls (the backward
    /// pass guarantees a higher-scoring target never has a larger q than a
    /// lower-scoring one).
    #[test]
    fn qvalues_monotone_in_score() {
        let entries = vec![
            entry(0, 0, 0.95, false),
            entry(1, 0, 0.90, true),
            entry(2, 0, 0.85, false),
            entry(3, 0, 0.80, false),
            entry(4, 0, 0.75, true),
            entry(5, 0, 0.70, false),
        ];
        let q = compute_qvalues(&entries);
        // targets in descending score order
        let qs: Vec<f64> = [(0, 0), (2, 0), (3, 0), (5, 0)].iter().map(|k| q[k]).collect();
        for w in qs.windows(2) {
            assert!(w[1] >= w[0] - 1e-12, "q not monotone non-decreasing: {qs:?}");
        }
    }
