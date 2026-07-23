    use super::*;

    #[test]
    fn empty_histogram_returns_none() {
        let h = MzDeltaHistogram::new(50.0, 0.01);
        assert_eq!(h.count_in_range(), 0);
        assert!(h.median().is_none());
        assert!(h.std_dev().is_none());
        assert!(h.calibrated_ppm(3.0).is_none());
    }

    #[test]
    fn single_sample_has_median_but_no_std() {
        let mut h = MzDeltaHistogram::new(50.0, 0.01);
        h.record(2.5);
        assert_eq!(h.count_in_range(), 1);
        assert!(h.median().is_some());
        assert!(h.std_dev().is_none());
    }

    #[test]
    fn median_recovers_constant() {
        let mut h = MzDeltaHistogram::new(50.0, 0.01);
        for _ in 0..1000 {
            h.record(2.34);
        }
        let m = h.median().unwrap();
        assert!((m - 2.34).abs() < 0.02, "median {m} not near 2.34");
    }

    #[test]
    fn std_zero_for_constant() {
        let mut h = MzDeltaHistogram::new(50.0, 0.01);
        for _ in 0..1000 {
            h.record(1.0);
        }
        let s = h.std_dev().unwrap();
        assert!(s < 0.02, "std {s} not near 0");
    }

    #[test]
    fn calibrated_value_combines_median_and_sigma() {
        // Skewed mix: 800 at 1.0 and 200 at 3.0. Lower-median falls at
        // 1.0, mean = 1.4, variance = 0.8·(1−1.4)² + 0.2·(3−1.4)² = 0.64,
        // σ = 0.8. So median + 3σ ≈ 1.0 + 2.4 = 3.4.
        let mut h = MzDeltaHistogram::new(50.0, 0.01);
        for _ in 0..800 {
            h.record(1.0);
        }
        for _ in 0..200 {
            h.record(3.0);
        }
        let cal = h.calibrated_ppm(3.0).unwrap();
        assert!(
            (cal - 3.4).abs() < 0.2,
            "expected calibrated ~3.4, got {cal}"
        );
    }

    #[test]
    fn negative_and_nan_dropped() {
        let mut h = MzDeltaHistogram::new(50.0, 0.01);
        h.record(-1.0);
        h.record(f64::NAN);
        h.record(f64::INFINITY);
        assert_eq!(h.count_in_range(), 0);
    }

    #[test]
    fn over_range_counted_separately() {
        let mut h = MzDeltaHistogram::new(10.0, 0.1);
        h.record(5.0);
        h.record(15.0);
        assert_eq!(h.count_in_range(), 1);
        assert_eq!(h.count_total(), 2);
    }
