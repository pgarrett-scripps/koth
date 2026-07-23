    use super::*;

    #[test]
    fn too_few_samples_returns_none() {
        let mut b = MzRecalBuilder::new(50.0);
        for _ in 0..5 {
            b.add(500.0, 10.0, 1.0);
        }
        assert!(b.finalize(4, 4, 50).is_none());
    }

    #[test]
    fn recovers_constant_offset() {
        let mut b = MzRecalBuilder::new(50.0);
        for i in 0..1000 {
            let mz = 300.0 + (i % 500) as f64;
            b.add(mz, (i % 30) as f64, 2.5);
        }
        let m = b.finalize(10, 5, 20).unwrap();
        assert!((m.predict(500.0, 15.0) - 2.5).abs() < 0.05);
        // Out-of-range clamps to edge bins, still ~2.5.
        assert!((m.predict(50.0, 100.0) - 2.5).abs() < 0.05);
    }

    #[test]
    fn recovers_mz_dependent_offset() {
        // Low m/z biased +1 ppm, high m/z biased -3 ppm.
        let mut b = MzRecalBuilder::new(50.0);
        for i in 0..2000 {
            let rt = (i % 20) as f64;
            b.add(400.0, rt, 1.0);
            b.add(1200.0, rt, -3.0);
        }
        let m = b.finalize(8, 4, 20).unwrap();
        assert!((m.predict(400.0, 5.0) - 1.0).abs() < 0.1, "low-mz {}", m.predict(400.0, 5.0));
        assert!((m.predict(1200.0, 5.0) + 3.0).abs() < 0.1, "high-mz {}", m.predict(1200.0, 5.0));
    }

    #[test]
    fn sparse_cell_falls_back_to_marginal_then_global() {
        let mut b = MzRecalBuilder::new(50.0);
        // Dense at rt≈0 with offset 4.0; a single sparse sample at rt≈100.
        for _ in 0..500 {
            b.add(600.0, 0.0, 4.0);
        }
        b.add(600.0, 100.0, -10.0);
        let m = b.finalize(4, 4, 50).unwrap();
        // The lone rt≈100 sample can't meet min_samples for its own cell, so
        // it falls back to the m/z-marginal median (dominated by the 4.0s).
        let p = m.predict(600.0, 100.0);
        assert!((p - 4.0).abs() < 0.5, "expected marginal fallback ~4.0, got {p}");
    }

    #[test]
    fn out_of_range_ppm_rejected() {
        let mut b = MzRecalBuilder::new(10.0);
        b.add(500.0, 1.0, 5.0);
        b.add(500.0, 1.0, 500.0); // rejected
        b.add(500.0, 1.0, f64::NAN); // rejected
        assert_eq!(b.len(), 1);
    }

    #[test]
    fn deterministic_across_builds() {
        let build = || {
            let mut b = MzRecalBuilder::new(50.0);
            for i in 0..777 {
                let mz = 300.0 + ((i * 37) % 900) as f64;
                let rt = ((i * 13) % 40) as f64;
                let resid = ((i % 7) as f64) - 3.0;
                b.add(mz, rt, resid);
            }
            b.finalize(9, 6, 15).unwrap()
        };
        let a = build();
        let c = build();
        for &(mz, rt) in &[(350.0, 5.0), (700.0, 20.0), (1100.0, 35.0)] {
            assert_eq!(a.predict(mz, rt), c.predict(mz, rt));
        }
    }
