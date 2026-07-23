    use super::*;
    use crate::config::{FileConfig, HillsConfig};
    use crate::models::{Peak, Spectrum};

    fn spec(scan_index: usize, rt: f64, mzs: &[f32]) -> Spectrum {
        Spectrum {
            scan_index,
            retention_time: rt,
            peaks: mzs
                .iter()
                .map(|&mz| Peak {
                    mz,
                    intensity: 1000.0,
                    ion_mobility: 0.0,
                })
                .collect(),
            ms_level: 1,
            isolation_window: None,
        }
    }

    fn default_cfgs(max_gap: usize) -> (HillsConfig, FileConfig) {
        let mut h = HillsConfig::default();
        h.max_gap = max_gap;
        h.split_hills = false;
        h.min_scans = 3;
        (h, FileConfig::default())
    }

    /// Regression: with `max_gap = 0`, a hill that matches the same m/z in N
    /// consecutive scans must be retained (previously a stale-cleanup off-by-one
    /// killed every hill before it could extend past 1 scan).
    #[test]
    fn max_gap_zero_retains_consecutive_hill() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..5 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 1, "expected 1 hill, got {}", hills.len());
        assert_eq!(hills[0].n_scans, 5);
        assert_eq!(hills[0].skipped_scans, 0);
    }

    /// With `max_gap = 0`, a one-scan gap must terminate the hill.
    #[test]
    fn max_gap_zero_disallows_gaps() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..3 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(3, 0.3, &[]));
        for i in 4..7 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 2, "expected 2 hills split across the gap, got {}", hills.len());
        for hl in &hills {
            assert_eq!(hl.skipped_scans, 0);
        }
    }

    /// With `max_gap = 1`, a single missed scan must be bridged.
    #[test]
    fn max_gap_one_bridges_single_miss() {
        let (h, f) = default_cfgs(1);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..2 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(2, 0.2, &[])); // miss
        for i in 3..5 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 1, "expected 1 bridged hill, got {}", hills.len());
        assert_eq!(hills[0].skipped_scans, 1);
    }

    /// Calibration recording: when `with_calibration_recording` is set,
    /// every accepted peak-to-hill match contributes one |ppm| delta to
    /// the histogram. A clean signal (same m/z every scan) lands at 0 ppm
    /// so the calibrated tolerance approaches 0 — the calibration loop
    /// has nothing to widen against.
    #[test]
    fn calibration_records_zero_ppm_for_constant_signal() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f).with_calibration_recording();
        for i in 0..10 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        // 9 accepted matches (peak in scans 1..9 each match the hill
        // started in scan 0). All are at exactly 500.0 m/z so |ppm|=0.
        assert_eq!(det.calibration_sample_count(), 9);
        let cal = det.calibrated_tolerance_ppm(3.0).unwrap();
        assert!(cal < 0.1, "expected calibrated ~0, got {cal}");
    }

    /// Calibration with a drift: m/z walks up by ~1 ppm per scan. The
    /// histogram median + 3σ should land somewhere comfortably above 0
    /// (the drift magnitude × sigma_mult).
    #[test]
    fn calibration_widens_for_drifting_signal() {
        let (h, f) = default_cfgs(0);
        let mut det = HillDetector::new(&h, &f).with_calibration_recording();
        // 1 ppm/scan walk on a 500 Da peak = +0.0005 Da/scan.
        for i in 0..20 {
            let mz = 500.0 + (i as f32) * 5e-4;
            det.process_scan(&spec(i, i as f64 * 0.1, &[mz]));
        }
        assert!(det.calibration_sample_count() > 0);
        let cal = det.calibrated_tolerance_ppm(3.0).unwrap();
        assert!(
            cal > 0.5,
            "expected calibrated > 0.5 ppm for drifting signal, got {cal}"
        );
    }

    /// With `max_gap = 1`, two consecutive missed scans must split the hill.
    #[test]
    fn max_gap_one_splits_at_two_misses() {
        let (h, f) = default_cfgs(1);
        let mut det = HillDetector::new(&h, &f);
        for i in 0..3 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        det.process_scan(&spec(3, 0.3, &[])); // miss
        det.process_scan(&spec(4, 0.4, &[])); // miss
        for i in 5..8 {
            det.process_scan(&spec(i, i as f64 * 0.1, &[500.0]));
        }
        let hills = det.finish();
        assert_eq!(hills.len(), 2, "expected 2 hills, got {}", hills.len());
    }
