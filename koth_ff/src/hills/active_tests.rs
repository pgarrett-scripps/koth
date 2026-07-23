    use super::*;
    use crate::hills::smooth;

    fn make_hill_with_gaps() -> ActiveHill {
        // 6 scans, gaps at indices 1, 2, 4. Real m/z values cluster within ~0.005 Da.
        let mut h = ActiveHill::new(524.270, 100.0, 0.0, 0.9, 0, true);
        h.add_gap();
        h.add_gap();
        h.mz_profile.push(524.275);
        h.intensity_profile.push(150.0);
        h.rt_profile.push(0.3);
        h.im_profile.push(0.91);
        h.n_real += 1;
        h.last_scan_seen = 3;
        h.add_gap();
        h.mz_profile.push(524.272);
        h.intensity_profile.push(120.0);
        h.rt_profile.push(0.5);
        h.im_profile.push(0.905);
        h.n_real += 1;
        h.last_scan_seen = 5;
        h
    }

    /// Kish SE for a single peak is 0 (no variance to estimate from).
    #[test]
    fn mz_kish_se_zero_for_singleton() {
        let h = ActiveHill::new(500.0, 1000.0, 0.0, 0.0, 0, true);
        assert_eq!(h.mz_kish_se(), 0.0);
    }

    /// Kish SE for uniform intensities equals `weighted_std / √n`. A 5-peak
    /// hill at constant intensity with m/z spread σ should give SE ≈ σ / √5.
    #[test]
    fn mz_kish_se_uniform_intensity_matches_classic_se() {
        // 5 peaks at constant intensity 1000.0, m/z [499.99, 499.995, 500.0, 500.005, 500.01]
        // weighted_std = sqrt(mean(d²)) where d = mz - 500.0 = [-0.01, -0.005, 0, 0.005, 0.01]
        // = sqrt((1e-4 + 0.25e-4 + 0 + 0.25e-4 + 1e-4) / 5) = sqrt(2.5e-4 / 5) ≈ 0.00707
        // neff = 5, SE = weighted_std / √5 ≈ 0.00316
        let mut h = ActiveHill::new(499.99, 1000.0, 0.0, 0.0, 0, true);
        for (i, mz) in [499.995f32, 500.0, 500.005, 500.01].iter().enumerate() {
            h.mz_profile.push(*mz);
            h.intensity_profile.push(1000.0);
            h.rt_profile.push((i + 1) as f32 * 0.1);
            h.im_profile.push(0.0);
            h.n_real += 1;
            h.last_scan_seen = i + 1;
        }
        let std = h.mz_weighted_std();
        let se = h.mz_kish_se();
        let expected_se = std / (5.0_f64).sqrt();
        assert!(
            (se - expected_se).abs() < 1e-6,
            "SE {se} not within 1e-6 of weighted_std/√n = {expected_se}",
        );
    }

    /// With skewed intensities (one peak dominates), `neff < n` and the SE is
    /// larger than the unweighted `std/√n` estimate. Captures the AlphaPept
    /// rationale: a hill whose weighted-mean is dominated by one peak has a
    /// less-confident m/z estimate than the simple count would suggest.
    #[test]
    fn mz_kish_se_skewed_intensity_inflates_se() {
        // 5 peaks, but peak 0 has 100× the intensity of the others.
        let mut h = ActiveHill::new(500.000, 100_000.0, 0.0, 0.0, 0, true);
        for (i, mz) in [500.005f32, 500.010, 500.015, 500.020].iter().enumerate() {
            h.mz_profile.push(*mz);
            h.intensity_profile.push(1000.0);
            h.rt_profile.push((i + 1) as f32 * 0.1);
            h.im_profile.push(0.0);
            h.n_real += 1;
            h.last_scan_seen = i + 1;
        }
        // neff = (100k + 4k)² / (100k² + 4·1k²) ≈ 104k² / 1.0004e10 ≈ 1.08
        // i.e. effective sample size barely > 1, so SE ≈ weighted_std
        let std = h.mz_weighted_std();
        let se = h.mz_kish_se();
        assert!(
            se > std * 0.5,
            "skewed-intensity SE {se} should be close to weighted_std {std} (neff ≈ 1)",
        );
        assert!(
            se >= std / (5.0_f64).sqrt(),
            "skewed-intensity SE {se} should exceed the unweighted √n estimate {}",
            std / (5.0_f64).sqrt(),
        );
    }

    /// Gaps (NaN m/z) are excluded from the Kish calculation, mirroring the
    /// behaviour of `mz_weighted_std`.
    #[test]
    fn mz_kish_se_ignores_gaps() {
        let h = make_hill_with_gaps();
        // 3 real peaks at intensities [100, 150, 120], so neff < 3 but > 1.
        let se = h.mz_kish_se();
        assert!(se > 0.0 && se.is_finite(), "expected non-zero SE, got {se}");
        // Should be much smaller than the m/z range (~0.005 Da).
        assert!(se < 0.005, "SE {se} seems implausibly large");
    }

    /// Without smoothing the std is fine — guarantee the new implementation matches.
    #[test]
    fn mz_std_without_smoothing_is_small() {
        let h = make_hill_with_gaps();
        let mean = h.mz_weighted_mean();
        let std = h.mz_weighted_std();

        assert!((mean - 524.273).abs() < 1e-3, "mean was {mean}");
        assert!(std < 0.01, "expected tight m/z std (<0.01 Da), got {std}");
    }

    /// Regression: with smoothing, `fill_gaps` makes gap-position intensities
    /// non-zero, but `mz_profile[gap] = NaN` is untouched. The old code summed
    /// ALL intensities in the denominator and filtered NaN m/z only in the
    /// numerator — producing a wildly biased mean (toward 0) and an exploding
    /// std. The fix filters NaN positions in both numerator and denominator.
    #[test]
    fn mz_std_unaffected_by_smoothing_gap_fill() {
        let mut h = make_hill_with_gaps();
        smooth::apply_intensity_filters(&mut h.intensity_profile, true, true, 1);

        let mean = h.mz_weighted_mean();
        let std = h.mz_weighted_std();

        // Mean must remain close to the real m/z cluster, NOT collapse to 0.
        assert!(
            (mean - 524.273).abs() < 0.01,
            "smoothing should not shift mean far from m/z cluster, got {mean}",
        );
        // Std must stay below the Da-scale range across the real points.
        assert!(
            std < 0.01,
            "smoothing should not inflate mz_std (gaps must not contribute); got {std}",
        );
    }

    /// Regression: after smoothing fills gap-scan intensities, the apex can
    /// land on a slot whose `rt_profile` entry is NaN. `apex_rt()` used to
    /// return 0.0 in that case (so hills.tsv showed `rt = 0.0` for hills with
    /// non-zero `rt_start`/`rt_end`). The fix interpolates from observed bounds.
    #[test]
    fn apex_rt_interpolates_when_apex_lands_on_gap() {
        let mut h = make_hill_with_gaps();
        // Boost the gap at index 1 so that after smoothing it becomes the apex.
        // rt_profile[1] is NaN — pre-fix, apex_rt() returned 0.0.
        h.intensity_profile[1] = 0.0; // gap, untouched
        h.intensity_profile[3] = 10_000.0; // anchor neighbor so smoothing peaks near idx 2-3
        smooth::apply_intensity_filters(&mut h.intensity_profile, true, true, 2);

        let rt = h.apex_rt();
        assert!(
            rt > 0.0 && rt.is_finite(),
            "apex_rt should not collapse to 0 when apex lands on a NaN slot; got {rt}",
        );
        let lo = h.min_rt();
        let hi = h.max_rt();
        assert!(
            rt >= lo && rt <= hi,
            "interpolated apex_rt {rt} must lie within [{lo}, {hi}]",
        );
    }

    /// Same regression for ion mobility.
    #[test]
    fn im_std_unaffected_by_smoothing_gap_fill() {
        let mut h = make_hill_with_gaps();
        smooth::apply_intensity_filters(&mut h.intensity_profile, true, true, 1);

        let mean = h.im_weighted_mean();
        let std = h.im_weighted_std();

        assert!((mean - 0.905).abs() < 0.01, "im mean drifted: {mean}");
        assert!(std < 0.01, "im_std inflated by gap-fill: {std}");
    }
