use super::*;
use crate::alignment::anchors::AnchorPair;

/// Anchor whose `ppm_error()` equals `ppm` at normalised RT `x`.
/// `ref_mz` is fixed at 1000 so `run_mz = 1000 * (1 + ppm/1e6)` reproduces `ppm`
/// exactly through `AnchorPair::ppm_error`.
fn anchor_ppm(x: f64, ppm: f64) -> AnchorPair {
    let ref_mz = 1000.0;
    AnchorPair {
        ref_rt_norm: x,
        run_rt_norm: x,
        ref_mz,
        run_mz: ref_mz * (1.0 + ppm / 1e6),
        ref_im: 0.0,
        run_im: 0.0,
    }
}

/// Anchor with IM present on both sides whose `im_delta()` equals `delta` at `x`.
fn anchor_im(x: f64, delta: f64) -> AnchorPair {
    AnchorPair {
        ref_rt_norm: x,
        run_rt_norm: x,
        ref_mz: 1000.0,
        run_mz: 1000.0,
        ref_im: 1.0,
        run_im: 1.0 + delta,
    }
}

#[test]
fn zero_fit_predicts_zero_everywhere() {
    let f = DriftFit::zero();
    assert_eq!(f.predict(0.0), 0.0);
    assert_eq!(f.predict(1.0), 0.0);
    assert_eq!(f.predict(-3.7), 0.0);
}

#[test]
fn predict_is_affine() {
    let f = DriftFit {
        intercept: 2.0,
        slope: 3.0,
    };
    assert!((f.predict(4.0) - 14.0).abs() < 1e-12);
    assert!((f.predict(0.0) - 2.0).abs() < 1e-12);
}

#[test]
fn ols_recovers_exact_line() {
    // y = 1.5 + 4.0 x, sampled without noise.
    let pairs: Vec<(f64, f64)> = (0..10)
        .map(|i| i as f64 / 10.0)
        .map(|x| (x, 1.5 + 4.0 * x))
        .collect();
    let fit = ols(&pairs);
    assert!(
        (fit.intercept - 1.5).abs() < 1e-9,
        "intercept {}",
        fit.intercept
    );
    assert!((fit.slope - 4.0).abs() < 1e-9, "slope {}", fit.slope);
}

#[test]
fn std_dev_matches_sample_formula() {
    // [2,4,6,8]: mean 5, sample variance (9+1+1+9)/3 = 20/3, sd = sqrt(6.6667).
    let sd = std_dev(&[2.0, 4.0, 6.0, 8.0]);
    assert!((sd - (20.0f64 / 3.0).sqrt()).abs() < 1e-9, "sd {sd}");
    assert_eq!(std_dev(&[42.0]), 0.0, "n<2 -> 0");
}

#[test]
fn mass_drift_recovers_clean_linear_trend() {
    // ppm = 2.0 + 5.0 * rt_norm over the unit interval.
    let anchors: Vec<AnchorPair> = (0..11)
        .map(|i| i as f64 / 10.0)
        .map(|x| anchor_ppm(x, 2.0 + 5.0 * x))
        .collect();
    let (fit, active) = fit_mass_drift(&anchors);
    assert!(
        (fit.intercept - 2.0).abs() < 1e-6,
        "intercept {}",
        fit.intercept
    );
    assert!((fit.slope - 5.0).abs() < 1e-6, "slope {}", fit.slope);
    assert!(
        active.iter().all(|&ok| ok),
        "clean anchors should all stay active"
    );
}

#[test]
fn mass_drift_sigma_clips_gross_outlier() {
    // 20 clean points on ppm = 5 x, plus one 500-ppm flyer at index 20.
    let mut anchors: Vec<AnchorPair> = (0..20)
        .map(|i| i as f64 / 20.0)
        .map(|x| anchor_ppm(x, 5.0 * x))
        .collect();
    anchors.push(anchor_ppm(0.5, 5.0 * 0.5 + 500.0)); // outlier
    let (fit, active) = fit_mass_drift(&anchors);

    assert!(!active[20], "the 500-ppm outlier must be clipped");
    assert!(
        active[..20].iter().all(|&ok| ok),
        "clean anchors must survive"
    );
    // With the outlier removed the fit collapses back onto the true line.
    assert!(
        (fit.intercept - 0.0).abs() < 1e-6,
        "intercept {}",
        fit.intercept
    );
    assert!((fit.slope - 5.0).abs() < 1e-6, "slope {}", fit.slope);
}

#[test]
fn mass_drift_too_few_pairs_is_zero() {
    let (fit, active) = fit_mass_drift(&[anchor_ppm(0.5, 3.0)]);
    assert_eq!(fit.slope, 0.0);
    assert_eq!(fit.intercept, 0.0);
    assert_eq!(active, vec![false]);
}

#[test]
fn im_drift_requires_four_pairs() {
    // Only 3 IM-bearing anchors -> zero fit, all-false mask.
    let anchors: Vec<AnchorPair> = (0..3).map(|i| anchor_im(i as f64 / 3.0, 0.1)).collect();
    let (fit, active) = fit_im_drift(&anchors);
    assert_eq!(fit.slope, 0.0);
    assert_eq!(fit.intercept, 0.0);
    assert!(active.iter().all(|&ok| !ok));
}

#[test]
fn im_drift_fits_and_masks_non_im_anchors() {
    // 5 IM anchors on delta = 0.1 + 0.3 x, interleaved with 2 IM-less anchors
    // (indices 2 and 5) that must be masked out even though they survive nothing.
    let mut anchors = Vec::new();
    let mut im_indices = Vec::new();
    let xs = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (k, &x) in xs.iter().enumerate() {
        if k == 2 {
            anchors.push(anchor_ppm(x, 0.0)); // no IM
        }
        im_indices.push(anchors.len());
        anchors.push(anchor_im(x, 0.1 + 0.3 * x));
        if k == 2 {
            anchors.push(anchor_ppm(x, 0.0)); // no IM
        }
    }
    let (fit, active) = fit_im_drift(&anchors);
    assert!(
        (fit.intercept - 0.1).abs() < 1e-6,
        "intercept {}",
        fit.intercept
    );
    assert!((fit.slope - 0.3).abs() < 1e-6, "slope {}", fit.slope);
    for (i, &ok) in active.iter().enumerate() {
        assert_eq!(ok, im_indices.contains(&i), "mask mismatch at {i}");
    }
}
