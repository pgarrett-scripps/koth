use super::*;
use crate::alignment::{AlignmentConfig, WarpKind};

#[test]
fn isotonic_pools_violators() {
    // The dip at index 2 must be pooled with index 1 into their mean (2.5).
    let out = isotonic(&[1.0, 3.0, 2.0, 4.0]);
    assert_eq!(out.len(), 4);
    assert!((out[0] - 1.0).abs() < 1e-9);
    assert!((out[1] - 2.5).abs() < 1e-9);
    assert!((out[2] - 2.5).abs() < 1e-9);
    assert!((out[3] - 4.0).abs() < 1e-9);
    // Result is non-decreasing.
    assert!(out.windows(2).all(|w| w[1] >= w[0] - 1e-12));
}

#[test]
fn isotonic_leaves_monotone_untouched() {
    let y = vec![0.0, 0.1, 0.1, 0.5, 0.9];
    assert_eq!(isotonic(&y), y);
}

fn anchor(run: f64, refn: f64) -> AnchorPair {
    AnchorPair {
        ref_rt_norm: refn,
        run_rt_norm: run,
        ref_mz: 500.0,
        run_mz: 500.0,
        ref_im: 0.0,
        run_im: 0.0,
    }
}

/// Adversarial anchors: the reference RT *decreases* as the run RT increases
/// over the middle of the gradient, which without monotonicity enforcement
/// folds the warp back on itself. The fitted warp must still be a
/// non-decreasing map of run RT → reference RT.
fn assert_warp_monotone(warp: &RtWarp) {
    let mut prev = f64::NEG_INFINITY;
    for i in 0..=200 {
        let run_norm = i as f64 / 200.0;
        let mapped = run_norm + piecewise_linear(&warp.knot_x, &warp.knot_y, run_norm);
        assert!(
            mapped >= prev - 1e-9,
            "warp inverted at run_norm={run_norm}: {mapped} < {prev}"
        );
        prev = mapped;
    }
}

#[test]
fn ransac_recovers_diagonal_under_contamination() {
    let mut config = AlignmentConfig::default();
    config.rt_warp_kind = WarpKind::Ransac;
    config.rt_warp_ransac_thresh = 0.01;

    // 200 true anchors on a mild-drift diagonal (ref = run + 0.02) with tiny
    // jitter, plus 100 uniform-background contaminants (~33% contamination).
    let mut anchors = Vec::new();
    let mut seed: u64 = 12345;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((seed >> 33) as f64) / (1u64 << 31) as f64 // [0,1)
    };
    let n_true = 200;
    for i in 0..n_true {
        let run = (i as f64 + 0.5) / n_true as f64;
        let jitter = (rnd() - 0.5) * 0.004;
        anchors.push(anchor(run, (run + 0.02 + jitter).clamp(0.0, 1.0)));
    }
    let n_contam = 100;
    for _ in 0..n_contam {
        anchors.push(anchor(rnd(), rnd()));
    }

    let (warp, mask) = fit_ransac_warp(&anchors, &config);
    assert_warp_monotone(&warp);

    // True anchors (first n_true) should mostly be inliers; contaminants mostly not.
    let true_inliers = mask[..n_true].iter().filter(|&&m| m).count();
    let contam_inliers = mask[n_true..].iter().filter(|&&m| m).count();
    assert!(
        true_inliers as f64 / n_true as f64 > 0.9,
        "RANSAC should keep >90% of true anchors, kept {true_inliers}/{n_true}"
    );
    // A uniform contaminant lands within the ±0.01 band of the diagonal only
    // ~2% of the time by chance, so the inlier set should be nearly pure.
    assert!(
        (contam_inliers as f64) / (n_contam as f64) < 0.15,
        "RANSAC should reject most contaminants, kept {contam_inliers}/{n_contam}"
    );

    // The recovered warp delta near mid-gradient should be ~+0.02.
    let delta_mid = piecewise_linear(&warp.knot_x, &warp.knot_y, 0.5);
    assert!(
        (delta_mid - 0.02).abs() < 0.01,
        "expected mid-gradient delta ≈ 0.02, got {delta_mid}"
    );
}

#[test]
fn warp_is_monotone_on_inverting_anchors() {
    let mut config = AlignmentConfig::default();
    config.rt_warp_kind = WarpKind::Ransac;
    config.rt_warp_ransac_thresh = 0.05;
    let mut anchors = Vec::new();
    // Clusters of 4 anchors per run-RT center. Reference RT rises, then
    // sharply inverts in the middle, then rises again.
    let centers = [
        (0.10, 0.10),
        (0.25, 0.30),
        (0.40, 0.65),
        (0.55, 0.25),
        (0.70, 0.30),
        (0.85, 0.90),
    ];
    for &(run, refn) in &centers {
        for k in 0..4 {
            let jitter = (k as f64 - 1.5) * 0.005;
            anchors.push(anchor(run + jitter, refn + jitter));
        }
    }
    let (warp, _) = fit_ransac_warp(&anchors, &config);
    assert_warp_monotone(&warp);
}
