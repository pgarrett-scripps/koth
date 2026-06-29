use serde::Serialize;

use super::{anchors::AnchorPair, AlignmentConfig};

/// Piecewise-linear RT warp in normalised [0, 1] space.
///
/// Knots are (run_rt_norm, delta) pairs where delta = ref_rt_norm - run_rt_norm.
/// apply() converts absolute run RT → absolute reference RT.
#[derive(Debug, Clone, Serialize)]
pub struct RtWarp {
    pub knot_x: Vec<f64>,
    pub knot_y: Vec<f64>,
}

impl RtWarp {
    /// Map an absolute run RT value into the reference RT space.
    pub fn apply(&self, run_rt: f64, run_rt_range: (f64, f64), ref_rt_range: (f64, f64)) -> f64 {
        let run_norm = normalize(run_rt, run_rt_range);
        let delta = piecewise_linear(&self.knot_x, &self.knot_y, run_norm);
        let ref_norm = (run_norm + delta).clamp(0.0, 1.0);
        denormalize(ref_norm, ref_rt_range)
    }

    /// Evaluate the warp delta at a position given in reference-normalised space.
    /// Used to approximately invert the warp: run_norm ≈ ref_norm - delta_at_ref_norm(ref_norm).
    /// This approximation holds well for the small-to-moderate shifts typical in LC-MS alignment.
    pub fn delta_at_ref_norm(&self, ref_norm: f64) -> f64 {
        piecewise_linear(&self.knot_x, &self.knot_y, ref_norm)
    }

    /// Evaluate the warp delta at a position given in run-normalised space.
    /// This is the value used directly by `apply()` — the knots are keyed on
    /// run-normalised RT — so it is the correct quantity for residual analysis.
    pub fn delta_at_run_norm(&self, run_norm: f64) -> f64 {
        piecewise_linear(&self.knot_x, &self.knot_y, run_norm)
    }
}

/// Identity warp — used when anchor count is too low.
pub fn identity_warp() -> RtWarp {
    RtWarp {
        knot_x: vec![0.0, 1.0],
        knot_y: vec![0.0, 0.0],
    }
}

/// Fit a *linear* RT warp by least squares on anchor pairs.
///
/// Models `ref_norm = a · run_norm + b` and returns the equivalent
/// `RtWarp` (delta = ref_norm − run_norm), expressed as two knots at
/// run_norm = 0 and run_norm = 1.
///
/// Compared to the piecewise-local-median warp, this is the lowest-noise
/// warp that still corrects for global RT drift between runs — useful when
/// the local-median warp over-fits anchor density and scatters projected
/// features beyond the consensus grouping tolerance. Iterative
/// sigma-clipping (same schedule as `fit_rt_warp`) removes outlier anchors
/// before the final fit.
pub fn fit_linear_warp(anchors: &[AnchorPair], config: &AlignmentConfig) -> (RtWarp, Vec<bool>) {
    let mut active = vec![true; anchors.len()];

    for _ in 0..config.rt_warp_clip_iters {
        let (a, b) = lstsq(anchors, &active);

        let residuals: Vec<f64> = anchors
            .iter()
            .zip(active.iter())
            .filter_map(|(ap, &ok)| {
                if !ok {
                    return None;
                }
                Some(ap.ref_rt_norm - (a * ap.run_rt_norm + b))
            })
            .collect();

        if residuals.is_empty() {
            break;
        }
        let std = std_dev(&residuals);
        if std < 1e-9 {
            break;
        }
        let threshold = config.rt_warp_sigma_clip * std;

        let mut res_iter = residuals.iter();
        for (_, ok) in anchors.iter().zip(active.iter_mut()) {
            if !*ok {
                continue;
            }
            if res_iter.next().unwrap().abs() > threshold {
                *ok = false;
            }
        }

        let n_active = active.iter().filter(|&&ok| ok).count();
        if n_active < config.min_anchor_count {
            active.iter_mut().for_each(|ok| *ok = true);
            break;
        }
    }

    let (a, b) = lstsq(anchors, &active);
    // Express ref_norm = a · run_norm + b as the delta(run_norm) form used by `apply`:
    //   delta(x) = ref_norm(x) − x = (a − 1) · x + b
    let delta_at_0 = b;
    let delta_at_1 = a - 1.0 + b;
    (
        RtWarp {
            knot_x: vec![0.0, 1.0],
            knot_y: vec![delta_at_0, delta_at_1],
        },
        active,
    )
}

/// Ordinary least squares `ref = a · run + b` over the active anchor mask.
/// Falls back to identity (`a=1, b=0`) on degenerate cases.
fn lstsq(anchors: &[AnchorPair], active: &[bool]) -> (f64, f64) {
    let n: usize = active.iter().filter(|&&ok| ok).count();
    if n < 2 {
        return (1.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (ap, &ok) in anchors.iter().zip(active.iter()) {
        if !ok {
            continue;
        }
        let x = ap.run_rt_norm;
        let y = ap.ref_rt_norm;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let n_f = n as f64;
    let denom = n_f * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return (1.0, 0.0);
    }
    let a = (n_f * sxy - sx * sy) / denom;
    let b = (sy - a * sx) / n_f;
    (a, b)
}

/// Number of random 2-point line samples drawn during RANSAC consensus search.
/// Fixed (not a config knob) — 4000 is far above saturation for the anchor
/// counts seen here, and keeping it internal is part of the point of this mode:
/// fewer settings than the sigma-clip schedule it replaces.
const RANSAC_ITERS: usize = 4000;

/// Maximum refinement passes after the initial RANSAC consensus: refit the
/// median PWL on the current inliers, re-collect inliers, repeat until the
/// inlier set stops growing.
const RANSAC_REFINE_PASSES: usize = 4;

/// Fit a RT warp by RANSAC inlier selection, then a median piecewise-linear fit
/// (with isotonic monotonicity, via [`build_knots`]) on the inliers.
///
/// Anchor matching is purely coordinate-based (charge + ppm + RT window), so a
/// sizeable fraction of anchors are wrong-peptide matches scattered uniformly
/// across the RT window. Sigma-clipping refits on that contaminated set and is
/// fragile in sparse windows; RANSAC instead defines the diagonal from a global
/// 2-point line consensus over *all* anchors, then refines with the robust
/// median PWL on the consensus inliers. Returns the warp and a per-anchor mask
/// (`true` = RANSAC inlier), aligned 1:1 with `anchors`.
///
/// Deterministic: the sample stream comes from a fixed-seed LCG, so repeated
/// runs on the same anchors produce identical warps.
pub fn fit_ransac_warp(anchors: &[AnchorPair], config: &AlignmentConfig) -> (RtWarp, Vec<bool>) {
    let n = anchors.len();
    if n < config.min_anchor_count {
        return (identity_warp(), vec![false; n]);
    }

    let xs: Vec<f64> = anchors.iter().map(|a| a.run_rt_norm).collect();
    let ys: Vec<f64> = anchors.iter().map(|a| a.ref_rt_norm).collect();
    let thresh = config.rt_warp_ransac_thresh;

    // Deterministic LCG (same constants as uno's calibrate::fit_ransac).
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = |m: usize| -> usize {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as usize) % m
    };

    // Initial consensus: best 2-point line `ref_norm = m·run_norm + b`.
    let mut best_mask = vec![true; n];
    let mut best_count = 0usize;
    for _ in 0..RANSAC_ITERS {
        let (i, j) = (next(n), next(n));
        if (xs[i] - xs[j]).abs() < 1e-9 {
            continue;
        }
        let m = (ys[j] - ys[i]) / (xs[j] - xs[i]);
        let b = ys[i] - m * xs[i];
        let mut mask = vec![false; n];
        let mut count = 0usize;
        for k in 0..n {
            if (ys[k] - (m * xs[k] + b)).abs() < thresh {
                mask[k] = true;
                count += 1;
            }
        }
        if count > best_count {
            best_count = count;
            best_mask = mask;
        }
    }

    // Refine: median PWL on inliers → re-collect inliers, until it stops growing.
    let mut active = best_mask;
    for _ in 0..RANSAC_REFINE_PASSES {
        let prev = active.iter().filter(|&&a| a).count();
        let (kx, ky) = build_knots(anchors, &active, config.rt_warp_bandwidth);
        let mut new = vec![false; n];
        let mut count = 0usize;
        for k in 0..n {
            let pred = xs[k] + piecewise_linear(&kx, &ky, xs[k]);
            if (ys[k] - pred).abs() < thresh {
                new[k] = true;
                count += 1;
            }
        }
        let grew = count > prev;
        active = new;
        if !grew {
            break;
        }
    }

    // Too few inliers ⇒ data too sparse/contaminated for a confident consensus;
    // fall back to the robust median fit over all anchors (still monotone).
    if active.iter().filter(|&&a| a).count() < config.min_anchor_count {
        active = vec![true; n];
    }

    let (kx, ky) = build_knots(anchors, &active, config.rt_warp_bandwidth);
    (RtWarp { knot_x: kx, knot_y: ky }, active)
}

/// Fit a RT warp from anchor pairs using sliding-window medians + sigma-clipping.
///
/// Returns the warp plus a per-anchor active mask (true = used by the final fit).
pub fn fit_rt_warp(anchors: &[AnchorPair], config: &AlignmentConfig) -> (RtWarp, Vec<bool>) {
    let mut active = vec![true; anchors.len()];

    for _ in 0..config.rt_warp_clip_iters {
        let (kx, ky) = build_knots(anchors, &active, config.rt_warp_bandwidth);

        let residuals: Vec<f64> = anchors
            .iter()
            .zip(active.iter())
            .filter_map(|(a, &ok)| {
                if !ok {
                    return None;
                }
                let predicted = piecewise_linear(&kx, &ky, a.run_rt_norm);
                let actual = a.ref_rt_norm - a.run_rt_norm;
                Some(actual - predicted)
            })
            .collect();

        if residuals.is_empty() {
            break;
        }
        let std = std_dev(&residuals);
        if std < 1e-9 {
            break;
        }
        let threshold = config.rt_warp_sigma_clip * std;

        let mut res_iter = residuals.iter();
        for (_, ok) in anchors.iter().zip(active.iter_mut()) {
            if !*ok {
                continue;
            }
            if res_iter.next().unwrap().abs() > threshold {
                *ok = false;
            }
        }

        let n_active = active.iter().filter(|&&ok| ok).count();
        if n_active < config.min_anchor_count {
            // Clipping too aggressive — restore and stop
            active.iter_mut().for_each(|ok| *ok = true);
            break;
        }
    }

    let (kx, ky) = build_knots(anchors, &active, config.rt_warp_bandwidth);
    (RtWarp { knot_x: kx, knot_y: ky }, active)
}

/// Pool-adjacent-violators (PAVA) isotonic regression: the least-squares
/// non-decreasing fit to `y` with equal weights (one per knot). Returns one
/// fitted value per input element. Used to project the warp's mapped knot
/// positions onto the nearest monotone sequence; ported from uno's
/// `calibrate::isotonic`, but tracking block lengths explicitly so the
/// expansion is correct regardless of pooled-weight magnitude.
fn isotonic(y: &[f64]) -> Vec<f64> {
    let mut vals: Vec<f64> = Vec::with_capacity(y.len());
    let mut wts: Vec<f64> = Vec::with_capacity(y.len());
    let mut len: Vec<usize> = Vec::with_capacity(y.len());
    for &yi in y {
        vals.push(yi);
        wts.push(1.0);
        len.push(1);
        // Merge the last block into its predecessor while it violates monotonicity.
        while vals.len() > 1 && vals[vals.len() - 2] > vals[vals.len() - 1] {
            let (v2, w2, l2) = (vals.pop().unwrap(), wts.pop().unwrap(), len.pop().unwrap());
            let (v1, w1, l1) = (vals.pop().unwrap(), wts.pop().unwrap(), len.pop().unwrap());
            let w = w1 + w2;
            vals.push((v1 * w1 + v2 * w2) / w);
            wts.push(w);
            len.push(l1 + l2);
        }
    }
    let mut out = Vec::with_capacity(y.len());
    for (&v, &l) in vals.iter().zip(&len) {
        for _ in 0..l {
            out.push(v);
        }
    }
    out
}

fn build_knots(anchors: &[AnchorPair], active: &[bool], bandwidth: f64) -> (Vec<f64>, Vec<f64>) {
    let active_anchors: Vec<&AnchorPair> = anchors
        .iter()
        .zip(active.iter())
        .filter_map(|(a, &ok)| if ok { Some(a) } else { None })
        .collect();

    if active_anchors.is_empty() {
        return (vec![0.0, 1.0], vec![0.0, 0.0]);
    }

    let half_bw = bandwidth / 2.0;
    let step = bandwidth / 2.0;
    let mut kx: Vec<f64> = Vec::new();
    let mut ky: Vec<f64> = Vec::new();

    let mut c = 0.0f64;
    while c <= 1.0 + step * 0.5 {
        let center = c.min(1.0);
        let deltas: Vec<f64> = active_anchors
            .iter()
            .filter(|a| (a.run_rt_norm - center).abs() <= half_bw)
            .map(|a| a.ref_rt_norm - a.run_rt_norm)
            .collect();

        if deltas.len() >= 3 {
            kx.push(center);
            ky.push(median(&deltas));
        }
        c += step;
        if c > 1.0 + step * 0.5 {
            break;
        }
    }

    if kx.is_empty() {
        // Fall back to global median delta
        let deltas: Vec<f64> = active_anchors
            .iter()
            .map(|a| a.ref_rt_norm - a.run_rt_norm)
            .collect();
        let m = median(&deltas);
        return (vec![0.0, 1.0], vec![m, m]);
    }

    // Pad boundaries with constant extrapolation
    if kx[0] > 1e-6 {
        kx.insert(0, 0.0);
        ky.insert(0, ky[0]);
    }
    if *kx.last().unwrap() < 1.0 - 1e-6 {
        kx.push(1.0);
        ky.push(*ky.last().unwrap());
    }

    // Enforce a monotone (non-decreasing) warp. `apply()` maps a knot to
    // `m_i = kx_i + ky_i` (ref_norm = run_norm + delta), so the warp folds the
    // RT axis back on itself wherever m_i+1 < m_i — which a noisy median delta
    // in a sparse-anchor window can produce, even though elution order cannot
    // invert between runs. Project the mapped knots onto the nearest monotone
    // sequence with PAVA isotonic regression, then convert back to delta form.
    let mapped: Vec<f64> = kx.iter().zip(ky.iter()).map(|(x, y)| x + y).collect();
    let mapped = isotonic(&mapped);
    let ky: Vec<f64> = kx.iter().zip(mapped.iter()).map(|(x, m)| m - x).collect();

    (kx, ky)
}

/// Piecewise-linear interpolation; clamps to boundary values outside range.
pub fn piecewise_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    if xs.is_empty() {
        return 0.0;
    }
    if xs.len() == 1 {
        return ys[0];
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= *xs.last().unwrap() {
        return *ys.last().unwrap();
    }
    let pos = xs.partition_point(|&xi| xi <= x);
    let i = pos.saturating_sub(1).min(xs.len() - 2);
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    ys[i] + t * (ys[i + 1] - ys[i])
}

fn normalize(rt: f64, range: (f64, f64)) -> f64 {
    let span = range.1 - range.0;
    if span.abs() < 1e-9 {
        return 0.5;
    }
    ((rt - range.0) / span).clamp(0.0, 1.0)
}

fn denormalize(norm: f64, range: (f64, f64)) -> f64 {
    norm * (range.1 - range.0) + range.0
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

fn std_dev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64;
    var.sqrt()
}

#[cfg(test)]
mod tests {
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
        let config = AlignmentConfig::default();
        let mut anchors = Vec::new();
        // Clusters of 4 anchors per run-RT center. Reference RT rises, then
        // sharply inverts in the middle, then rises again.
        let centers = [
            (0.10, 0.10), (0.25, 0.30), (0.40, 0.65),
            (0.55, 0.25), (0.70, 0.30), (0.85, 0.90),
        ];
        for &(run, refn) in &centers {
            for k in 0..4 {
                let jitter = (k as f64 - 1.5) * 0.005;
                anchors.push(anchor(run + jitter, refn + jitter));
            }
        }
        let (warp, _) = fit_rt_warp(&anchors, &config);
        assert_warp_monotone(&warp);
    }
}
