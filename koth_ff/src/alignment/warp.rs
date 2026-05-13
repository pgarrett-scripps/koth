use super::{anchors::AnchorPair, AlignmentConfig};

/// Piecewise-linear RT warp in normalised [0, 1] space.
///
/// Knots are (run_rt_norm, delta) pairs where delta = ref_rt_norm - run_rt_norm.
/// apply() converts absolute run RT → absolute reference RT.
#[derive(Debug, Clone)]
pub struct RtWarp {
    knot_x: Vec<f64>,
    knot_y: Vec<f64>,
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
}

/// Identity warp — used when anchor count is too low.
pub fn identity_warp() -> RtWarp {
    RtWarp {
        knot_x: vec![0.0, 1.0],
        knot_y: vec![0.0, 0.0],
    }
}

/// Fit a RT warp from anchor pairs using sliding-window medians + sigma-clipping.
pub fn fit_rt_warp(anchors: &[AnchorPair], config: &AlignmentConfig) -> RtWarp {
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
    RtWarp { knot_x: kx, knot_y: ky }
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
