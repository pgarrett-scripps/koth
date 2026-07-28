use serde::Serialize;

use super::anchors::AnchorPair;

#[cfg(test)]
#[path = "drift_tests.rs"]
mod tests;

/// Sigma-clip schedule for the mass / IM drift fits. Formerly shared the RT
/// warp's `rt_warp_clip_iters` / `rt_warp_sigma_clip` config knobs; those were
/// removed with the legacy piecewise/linear warps (the shipped configs and the
/// old defaults all used these exact values), so the schedule is now fixed.
const DRIFT_CLIP_ITERS: usize = 5;
const DRIFT_SIGMA_CLIP: f64 = 3.0;

/// A simple linear drift model: y = intercept + slope * x.
#[derive(Debug, Clone, Serialize)]
pub struct DriftFit {
    pub intercept: f64,
    pub slope: f64,
}

impl DriftFit {
    pub fn zero() -> Self {
        Self { intercept: 0.0, slope: 0.0 }
    }

    /// Predict the drift value at a given (normalised) RT.
    pub fn predict(&self, rt_norm: f64) -> f64 {
        self.intercept + self.slope * rt_norm
    }
}

/// Fit mass PPM drift as a linear function of normalised reference RT.
///
/// Returns the fit plus a per-anchor active mask (true = used by the final fit,
/// false = sigma-clipped). The mask is aligned 1:1 with the input `anchors` slice.
pub fn fit_mass_drift(anchors: &[AnchorPair]) -> (DriftFit, Vec<bool>) {
    let pairs: Vec<(f64, f64)> = anchors
        .iter()
        .map(|a| (a.ref_rt_norm, a.ppm_error()))
        .collect();
    fit_with_sigma_clip(&pairs)
}

/// Fit ion-mobility drift as a linear function of normalised reference RT.
/// Returns zero drift if fewer than 4 IM anchor pairs are available.
///
/// The returned mask is aligned 1:1 with the input `anchors` slice: an entry is
/// `true` only if the anchor had IM on both sides AND survived sigma-clipping.
pub fn fit_im_drift(anchors: &[AnchorPair]) -> (DriftFit, Vec<bool>) {
    let mut full_active = vec![false; anchors.len()];

    let mut idx_map: Vec<usize> = Vec::with_capacity(anchors.len());
    let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(anchors.len());
    for (i, a) in anchors.iter().enumerate() {
        if a.ref_im != 0.0 && a.run_im != 0.0 {
            idx_map.push(i);
            pairs.push((a.ref_rt_norm, a.im_delta()));
        }
    }

    if pairs.len() < 4 {
        return (DriftFit::zero(), full_active);
    }

    let (fit, sub_active) = fit_with_sigma_clip(&pairs);
    for (j, &full_idx) in idx_map.iter().enumerate() {
        full_active[full_idx] = sub_active[j];
    }
    (fit, full_active)
}

fn fit_with_sigma_clip(pairs: &[(f64, f64)]) -> (DriftFit, Vec<bool>) {
    if pairs.len() < 2 {
        return (DriftFit::zero(), vec![false; pairs.len()]);
    }

    let mut active = vec![true; pairs.len()];
    let mut fit = ols_active(pairs, &active);

    for _ in 0..DRIFT_CLIP_ITERS {
        let residuals: Vec<f64> = pairs
            .iter()
            .zip(active.iter())
            .filter_map(|((x, y), &ok)| {
                if ok {
                    Some(y - fit.predict(*x))
                } else {
                    None
                }
            })
            .collect();

        if residuals.is_empty() {
            break;
        }
        let std = std_dev(&residuals);
        if std < 1e-9 {
            break;
        }
        let threshold = DRIFT_SIGMA_CLIP * std;

        let mut res_iter = residuals.iter();
        for ((x, y), ok) in pairs.iter().zip(active.iter_mut()) {
            let _ = (x, y);
            if !*ok {
                continue;
            }
            if res_iter.next().unwrap().abs() > threshold {
                *ok = false;
            }
        }

        let n_active = active.iter().filter(|&&ok| ok).count();
        if n_active < 2 {
            active.iter_mut().for_each(|ok| *ok = true);
            break;
        }
        fit = ols_active(pairs, &active);
    }

    (fit, active)
}

fn ols_active(pairs: &[(f64, f64)], active: &[bool]) -> DriftFit {
    let active_pairs: Vec<(f64, f64)> = pairs
        .iter()
        .zip(active.iter())
        .filter_map(|(&(x, y), &ok)| if ok { Some((x, y)) } else { None })
        .collect();
    ols(&active_pairs)
}

fn ols(pairs: &[(f64, f64)]) -> DriftFit {
    let n = pairs.len() as f64;
    if n < 2.0 {
        return DriftFit::zero();
    }
    let sx: f64 = pairs.iter().map(|(x, _)| x).sum();
    let sy: f64 = pairs.iter().map(|(_, y)| y).sum();
    let sxx: f64 = pairs.iter().map(|(x, _)| x * x).sum();
    let sxy: f64 = pairs.iter().map(|(x, y)| x * y).sum();

    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return DriftFit { intercept: sy / n, slope: 0.0 };
    }
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    DriftFit { intercept, slope }
}

fn std_dev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64;
    var.sqrt()
}
