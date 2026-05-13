use super::{anchors::AnchorPair, AlignmentConfig};

/// A simple linear drift model: y = intercept + slope * x.
#[derive(Debug, Clone)]
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
pub fn fit_mass_drift(anchors: &[AnchorPair], config: &AlignmentConfig) -> DriftFit {
    let pairs: Vec<(f64, f64)> = anchors
        .iter()
        .map(|a| (a.ref_rt_norm, a.ppm_error()))
        .collect();
    fit_with_sigma_clip(&pairs, config)
}

/// Fit ion-mobility drift as a linear function of normalised reference RT.
/// Returns zero drift if fewer than 4 IM anchor pairs are available.
pub fn fit_im_drift(anchors: &[AnchorPair], config: &AlignmentConfig) -> DriftFit {
    let pairs: Vec<(f64, f64)> = anchors
        .iter()
        .filter(|a| a.ref_im != 0.0 && a.run_im != 0.0)
        .map(|a| (a.ref_rt_norm, a.im_delta()))
        .collect();
    if pairs.len() < 4 {
        return DriftFit::zero();
    }
    fit_with_sigma_clip(&pairs, config)
}

fn fit_with_sigma_clip(pairs: &[(f64, f64)], config: &AlignmentConfig) -> DriftFit {
    if pairs.len() < 2 {
        return DriftFit::zero();
    }

    let mut active = vec![true; pairs.len()];
    let mut fit = ols_active(pairs, &active);

    for _ in 0..config.rt_warp_clip_iters {
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
        let threshold = config.rt_warp_sigma_clip * std;

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

    fit
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
