use crate::models::Spectrum;

/// Remove per-scan noise using iterative sigma clipping.
///
/// Works by repeatedly discarding signal peaks (the tall ones) until only the
/// noise floor is left, then setting a threshold at `median + sigma * MAD_sigma`
/// of that noise floor.  Because real LC-MS signal is orders of magnitude above
/// the noise, convergence is fast (typically 2–5 iterations).
///
/// A `sigma` of 3.0 is a reasonable starting point; lower values cut more
/// aggressively.
pub fn filter_spectrum(spectrum: &mut Spectrum, sigma: f64) {
    if spectrum.peaks.len() < 3 {
        return;
    }

    let mut noise: Vec<f32> = spectrum.peaks.iter().map(|p| p.intensity).collect();
    noise.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    loop {
        let threshold = noise_threshold(&noise, sigma);
        let new_end = noise.partition_point(|&x| x <= threshold);
        if new_end == noise.len() {
            break;
        }
        noise.truncate(new_end);
        if noise.len() < 3 {
            break;
        }
    }

    let threshold = noise_threshold(&noise, sigma);
    spectrum.peaks.retain(|p| p.intensity > threshold);
}

/// Compute `median + sigma * (1.4826 * MAD)` from a sorted intensity slice.
///
/// The 1.4826 factor makes MAD a consistent estimator of σ for Gaussian noise.
fn noise_threshold(sorted: &[f32], sigma: f64) -> f32 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }

    let median = if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    };

    let mut deviations: Vec<f32> = sorted.iter().map(|&x| (x - median).abs()).collect();
    deviations.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = if n % 2 == 0 {
        (deviations[n / 2 - 1] + deviations[n / 2]) / 2.0
    } else {
        deviations[n / 2]
    };

    median + (sigma * 1.4826 * mad as f64) as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{Peak, Spectrum};

    fn peak(intensity: f32) -> Peak {
        Peak { mz: 500.0, intensity, ion_mobility: 0.0 }
    }

    fn spectrum(intensities: &[f32]) -> Spectrum {
        Spectrum {
            scan_index: 0,
            retention_time: 0.0,
            peaks: intensities.iter().map(|&i| peak(i)).collect(),
            ms_level: 1,
            isolation_window: None,
        }
    }

    #[test]
    fn threshold_empty_is_zero() {
        assert_eq!(noise_threshold(&[], 3.0), 0.0);
    }

    #[test]
    fn threshold_is_median_plus_scaled_mad() {
        // sorted [1,2,3,4,5]: median 3, |x-3| = [2,1,0,1,2] -> MAD 1.
        // threshold = 3 + 2 * 1.4826 * 1 = 5.9652.
        let t = noise_threshold(&[1.0, 2.0, 3.0, 4.0, 5.0], 2.0);
        assert!((t - (3.0 + 2.0 * 1.4826)).abs() < 1e-4, "threshold {t}");
    }

    #[test]
    fn zero_mad_floor_gives_median() {
        // A flat noise floor has MAD 0, so the threshold is exactly the median.
        assert_eq!(noise_threshold(&[10.0; 5], 3.0), 10.0);
    }

    #[test]
    fn filter_removes_flat_floor_keeps_spikes() {
        // Eight equal noise peaks + two tall spikes. Clipping strips the spikes
        // from the noise set, leaving a flat floor with threshold == 1.0; the
        // spectrum then retains only intensities strictly above it.
        let mut s = spectrum(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1000.0, 1000.0]);
        filter_spectrum(&mut s, 3.0);
        assert_eq!(s.peaks.len(), 2, "only the two spikes survive");
        assert!(s.peaks.iter().all(|p| p.intensity == 1000.0));
    }

    #[test]
    fn filter_is_noop_below_three_peaks() {
        let mut s = spectrum(&[1.0, 9999.0]);
        filter_spectrum(&mut s, 3.0);
        assert_eq!(s.peaks.len(), 2, "fewer than 3 peaks are left untouched");
    }
}
