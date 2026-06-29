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

/// Keep only the `n` most intense peaks in the spectrum, preserving m/z order.
///
/// Mirrors AlphaPept's `n_most_abundant` per-scan peak cap. A no-op when the
/// spectrum already has `<= n` peaks.
pub fn keep_most_abundant(spectrum: &mut Spectrum, n: usize) {
    if spectrum.peaks.len() <= n {
        return;
    }
    // Partition so the n highest-intensity peaks are first, then keep them.
    let mut idx: Vec<usize> = (0..spectrum.peaks.len()).collect();
    idx.select_nth_unstable_by(n, |&a, &b| {
        spectrum.peaks[b]
            .intensity
            .partial_cmp(&spectrum.peaks[a].intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    idx.truncate(n);
    idx.sort_unstable(); // restore m/z (index) order
    spectrum.peaks = idx.into_iter().map(|i| spectrum.peaks[i].clone()).collect();
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
