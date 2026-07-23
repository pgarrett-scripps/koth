//! Empirical m/z tolerance calibration via the absolute-ppm delta histogram
//! collected from accepted peak-to-hill matches during a wide-tolerance
//! "pass 1" hill detection.
//!
//! Inspired by AlphaPept's two-pass calibration: pass 1 runs at a permissive
//! ppm; the actual ppm-delta distribution of accepted matches is recorded;
//! pass 2 uses `median + N·σ` as the calibrated tolerance.

/// Fixed-resolution histogram over absolute ppm deltas. Memory is bounded
/// by `n_bins` regardless of how many samples are recorded, so a 95-min
/// Orbitrap or a long Bruker .d run stay flat.
#[derive(Debug, Clone)]
pub struct MzDeltaHistogram {
    bins: Vec<u32>,
    bin_width_ppm: f64,
    max_ppm: f64,
    /// Total samples landed in `bins` (excludes overflow).
    in_range: u64,
    /// Samples that exceeded `max_ppm` (still counted for total, not used
    /// for median/std).
    over_count: u64,
}

impl MzDeltaHistogram {
    /// `max_ppm` is the upper edge of the histogram; samples beyond it
    /// fall into the overflow bucket. `bin_width_ppm` is the bin
    /// resolution (e.g. 0.01 ppm).
    pub fn new(max_ppm: f64, bin_width_ppm: f64) -> Self {
        assert!(max_ppm > 0.0 && bin_width_ppm > 0.0 && bin_width_ppm < max_ppm);
        let n_bins = (max_ppm / bin_width_ppm).ceil() as usize;
        Self {
            bins: vec![0; n_bins],
            bin_width_ppm,
            max_ppm,
            in_range: 0,
            over_count: 0,
        }
    }

    /// Record one sample. Negatives, NaN, and infinities are silently
    /// dropped — only finite, non-negative ppm deltas are counted.
    #[inline]
    pub fn record(&mut self, ppm: f64) {
        if !ppm.is_finite() || ppm < 0.0 {
            return;
        }
        let bin = (ppm / self.bin_width_ppm) as usize;
        if bin >= self.bins.len() {
            self.over_count += 1;
        } else {
            self.bins[bin] = self.bins[bin].saturating_add(1);
            self.in_range += 1;
        }
    }

    pub fn count_in_range(&self) -> u64 {
        self.in_range
    }

    pub fn count_total(&self) -> u64 {
        self.in_range + self.over_count
    }

    /// Median of the in-range samples, computed by walking the cumulative
    /// histogram. Returns `None` if no in-range samples were recorded.
    pub fn median(&self) -> Option<f64> {
        if self.in_range == 0 {
            return None;
        }
        // Use the lower-median convention: first bin whose cumulative
        // count reaches `ceil(n/2)`.
        let target = (self.in_range + 1) / 2;
        let mut cum: u64 = 0;
        for (i, &c) in self.bins.iter().enumerate() {
            cum += c as u64;
            if cum >= target {
                return Some((i as f64 + 0.5) * self.bin_width_ppm);
            }
        }
        Some(self.max_ppm)
    }

    /// Standard deviation of the in-range samples, computed from bin
    /// centers. Two-pass for numerical stability isn't necessary here —
    /// values are bounded in [0, max_ppm], so the naive single-pass
    /// `E[x²] - E[x]²` is fine.
    pub fn std_dev(&self) -> Option<f64> {
        if self.in_range < 2 {
            return None;
        }
        let n = self.in_range as f64;
        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        for (i, &c) in self.bins.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let center = (i as f64 + 0.5) * self.bin_width_ppm;
            let count = c as f64;
            sum += center * count;
            sum_sq += center * center * count;
        }
        let mean = sum / n;
        let var = (sum_sq / n - mean * mean).max(0.0);
        Some(var.sqrt())
    }

    /// `median + sigma_mult × std_dev`. Returns `None` when fewer than 2
    /// in-range samples were recorded (calibration is not statistically
    /// meaningful).
    pub fn calibrated_ppm(&self, sigma_mult: f64) -> Option<f64> {
        let m = self.median()?;
        let s = self.std_dev()?;
        Some(m + sigma_mult * s)
    }
}

#[cfg(test)]
#[path = "calibration_tests.rs"]
mod tests;
