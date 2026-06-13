/// An active (still-being-built) chromatographic hill.
///
/// Profiles use f32 with f32::NAN to mark gap positions (scans where no
/// peak was matched). Intensity at gaps is always 0.0. Using f32 cuts
/// per-element storage from 16 bytes (Option<f64>) to 4 bytes.
#[derive(Debug)]
pub struct ActiveHill {
    /// Running mean mz (for fast distance comparison, updated online)
    pub running_mz_mean: f64,
    /// Running mean ion mobility (0.0 if IM not in use)
    pub running_im_mean: f64,
    /// Absolute scan index when this hill started
    pub scan_start: usize,
    /// Last scan index at which a real peak was added
    pub last_scan_seen: usize,
    /// m/z at each position (f32::NAN at gap positions)
    pub mz_profile: Vec<f32>,
    /// Intensity at each position (0.0 at gaps)
    pub intensity_profile: Vec<f32>,
    /// Retention time at each position (f32::NAN at gap positions)
    pub rt_profile: Vec<f32>,
    /// Ion mobility at each position (f32::NAN at gap positions)
    pub im_profile: Vec<f32>,
    /// Monotonically increasing count of real peaks added (for running mean update)
    pub n_real: usize,
}

impl ActiveHill {
    pub fn new(mz: f64, intensity: f64, rt: f64, im: f64, scan_index: usize) -> Self {
        Self {
            running_mz_mean: mz,
            running_im_mean: im,
            scan_start: scan_index,
            last_scan_seen: scan_index,
            mz_profile: vec![mz as f32],
            intensity_profile: vec![intensity as f32],
            rt_profile: vec![rt as f32],
            im_profile: vec![im as f32],
            n_real: 1,
        }
    }

    pub fn len(&self) -> usize {
        self.intensity_profile.len()
    }

    pub fn add_point(&mut self, mz: f64, intensity: f64, rt: f64, im: f64, scan_index: usize) {
        // Fill any gaps between last_scan_seen and scan_index
        let gap = scan_index - self.last_scan_seen - 1;
        for _ in 0..gap {
            self.add_gap();
        }

        self.mz_profile.push(mz as f32);
        self.intensity_profile.push(intensity as f32);
        self.rt_profile.push(rt as f32);
        self.im_profile.push(im as f32);
        self.last_scan_seen = scan_index;

        // Update running means (online Welford)
        self.n_real += 1;
        self.running_mz_mean += (mz - self.running_mz_mean) / self.n_real as f64;
        if im != 0.0 {
            self.running_im_mean += (im - self.running_im_mean) / self.n_real as f64;
        }
    }

    pub fn add_gap(&mut self) {
        self.mz_profile.push(f32::NAN);
        self.intensity_profile.push(0.0);
        self.rt_profile.push(f32::NAN);
        self.im_profile.push(f32::NAN);
    }

    pub fn end_scan(&self) -> usize {
        self.scan_start + self.intensity_profile.len() - 1
    }

    pub fn last_real_intensity(&self) -> Option<f64> {
        self.intensity_profile
            .iter()
            .rev()
            .find(|&&x| x > 0.0)
            .map(|&x| x as f64)
    }

    /// Trim leading and trailing zeros from intensity profile (and other profiles).
    pub fn trim(&mut self) {
        let start = self
            .intensity_profile
            .iter()
            .position(|&x| x > 0.0)
            .unwrap_or(0);
        let end = self
            .intensity_profile
            .iter()
            .rposition(|&x| x > 0.0)
            .map(|i| i + 1)
            .unwrap_or(0);

        if start > 0 || end < self.intensity_profile.len() {
            self.scan_start += start;
            self.mz_profile = self.mz_profile[start..end].to_vec();
            self.intensity_profile = self.intensity_profile[start..end].to_vec();
            self.rt_profile = self.rt_profile[start..end].to_vec();
            self.im_profile = self.im_profile[start..end].to_vec();
        }
    }

    /// Trim the hill so that the retained scans account for at least `coverage`
    /// fraction of the total intensity.
    ///
    /// Scans are included by expanding outward from the apex, at each step
    /// choosing the neighbour (left or right) with the higher intensity, until
    /// the accumulated total meets the target.  Leading/trailing zeros left
    /// after this operation are removed by a final `trim()` call.
    pub fn trim_to_coverage(&mut self, coverage: f64) {
        if coverage >= 1.0 {
            return;
        }
        let n = self.intensity_profile.len();
        if n == 0 {
            return;
        }

        let total: f64 = self.intensity_profile.iter().map(|&x| x as f64).sum();
        if total == 0.0 {
            return;
        }
        let target = total * coverage.clamp(0.0, 1.0);

        let apex = self.apex_index();
        let mut left = apex;
        let mut right = apex;
        let mut accumulated = self.intensity_profile[apex] as f64;

        while accumulated < target {
            let left_val = if left > 0 { self.intensity_profile[left - 1] as f64 } else { 0.0 };
            let right_val =
                if right + 1 < n { self.intensity_profile[right + 1] as f64 } else { 0.0 };

            if left_val == 0.0 && right_val == 0.0 {
                // Both neighbours are gaps; still extend to capture any intensity
                // that lies beyond consecutive zeros, but only if there is room.
                if left > 0 { left -= 1; } else if right + 1 < n { right += 1; } else { break; }
            } else if left_val >= right_val && left > 0 {
                left -= 1;
                accumulated += left_val;
            } else if right + 1 < n {
                right += 1;
                accumulated += right_val;
            } else if left > 0 {
                left -= 1;
                accumulated += left_val;
            } else {
                break;
            }
        }

        if left > 0 || right + 1 < n {
            self.scan_start += left;
            self.mz_profile = self.mz_profile[left..=right].to_vec();
            self.intensity_profile = self.intensity_profile[left..=right].to_vec();
            self.rt_profile = self.rt_profile[left..=right].to_vec();
            self.im_profile = self.im_profile[left..=right].to_vec();
        }

        // Remove any zero-intensity prefix/suffix exposed by the slice.
        self.trim();
    }

    // ---- Statistics ----

    /// Intensity-weighted mean m/z over real (non-gap) scans only.
    ///
    /// Both numerator and denominator restrict to positions with non-NaN m/z.
    /// Including gap positions in the denominator (as a naïve `intensity_profile.iter().sum()`
    /// would) breaks once smoothing fills gap intensities, because the numerator
    /// still skips NaN m/z — the resulting mean is biased toward zero.
    pub fn mz_weighted_mean(&self) -> f64 {
        let mut sum_iw = 0.0f64;
        let mut sum_w = 0.0f64;
        for (&mz, &i) in self.mz_profile.iter().zip(&self.intensity_profile) {
            if !mz.is_nan() {
                let w = i as f64;
                sum_iw += mz as f64 * w;
                sum_w += w;
            }
        }
        if sum_w == 0.0 {
            0.0
        } else {
            sum_iw / sum_w
        }
    }

    /// Intensity-weighted std of m/z over real (non-gap) scans only.
    ///
    /// Gap positions are excluded from both the variance sum and the weight
    /// total, so the result is unaffected by gap count or by gap-fill smoothing.
    pub fn mz_weighted_std(&self) -> f64 {
        let mean = self.mz_weighted_mean();
        let mut sum_w = 0.0f64;
        let mut sum_var = 0.0f64;
        for (&mz, &i) in self.mz_profile.iter().zip(&self.intensity_profile) {
            if !mz.is_nan() {
                let w = i as f64;
                let d = mz as f64 - mean;
                sum_var += w * d * d;
                sum_w += w;
            }
        }
        if sum_w == 0.0 {
            0.0
        } else {
            (sum_var / sum_w).sqrt()
        }
    }

    /// Standard error of the intensity-weighted m/z mean, computed via Kish's
    /// effective sample size:
    ///
    /// ```text
    /// neff   = (Σ wᵢ)² / Σ wᵢ²
    /// mz_se  = mz_weighted_std / √neff
    /// ```
    ///
    /// where the sums run over real (non-gap) peaks with `wᵢ = intensityᵢ`.
    /// This is the closed-form analytical analog of the AlphaPept bootstrap
    /// (B=150 resamples of `(mz, intensity)` pairs, std of weighted-mean
    /// estimates). Equivalent for unimodal m/z distributions; cheaper by
    /// 2-3 orders of magnitude. Returns 0.0 when the hill has fewer than
    /// two real peaks or zero total weight (the instrument tolerance floor
    /// alone handles those cases downstream).
    pub fn mz_kish_se(&self) -> f64 {
        let mut sum_w = 0.0f64;
        let mut sum_w2 = 0.0f64;
        let mut n_real = 0usize;
        for (&mz, &i) in self.mz_profile.iter().zip(&self.intensity_profile) {
            if !mz.is_nan() {
                let w = i as f64;
                sum_w += w;
                sum_w2 += w * w;
                n_real += 1;
            }
        }
        if n_real < 2 || sum_w == 0.0 || sum_w2 == 0.0 {
            return 0.0;
        }
        // neff is bounded in [1, n_real]; clamp to ≥ 1 for the rare
        // numerical edge where round-off pushes it slightly below 1.
        let neff = (sum_w * sum_w / sum_w2).max(1.0);
        self.mz_weighted_std() / neff.sqrt()
    }

    /// Intensity-weighted mean ion mobility over real (non-gap) scans only.
    /// See [`mz_weighted_mean`] for why the denominator must filter NaN.
    pub fn im_weighted_mean(&self) -> f64 {
        let mut sum_iw = 0.0f64;
        let mut sum_w = 0.0f64;
        for (&im, &i) in self.im_profile.iter().zip(&self.intensity_profile) {
            if !im.is_nan() {
                let w = i as f64;
                sum_iw += im as f64 * w;
                sum_w += w;
            }
        }
        if sum_w == 0.0 {
            0.0
        } else {
            sum_iw / sum_w
        }
    }

    /// Intensity-weighted std of ion mobility over real (non-gap) scans only.
    pub fn im_weighted_std(&self) -> f64 {
        let mean = self.im_weighted_mean();
        let mut sum_w = 0.0f64;
        let mut sum_var = 0.0f64;
        for (&im, &i) in self.im_profile.iter().zip(&self.intensity_profile) {
            if !im.is_nan() {
                let w = i as f64;
                let d = im as f64 - mean;
                sum_var += w * d * d;
                sum_w += w;
            }
        }
        if sum_w == 0.0 {
            0.0
        } else {
            (sum_var / sum_w).sqrt()
        }
    }

    pub fn min_rt(&self) -> f64 {
        self.rt_profile
            .iter()
            .filter(|x| !x.is_nan())
            .map(|&x| x as f64)
            .fold(f64::INFINITY, f64::min)
    }

    pub fn max_rt(&self) -> f64 {
        self.rt_profile
            .iter()
            .filter(|x| !x.is_nan())
            .map(|&x| x as f64)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    pub fn apex_index(&self) -> usize {
        self.intensity_profile
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    pub fn apex_rt(&self) -> f64 {
        let idx = self.apex_index();
        let v = self.rt_profile[idx];
        if !v.is_nan() {
            return v as f64;
        }
        // Apex landed on a gap scan (NaN RT). This happens when smoothing
        // pushes the apex onto a previously-zero slot. Interpolate from the
        // observed RT bounds using the apex's position within the profile.
        let start = self.min_rt();
        let end = self.max_rt();
        let n = self.intensity_profile.len();
        if !start.is_finite() {
            return 0.0;
        }
        if !end.is_finite() || n <= 1 {
            return start;
        }
        start + (end - start) * (idx as f64) / ((n - 1) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hills::smooth;

    fn make_hill_with_gaps() -> ActiveHill {
        // 6 scans, gaps at indices 1, 2, 4. Real m/z values cluster within ~0.005 Da.
        let mut h = ActiveHill::new(524.270, 100.0, 0.0, 0.9, 0);
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
        let h = ActiveHill::new(500.0, 1000.0, 0.0, 0.0, 0);
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
        let mut h = ActiveHill::new(499.99, 1000.0, 0.0, 0.0, 0);
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
        let mut h = ActiveHill::new(500.000, 100_000.0, 0.0, 0.0, 0);
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
        smooth::smooth_profile(&mut h.intensity_profile, 1);

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
        smooth::smooth_profile(&mut h.intensity_profile, 2);

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
        smooth::smooth_profile(&mut h.intensity_profile, 1);

        let mean = h.im_weighted_mean();
        let std = h.im_weighted_std();

        assert!((mean - 0.905).abs() < 0.01, "im mean drifted: {mean}");
        assert!(std < 0.01, "im_std inflated by gap-fill: {std}");
    }
}
