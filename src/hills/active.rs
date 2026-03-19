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

    pub fn mz_weighted_mean(&self) -> f64 {
        let total: f64 = self.intensity_profile.iter().map(|&x| x as f64).sum();
        if total == 0.0 {
            return 0.0;
        }
        self.mz_profile
            .iter()
            .zip(&self.intensity_profile)
            .filter(|(mz, _)| !mz.is_nan())
            .map(|(&mz, &i)| mz as f64 * i as f64)
            .sum::<f64>()
            / total
    }

    pub fn mz_weighted_std(&self) -> f64 {
        let mean = self.mz_weighted_mean();
        let total: f64 = self.intensity_profile.iter().map(|&x| x as f64).sum();
        if total == 0.0 {
            return 0.0;
        }
        let var: f64 = self
            .mz_profile
            .iter()
            .zip(&self.intensity_profile)
            .filter(|(mz, _)| !mz.is_nan())
            .map(|(&mz, &i)| i as f64 * (mz as f64 - mean).powi(2))
            .sum::<f64>()
            / total;
        var.sqrt()
    }

    pub fn im_weighted_mean(&self) -> f64 {
        let total: f64 = self.intensity_profile.iter().map(|&x| x as f64).sum();
        if total == 0.0 {
            return 0.0;
        }
        self.im_profile
            .iter()
            .zip(&self.intensity_profile)
            .filter(|(im, _)| !im.is_nan())
            .map(|(&im, &i)| im as f64 * i as f64)
            .sum::<f64>()
            / total
    }

    pub fn im_weighted_std(&self) -> f64 {
        let mean = self.im_weighted_mean();
        let total: f64 = self.intensity_profile.iter().map(|&x| x as f64).sum();
        if total == 0.0 {
            return 0.0;
        }
        let var: f64 = self
            .im_profile
            .iter()
            .zip(&self.intensity_profile)
            .filter(|(im, _)| !im.is_nan())
            .map(|(&im, &i)| i as f64 * (im as f64 - mean).powi(2))
            .sum::<f64>()
            / total;
        var.sqrt()
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
        if v.is_nan() { 0.0 } else { v as f64 }
    }
}
