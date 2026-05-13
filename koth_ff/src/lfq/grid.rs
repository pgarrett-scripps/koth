use crate::{alignment::RunInput, models::Hill};

use super::LfqConfig;

pub const C13_NEUTRON: f64 = 1.003_354_835;

/// An extracted-ion-chromatogram grid for a single feature in a single run.
///
/// Layout: `intensities[isotopologue][rt_bin]`
/// - Rows 0..n_isotopes correspond to M, M+1, M+2, ...
/// - Columns 0..grid_cols span [rt_min, rt_max] in equal-width bins.
#[derive(Debug, Clone)]
pub struct XicGrid {
    /// `[isotopologue][rt_bin]` intensities (summed from matching hills)
    pub intensities: Vec<Vec<f32>>,
    pub rt_min: f64,
    pub rt_max: f64,
    /// How many of the n_isotopes rows have any non-zero signal
    pub n_slots_filled: u8,
}

impl XicGrid {
    pub fn empty(n_rows: usize, n_cols: usize, rt_min: f64, rt_max: f64) -> Self {
        Self {
            intensities: vec![vec![0.0f32; n_cols]; n_rows],
            rt_min,
            rt_max,
            n_slots_filled: 0,
        }
    }

    /// Zero all intensity rows and update RT bounds for reuse across (feature, run) pairs.
    pub fn reset(&mut self, rt_min: f64, rt_max: f64) {
        for row in &mut self.intensities {
            row.fill(0.0);
        }
        self.rt_min = rt_min;
        self.rt_max = rt_max;
        self.n_slots_filled = 0;
    }

    pub fn n_rows(&self) -> usize {
        self.intensities.len()
    }

    pub fn n_cols(&self) -> usize {
        self.intensities.first().map_or(0, |r| r.len())
    }

    pub fn is_empty(&self) -> bool {
        self.n_slots_filled == 0
    }
}

/// Build an XIC grid for a single feature in a single run.
///
/// `grid` is reset (zeroed) at entry and filled in-place; callers should
/// pre-allocate it once and pass the same instance across (feature, run) pairs.
/// `target_mz` / `target_rt` / `target_im` are already expressed in the
/// reference coordinate space (alignment-corrected).  The function searches
/// `run.hills` (via the pre-sorted `sorted_hill_idx`) for hills matching each
/// isotopologue and fills the grid from their intensity profiles.
pub fn build_grid(
    grid: &mut XicGrid,
    run: &RunInput,
    sorted_hill_idx: &[usize],
    target_mz: f64,
    charge: u8,
    target_rt: f64,
    target_im: f64,
    half_window: f64,
    config: &LfqConfig,
) {
    let n_rows = config.n_isotopes;
    let n_cols = config.grid_cols;

    let rt_min = target_rt - half_window;
    let rt_max = target_rt + half_window;

    grid.reset(rt_min, rt_max);
    let has_im = target_im != 0.0;

    for iso in 0..n_rows {
        let expected_mz = target_mz + iso as f64 * C13_NEUTRON / charge as f64;
        let mz_tol = expected_mz * config.mz_ppm / 1e6;

        // Binary search for the mz window within the sorted hill index
        let lo = sorted_hill_idx
            .partition_point(|&i| run.hills[i].mz < expected_mz - mz_tol);
        let hi = sorted_hill_idx
            .partition_point(|&i| run.hills[i].mz <= expected_mz + mz_tol);

        // Among all candidates in the mz band, pick the one whose apex is
        // closest to the target RT and whose RT range overlaps the window.
        let best = sorted_hill_idx[lo..hi]
            .iter()
            .filter_map(|&hill_idx| {
                let hill = &run.hills[hill_idx];
                if hill.rt_end < rt_min || hill.rt_start > rt_max {
                    return None;
                }
                if has_im && hill.im != 0.0 && (hill.im - target_im).abs() > config.im_tolerance {
                    return None;
                }
                let rt_dist = (hill.rt - target_rt).abs();
                Some((hill_idx, rt_dist))
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        if let Some((hill_idx, _)) = best {
            let hill = &run.hills[hill_idx];
            let added =
                fill_row_from_hill(hill, &mut grid.intensities[iso], rt_min, rt_max, n_cols, &run.scan_times);
            if added {
                grid.n_slots_filled += 1;
            }
        }
    }
}

/// Splat a hill's per-scan intensities into a single grid row.
/// Returns true if at least one sample fell within the window.
fn fill_row_from_hill(
    hill: &Hill,
    row: &mut Vec<f32>,
    rt_min: f64,
    rt_max: f64,
    n_cols: usize,
    scan_times: &[f64],
) -> bool {
    let rt_span = rt_max - rt_min;
    if rt_span <= 0.0 || hill.intensity_profile.is_empty() {
        return false;
    }

    let mut any = false;

    for (i, &intensity) in hill.intensity_profile.iter().enumerate() {
        if intensity <= 0.0 {
            continue;
        }

        let rt = scan_rt(hill, i, scan_times);
        let t = (rt - rt_min) / rt_span;
        if t < 0.0 || t >= 1.0 {
            continue;
        }

        let col = (t * n_cols as f64) as usize;
        let col = col.min(n_cols - 1);
        row[col] += intensity;
        any = true;
    }

    any
}

/// Return the retention time for the i-th sample of a hill's intensity profile.
/// Uses the run's scan_times array when available; falls back to linear
/// interpolation between hill.rt_start and hill.rt_end otherwise.
fn scan_rt(hill: &Hill, profile_index: usize, scan_times: &[f64]) -> f64 {
    if !scan_times.is_empty() {
        let scan_abs = hill.scan_start + profile_index;
        if scan_abs < scan_times.len() {
            return scan_times[scan_abs];
        }
    }
    // Linear interpolation fallback
    let n = hill.intensity_profile.len();
    if n <= 1 {
        return hill.rt;
    }
    hill.rt_start + (hill.rt_end - hill.rt_start) * profile_index as f64 / (n - 1) as f64
}
