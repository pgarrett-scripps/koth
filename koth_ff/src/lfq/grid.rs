use crate::models::Hill;

use super::LfqConfig;

#[cfg(test)]
#[path = "grid_tests.rs"]
mod tests;

pub const C13_NEUTRON: f64 = 1.003_354_835;

/// Cache-friendly, mz-sorted lookup record for `build_grid`.
///
/// Holds only the fields the hot loop reads (mz / rt / rt_start / rt_end / im)
/// plus an index back into `run.hills` for the eventual `fill_row_from_hill`
/// call. Packed to 48 bytes so each candidate access during the in-band scan
/// is a single cache-line fetch — versus the 200-byte `Hill` struct, where
/// the same five fields are scattered across multiple cache lines and 75% of
/// each loaded line is unused.
#[derive(Clone, Copy)]
#[repr(C)]
pub struct HillKey {
    pub mz: f64,
    pub rt: f64,
    pub rt_start: f64,
    pub rt_end: f64,
    pub im: f64,
    /// Index back into `run.hills` — only used for `fill_row_from_hill`.
    pub hill_idx: u32,
    _pad: u32,
}

/// `Vec<HillKey>` sorted ascending by mz, plus a stable name for the SoA-style
/// callsite signature.
pub struct SortedHills {
    pub keys: Vec<HillKey>,
}

impl SortedHills {
    /// Build the sorted view from a run's hills. Records are sorted ascending by mz.
    pub fn from_hills(hills: &[Hill]) -> Self {
        let mut keys: Vec<HillKey> = hills
            .iter()
            .enumerate()
            .map(|(i, h)| HillKey {
                mz: h.mz,
                rt: h.rt,
                rt_start: h.rt_start,
                rt_end: h.rt_end,
                im: h.im,
                hill_idx: i as u32,
                _pad: 0,
            })
            .collect();
        keys.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
        SortedHills { keys }
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }
}

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
    /// Centroid m/z of the winning hill per isotopologue row.
    /// NaN if no hill matched that row. Length = n_rows.
    pub winner_mz: Vec<f64>,
    /// Ion mobility of the winning hill per isotopologue row.
    /// NaN if no hill matched. Length = n_rows.
    pub winner_im: Vec<f64>,
}

impl XicGrid {
    pub fn empty(n_rows: usize, n_cols: usize, rt_min: f64, rt_max: f64) -> Self {
        Self {
            intensities: vec![vec![0.0f32; n_cols]; n_rows],
            rt_min,
            rt_max,
            n_slots_filled: 0,
            winner_mz: vec![f64::NAN; n_rows],
            winner_im: vec![f64::NAN; n_rows],
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
        for v in &mut self.winner_mz {
            *v = f64::NAN;
        }
        for v in &mut self.winner_im {
            *v = f64::NAN;
        }
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
/// `hills` (a mz-sorted `SortedHills` SoA view of the run) is binary-searched
/// for hills matching each isotopologue; the winner is then materialised
/// from `hills_data` (the run's full `Vec<Hill>`) to read its `intensity_profile`.
pub fn build_grid(
    grid: &mut XicGrid,
    hills_data: &[Hill],
    hills: &SortedHills,
    scan_times: &[f64],
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
        let mz_lo = expected_mz - mz_tol;
        let mz_hi = expected_mz + mz_tol;

        // Binary search over the packed mz field (stride-48 reads, but the
        // candidate slice below is contiguous).
        let lo = hills.keys.partition_point(|k| k.mz < mz_lo);
        let hi = hills.keys.partition_point(|k| k.mz <= mz_hi);

        // Sum *every* in-box hill into this isotope row. Each hill emits
        // per-scan intensities (its `intensity_profile`) which are binned
        // into the grid's RT columns; calling `fill_row_from_hill` more
        // than once accumulates because the inner loop does `row[col] +=`.
        // The net effect is a reconstruction of the raw MS1 XIC at this
        // m/z across the grid's RT window, built from the hill database
        // without going back to the mzML.
        //
        // The "winner" hill (highest `intensity_sum` among contributors)
        // is tracked only for the diagnostic `observed_mz` / `observed_im`
        // fields written into LfqEntry — it does not change which signal
        // is integrated.
        let mut any_added = false;
        let mut winner_intensity: f32 = 0.0;
        for k_idx in lo..hi {
            let k = &hills.keys[k_idx];
            if k.rt_end < rt_min || k.rt_start > rt_max {
                continue;
            }
            if has_im && k.im != 0.0 && (k.im - target_im).abs() > config.im_tolerance {
                continue;
            }
            let hill = &hills_data[k.hill_idx as usize];
            let added = fill_row_from_hill(
                hill,
                &mut grid.intensities[iso],
                rt_min,
                rt_max,
                n_cols,
                scan_times,
            );
            if added {
                any_added = true;
                let intensity_sum = hill.intensity_sum as f32;
                if intensity_sum > winner_intensity {
                    winner_intensity = intensity_sum;
                    grid.winner_mz[iso] = hill.mz;
                    grid.winner_im[iso] = hill.im;
                }
            }
        }
        if any_added {
            grid.n_slots_filled += 1;
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
