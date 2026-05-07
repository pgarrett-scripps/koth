use std::sync::Arc;

use crate::models::Hill;
use super::smooth;

/// Split a hill that contains multiple co-eluting peaks.
///
/// Algorithm:
/// 1. Find local maxima above `min_peak_height` with `min_peak_distance` spacing.
/// 2. Compute the prominence of each candidate peak and discard those below
///    `min_prominence * max_intensity`.  Prominence = peak height minus the
///    highest valley between the peak and the nearest taller peak on either side
///    (or the array boundary).  This is noise-robust: a wiggle on the flank of a
///    large peak has near-zero prominence regardless of its local neighbourhood.
/// 3. Split at the valley (minimum) between adjacent surviving peaks, provided
///    both segments have >= min_scans valid data points.
///
/// Returns the original hill unchanged if no split passes all criteria.
pub fn split_hill(
    hill: &Hill,
    min_peak_distance: usize,
    min_peak_height: f64,
    min_scans: usize,
    min_prominence: f64,
) -> Vec<Hill> {
    let profile = &hill.intensity_profile;
    if profile.len() < min_scans * 2 {
        return vec![hill.clone()];
    }

    let max_intensity = profile.iter().cloned().fold(0.0f32, f32::max);
    if max_intensity == 0.0 {
        return vec![hill.clone()];
    }

    let min_height = max_intensity * min_peak_height as f32;
    let prom_threshold = max_intensity * min_prominence as f32;

    let peaks = find_prominent_peaks(profile, min_height, min_peak_distance, prom_threshold);
    if peaks.len() <= 1 {
        return vec![hill.clone()];
    }

    // Split at the valley (minimum) between each pair of adjacent prominent peaks.
    let mut split_points: Vec<usize> = Vec::new();
    for window in peaks.windows(2) {
        let (left_peak, right_peak) = (window[0], window[1]);
        let valley_idx = (left_peak..=right_peak)
            .min_by(|&a, &b| profile[a].partial_cmp(&profile[b]).unwrap())
            .unwrap_or(left_peak);

        let left_size = valley_idx;
        let right_size = profile.len() - valley_idx;
        if left_size >= min_scans && right_size >= min_scans {
            split_points.push(valley_idx);
        }
    }

    if split_points.is_empty() {
        return vec![hill.clone()];
    }

    let boundaries: Vec<usize> = std::iter::once(0)
        .chain(split_points.iter().copied())
        .chain(std::iter::once(profile.len()))
        .collect();

    let rt_per_scan = if hill.scan_end > hill.scan_start {
        (hill.rt_end - hill.rt_start) / (hill.scan_end - hill.scan_start) as f64
    } else {
        0.0
    };

    let mut result = Vec::new();
    for window in boundaries.windows(2) {
        let seg_start = window[0];
        let seg_end = window[1];
        let segment_profile: Vec<f32> = profile[seg_start..seg_end].to_vec();

        if segment_profile.len() < min_scans {
            continue;
        }

        let apex_relative = segment_profile
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        let abs_start = hill.scan_start + seg_start;
        let abs_end = hill.scan_start + seg_end - 1;
        let abs_apex = abs_start + apex_relative;

        let rt_start = hill.rt_start + (abs_start - hill.scan_start) as f64 * rt_per_scan;
        let rt_end = hill.rt_start + (abs_end - hill.scan_start) as f64 * rt_per_scan;
        let rt_apex = hill.rt_start + (abs_apex - hill.scan_start) as f64 * rt_per_scan;

        let intensity_sum: f64 = segment_profile.iter().map(|&x| x as f64).sum();
        let intensity_max = segment_profile.iter().map(|&x| x as f64).fold(0.0f64, f64::max);
        let skipped = segment_profile.iter().filter(|&&x| x == 0.0).count();
        let hill_score = smooth::compute_hill_score(&segment_profile);

        result.push(Hill {
            mz: hill.mz,
            mz_std: hill.mz_std,
            rt: rt_apex,
            rt_start,
            rt_end,
            rt_width: rt_end - rt_start,
            im: hill.im,
            im_std: hill.im_std,
            scan_start: abs_start,
            scan_apex: abs_apex,
            scan_end: abs_end,
            n_scans: segment_profile.len(),
            skipped_scans: skipped,
            intensity_sum,
            intensity_max,
            hill_score,
            intensity_profile: Arc::from(segment_profile.as_slice()),
        });
    }

    if result.is_empty() { vec![hill.clone()] } else { result }
}

/// Compute the prominence of a peak at `peak_idx`.
///
/// Walks left and right from the peak until a taller value is encountered (or
/// the array boundary), recording the minimum along each path.  Prominence is
/// the peak height minus the higher of the two path minima.
///
/// A noise wiggle on the flank of a large peak will have a high path-minimum on
/// the side facing the large peak (because it never descends far before hitting
/// taller terrain), so its prominence is near zero.
fn compute_prominence(profile: &[f32], peak_idx: usize) -> f32 {
    let h = profile[peak_idx];

    let mut left_min = h;
    for i in (0..peak_idx).rev() {
        if profile[i] > h {
            break;
        }
        if profile[i] < left_min {
            left_min = profile[i];
        }
    }

    let mut right_min = h;
    for i in (peak_idx + 1)..profile.len() {
        if profile[i] > h {
            break;
        }
        if profile[i] < right_min {
            right_min = profile[i];
        }
    }

    h - left_min.max(right_min)
}

/// Find local maxima above `min_height` with `min_distance` spacing, then
/// filter to those whose prominence is >= `min_prominence`.
fn find_prominent_peaks(
    profile: &[f32],
    min_height: f32,
    min_distance: usize,
    min_prominence: f32,
) -> Vec<usize> {
    let n = profile.len();
    if n < 3 {
        return Vec::new();
    }

    // Candidate local maxima: >= both neighbours, above min_height.
    let mut candidates: Vec<(usize, f32)> = Vec::new();
    for i in 0..n {
        let left_ok  = i == 0     || profile[i] >= profile[i - 1];
        let right_ok = i == n - 1 || profile[i] >= profile[i + 1];
        if profile[i] >= min_height && left_ok && right_ok {
            candidates.push((i, profile[i]));
        }
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    // Enforce min_distance: within each window keep the tallest.
    let mut spaced: Vec<usize> = Vec::new();
    'outer: for &(idx, height) in &candidates {
        let mut to_remove: Option<usize> = None;
        for (ri, &existing) in spaced.iter().enumerate() {
            let dist = idx.abs_diff(existing);
            if dist < min_distance {
                if height <= profile[existing] {
                    continue 'outer;
                } else {
                    to_remove = Some(ri);
                    break;
                }
            }
        }
        if let Some(ri) = to_remove {
            spaced.remove(ri);
        }
        spaced.push(idx);
    }
    spaced.sort_unstable();

    // Prominence filter.
    spaced
        .into_iter()
        .filter(|&idx| compute_prominence(profile, idx) >= min_prominence)
        .collect()
}

/// Apply co-elution splitting to a list of hills.
pub fn split_coeluting(
    hills: Vec<Hill>,
    min_peak_distance: usize,
    min_peak_height: f64,
    min_scans: usize,
    min_prominence: f64,
) -> Vec<Hill> {
    let mut result = Vec::with_capacity(hills.len());
    let mut split_count = 0usize;

    for hill in &hills {
        let parts = split_hill(hill, min_peak_distance, min_peak_height, min_scans, min_prominence);
        if parts.len() > 1 {
            split_count += 1;
        }
        result.extend(parts);
    }

    log::debug!("Split {} co-eluting hills; total hills now: {}", split_count, result.len());
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prominence_isolated_peak() {
        // A lone peak rising from zero: prominence = peak height.
        let p = vec![0.0f32, 0.0, 100.0, 0.0, 0.0];
        assert!((compute_prominence(&p, 2) - 100.0).abs() < 1e-4);
    }

    #[test]
    fn prominence_shoulder_peak() {
        // Small bump on the flank of a large peak: path toward the large peak
        // never descends below the bump, so prominence is small.
        // Profile: [0, 50, 80, 60, 100, 0]  — bump at 2, large peak at 4.
        let p = vec![0.0f32, 50.0, 80.0, 60.0, 100.0, 0.0];
        let prom = compute_prominence(&p, 2); // peak=80, right path min=60 before hitting 100
        // left_min = 0 (walks to edge), right_min = 60 (stops at 100 > 80)
        // reference = max(0, 60) = 60, prominence = 80 - 60 = 20
        assert!((prom - 20.0).abs() < 1e-4, "got {prom}");
    }

    #[test]
    fn prominence_two_equal_peaks() {
        // Two equal peaks with a deep valley: both should have high prominence.
        let p = vec![0.0f32, 100.0, 0.0, 20.0, 0.0, 100.0, 0.0];
        let p0 = compute_prominence(&p, 1); // left path min=0, right path min=0 (stops at equal peak)
        let p1 = compute_prominence(&p, 5);
        // For peak at 1: walks right, hits 100 at idx 5 (equal, not strictly greater — keeps walking)
        // Actually profile[5]=100 is NOT > profile[1]=100, so we keep walking to edge. right_min=0.
        // prominence = 100 - max(0, 0) = 100
        assert!(p0 > 50.0, "p0={p0}");
        assert!(p1 > 50.0, "p1={p1}");
    }
}
