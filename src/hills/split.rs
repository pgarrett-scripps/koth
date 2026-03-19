use std::sync::Arc;

use crate::models::Hill;

/// Split a hill that contains multiple co-eluting peaks.
///
/// Algorithm:
/// 1. Find local maxima in the intensity profile
/// 2. Find valleys between adjacent peak pairs
/// 3. Split at valleys if valley_intensity <= min(peaks) * min_valley_ratio
///    and both resulting segments have >= min_scans peaks
///
/// Returns the original hill (as a vec of 1) if no split is performed,
/// or multiple hills if splits are made.
pub fn split_hill(
    hill: &Hill,
    min_peak_distance: usize,
    min_peak_height: f64,
    min_scans: usize,
    min_valley_ratio: f64,
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

    // Find local maxima
    let peaks = find_local_maxima(profile, min_height, min_peak_distance);
    if peaks.len() <= 1 {
        return vec![hill.clone()];
    }

    // Find split points (valleys between consecutive peaks)
    let mut split_points: Vec<usize> = Vec::new();
    for window in peaks.windows(2) {
        let left_peak = window[0];
        let right_peak = window[1];

        // Find minimum between the two peaks
        let valley_idx = (left_peak..=right_peak)
            .min_by(|&a, &b| profile[a].partial_cmp(&profile[b]).unwrap())
            .unwrap_or(left_peak);
        let valley_intensity = profile[valley_idx];

        let min_peak_intensity = profile[left_peak].min(profile[right_peak]);
        if valley_intensity <= min_peak_intensity * min_valley_ratio as f32 {
            // Both sides must have enough scans
            let left_size = valley_idx; // from 0 to valley_idx (exclusive)
            let right_size = profile.len() - valley_idx;
            if left_size >= min_scans && right_size >= min_scans {
                split_points.push(valley_idx);
            }
        }
    }

    if split_points.is_empty() {
        return vec![hill.clone()];
    }

    // Build segments
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
            intensity_profile: Arc::from(segment_profile.as_slice()),
        });
    }

    if result.is_empty() {
        vec![hill.clone()]
    } else {
        result
    }
}

/// Find indices of local maxima in profile.
fn find_local_maxima(profile: &[f32], min_height: f32, min_distance: usize) -> Vec<usize> {
    let n = profile.len();
    if n < 3 {
        return Vec::new();
    }

    // Find all candidate local maxima
    let mut candidates: Vec<(usize, f32)> = Vec::new();
    for i in 1..n - 1 {
        if profile[i] >= min_height && profile[i] >= profile[i - 1] && profile[i] >= profile[i + 1] {
            candidates.push((i, profile[i]));
        }
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    // Enforce min_distance: keep highest peak within each distance window
    let mut result: Vec<usize> = Vec::new();
    'outer: for &(idx, height) in &candidates {
        // Check against already selected peaks
        let mut to_remove: Option<usize> = None;
        for (ri, &existing) in result.iter().enumerate() {
            let dist = if idx > existing { idx - existing } else { existing - idx };
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
            result.remove(ri);
        }
        result.push(idx);
    }
    result.sort_unstable();
    result
}

/// Apply co-elution splitting to a list of hills.
pub fn split_coeluting(
    hills: Vec<Hill>,
    min_peak_distance: usize,
    min_peak_height: f64,
    min_scans: usize,
    min_valley_ratio: f64,
) -> Vec<Hill> {
    let mut result = Vec::with_capacity(hills.len());
    let mut split_count = 0usize;

    for hill in &hills {
        let parts = split_hill(hill, min_peak_distance, min_peak_height, min_scans, min_valley_ratio);
        if parts.len() > 1 {
            split_count += 1;
        }
        result.extend(parts);
    }

    log::debug!("Split {} co-eluting hills; total hills now: {}", split_count, result.len());
    result
}
