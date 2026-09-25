use std::sync::Arc;

use super::smooth;
use crate::config::HillsConfig;
use crate::models::Hill;

/// Build sub-hills by cutting `hill.intensity_profile` at the given (sorted,
/// interior) split-point indices. Each segment's statistics are recomputed from
/// the RAW segment intensities. Segments shorter than `min_scans` are dropped;
/// the original hill is returned unchanged if nothing valid remains.
fn build_sub_hills(hill: &Hill, split_points: &[usize], min_scans: usize) -> Vec<Hill> {
    let profile = &hill.intensity_profile;
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
        let intensity_max = segment_profile
            .iter()
            .map(|&x| x as f64)
            .fold(0.0f64, f64::max);
        let skipped = segment_profile.iter().filter(|&&x| x == 0.0).count();
        let hill_score = smooth::compute_hill_score(&segment_profile);

        result.push(Hill {
            hill_id: 0, // re-assigned by `assign_hill_ids` after splitting
            mz: hill.mz,
            mz_std: hill.mz_std,
            // SE inherited from the parent hill. A split sub-hill has fewer
            // contributing scans than the parent, so its true SE would be
            // larger by √(n_parent / n_sub), but the raw peak data isn't
            // retained past finalization. Inheritance is a conservative
            // under-estimate that keeps Kish-aware chain extension functional.
            mz_se: hill.mz_se,
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
            isolation_window: hill.isolation_window,
            faims_cv: hill.faims_cv,
        });
    }

    if result.is_empty() {
        vec![hill.clone()]
    } else {
        result
    }
}

/// Split points = the valley (minimum) index between each adjacent peak pair,
/// retained only when both resulting segments have >= `min_scans` points.
/// Valleys are located on whatever profile `prof` is passed (raw for the
/// prominence splitter, smoothed for the persistence splitter).
fn valleys_between(prof: &[f32], peaks: &[usize], min_scans: usize) -> Vec<usize> {
    let mut split_points = Vec::new();
    for window in peaks.windows(2) {
        let (left_peak, right_peak) = (window[0], window[1]);
        let valley_idx = (left_peak..=right_peak)
            .min_by(|&a, &b| prof[a].partial_cmp(&prof[b]).unwrap())
            .unwrap_or(left_peak);
        if valley_idx >= min_scans && prof.len() - valley_idx >= min_scans {
            split_points.push(valley_idx);
        }
    }
    split_points
}

/// Split a hill containing multiple co-eluting peaks (the "persistence"
/// algorithm). Designed to be robust to noisy scans while separating
/// visually-obvious dual peaks.
///
/// 1. Despike (median-3) and smooth (moving average) a WORKING COPY — all peak
///    finding runs on this; sub-hill intensities still come from the raw trace.
/// 2. Candidate maxima above `height_frac * robust_max` (95th-percentile max,
///    so a single noise spike can't inflate the threshold). No hard min-distance
///    eviction — close real peaks survive.
/// 3. Agglomeratively merge adjacent peaks unless the valley between them is
///    BOTH relatively deep (`valley / smaller_peak <= valley_ratio`) AND
///    absolutely deep (`smaller_peak - valley >= sigma_mult * noise_sigma`).
/// 4. Width gate: drop surviving 1–2 scan spikes median-3 couldn't remove.
/// 5. Split at the (smoothed) valley between adjacent survivors.
pub fn split_hill_persistence(
    hill: &Hill,
    height_frac: f64,
    valley_ratio: f64,
    sigma_mult: f64,
    min_scans: usize,
) -> Vec<Hill> {
    let profile = &hill.intensity_profile;
    if profile.len() < min_scans * 2 {
        return vec![hill.clone()];
    }

    // Working copy: nonlinear despike, then linear smooth (half-window 2).
    let despiked = median3(profile);
    let mut s = despiked.clone();
    smooth::running_average(&mut s, 2);

    let sigma = noise_sigma(&despiked).max(1e-6);
    let floor = sigma_mult as f32 * sigma;
    let robust_max = percentile(&s, 0.95).max(1e-6);
    let min_h = (height_frac as f32 * robust_max).max(floor);

    let mut peaks = local_maxima(&s, min_h); // deliberately NO spacing eviction
    if peaks.len() <= 1 {
        return vec![hill.clone()];
    }

    // Agglomeratively drop the single worst-separated peak until every
    // surviving notch passes both the relative and absolute depth tests.
    let ratio = valley_ratio as f32;
    loop {
        if peaks.len() <= 1 {
            break;
        }
        let mut merge_at: Option<(usize, f32)> = None; // (index, badness)
        for (idx, w) in peaks.windows(2).enumerate() {
            let v = valley_min(&s, w[0], w[1]);
            let smaller = s[w[0]].min(s[w[1]]).max(1e-6);
            let r = v / smaller; // higher = shallower notch
            let depth = smaller - v;
            if r > ratio || depth < floor {
                let badness = (r - ratio).max(0.0) + (floor - depth).max(0.0) / floor;
                if merge_at.is_none_or(|(_, b)| badness > b) {
                    merge_at = Some((idx, badness));
                }
            }
        }
        let Some((idx, _)) = merge_at else { break };
        let (a, b) = (peaks[idx], peaks[idx + 1]);
        let drop = if s[a] <= s[b] { idx } else { idx + 1 };
        peaks.remove(drop);
    }

    // Width gate: a real chromatographic peak is several scans wide at half-max.
    const MIN_WIDTH: usize = 4;
    peaks.retain(|&pk| {
        let half = 0.5 * s[pk];
        let lo = pk.saturating_sub(8);
        let hi = (pk + 9).min(s.len());
        (lo..hi).filter(|&i| s[i] >= half).count() >= MIN_WIDTH
    });
    if peaks.len() <= 1 {
        return vec![hill.clone()];
    }

    let split_points = valleys_between(&s, &peaks, min_scans);
    build_sub_hills(hill, &split_points, min_scans)
}

// ---- Persistence-splitter primitives ----

/// Width-3 median filter (nonlinear despike). Removes lone 1-sample spikes
/// completely while preserving genuine peak edges — something no linear
/// (moving-average) smoother can do. Endpoints are left untouched.
fn median3(p: &[f32]) -> Vec<f32> {
    let n = p.len();
    if n < 3 {
        return p.to_vec();
    }
    let mut out = p.to_vec();
    for i in 1..n - 1 {
        let mut t = [p[i - 1], p[i], p[i + 1]];
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out[i] = t[1];
    }
    out
}

/// `q`-quantile (0..=1) of a copy of `p`. Used for a robust max that a single
/// noise spike cannot inflate.
fn percentile(p: &[f32], q: f32) -> f32 {
    if p.is_empty() {
        return 0.0;
    }
    let mut s = p.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((s.len() as f32 - 1.0) * q).round() as usize;
    s[idx]
}

/// Robust high-frequency noise σ via the median absolute first-difference.
/// Peaks contribute a few large diffs; the median ignores them. The √2 divisor
/// undoes the variance inflation from differencing; 1.4826 converts MAD→σ.
fn noise_sigma(p: &[f32]) -> f32 {
    if p.len() < 3 {
        return 0.0;
    }
    let mut d: Vec<f32> = p.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = d[d.len() / 2];
    med * 1.4826 / std::f32::consts::SQRT_2
}

/// Indices that are >= both neighbours and above `min_height`.
fn local_maxima(p: &[f32], min_height: f32) -> Vec<usize> {
    let n = p.len();
    let mut out = Vec::new();
    for i in 0..n {
        let l = i == 0 || p[i] >= p[i - 1];
        let r = i == n - 1 || p[i] >= p[i + 1];
        if p[i] >= min_height && l && r {
            out.push(i);
        }
    }
    out
}

/// Lowest value in the inclusive range `[a, b]`.
fn valley_min(p: &[f32], a: usize, b: usize) -> f32 {
    p[a..=b].iter().cloned().fold(f32::INFINITY, f32::min)
}

/// Apply co-elution splitting to a list of hills using the persistence
/// splitter (despike + smooth + robust-max threshold + valley merge).
pub fn split_coeluting(hills: Vec<Hill>, config: &HillsConfig) -> Vec<Hill> {
    let mut result = Vec::with_capacity(hills.len());
    let mut split_count = 0usize;

    for hill in &hills {
        let parts = split_hill_persistence(
            hill,
            config.split_height_frac,
            config.split_valley_ratio,
            config.split_sigma_mult,
            config.min_scans,
        );
        if parts.len() > 1 {
            split_count += 1;
        }
        result.extend(parts);
    }

    log::debug!(
        "Split {} co-eluting hills (persistence algo); total hills now: {}",
        split_count,
        result.len()
    );
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hill_from_profile(profile: Vec<f32>) -> Hill {
        let n = profile.len();
        let max = profile.iter().cloned().fold(0.0f32, f32::max) as f64;
        Hill {
            hill_id: 0,
            mz: 500.0,
            mz_std: 0.001,
            mz_se: 0.0005,
            rt: 0.0,
            rt_start: 0.0,
            rt_end: (n - 1) as f64,
            rt_width: (n - 1) as f64,
            im: 0.0,
            im_std: 0.0,
            scan_start: 0,
            scan_apex: 0,
            scan_end: n - 1,
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: profile.iter().map(|&x| x as f64).sum(),
            intensity_max: max,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
            faims_cv: None,
        }
    }

    fn gaussian(n: usize, center: f32, sigma: f32, amp: f32) -> Vec<f32> {
        (0..n)
            .map(|i| {
                let d = (i as f32 - center) / sigma;
                amp * (-0.5 * d * d).exp()
            })
            .collect()
    }

    #[test]
    fn median3_kills_lone_spike() {
        let p = vec![10.0f32, 9.0, 500.0, 11.0, 10.0];
        let m = median3(&p);
        assert!(
            m[2] < 12.0,
            "spike should be replaced by a neighbour: {m:?}"
        );
    }

    #[test]
    fn persistence_splits_clean_dual() {
        let mut p = gaussian(100, 35.0, 5.0, 1000.0);
        for (i, v) in gaussian(100, 65.0, 5.0, 1000.0).into_iter().enumerate() {
            p[i] += v;
        }
        let hill = hill_from_profile(p);
        let parts = split_hill_persistence(&hill, 0.10, 0.70, 4.0, 3);
        assert_eq!(parts.len(), 2, "clean dual should split into 2");
    }

    #[test]
    fn persistence_keeps_single_peak() {
        let p = gaussian(80, 40.0, 6.0, 1000.0);
        let hill = hill_from_profile(p);
        let parts = split_hill_persistence(&hill, 0.10, 0.70, 4.0, 3);
        assert_eq!(parts.len(), 1, "single peak must not split");
    }

    #[test]
    fn persistence_ignores_lone_spike() {
        // One real peak plus a 1-sample spike taller than it: must stay 1 hill.
        let mut p = gaussian(80, 40.0, 6.0, 1000.0);
        p[12] += 1800.0;
        let hill = hill_from_profile(p);
        let parts = split_hill_persistence(&hill, 0.10, 0.70, 4.0, 3);
        assert_eq!(parts.len(), 1, "lone spike must not create a second hill");
    }
}
