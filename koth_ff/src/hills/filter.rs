//! Drop "large baseline-like" hills before isotope chain assembly.
//!
//! Ported from AlphaPept's `filter_hills` / `check_large_hills`
//! ([feature_finding.py:520-599](https://github.com/MannLabs/alphapept)).
//!
//! Rationale: koth's hill detector clusters every chromatographically-
//! coherent m/z trace into a hill, regardless of shape. That captures real
//! peptide elution profiles, but it also captures:
//!
//! - **Column bleed / siloxanes / plasticizers** — contaminants that elute
//!   weakly across most of the LC gradient at consistent m/z.
//! - **Solvent ion adducts** — present in nearly every scan.
//! - **Baseline drift** — slowly-varying background the centroider
//!   clusters as a single hill because m/z stays inside tolerance.
//!
//! These are real "hills" in the m/z sense but not chromatographic peaks
//! of analytes. They pollute downstream chain assembly (every contaminant
//! pair within a neutron offset can produce a spurious feature).
//!
//! The filter keeps a large hill only when its unsmoothed intensity
//! maximum is at least `peak_factor` (default 2.0) times the smoothed
//! intensity at *both* endpoints. A clear chromatographic peak passes
//! easily; a flat baseline trace or a slow ramp fails.

use crate::models::Hill;

/// Smoothing half-window used internally by the endpoint median+mean
/// filter. AlphaPept uses 1 (= 3-point neighbourhoods); we match it.
const SMOOTH_WINDOW: usize = 1;

/// Drop hills that span at least `min_scans` scans but don't show a clear
/// chromatographic peak shape — the unsmoothed intensity maximum has to
/// exceed the smoothed endpoint intensity by at least `peak_factor` on
/// both sides. Small hills (< min_scans) pass through untouched.
pub fn filter_large_baseline_hills(
    hills: Vec<Hill>,
    min_scans: usize,
    peak_factor: f64,
) -> Vec<Hill> {
    let n_before = hills.len();
    let filtered: Vec<Hill> = hills
        .into_iter()
        .filter(|h| keep_hill(h, min_scans, peak_factor))
        .collect();
    let n_after = filtered.len();
    log::info!(
        "filter_large_baseline_hills: kept {}/{} hills ({} dropped as baseline-like)",
        n_after,
        n_before,
        n_before - n_after,
    );
    filtered
}

/// Returns `true` if the hill survives the filter.
///
/// Conservative defaults: hills shorter than `min_scans`, hills with a
/// zero smoothed endpoint (degenerate; should never happen for real
/// signal but possible on edge cases), and hills with a non-positive max
/// are all kept.
fn keep_hill(hill: &Hill, min_scans: usize, peak_factor: f64) -> bool {
    let profile = &hill.intensity_profile;
    let n = profile.len();
    if n < min_scans {
        return true;
    }

    let smoothed = median_then_mean_smooth(profile, SMOOTH_WINDOW);
    let max_unsmoothed: f32 = profile.iter().copied().fold(0.0f32, f32::max);
    if max_unsmoothed <= 0.0 {
        return true;
    }

    let edge_start = smoothed[0];
    let edge_end = smoothed[n - 1];
    if edge_start <= 0.0 || edge_end <= 0.0 {
        return true;
    }

    let ratio_start = max_unsmoothed as f64 / edge_start as f64;
    let ratio_end = max_unsmoothed as f64 / edge_end as f64;

    ratio_start > peak_factor && ratio_end > peak_factor
}

/// Sequential in-place median filter followed by sequential in-place mean
/// filter, both with half-window `window`. Matches AlphaPept's behaviour
/// at [feature_finding.py:544-552](https://github.com/MannLabs/alphapept):
/// each output cell is written before the next cell reads, so a later
/// cell sees the already-smoothed earlier cells. The result differs from
/// a parallel filter in detail but the qualitative endpoint behaviour
/// (suppress isolated spikes at the trace ends) is what the peak-shape
/// test relies on.
fn median_then_mean_smooth(profile: &[f32], window: usize) -> Vec<f32> {
    let n = profile.len();
    let mut buf: Vec<f32> = profile.to_vec();
    if n == 0 {
        return buf;
    }
    let mut scratch: Vec<f32> = Vec::with_capacity(2 * window + 1);

    for i in 0..n {
        let lo = i.saturating_sub(window);
        let hi = (i + window + 1).min(n);
        scratch.clear();
        scratch.extend_from_slice(&buf[lo..hi]);
        scratch.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let m = scratch.len();
        buf[i] = if m % 2 == 1 {
            scratch[m / 2]
        } else {
            0.5 * (scratch[m / 2 - 1] + scratch[m / 2])
        };
    }

    for i in 0..n {
        let lo = i.saturating_sub(window);
        let hi = (i + window + 1).min(n);
        let sum: f32 = buf[lo..hi].iter().sum();
        buf[i] = sum / (hi - lo) as f32;
    }

    buf
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::Hill;
    use std::sync::Arc;

    fn make_hill(profile: Vec<f32>) -> Hill {
        let n = profile.len();
        let int_sum: f64 = profile.iter().map(|&x| x as f64).sum();
        let int_max: f64 = profile.iter().copied().fold(0.0f32, f32::max) as f64;
        Hill {
            hill_id: 0,
            mz: 500.0,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: 0.5,
            rt_start: 0.0,
            rt_end: 1.0,
            rt_width: 1.0,
            im: 0.0,
            im_std: 0.0,
            scan_start: 0,
            scan_apex: n / 2,
            scan_end: n.saturating_sub(1),
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: int_sum,
            intensity_max: int_max,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
        }
    }

    /// Small hills (below threshold) are always kept regardless of shape.
    #[test]
    fn small_hill_always_kept() {
        let flat = make_hill(vec![1000.0; 20]); // 20 < 40
        assert!(keep_hill(&flat, 40, 2.0));
    }

    /// A large hill with a clear chromatographic apex (max in the middle,
    /// low intensity at endpoints) is kept.
    #[test]
    fn large_peaked_hill_kept() {
        // 50-scan Gaussian-ish profile peaking at center 25.
        let n = 50;
        let profile: Vec<f32> = (0..n)
            .map(|i| {
                let d = (i as f32 - 25.0) / 8.0;
                10_000.0 * (-d * d).exp() + 50.0 // small baseline at edges
            })
            .collect();
        let h = make_hill(profile);
        assert!(keep_hill(&h, 40, 2.0));
    }

    /// A large flat hill (constant intensity across all scans) — max ≈
    /// endpoints — is dropped. This is the "column-bleed contaminant"
    /// case the filter is designed to catch.
    #[test]
    fn large_flat_hill_dropped() {
        let flat = make_hill(vec![1000.0; 50]);
        assert!(!keep_hill(&flat, 40, 2.0));
    }

    /// A large hill with a monotonic ramp shape (max at one end, low at
    /// the other) is dropped — the apex-ratio test requires both
    /// endpoints to be low.
    #[test]
    fn large_ramp_hill_dropped() {
        let n = 50;
        let profile: Vec<f32> = (0..n).map(|i| 100.0 + i as f32 * 200.0).collect();
        let h = make_hill(profile);
        assert!(!keep_hill(&h, 40, 2.0));
    }

    /// Boundary check: a hill with exactly `min_scans` scans is evaluated
    /// (not skipped as "small").
    #[test]
    fn boundary_at_min_scans_is_evaluated() {
        let flat = make_hill(vec![1000.0; 40]);
        assert!(!keep_hill(&flat, 40, 2.0));
    }

    /// Degenerate case: smoothed endpoint is zero (e.g., heavy gap fill).
    /// The filter conservatively keeps such hills rather than dividing by
    /// zero.
    #[test]
    fn zero_endpoint_hill_kept_conservatively() {
        // Heavy zeros at both ends → smoothed endpoints will be zero.
        let mut profile = vec![0.0f32; 50];
        for i in 20..30 {
            profile[i] = 1000.0;
        }
        let h = make_hill(profile);
        assert!(keep_hill(&h, 40, 2.0));
    }

    /// End-to-end: a mixed input of 3 peaked + 2 flat large hills + 1
    /// small hill is filtered down to 4 (3 peaked + 1 small).
    #[test]
    fn filter_drops_flat_keeps_peaked_and_small() {
        let peaked = |n: usize| {
            (0..n)
                .map(|i| {
                    let d = (i as f32 - n as f32 / 2.0) / (n as f32 / 6.0);
                    10_000.0 * (-d * d).exp() + 50.0
                })
                .collect::<Vec<f32>>()
        };
        let hills = vec![
            make_hill(peaked(45)),
            make_hill(peaked(50)),
            make_hill(peaked(60)),
            make_hill(vec![1000.0; 45]),
            make_hill(vec![1000.0; 50]),
            make_hill(vec![1000.0; 20]), // small — pass-through
        ];
        let filtered = filter_large_baseline_hills(hills, 40, 2.0);
        assert_eq!(filtered.len(), 4);
    }
}
