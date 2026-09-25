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
#[path = "filter_tests.rs"]
mod tests;
