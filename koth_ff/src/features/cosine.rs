use crate::models::Hill;

/// Minimum number of scans where BOTH hills have non-gap signal,
/// required before any cosine similarity is computed. Below this, the
/// function returns 0 outright. Inspired by AlphaPept's `correlate`
/// ([feature_finding.py:781-784]) — AlphaPept's comment says "at least
/// an overlap of 3 elements" but their `+3 >` formula actually requires
/// 4 scans. We use the literal intent here: overlap ≥ 3.
///
/// Strict cutoff > smooth attenuation when downstream thresholds
/// (0.5–0.6) might otherwise let edge-overlap pairs slip through.
const MIN_MUTUAL_OVERLAP_SCANS: usize = 3;

/// Compute cosine similarity between two hill intensity profiles.
///
/// Profiles are aligned by absolute scan index and padded with zeros
/// where one hill has no data. Standard L2-normalized cosine is then
/// computed — no baseline subtraction, so two hills peaking at different
/// retention times correctly score low.
///
/// Returns 0.0 when the two hills do not mutually overlap by at least
/// `MIN_MUTUAL_OVERLAP_SCANS` scans — i.e., each hill's range must
/// extend at least that far into the other's range.
pub fn cosine_similarity(hill1: &Hill, hill2: &Hill) -> f64 {
    // Mutual-overlap gate (AlphaPept-style). The overlap region is
    // `[max(start_a, start_b), min(end_a, end_b)]`, size
    // `min(end_a, end_b) - max(start_a, start_b) + 1` (or zero if
    // disjoint). Reject when size < MIN_MUTUAL_OVERLAP_SCANS.
    //
    // Equivalent (and faster than computing max/min twice) to:
    //   start_a + (MIN - 1) > end_b  ||  start_b + (MIN - 1) > end_a
    if hill1.scan_start + MIN_MUTUAL_OVERLAP_SCANS - 1 > hill2.scan_end
        || hill2.scan_start + MIN_MUTUAL_OVERLAP_SCANS - 1 > hill1.scan_end
    {
        return 0.0;
    }

    let min_scan = hill1.scan_start.min(hill2.scan_start);
    let max_scan = hill1.scan_end.max(hill2.scan_end);
    let n = max_scan - min_scan + 1;

    let mut p1 = vec![0.0f64; n];
    let mut p2 = vec![0.0f64; n];

    for (i, &v) in hill1.intensity_profile.iter().enumerate() {
        let idx = hill1.scan_start + i - min_scan;
        if idx < n {
            p1[idx] = v as f64;
        }
    }
    for (i, &v) in hill2.intensity_profile.iter().enumerate() {
        let idx = hill2.scan_start + i - min_scan;
        if idx < n {
            p2[idx] = v as f64;
        }
    }

    let dot: f64 = p1.iter().zip(p2.iter()).map(|(a, b)| a * b).sum();
    let norm1: f64 = p1.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm2: f64 = p2.iter().map(|x| x * x).sum::<f64>().sqrt();

    if norm1 == 0.0 || norm2 == 0.0 {
        return 0.0;
    }

    (dot / (norm1 * norm2)).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn h(scan_start: usize, profile: Vec<f32>) -> Hill {
        let n = profile.len();
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
            scan_start,
            scan_apex: scan_start + n / 2,
            scan_end: scan_start + n.saturating_sub(1),
            n_scans: n,
            skipped_scans: 0,
            intensity_sum: profile.iter().map(|&x| x as f64).sum(),
            intensity_max: profile.iter().copied().fold(0.0f32, f32::max) as f64,
            hill_score: 1.0,
            intensity_profile: Arc::from(profile.as_slice()),
            isolation_window: None,
        }
    }

    /// Identical co-eluting hills score ~1.
    #[test]
    fn identical_profiles_score_one() {
        let p = vec![1.0, 5.0, 10.0, 5.0, 1.0];
        let cos = cosine_similarity(&h(10, p.clone()), &h(10, p));
        assert!((cos - 1.0).abs() < 1e-9, "expected ~1.0, got {cos}");
    }

    /// Hills with NO scan overlap return 0.
    #[test]
    fn disjoint_hills_score_zero() {
        let h1 = h(10, vec![1.0, 5.0, 1.0]);
        let h2 = h(20, vec![1.0, 5.0, 1.0]);
        assert_eq!(cosine_similarity(&h1, &h2), 0.0);
    }

    /// Hills with only 2-scan mutual overlap fail the gate (returns 0)
    /// even though the union-padded cosine would yield a small value.
    #[test]
    fn two_scan_overlap_fails_gate() {
        // h1 scans 10..=14, h2 scans 13..=17. Mutual overlap = scans 13, 14
        // → 2 scans. Below the 3-scan threshold → gate returns 0.
        let h1 = h(10, vec![1.0, 2.0, 10.0, 5.0, 1.0]);
        let h2 = h(13, vec![5.0, 10.0, 2.0, 1.0, 0.5]);
        assert_eq!(cosine_similarity(&h1, &h2), 0.0);
    }

    /// Hills with exactly 3-scan mutual overlap pass the gate and
    /// produce a non-zero similarity.
    #[test]
    fn three_scan_overlap_passes_gate() {
        // h1 scans 10..=14, h2 scans 12..=16. Mutual overlap = scans 12, 13, 14
        // → 3 scans. Gate passes (overlap >= MIN_MUTUAL_OVERLAP_SCANS = 3).
        let h1 = h(10, vec![1.0, 2.0, 10.0, 5.0, 1.0]);
        let h2 = h(12, vec![10.0, 5.0, 1.0, 0.5, 0.1]);
        let cos = cosine_similarity(&h1, &h2);
        assert!(cos > 0.0, "expected non-zero cosine, got {cos}");
        assert!(cos <= 1.0);
    }

    /// Two hills peaking at different RTs (offset by half the profile)
    /// produce a low cosine even when they overlap.
    #[test]
    fn offset_peaks_score_low() {
        // h1: peak at scan 12. h2: peak at scan 17.
        let h1 = h(10, vec![1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let h2 = h(10, vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0]);
        let cos = cosine_similarity(&h1, &h2);
        assert!(cos < 0.6, "expected low cosine for offset peaks, got {cos}");
    }
}
