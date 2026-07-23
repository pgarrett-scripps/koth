use crate::models::Hill;

/// Compute cosine similarity between two hill intensity profiles.
///
/// Profiles are aligned by absolute scan index and padded with zeros
/// where one hill has no data. Standard L2-normalized cosine is then
/// computed — no baseline subtraction, so two hills peaking at different
/// retention times correctly score low.
///
/// Returns 0.0 when the two hills do not mutually overlap by at least
/// `min_overlap` scans — i.e., their scan ranges must intersect in at
/// least that many scans. This is the `features.min_scan_overlap` config
/// value (default 3); inspired by AlphaPept's `correlate` overlap cutoff.
/// A strict cutoff beats smooth attenuation when downstream thresholds
/// (0.5–0.6) might otherwise let edge-overlap pairs slip through — but on
/// fast gradients (3–5 scan hills) a value of 2 recovers real isotope
/// pairs the 3-scan gate rejects (see `min_scan_overlap` docs).
///
/// Normalisation is **union-padded**: each profile's L2 norm is taken over its
/// full length while the dot product covers only the scan overlap, so two hills
/// peaking at different RTs (or with very different lengths) correctly score
/// low. This mismatched-length penalty acts as an implicit noise filter and is
/// what the PXD003881 LFQ quant tuning was validated against.
pub fn cosine_similarity(
    hill1: &Hill,
    hill2: &Hill,
    min_overlap: usize,
) -> f64 {
    // Mutual-overlap gate. Overlap size = `min(end_a, end_b) - max(start_a,
    // start_b) + 1` (or 0 if disjoint). Reject when it is below `min_overlap`.
    // Computed directly on scan_start/scan_end (underflow-safe for small
    // `min_overlap`, unlike the old `start + MIN - 1 > end` form).
    let gate_lo = hill1.scan_start.max(hill2.scan_start);
    let gate_hi = hill1.scan_end.min(hill2.scan_end);
    let overlap = if gate_hi >= gate_lo {
        gate_hi - gate_lo + 1
    } else {
        0
    };
    if overlap < min_overlap {
        return 0.0;
    }

    // Overlap is bounded by each hill's *actual profile extent*
    // (`[scan_start, scan_start + len - 1]`) rather than `scan_end`, so the
    // direct indexing below can never run off the end of a profile even if a
    // hill's `scan_end` were ever out of sync with its profile length.
    let end1 = hill1.scan_start + hill1.intensity_profile.len().saturating_sub(1);
    let end2 = hill2.scan_start + hill2.intensity_profile.len().saturating_sub(1);
    let ov_start = hill1.scan_start.max(hill2.scan_start);
    let ov_end = end1.min(end2);
    if ov_end < ov_start {
        return 0.0;
    }

    // Union-padded. Allocation-free equivalent of zero-padding both
    // profiles onto a common scan axis and taking the L2-normalized cosine:
    //
    //   * Each profile's L2 norm is invariant to zero-padding — `√(Σ xᵢ²)` is
    //     the same with or without padding — so we sum each hill's own squares.
    //   * The dot product is nonzero only over the scan overlap; everywhere else
    //     one side is a padding zero. So we walk just the overlap.
    let norm1: f64 = hill1
        .intensity_profile
        .iter()
        .map(|&x| {
            let x = x as f64;
            x * x
        })
        .sum::<f64>()
        .sqrt();
    let norm2: f64 = hill2
        .intensity_profile
        .iter()
        .map(|&x| {
            let x = x as f64;
            x * x
        })
        .sum::<f64>()
        .sqrt();

    if norm1 == 0.0 || norm2 == 0.0 {
        return 0.0;
    }

    let mut dot = 0.0f64;
    for scan in ov_start..=ov_end {
        let a = hill1.intensity_profile[scan - hill1.scan_start] as f64;
        let b = hill2.intensity_profile[scan - hill2.scan_start] as f64;
        dot += a * b;
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

    // Default overlap gate used by the pre-existing tests (config default).
    const OV: usize = 3;

    /// Identical co-eluting hills score ~1.
    #[test]
    fn identical_profiles_score_one() {
        let p = vec![1.0, 5.0, 10.0, 5.0, 1.0];
        let cos = cosine_similarity(&h(10, p.clone()), &h(10, p.clone()), OV);
        assert!((cos - 1.0).abs() < 1e-9, "expected ~1.0, got {cos}");
    }

    /// Hills with NO scan overlap return 0.
    #[test]
    fn disjoint_hills_score_zero() {
        let h1 = h(10, vec![1.0, 5.0, 1.0]);
        let h2 = h(20, vec![1.0, 5.0, 1.0]);
        assert_eq!(cosine_similarity(&h1, &h2, OV), 0.0);
    }

    /// Hills with only 2-scan mutual overlap fail the default 3-scan gate
    /// (returns 0) even though the union-padded cosine would yield a value.
    #[test]
    fn two_scan_overlap_fails_default_gate() {
        // h1 scans 10..=14, h2 scans 13..=17. Mutual overlap = scans 13, 14
        // → 2 scans. Below min_overlap=3 → gate returns 0.
        let h1 = h(10, vec![1.0, 2.0, 10.0, 5.0, 1.0]);
        let h2 = h(13, vec![5.0, 10.0, 2.0, 1.0, 0.5]);
        assert_eq!(cosine_similarity(&h1, &h2, 3), 0.0);
    }

    /// The SAME 2-scan-overlap pair passes once `min_overlap` is lowered to 2,
    /// producing a non-zero similarity. This is the configurable rescue for
    /// fast-gradient data where isotopes overlap the mono by only ~2 scans.
    #[test]
    fn two_scan_overlap_passes_when_gate_lowered() {
        let h1 = h(10, vec![1.0, 2.0, 10.0, 5.0, 1.0]);
        let h2 = h(13, vec![5.0, 10.0, 2.0, 1.0, 0.5]);
        let cos = cosine_similarity(&h1, &h2, 2);
        assert!(cos > 0.0, "expected non-zero cosine at min_overlap=2, got {cos}");
        assert!(cos <= 1.0);
    }

    /// Hills with exactly 3-scan mutual overlap pass the default gate.
    #[test]
    fn three_scan_overlap_passes_gate() {
        // h1 scans 10..=14, h2 scans 12..=16. Mutual overlap = scans 12, 13, 14
        // → 3 scans. Gate passes (overlap >= min_overlap = 3).
        let h1 = h(10, vec![1.0, 2.0, 10.0, 5.0, 1.0]);
        let h2 = h(12, vec![10.0, 5.0, 1.0, 0.5, 0.1]);
        let cos = cosine_similarity(&h1, &h2, OV);
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
        let cos = cosine_similarity(&h1, &h2, OV);
        assert!(cos < 0.6, "expected low cosine for offset peaks, got {cos}");
    }
}
