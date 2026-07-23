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
#[path = "cosine_tests.rs"]
mod tests;
