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
    assert!(
        cos > 0.0,
        "expected non-zero cosine at min_overlap=2, got {cos}"
    );
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
