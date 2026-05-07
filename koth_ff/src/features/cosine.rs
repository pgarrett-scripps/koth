use crate::models::Hill;

/// Compute cosine similarity between two hill intensity profiles.
///
/// Profiles are aligned by absolute scan index and padded with zeros
/// where one hill has no data. Standard L2-normalized cosine is then
/// computed — no baseline subtraction, so two hills peaking at different
/// retention times correctly score low.
pub fn cosine_similarity(hill1: &Hill, hill2: &Hill) -> f64 {
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
