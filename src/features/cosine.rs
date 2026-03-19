use crate::models::Hill;

/// Compute cosine similarity between two hill intensity profiles.
///
/// Profiles are aligned by absolute scan index, padded with zeros
/// where one hill has no data. Each profile is independently
/// min-max normalized before computing the cosine.
///
/// Matches the Python zenith_feature_finder implementation exactly.
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

    // Min-max normalize each profile independently
    min_max_normalize(&mut p1);
    min_max_normalize(&mut p2);

    // Cosine similarity
    let dot: f64 = p1.iter().zip(p2.iter()).map(|(a, b)| a * b).sum();
    let norm1: f64 = p1.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm2: f64 = p2.iter().map(|x| x * x).sum::<f64>().sqrt();

    if norm1 == 0.0 || norm2 == 0.0 {
        return 0.0;
    }

    (dot / (norm1 * norm2)).clamp(0.0, 1.0)
}

fn min_max_normalize(v: &mut [f64]) {
    let max = v.iter().cloned().fold(0.0f64, f64::max);
    if max == 0.0 {
        return;
    }
    let min = v.iter().cloned().fold(max, f64::min);
    let range = max - min;
    if range == 0.0 {
        for x in v.iter_mut() {
            *x = 1.0;
        }
        return;
    }
    for x in v.iter_mut() {
        *x = (*x - min) / range;
    }
}
