/// Fill internal zero-runs in an intensity profile using linear interpolation.
///
/// After `ActiveHill::trim()` all leading/trailing zeros are gone, so every
/// zero is an internal gap. Each gap run is filled by interpolating between the
/// nearest non-zero neighbour on each side.
pub fn fill_gaps(profile: &mut Vec<f32>) {
    let n = profile.len();
    let mut i = 0;
    while i < n {
        if profile[i] == 0.0 {
            // Find the end of this zero-run.
            let gap_start = i;
            while i < n && profile[i] == 0.0 {
                i += 1;
            }
            let gap_end = i; // exclusive

            let left_val = if gap_start > 0 { profile[gap_start - 1] } else { 0.0 };
            let right_val = if gap_end < n { profile[gap_end] } else { 0.0 };
            let span = (gap_end - gap_start + 1) as f32; // +1 to include the right anchor

            for (offset, slot) in profile[gap_start..gap_end].iter_mut().enumerate() {
                let t = (offset + 1) as f32 / span;
                *slot = left_val + (right_val - left_val) * t;
            }
        } else {
            i += 1;
        }
    }
}

/// Apply a symmetric box-filter running average with the given half-window.
///
/// At each position `i` the output is the mean of
/// `profile[i-half_window .. i+half_window]` (clamped to array bounds).
/// `half_window == 0` is a no-op.
pub fn running_average(profile: &mut Vec<f32>, half_window: usize) {
    if half_window == 0 || profile.len() < 2 {
        return;
    }
    let n = profile.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = i.saturating_sub(half_window);
        let hi = (i + half_window + 1).min(n);
        let count = (hi - lo) as f32;
        let sum: f32 = profile[lo..hi].iter().sum();
        out.push(sum / count);
    }
    *profile = out;
}

/// Fill gaps then apply running average — the single entry point used by the detector.
pub fn smooth_profile(profile: &mut Vec<f32>, half_window: usize) {
    fill_gaps(profile);
    running_average(profile, half_window);
}

/// Monotonicity-fraction hill score (0–1, higher = more hill-like).
///
/// Counts transitions on the left side that increase toward the apex and
/// transitions on the right side that decrease away from it. Ties count as
/// correct. Returns 1.0 for profiles shorter than 2 scans.
pub fn compute_hill_score(profile: &[f32]) -> f64 {
    let n = profile.len();
    if n < 2 {
        return 1.0;
    }

    let apex = profile
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let mut n_correct = 0usize;
    let mut n_total = 0usize;

    for i in 0..apex {
        n_total += 1;
        if profile[i] <= profile[i + 1] {
            n_correct += 1;
        }
    }

    for i in apex..n - 1 {
        n_total += 1;
        if profile[i] >= profile[i + 1] {
            n_correct += 1;
        }
    }

    if n_total == 0 {
        return 1.0;
    }
    n_correct as f64 / n_total as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fill_single_gap() {
        let mut p = vec![100.0f32, 0.0, 200.0];
        fill_gaps(&mut p);
        assert!((p[1] - 150.0).abs() < 1e-4, "single gap should be midpoint");
    }

    #[test]
    fn fill_multi_gap() {
        let mut p = vec![0.0f32, 0.0, 100.0, 0.0, 0.0, 0.0, 200.0, 0.0, 0.0];
        // After trim the caller guarantees no leading/trailing zeros but we test raw here.
        fill_gaps(&mut p);
        // positions 3,4,5 between 100 and 200: span=4, t=1/4,2/4,3/4
        assert!((p[3] - 125.0).abs() < 1e-3);
        assert!((p[4] - 150.0).abs() < 1e-3);
        assert!((p[5] - 175.0).abs() < 1e-3);
    }

    #[test]
    fn running_average_half1() {
        let mut p = vec![0.0f32, 2.0, 4.0, 6.0, 8.0];
        running_average(&mut p, 1);
        // edge: (0+2)/2=1, interior: (0+2+4)/3=2, (2+4+6)/3=4, (4+6+8)/3=6, edge: (6+8)/2=7
        assert!((p[0] - 1.0).abs() < 1e-4);
        assert!((p[2] - 4.0).abs() < 1e-4);
        assert!((p[4] - 7.0).abs() < 1e-4);
    }

    #[test]
    fn noop_on_zero_window() {
        let mut p = vec![1.0f32, 2.0, 3.0];
        running_average(&mut p, 0);
        assert_eq!(p, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn score_perfect_hill() {
        let p = vec![10.0f32, 50.0, 100.0, 80.0, 20.0];
        assert!((compute_hill_score(&p) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn score_noisy_hill() {
        // Left: (10→80)✓ (80→50)✗ (50→100)✓ → 2/3; Right: (100→60)✓ (60→30)✓ (30→40)✗ (40→5)✓ → 3/4; total 5/7
        let p = vec![10.0f32, 80.0, 50.0, 100.0, 60.0, 30.0, 40.0, 5.0];
        let s = compute_hill_score(&p);
        assert!((s - 5.0 / 7.0).abs() < 1e-6, "got {s}");
    }

    #[test]
    fn score_single_scan() {
        assert!((compute_hill_score(&[42.0f32]) - 1.0).abs() < 1e-6);
    }
}
