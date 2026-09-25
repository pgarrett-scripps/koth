//! Small shared numeric helpers.
//!
//! These consolidate copy-pasted implementations that had *identical*
//! semantics. Where existing local copies differ subtly (e.g. returning `0.0`
//! rather than `NAN` for degenerate input, or computing a lower median instead
//! of averaging the two central values) they are intentionally left in place —
//! see the Phase-3b dedup notes.

/// Sample standard deviation (divides by `n − 1`).
///
/// Returns `f64::NAN` when fewer than two values are supplied — matching the
/// report/align-report summary helpers this replaces.
pub fn std_dev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return f64::NAN;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64;
    var.sqrt()
}

/// Median of `v`, sorting the slice in place.
///
/// Returns `0.0` for an empty slice; for even lengths returns the average of
/// the two central values. Matches the median-of-ratios helper this replaces.
pub fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn std_dev_nan_below_two() {
        assert!(std_dev(&[]).is_nan());
        assert!(std_dev(&[1.0]).is_nan());
    }

    #[test]
    fn std_dev_sample_variance() {
        // Sample std of [2, 4, 4, 4, 5, 5, 7, 9]: mean 5, var = 32/7, sqrt ≈ 2.1381
        let s = std_dev(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);
        assert!((s - (32.0f64 / 7.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn median_empty_is_zero() {
        let mut v: Vec<f64> = vec![];
        assert_eq!(median(&mut v), 0.0);
    }

    #[test]
    fn median_odd_and_even() {
        let mut odd = [3.0, 1.0, 2.0];
        assert_eq!(median(&mut odd), 2.0);
        let mut even = [4.0, 1.0, 3.0, 2.0];
        assert_eq!(median(&mut even), 2.5);
    }

    #[test]
    fn median_sorts_in_place() {
        let mut v = [3.0, 1.0, 2.0];
        median(&mut v);
        assert_eq!(v, [1.0, 2.0, 3.0]);
    }
}
