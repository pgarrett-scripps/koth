//! Unit tests for `averagine.rs`. Wired in via `#[path]`, so this is a child
//! module of the production module and sees its private items.

use super::*;

use crate::config::FeaturesConfig;

#[test]
fn averagine_counts_scales_linearly() {
    // 1500 Da peptide → averagine integer composition
    let (c, h, n, o, s) = averagine_counts(1500.0);
    let scale = 1500.0 / AVG_RESIDUE_MASS;
    assert!((c as f64 - 4.9384 * scale).abs() <= 1.0);
    assert!((h as f64 - 7.7583 * scale).abs() <= 1.0);
    assert!((n as f64 - 1.3577 * scale).abs() <= 1.0);
    assert!((o as f64 - 1.4773 * scale).abs() <= 1.0);
    assert!((s as f64 - 0.0417 * scale).abs() <= 1.0);
}

#[test]
fn averagine_distribution_sums_to_one() {
    for mass in [500.0, 1500.0, 3000.0, 5000.0] {
        let d = averagine_distribution(mass);
        let sum: f64 = d.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "{mass}: sum={sum}");
    }
}

/// Well-aligned 5-hill chain — BC should be high.
#[test]
fn well_aligned_chain_scores_high() {
    let template: [f64; K_PATTERN] = [0.50, 0.30, 0.15, 0.04, 0.008, 0.001, 0.0005, 0.0, 0.0, 0.0];
    let obs = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
    let bc = bhattacharyya_score(&obs, &template);
    assert!(
        bc > 0.95,
        "BC for well-aligned chain should be > 0.95, got {bc}"
    );
}

/// Over-extended chain — same well-aligned first 5 positions plus 3 noise hills
/// where the template predicts essentially nothing. New BC must drop noticeably.
#[test]
fn over_extended_chain_drops_score() {
    let template: [f64; K_PATTERN] = [
        0.50, 0.30, 0.15, 0.04, 0.008, 0.001, 0.0005, 0.0001, 0.0, 0.0,
    ];
    let obs_real = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
    let bc_real = bhattacharyya_score(&obs_real, &template);
    let obs_overextended = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5, 5.0e6, 5.0e6, 5.0e6];
    let bc_over = bhattacharyya_score(&obs_overextended, &template);
    assert!(
        bc_over < bc_real - 0.05,
        "over-extended BC ({bc_over:.3}) should drop at least 0.05 below clean BC ({bc_real:.3})"
    );
}

/// Empty / zero input safe.
#[test]
fn empty_input_returns_zero() {
    let template: [f64; K_PATTERN] = [0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert_eq!(bhattacharyya_score(&[], &template), 0.0);
    assert_eq!(bhattacharyya_score(&[0.0, 0.0], &template), 0.0);
}

/// The default (spanning) variant set is unchanged — this is the byte-safety
/// The declared offsets resolve against `ceil(expected n_S)` per mass band.
#[test]
fn offsets_resolve_against_ceil_bands() {
    // s_float = 3.7524e-4 * M -> ceil steps at ~2665 / 5330 / 7995 Da.
    let d = FeaturesConfig::default().sulfur_offsets;
    for mass in [500.0, 1500.0, 2600.0] {
        assert_eq!(averagine_sulfur_ceil(mass), 1, "mass {mass}");
        assert_eq!(resolve_sulfur_counts(mass, &d), vec![0, 1, 2]);
    }
    for mass in [2700.0, 4000.0, 5300.0] {
        assert_eq!(averagine_sulfur_ceil(mass), 2, "mass {mass}");
        assert_eq!(resolve_sulfur_counts(mass, &d), vec![1, 2, 3]);
    }
    for mass in [5400.0, 7000.0, 7900.0] {
        assert_eq!(averagine_sulfur_ceil(mass), 3, "mass {mass}");
        assert_eq!(resolve_sulfur_counts(mass, &d), vec![2, 3, 4]);
    }
}

/// Negative results saturate at 0 and duplicates collapse, so a wide list on a
/// small peptide does not score the same template repeatedly.
#[test]
fn offsets_saturate_and_dedup() {
    assert_eq!(
        resolve_sulfur_counts(1500.0, &[-3, -2, -1, 0, 1]),
        vec![0, 1, 2]
    );
    assert_eq!(resolve_sulfur_counts(1500.0, &[0, 0, 0]), vec![1]);
    assert_eq!(resolve_sulfur_counts(1500.0, &[]), Vec::<u32>::new());
    // Order is preserved as declared, not sorted.
    assert_eq!(resolve_sulfur_counts(4000.0, &[1, -1, 0]), vec![3, 1, 2]);
}

/// The default list reaches n_S = 2 across the tryptic band — the count the
/// old spanning set `{0, 1, 3, 5}` skipped over 1332-3997 Da.
#[test]
fn default_offsets_reach_two_over_tryptic_band() {
    let d = FeaturesConfig::default().sulfur_offsets;
    for mass in [800.0, 1332.0, 2000.0, 3000.0, 3997.0] {
        let v = resolve_sulfur_counts(mass, &d);
        assert!(v.contains(&2), "mass {mass} set {v:?} skips n_S = 2");
    }
    // Widening by one recovers the no-sulfur template above 2665 Da.
    assert!(!resolve_sulfur_counts(3000.0, &d).contains(&0));
    assert!(resolve_sulfur_counts(3000.0, &[-2, -1, 0, 1]).contains(&0));
}

/// Scoring equals the max over the resolved set.
#[test]
fn scoring_is_max_over_resolved_set() {
    let obs = [4.0e7, 2.6e7, 1.6e7, 5.0e6, 1.0e6];
    let d = FeaturesConfig::default().sulfur_offsets;
    for mass in [900.0, 2000.0, 3000.0, 6000.0] {
        let (c, h, n, o, _) = averagine_counts(mass);
        let manual = resolve_sulfur_counts(mass, &d)
            .iter()
            .map(|&s| bhattacharyya_score(&obs, &cache().distribution(c, h, n, o, s)))
            .fold(0.0f64, f64::max);
        let got = bhattacharyya_score_best_sulfur(&obs, mass, &d).0;
        assert!((got - manual).abs() < 1e-12, "{mass}: {got} vs {manual}");
    }
}

/// An empty offset list disables sulfur awareness: the plain averagine template.
#[test]
fn empty_offsets_match_plain_template() {
    let obs = [5.0e7, 3.0e7, 1.5e7, 4.0e6, 8.0e5];
    for mass in [900.0, 1500.0, 3000.0, 6000.0] {
        let plain = bhattacharyya_score(&obs, &lookup_template(mass));
        let got = bhattacharyya_score_best_sulfur(&obs, mass, &[]).0;
        assert!((got - plain).abs() < 1e-12, "{mass}: {got} vs {plain}");
    }
}

/// The shipped default is `[-1, 0, 1]`.
#[test]
fn shipped_default_offsets() {
    assert_eq!(FeaturesConfig::default().sulfur_offsets, vec![-1, 0, 1]);
}

/// Sanity: sulfur-override path returns a different template than the default.
#[test]
fn sulfur_override_changes_pattern() {
    let m = 1500.0;
    let avg = averagine_distribution(m);
    let high_s = averagine_distribution_with_sulfur(m, 3);
    assert!(
        (high_s[2] - avg[2]).abs() > 0.005,
        "3-S override should noticeably shift M+2 (avg={}, hi-S={})",
        avg[2],
        high_s[2]
    );
}
