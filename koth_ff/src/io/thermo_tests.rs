//! Unit tests for the pure Thermo DIA MS2 logic: isolation-window derivation and
//! DIA-vs-DDA schedule classification. These exercise no `.raw` file and no .NET
//! runtime, so they run under a plain `cargo test -p koth_ff --features thermo`.
//! End-to-end reading of a real `.raw` is covered separately, `#[ignore]`-gated
//! behind `KOTH_DIA_RAW` in `tests/thermo_dia_ms2.rs`.

use super::*;

/// Build a set of per-MS2-scan windows for a fixed DIA schedule of `n_windows`
/// evenly-tiled `width`-Th windows starting at `start`, repeated for `cycles`
/// (i.e. `n_windows * cycles` MS2 scans, each window recurring `cycles` times).
fn dia_schedule(start: f64, width: f64, n_windows: usize, cycles: usize) -> Vec<IsolationWindow> {
    let mut out = Vec::new();
    for _ in 0..cycles {
        for w in 0..n_windows {
            let lower = start + w as f64 * width;
            let upper = lower + width;
            out.push(IsolationWindow {
                target: (lower + upper) / 2.0,
                lower,
                upper,
            });
        }
    }
    out
}

/// Build DDA-like windows: `n` distinct narrow precursors, each selected once.
fn dda_windows(n: usize) -> Vec<IsolationWindow> {
    (0..n)
        .map(|i| {
            let center = 400.0 + i as f64 * 0.37; // scattered, unique per scan
            IsolationWindow {
                target: center,
                lower: center - 0.35,
                upper: center + 0.35,
            }
        })
        .collect()
}

// -------------------- derive_isolation_window --------------------

#[test]
fn derive_window_uses_recorded_bounds() {
    // Normal DIA scan: 25-Th window centered at 512.5.
    let w = derive_isolation_window(512.5, 500.0, 525.0, 512.4).expect("some window");
    assert!((w.target - 512.5).abs() < 1e-9);
    assert!((w.lower - 500.0).abs() < 1e-9);
    assert!((w.upper - 525.0).abs() < 1e-9);
}

#[test]
fn derive_window_falls_back_to_precursor_mz_for_center() {
    // Some files leave the isolation target unset (0.0); use the precursor m/z.
    let w = derive_isolation_window(0.0, 500.0, 525.0, 512.5).expect("some window");
    assert!((w.target - 512.5).abs() < 1e-9);
    assert!((w.lower - 500.0).abs() < 1e-9);
    assert!((w.upper - 525.0).abs() < 1e-9);
}

#[test]
fn derive_window_collapses_to_point_when_bounds_degenerate() {
    // Unset / degenerate bounds (0/0, or upper <= lower) collapse to the center.
    let w = derive_isolation_window(600.0, 0.0, 0.0, 600.0).expect("some window");
    assert!((w.target - 600.0).abs() < 1e-9);
    assert!((w.lower - 600.0).abs() < 1e-9);
    assert!((w.upper - 600.0).abs() < 1e-9);

    let w2 = derive_isolation_window(600.0, 610.0, 605.0, 600.0).expect("some window");
    assert!((w2.lower - 600.0).abs() < 1e-9 && (w2.upper - 600.0).abs() < 1e-9);
}

#[test]
fn derive_window_none_when_no_usable_center() {
    assert!(derive_isolation_window(0.0, 0.0, 0.0, 0.0).is_none());
    assert!(derive_isolation_window(-1.0, 0.0, 0.0, -5.0).is_none());
    assert!(derive_isolation_window(f64::NAN, 1.0, 2.0, f64::NAN).is_none());
}

// -------------------- is_dia_schedule --------------------

#[test]
fn classifies_fixed_dia_schedule_as_dia() {
    // 40 windows tiling 400–1000 at 15 Th, repeated over 100 cycles.
    let windows = dia_schedule(400.0, 15.0, 40, 100);
    let s = schedule_stats(&windows);
    assert_eq!(s.distinct, 40);
    assert!((s.recurrence - 100.0).abs() < 1e-9);
    assert!(is_dia_schedule(&windows));
}

#[test]
fn classifies_dda_precursor_stream_as_not_dia() {
    // 5000 distinct data-dependent precursors, each once: recurrence ≈ 1.
    let windows = dda_windows(5000);
    let s = schedule_stats(&windows);
    assert_eq!(s.distinct, 5000);
    assert!(s.recurrence < 1.5);
    assert!(!is_dia_schedule(&windows));
}

#[test]
fn single_wide_window_all_ion_is_dia() {
    // All-ion / MSᴱ style: one wide window repeating every cycle.
    let windows = dia_schedule(100.0, 900.0, 1, 500);
    let s = schedule_stats(&windows);
    assert_eq!(s.distinct, 1);
    assert!(is_dia_schedule(&windows));
}

#[test]
fn too_few_ms2_scans_is_not_dia() {
    // Below MIN_MS2_FOR_CLASSIFICATION we refuse to classify (emit nothing).
    let windows = dia_schedule(400.0, 25.0, 3, 2); // 6 scans < 8
    assert!(windows.len() < MIN_MS2_FOR_CLASSIFICATION);
    assert!(!is_dia_schedule(&windows));
}

#[test]
fn recurrence_at_threshold_boundary() {
    // Exactly 3x recurrence passes; just under does not.
    let at = dia_schedule(400.0, 25.0, 4, 3); // 12 scans / 4 windows = 3.0
    assert!((schedule_stats(&at).recurrence - MIN_RECURRENCE).abs() < 1e-9);
    assert!(is_dia_schedule(&at));

    // 4 windows x 2 cycles + 1 extra distinct = 9 scans / 5 windows = 1.8x.
    let mut under = dia_schedule(400.0, 25.0, 4, 2);
    under.push(IsolationWindow {
        target: 999.0,
        lower: 990.0,
        upper: 1008.0,
    });
    let s = schedule_stats(&under);
    assert!(s.recurrence < MIN_RECURRENCE);
    assert!(!is_dia_schedule(&under));
}

#[test]
fn too_many_distinct_windows_is_not_dia() {
    // Above MAX_DIA_WINDOWS distinct windows, even recurring, is treated as DDA.
    let windows = dia_schedule(200.0, 0.5, MAX_DIA_WINDOWS + 1, 4);
    assert!(schedule_stats(&windows).distinct > MAX_DIA_WINDOWS);
    assert!(!is_dia_schedule(&windows));
}

#[test]
fn empty_windows_is_not_dia() {
    assert!(!is_dia_schedule(&[]));
    assert_eq!(schedule_stats(&[]).recurrence, 0.0);
}
