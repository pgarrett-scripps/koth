use super::*;
use crate::models::Hill;
use std::sync::Arc;

/// Minimal hill for grid tests: only the fields `build_grid` / `scan_rt` /
/// `fill_row_from_hill` read are meaningful; the rest are inert.
fn hill(
    mz: f64,
    rt_start: f64,
    rt_end: f64,
    im: f64,
    scan_start: usize,
    profile: Vec<f32>,
) -> Hill {
    let n = profile.len();
    let int_sum: f64 = profile.iter().map(|&x| x as f64).sum();
    let int_max: f64 = profile.iter().copied().fold(0.0f32, f32::max) as f64;
    Hill {
        hill_id: 0,
        mz,
        mz_std: 0.0,
        mz_se: 0.0,
        rt: (rt_start + rt_end) / 2.0,
        rt_start,
        rt_end,
        rt_width: rt_end - rt_start,
        im,
        im_std: 0.0,
        scan_start,
        scan_apex: scan_start + n / 2,
        scan_end: scan_start + n.saturating_sub(1),
        n_scans: n,
        skipped_scans: 0,
        intensity_sum: int_sum,
        intensity_max: int_max,
        hill_score: 1.0,
        intensity_profile: Arc::from(profile.as_slice()),
        isolation_window: None,
        faims_cv: None,
    }
}

#[test]
fn sorted_hills_orders_by_mz_and_keeps_index() {
    let hills = vec![
        hill(500.0, 0.0, 1.0, 0.0, 0, vec![1.0]),
        hill(300.0, 0.0, 1.0, 0.0, 0, vec![1.0]),
        hill(400.0, 0.0, 1.0, 0.0, 0, vec![1.0]),
    ];
    let sorted = SortedHills::from_hills(&hills);
    assert_eq!(sorted.len(), 3);
    let mzs: Vec<f64> = sorted.keys.iter().map(|k| k.mz).collect();
    assert_eq!(mzs, vec![300.0, 400.0, 500.0]);
    // hill_idx points back to the original (unsorted) position.
    assert_eq!(sorted.keys[0].hill_idx, 1); // mz 300 was index 1
    assert_eq!(sorted.keys[2].hill_idx, 0); // mz 500 was index 0
}

#[test]
fn scan_rt_prefers_scan_times() {
    let h = hill(500.0, 0.0, 1.0, 0.0, 2, vec![1.0, 1.0, 1.0]);
    let scan_times = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
    // profile index 1 -> absolute scan 2 + 1 = 3 -> scan_times[3]
    assert_eq!(scan_rt(&h, 1, &scan_times), 0.3);
}

#[test]
fn scan_rt_interpolates_without_scan_times() {
    let h = hill(500.0, 1.0, 2.0, 0.0, 0, vec![0.0; 5]);
    // rt_start + (rt_end - rt_start) * i / (n - 1); i=2, n=5 -> 1.0 + 1.0*2/4 = 1.5
    assert_eq!(scan_rt(&h, 2, &[]), 1.5);
    // Single-sample profile has no span -> falls back to hill.rt.
    let single = hill(500.0, 1.0, 2.0, 0.0, 0, vec![7.0]);
    assert_eq!(scan_rt(&single, 0, &[]), single.rt);
}

#[test]
fn fill_row_bins_intensities_into_columns() {
    let h = hill(500.0, 0.0, 1.0, 0.0, 0, vec![10.0, 20.0, 30.0]);
    let scan_times = vec![0.0, 0.25, 0.75];
    let mut row = vec![0.0f32; 4];
    let filled = fill_row_from_hill(&h, &mut row, 0.0, 1.0, 4, &scan_times);
    assert!(filled);
    // t = 0.0 -> col 0; t = 0.25 -> col 1; t = 0.75 -> col 3.
    assert_eq!(row, vec![10.0, 20.0, 0.0, 30.0]);
}

#[test]
fn fill_row_rejects_empty_profile_and_degenerate_window() {
    let empty = hill(500.0, 0.0, 1.0, 0.0, 0, vec![]);
    let mut row = vec![0.0f32; 4];
    assert!(!fill_row_from_hill(&empty, &mut row, 0.0, 1.0, 4, &[]));

    let h = hill(500.0, 0.0, 1.0, 0.0, 0, vec![5.0]);
    // rt_max <= rt_min -> zero span -> nothing filled.
    assert!(!fill_row_from_hill(&h, &mut row, 1.0, 1.0, 4, &[]));
    assert_eq!(row, vec![0.0f32; 4]);
}

#[test]
fn fill_row_skips_nonpositive_and_out_of_window_samples() {
    // Middle sample is zero; last sample lands at t=1.0 (excluded, half-open).
    let h = hill(500.0, 0.0, 1.0, 0.0, 0, vec![10.0, 0.0, 40.0]);
    let scan_times = vec![0.0, 0.5, 1.0];
    let mut row = vec![0.0f32; 2];
    let filled = fill_row_from_hill(&h, &mut row, 0.0, 1.0, 2, &scan_times);
    assert!(filled);
    assert_eq!(
        row,
        vec![10.0, 0.0],
        "zero sample skipped, t=1.0 sample excluded"
    );
}

fn grid_config() -> LfqConfig {
    let mut cfg = LfqConfig::default();
    cfg.n_isotopes = 2;
    cfg.grid_cols = 4;
    cfg
}

#[test]
fn build_grid_matches_monoisotope_and_m1_rows() {
    let cfg = grid_config();
    let mono = 500.0;
    let m1 = mono + C13_NEUTRON; // charge 1
    let scan_times = vec![0.0, 0.5, 1.0];
    let hills = vec![
        hill(mono, 0.4, 0.6, 0.0, 1, vec![100.0]),   // -> iso 0
        hill(m1, 0.4, 0.6, 0.0, 1, vec![80.0]),      // -> iso 1
        hill(500.02, 0.4, 0.6, 0.0, 1, vec![999.0]), // outside 10 ppm -> ignored
    ];
    let sorted = SortedHills::from_hills(&hills);
    let mut grid = XicGrid::empty(cfg.n_isotopes, cfg.grid_cols, 0.0, 1.0);

    build_grid(
        &mut grid,
        &hills,
        &sorted,
        &scan_times,
        mono,
        1,
        0.5,
        0.0,
        0.5,
        &cfg,
    );

    assert_eq!(
        grid.n_slots_filled, 2,
        "mono + M+1 rows filled, off-tol hill ignored"
    );
    // scan_start 1 -> scan_times[1] = 0.5 -> t=0.5 -> col 2 of 4.
    assert_eq!(grid.intensities[0][2], 100.0);
    assert_eq!(grid.intensities[1][2], 80.0);
    assert!((grid.winner_mz[0] - mono).abs() < 1e-9);
    assert!((grid.winner_mz[1] - m1).abs() < 1e-9);
}

#[test]
fn build_grid_excludes_hills_outside_rt_window_and_im_tolerance() {
    let cfg = grid_config();
    let mono = 500.0;
    let scan_times = vec![0.0, 0.5, 1.0];
    let hills = vec![
        // Correct m/z and IM but entirely outside the RT window [0,1].
        hill(mono, 5.0, 6.0, 1.2, 5, vec![100.0]),
        // Correct m/z, in the window, but IM far from target 1.2 (> 0.05).
        hill(mono, 0.4, 0.6, 2.0, 1, vec![100.0]),
    ];
    let sorted = SortedHills::from_hills(&hills);
    let mut grid = XicGrid::empty(cfg.n_isotopes, cfg.grid_cols, 0.0, 1.0);

    build_grid(
        &mut grid,
        &hills,
        &sorted,
        &scan_times,
        mono,
        1,
        0.5,
        1.2,
        0.5,
        &cfg,
    );

    assert_eq!(
        grid.n_slots_filled, 0,
        "RT-out and IM-mismatch hills both rejected"
    );
    assert!(grid.is_empty());
}
