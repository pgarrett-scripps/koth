//! End-to-end test of the native Thermo DIA `.raw` MS2 reader.
//!
//! Needs the default `thermo` feature (pure-Rust `opentfraw` `.raw` reader) and
//! a real DIA `.raw` fixture, which CI does not have, so it is `#[ignore]`-gated (mirroring `bruker_diapasef_ms2`). Point it at
//! a DIA `.raw` with `KOTH_DIA_RAW`; it skips (does not fail) when the file is
//! absent. Run explicitly with:
//!
//! ```text
//! KOTH_DIA_RAW=/path/to/dia.raw \
//!   cargo test -p koth_ff --test thermo_dia_ms2 -- --ignored --nocapture
//! ```

#![cfg(feature = "thermo")]

use std::path::PathBuf;

use koth_ff::{
    config::KothConfig, group_ms2_hills_by_window, run_pipeline, run_pipeline_with_ms2,
    PipelineOptions,
};

/// DIA `.raw` fixture path from `KOTH_DIA_RAW`. No baked-in default — this data
/// is not on any known machine here — so the test skips unless the env var points
/// at an existing file.
fn dia_raw_fixture() -> Option<PathBuf> {
    std::env::var("KOTH_DIA_RAW").ok().map(PathBuf::from)
}

#[test]
#[ignore = "needs a DIA .raw via KOTH_DIA_RAW"]
fn thermo_dia_raw_ms2() {
    let Some(path) = dia_raw_fixture() else {
        eprintln!("SKIP thermo_dia_raw_ms2: set KOTH_DIA_RAW to a DIA .raw file");
        return;
    };
    if !path.is_file() {
        eprintln!("SKIP thermo_dia_raw_ms2: no .raw at {}", path.display());
        return;
    }

    let config = KothConfig::default();

    let baseline = run_pipeline(&path, &config, &PipelineOptions::default()).expect("baseline");
    let with_ms2 = run_pipeline_with_ms2(&path, &config, &PipelineOptions::default())
        .expect("with_ms2 on DIA .raw must not error");

    // Turning MS2 on must not perturb the MS1 side.
    assert_eq!(
        baseline.hills.len(),
        with_ms2.hills.len(),
        "MS1 hill count unaffected by emit_ms2 on DIA .raw"
    );

    // MS2 hills exist and every one carries a valid isolation window.
    assert!(
        !with_ms2.ms2_hills.is_empty(),
        "DIA .raw must yield MS2 hills"
    );
    assert!(
        with_ms2
            .ms2_hills
            .iter()
            .all(|h| h.isolation_window.is_some()),
        "every MS2 hill must carry its isolation window"
    );

    let windows = group_ms2_hills_by_window(with_ms2.ms2_hills.clone());
    eprintln!(
        "Thermo DIA MS2: {} hills across {} isolation windows",
        with_ms2.ms2_hills.len(),
        windows.len()
    );
    assert!(
        windows.len() >= 2,
        "expected multiple DIA windows, got {}",
        windows.len()
    );
    for (w, _) in &windows {
        assert!(w.lower <= w.target && w.target <= w.upper && w.lower <= w.upper);
    }
}
