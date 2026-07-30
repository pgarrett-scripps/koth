//! Path-based parity: the in-process `run_pipeline` API must yield exactly what
//! the `koth_ff` binary's file-writing path would serialize, on a *real* input.
//!
//! Uses the Bruker `.d` fixture, so it needs the `tdf` feature. Marked
//! `#[ignore]` (like `bruker_streaming`) because reading the `.d` is slow in
//! debug — run it explicitly:
//!
//! ```sh
//! cargo test -p koth_ff --features tdf --test pipeline_parity -- --ignored --nocapture
//! ```
#![cfg(feature = "tdf")]

use std::path::PathBuf;

use koth_ff::config::KothConfig;
use koth_ff::output::{write_features_tsv, write_hills_tsv};
use koth_ff::{
    run_features, run_hills_streaming, run_pipeline, run_pipeline_with_ms2, run_scoring,
    PipelineOptions,
};

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/data/example_dda.d")
}

fn hills_tsv(hills: &[koth_ff::Hill]) -> String {
    let p = std::env::temp_dir().join(format!("koth_pp_h_{}_{}.tsv", std::process::id(), rand_tag()));
    write_hills_tsv(hills, &p).expect("write hills tsv");
    let s = std::fs::read_to_string(&p).expect("read hills tsv");
    let _ = std::fs::remove_file(&p);
    s
}

fn features_tsv(features: &[koth_ff::ScoredFeature]) -> String {
    let p = std::env::temp_dir().join(format!("koth_pp_f_{}_{}.tsv", std::process::id(), rand_tag()));
    write_features_tsv(features, &p).expect("write features tsv");
    let s = std::fs::read_to_string(&p).expect("read features tsv");
    let _ = std::fs::remove_file(&p);
    s
}

fn rand_tag() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos() as u64
}

#[test]
#[ignore = "reads a large .d fixture; run explicitly with --features tdf"]
fn run_pipeline_matches_binary_staged_path() {
    let path = fixture();
    assert!(path.is_dir(), "fixture missing: {}", path.display());

    let config = KothConfig::default();

    // (a) The exact staged sequence main.rs runs.
    let hills = run_hills_streaming(&path, &config.hills, &config.file).expect("hills");
    let features = run_features(&hills, &config.features, &config.file).expect("features");
    let scored = run_scoring(&features, &config.scoring, &config.features);

    // (b) The in-process API.
    let out = run_pipeline(&path, &config, &PipelineOptions::default()).expect("pipeline");

    assert_eq!(out.hills.len(), hills.len(), "hill count");
    assert_eq!(out.features.len(), scored.len(), "feature count");
    assert_eq!(hills_tsv(&out.hills), hills_tsv(&hills), "hills TSV parity on real .d");
    assert_eq!(
        features_tsv(&out.features),
        features_tsv(&scored),
        "features TSV parity on real .d"
    );

    eprintln!(
        "run_pipeline == staged path: {} hills, {} features",
        out.hills.len(),
        out.features.len()
    );
}

#[test]
#[ignore = "reads a large .d fixture; run explicitly with --features tdf"]
fn emit_ms2_degrades_gracefully_on_bruker_d() {
    // MS2 on Bruker `.d` is diaPASEF-only. `example_dda.d` is ddaPASEF, so
    // requesting MS2 must NOT error and must NOT disturb the MS1 side — it just
    // yields an empty MS2 set (with a warning). (A diaPASEF `.d` instead recovers
    // per-window MS2 hills — see `bruker_diapasef_ms2` below.)
    let path = fixture();
    assert!(path.is_dir(), "fixture missing: {}", path.display());

    let config = KothConfig::default();

    let baseline = run_pipeline(&path, &config, &PipelineOptions::default()).expect("baseline");
    let with_ms2 = run_pipeline_with_ms2(&path, &config, &PipelineOptions::default())
        .expect("with_ms2 on .d must not error");

    assert!(with_ms2.ms2_hills.is_empty(), "no MS2 hills from a Bruker .d");
    assert_eq!(
        hills_tsv(&with_ms2.hills),
        hills_tsv(&baseline.hills),
        "MS1 hills unaffected by emit_ms2 on .d"
    );
    assert_eq!(
        features_tsv(&with_ms2.features),
        features_tsv(&baseline.features),
        "MS1 features unaffected by emit_ms2 on .d"
    );
}

/// Smallest diaPASEF `.d` on this machine (not checked into the repo). Override
/// with `KOTH_DIAPASEF_D`. The test skips (does not fail) when it is absent, so
/// CI without the raw data is unaffected.
fn diapasef_fixture() -> PathBuf {
    std::env::var("KOTH_DIAPASEF_D").map(PathBuf::from).unwrap_or_else(|_| {
        PathBuf::from(
            "/home/patrick-garrett/Repos/d_noise/benchmark/data/dia_5min/raw/\
             LFQ_Ultra2_diaPASEF_5min_50ng_Condition_A_REP1.d",
        )
    })
}

#[test]
#[ignore = "reads a large diaPASEF .d fixture; run explicitly with --features tdf"]
fn bruker_diapasef_ms2() {
    // End-to-end: on a real diaPASEF `.d`, `emit_ms2` must recover the fixed
    // isolation windows and stamp every MS2 hill with the window it came from,
    // while leaving the MS1 side byte-identical to a plain MS1 run.
    let path = diapasef_fixture();
    if !path.is_dir() {
        eprintln!("SKIP bruker_diapasef_ms2: no diaPASEF fixture at {}", path.display());
        return;
    }

    let config = KothConfig::default();

    let baseline = run_pipeline(&path, &config, &PipelineOptions::default()).expect("baseline");
    let with_ms2 = run_pipeline_with_ms2(&path, &config, &PipelineOptions::default())
        .expect("with_ms2 on diaPASEF .d must not error");

    // MS1 side is untouched by turning MS2 on.
    assert_eq!(
        hills_tsv(&with_ms2.hills),
        hills_tsv(&baseline.hills),
        "MS1 hills unaffected by emit_ms2 on diaPASEF .d"
    );

    // MS2 hills exist and every one carries a valid isolation window.
    assert!(!with_ms2.ms2_hills.is_empty(), "diaPASEF .d must yield MS2 hills");
    assert!(
        with_ms2.ms2_hills.iter().all(|h| h.isolation_window.is_some()),
        "every MS2 hill must carry its isolation window"
    );

    // Windows recovered (this 5-min diaPASEF method has 8 groups x 3 = 24).
    let windows = koth_ff::group_ms2_hills_by_window(with_ms2.ms2_hills.clone());
    eprintln!(
        "diaPASEF MS2: {} hills across {} isolation windows",
        with_ms2.ms2_hills.len(),
        windows.len()
    );
    assert!(windows.len() >= 2, "expected multiple DIA windows, got {}", windows.len());
    for (w, _) in &windows {
        assert!(w.lower < w.upper && w.lower <= w.target && w.target <= w.upper);
    }
}
