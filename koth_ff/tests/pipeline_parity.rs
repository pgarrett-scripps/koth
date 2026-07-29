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
use koth_ff::{run_features, run_hills_streaming, run_pipeline, run_scoring, PipelineOptions};

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
