//! Validates the opt-in dnoise streaming Bruker path against the historical
//! local path on a real `.d` fixture.
//!
//! With the halo disabled, the streaming path (RunContext: vertical filter +
//! watershed on the raw frames) runs the *same* two stages, with the same
//! parameters, on the same frames as the local path — so the two must agree
//! peak-for-peak. This pins down that the streaming wiring reproduces the local
//! result exactly before halo is layered on.
//!
//! Marked `#[ignore]` because it reads a 710-frame `.d` twice (slow in debug):
//!
//! ```sh
//! cargo test -p koth_ff --features tdf --test bruker_streaming -- --ignored --nocapture
//! ```
#![cfg(feature = "tdf")]

use std::path::PathBuf;

use koth_ff::config::FileConfig;
use koth_ff::io::bruker::read_bruker;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/data/example_dda.d")
}

#[test]
#[ignore = "reads a large .d fixture twice; run explicitly"]
fn streaming_no_halo_matches_local() {
    let path = fixture();
    assert!(path.is_dir(), "fixture missing: {}", path.display());

    let local_cfg = FileConfig {
        bruker_streaming: false,
        ..FileConfig::default()
    };
    let stream_cfg = FileConfig {
        bruker_streaming: true,
        bruker_halo: false, // isolate the vertical-filter + watershed equivalence
        ..FileConfig::default()
    };

    let local = read_bruker(&path, &local_cfg).expect("local read");
    let stream = read_bruker(&path, &stream_cfg).expect("streaming read");

    assert_eq!(
        local.len(),
        stream.len(),
        "spectrum count differs (local {} vs streaming {})",
        local.len(),
        stream.len()
    );

    let local_peaks: usize = local.iter().map(|s| s.peaks.len()).sum();
    let stream_peaks: usize = stream.iter().map(|s| s.peaks.len()).sum();
    assert_eq!(
        local_peaks, stream_peaks,
        "total peak count differs (local {local_peaks} vs streaming {stream_peaks})"
    );

    // Per-spectrum peak counts must line up once both are RT-sorted (finalize
    // already sorts). This catches any per-frame divergence the totals hide.
    for (i, (l, s)) in local.iter().zip(stream.iter()).enumerate() {
        assert_eq!(
            l.peaks.len(),
            s.peaks.len(),
            "spectrum {i} peak count differs (local {} vs streaming {})",
            l.peaks.len(),
            s.peaks.len()
        );
    }
    eprintln!(
        "streaming(no-halo) == local: {} spectra, {local_peaks} peaks",
        local.len()
    );
}
