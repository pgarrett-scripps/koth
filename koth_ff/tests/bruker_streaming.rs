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
//! cargo test -p koth-ms --features tdf --test bruker_streaming -- --ignored --nocapture
//! ```
#![cfg(feature = "tdf")]

use std::path::PathBuf;

use koth_ms::config::FileConfig;
use koth_ms::io::bruker::read_bruker;

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

/// Fingerprints captured from the collecting readers and detector at 0908b50.
/// Compare every spectrum field/peak and every hill field/profile, not just counts.
#[test]
#[ignore = "reads the real .d fixture six times; run in release mode"]
fn bounded_ms1_matches_frozen_batch_outputs() {
    use koth_ms::{config::KothConfig, io::stream_spectra, run_hills_streaming};
    use sha2::{Digest, Sha256};
    for (streaming, halo, spectrum_hash, hill_hash) in [
        (
            false,
            false,
            "1a09e43618b74d1d84cba031f13d0251a40747d602314645b1955f350c99483d",
            "3fd083ac9866751ab01de0720e9a7e762c2bedd08c3377898c740f58a49477c4",
        ),
        (
            true,
            false,
            "1a09e43618b74d1d84cba031f13d0251a40747d602314645b1955f350c99483d",
            "3fd083ac9866751ab01de0720e9a7e762c2bedd08c3377898c740f58a49477c4",
        ),
        (
            true,
            true,
            "972d279ae68d29ea83bfe092f35322e852c0a87203fb97f4f3a501f995670a53",
            "2fc7b7d85fadb4fb4eda37550e2925a49dbc85ad8aadead6d42c4a8da31e9fe4",
        ),
    ] {
        let mut cfg = KothConfig::default();
        cfg.file.bruker_streaming = streaming;
        cfg.file.bruker_halo = halo;
        let spectra = stream_spectra(&fixture(), &cfg.file)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(
            format!("{:x}", Sha256::digest(format!("{spectra:?}").as_bytes())),
            spectrum_hash
        );
        drop(spectra);
        let hills = run_hills_streaming(&fixture(), &cfg.hills, &cfg.file).unwrap();
        // HashMap channel iteration need not preserve hill order. Normalize IDs
        // and sort complete records; all numerical values must remain identical.
        let mut rows: Vec<_> = hills
            .into_iter()
            .map(|mut h| {
                h.hill_id = 0;
                format!("{h:?}")
            })
            .collect();
        rows.sort();
        assert_eq!(
            format!("{:x}", Sha256::digest(format!("{rows:?}").as_bytes())),
            hill_hash
        );
    }
}
