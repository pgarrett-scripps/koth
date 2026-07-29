//! Parity tests for the in-process / streaming API.
//!
//! Strategy: run a small **synthetic** run (a handful of co-eluting isotope
//! envelopes, no files, no `tdf`/`thermo` feature needed) through
//!   (a) the manual staged path the binary uses
//!       (`detect_hills_from_iter` → `run_features` → `run_scoring`), and
//!   (b) the new `pipeline` API (both the collect-all and the streaming-sink
//!       shapes),
//! then assert they are identical by serializing each result through the *real*
//! output writers and comparing the emitted TSV byte-for-byte. Serializing
//! through the writers proves the API yields exactly what the file path would
//! have written.

use super::*;

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::config::KothConfig;
use crate::models::{Peak, Spectrum};
use crate::output::{write_features_tsv, write_hills_tsv};

/// One synthetic isotope envelope: `n_isotopes` peaks spaced by neutron/`charge`
/// starting at `mono_mz`, each following the shared elution shape across scans.
struct Envelope {
    mono_mz: f32,
    charge: u8,
    n_isotopes: usize,
    base: f32,
}

/// Build a deterministic synthetic run: several co-eluting isotope envelopes
/// over `n_scans` scans. Isotope intensities decay geometrically and every
/// isotope of an envelope shares the same triangular elution shape, so chains
/// co-elute (high chromatographic cosine) and are charge-resolvable.
fn synthetic_run(n_scans: usize) -> Vec<Spectrum> {
    const NEUTRON: f32 = 1.003_354_8;
    let envelopes = [
        Envelope { mono_mz: 500.0, charge: 2, n_isotopes: 4, base: 1.0e5 },
        Envelope { mono_mz: 700.0, charge: 3, n_isotopes: 3, base: 6.0e4 },
        Envelope { mono_mz: 900.5, charge: 1, n_isotopes: 3, base: 4.0e4 },
    ];

    // Triangular elution shape peaking at the middle scan, in [0.05, 1.0].
    let center = (n_scans as f32 - 1.0) / 2.0;
    let shape = |i: usize| -> f32 {
        let d = (i as f32 - center).abs();
        (1.0 - d / (center + 1.0)).max(0.05)
    };

    (0..n_scans)
        .map(|i| {
            let s = shape(i);
            let mut peaks = Vec::new();
            for e in &envelopes {
                let step = NEUTRON / e.charge as f32;
                for j in 0..e.n_isotopes {
                    let iso_decay = 0.6f32.powi(j as i32);
                    peaks.push(Peak {
                        mz: e.mono_mz + j as f32 * step,
                        intensity: e.base * iso_decay * s,
                        ion_mobility: 0.0,
                    });
                }
            }
            peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
            Spectrum {
                scan_index: i,
                retention_time: i as f64 * 0.1,
                peaks,
                ms_level: 1,
                isolation_window: None,
            }
        })
        .collect()
}

/// Unique temp path so parallel test threads never collide.
fn tmp_path(tag: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let n = N.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!(
        "koth_pipeline_parity_{}_{}_{}.tsv",
        std::process::id(),
        n,
        tag
    ))
}

/// Serialize hills through the real TSV writer and return the file contents.
fn hills_tsv(hills: &[Hill]) -> String {
    let p = tmp_path("hills");
    write_hills_tsv(hills, &p).expect("write hills tsv");
    let s = std::fs::read_to_string(&p).expect("read hills tsv");
    let _ = std::fs::remove_file(&p);
    s
}

/// Serialize scored features through the real TSV writer and return contents.
fn features_tsv(features: &[ScoredFeature]) -> String {
    let p = tmp_path("features");
    write_features_tsv(features, &p).expect("write features tsv");
    let s = std::fs::read_to_string(&p).expect("read features tsv");
    let _ = std::fs::remove_file(&p);
    s
}

/// The manual staged path the binary runs (default config, scoring on).
fn manual_staged(config: &KothConfig, spectra: &[Spectrum]) -> (Vec<Hill>, Vec<ScoredFeature>) {
    let hills =
        crate::hills::detect_hills_from_iter(spectra.iter().cloned(), &config.hills, &config.file);
    let features = crate::run_features(&hills, &config.features, &config.file).expect("features");
    let scored = crate::run_scoring(&features, &config.scoring, &config.features);
    (hills, scored)
}

#[test]
fn collect_all_matches_manual_staged_path() {
    let config = KothConfig::default();
    let spectra = synthetic_run(11);

    // Sanity: the fixture must actually exercise charged, scored features,
    // otherwise the parity check would be vacuous.
    let (manual_hills, manual_scored) = manual_staged(&config, &spectra);
    assert!(!manual_hills.is_empty(), "fixture produced no hills");
    assert!(
        manual_scored.iter().any(|f| f.feature.charge > 0),
        "fixture produced no charged features"
    );

    let out = run_pipeline_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions::default(),
    )
    .expect("pipeline");

    // Same counts...
    assert_eq!(out.hills.len(), manual_hills.len(), "hill count");
    assert_eq!(out.features.len(), manual_scored.len(), "feature count");

    // ...and byte-identical once serialized through the real output writers.
    assert_eq!(hills_tsv(&out.hills), hills_tsv(&manual_hills), "hills TSV parity");
    assert_eq!(
        features_tsv(&out.features),
        features_tsv(&manual_scored),
        "features TSV parity"
    );
}

/// A minimal `koth_tracer`-style sink: takes ownership of the hill set, then
/// folds features incrementally without ever holding the whole feature vector.
#[derive(Default)]
struct CountingSink {
    hills: Vec<Hill>,
    on_hills_calls: usize,
    features: Vec<ScoredFeature>,
}

impl PipelineSink for CountingSink {
    fn on_hills(&mut self, hills: Vec<Hill>) {
        self.on_hills_calls += 1;
        self.hills = hills;
    }
    fn on_feature(&mut self, feature: ScoredFeature) {
        self.features.push(feature);
    }
}

#[test]
fn streaming_sink_matches_collect_all() {
    let config = KothConfig::default();
    let spectra = synthetic_run(11);

    let collected = run_pipeline_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions::default(),
    )
    .expect("collect");

    let mut sink = CountingSink::default();
    run_pipeline_streaming_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions::default(),
        &mut sink,
    )
    .expect("stream");

    // on_hills fires exactly once (whole finalized set), then one on_feature
    // per feature — same set/order as collect-all.
    assert_eq!(sink.on_hills_calls, 1, "on_hills must fire exactly once");
    assert_eq!(sink.hills.len(), collected.hills.len());
    assert_eq!(sink.features.len(), collected.features.len());
    assert_eq!(hills_tsv(&sink.hills), hills_tsv(&collected.hills));
    assert_eq!(features_tsv(&sink.features), features_tsv(&collected.features));
}

#[test]
fn no_scoring_option_matches_binary_no_scoring_wrap() {
    let config = KothConfig::default();
    let spectra = synthetic_run(11);

    // Manual --no-scoring: features wrapped with zeroed scores (mirrors main.rs).
    let hills =
        crate::hills::detect_hills_from_iter(spectra.iter().cloned(), &config.hills, &config.file);
    let features = crate::run_features(&hills, &config.features, &config.file).expect("features");
    let manual_wrapped: Vec<ScoredFeature> =
        features.into_iter().map(super::wrap_unscored).collect();

    let out = run_pipeline_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions { scoring: false },
    )
    .expect("pipeline no-scoring");

    assert_eq!(out.features.len(), manual_wrapped.len());
    assert_eq!(
        features_tsv(&out.features),
        features_tsv(&manual_wrapped),
        "no-scoring features TSV parity"
    );
}
