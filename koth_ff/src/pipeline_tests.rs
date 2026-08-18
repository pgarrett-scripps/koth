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
use crate::models::{IsolationWindow, Peak, Spectrum};
use crate::output::{write_features_tsv, write_hills_tsv, write_ms2_hills_tsv};

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
        Envelope {
            mono_mz: 500.0,
            charge: 2,
            n_isotopes: 4,
            base: 1.0e5,
        },
        Envelope {
            mono_mz: 700.0,
            charge: 3,
            n_isotopes: 3,
            base: 6.0e4,
        },
        Envelope {
            mono_mz: 900.5,
            charge: 1,
            n_isotopes: 3,
            base: 4.0e4,
        },
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

/// Build a synthetic DIA run: MS1 envelopes (from [`synthetic_run`]) interleaved
/// with MS2 fragment spectra in two fixed isolation windows. Each cycle emits
/// the MS1 spectrum, then one MS2 spectrum per window; every window carries a
/// couple of fragment traces following the same triangular elution shape so MS2
/// hills form. This is the mixed-level stream the file path splits into two
/// readers; the from-spectra API splits it internally.
fn synthetic_dia_run(n_scans: usize) -> Vec<Spectrum> {
    const NEUTRON: f32 = 1.003_354_8;
    let ms1 = synthetic_run(n_scans);

    // Two DIA windows; each isolates a precursor m/z band and shows fragments.
    let windows = [
        IsolationWindow {
            target: 500.0,
            lower: 495.0,
            upper: 505.0,
        },
        IsolationWindow {
            target: 700.0,
            lower: 695.0,
            upper: 705.0,
        },
    ];
    // Fragment mono m/z per window (2 traces each, one charge-1 doublet).
    let frags = [[300.0f32, 450.0], [280.0f32, 610.0]];

    let center = (n_scans as f32 - 1.0) / 2.0;
    let shape = |i: usize| -> f32 {
        let d = (i as f32 - center).abs();
        (1.0 - d / (center + 1.0)).max(0.05)
    };

    let mut out = Vec::with_capacity(n_scans * 3);
    for (i, ms1_spec) in ms1.into_iter().enumerate() {
        out.push(ms1_spec);
        let s = shape(i);
        for (w, iw) in windows.iter().enumerate() {
            let mut peaks = Vec::new();
            for &mono in &frags[w] {
                for j in 0..2usize {
                    let iso_decay = 0.6f32.powi(j as i32);
                    peaks.push(Peak {
                        mz: mono + j as f32 * NEUTRON,
                        intensity: 5.0e4 * iso_decay * s,
                        ion_mobility: 0.0,
                    });
                }
            }
            peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
            out.push(Spectrum {
                scan_index: out.len(),
                retention_time: i as f64 * 0.1,
                peaks,
                ms_level: 2,
                isolation_window: Some(*iw),
            });
        }
    }
    out
}

/// Serialize MS2 hills through the real TSV writer and return contents.
fn ms2_hills_tsv(hills: &[Hill]) -> String {
    let p = tmp_path("ms2_hills");
    write_ms2_hills_tsv(hills, &p).expect("write ms2 hills tsv");
    let s = std::fs::read_to_string(&p).expect("read ms2 hills tsv");
    let _ = std::fs::remove_file(&p);
    s
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
    assert_eq!(
        hills_tsv(&out.hills),
        hills_tsv(&manual_hills),
        "hills TSV parity"
    );
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
    assert_eq!(
        features_tsv(&sink.features),
        features_tsv(&collected.features)
    );
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
        &PipelineOptions {
            scoring: false,
            ..PipelineOptions::default()
        },
    )
    .expect("pipeline no-scoring");

    assert_eq!(out.features.len(), manual_wrapped.len());
    assert_eq!(
        features_tsv(&out.features),
        features_tsv(&manual_wrapped),
        "no-scoring features TSV parity"
    );
}

/// A sink that records both MS1 and MS2 deliveries, so we can assert the call
/// contract (`on_ms2_hills` fires exactly once, after the MS1 side).
#[derive(Default)]
struct Ms2Sink {
    hills: Vec<Hill>,
    features: Vec<ScoredFeature>,
    ms2_hills: Vec<Hill>,
    on_ms2_calls: usize,
    ms2_after_features: bool,
    saw_feature: bool,
}

impl PipelineSink for Ms2Sink {
    fn on_hills(&mut self, hills: Vec<Hill>) {
        self.hills = hills;
    }
    fn on_feature(&mut self, feature: ScoredFeature) {
        self.saw_feature = true;
        self.features.push(feature);
    }
    fn on_ms2_hills(&mut self, hills: Vec<Hill>) {
        self.on_ms2_calls += 1;
        // MS2 must arrive only after the MS1 features have all been delivered.
        self.ms2_after_features = self.saw_feature;
        self.ms2_hills = hills;
    }
}

#[test]
fn emit_ms2_yields_same_hills_as_detect_ms2_hills_from_iter() {
    let config = KothConfig::default();
    let spectra = synthetic_dia_run(11);

    // Ground truth: the MS2 detector run directly over the same spectra.
    let expected = crate::hills::detect_ms2_hills_from_iter(
        spectra.iter().cloned(),
        &config.hills,
        &config.file,
    );
    assert!(!expected.is_empty(), "fixture produced no MS2 hills");
    assert!(
        expected.iter().all(|h| h.isolation_window.is_some()),
        "every MS2 hill must carry its isolation window"
    );

    // Through the collect-all API with MS2 on.
    let out = run_pipeline_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions {
            emit_ms2: true,
            ..PipelineOptions::default()
        },
    )
    .expect("pipeline emit_ms2");

    assert_eq!(out.ms2_hills.len(), expected.len(), "MS2 hill count");
    assert_eq!(
        ms2_hills_tsv(&out.ms2_hills),
        ms2_hills_tsv(&expected),
        "MS2 hills TSV parity with detect_ms2_hills_from_iter"
    );

    // The whole run also carries MS1 features (the two sides coexist).
    assert!(
        !out.features.is_empty(),
        "MS1 features still produced with MS2 on"
    );
}

#[test]
fn emit_ms2_isolates_ms1_and_off_never_emits_ms2() {
    let config = KothConfig::default();
    let spectra = synthetic_dia_run(11);

    // Baseline: the *clean* MS1 result — the MS1-only spectra through the default
    // (MS2-off) pipeline, exactly what an MS1-only caller has always gotten.
    let ms1_only: Vec<Spectrum> = spectra
        .iter()
        .filter(|s| s.ms_level != 2)
        .cloned()
        .collect();
    let baseline =
        run_pipeline_from_spectra(ms1_only.into_iter(), &config, &PipelineOptions::default())
            .expect("baseline");
    assert!(!baseline.hills.is_empty(), "baseline produced MS1 hills");

    // With MS2 *on* over the mixed stream, the MS1 side must be byte-identical to
    // that clean baseline: the fragment (MS2) scans are internally routed to the
    // MS2 detector and never pollute the MS1 hills/features.
    let on = run_pipeline_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions {
            emit_ms2: true,
            ..PipelineOptions::default()
        },
    )
    .expect("ms2 on");
    assert_eq!(
        hills_tsv(&on.hills),
        hills_tsv(&baseline.hills),
        "MS2-on MS1 hills == clean MS1-only baseline"
    );
    assert_eq!(
        features_tsv(&on.features),
        features_tsv(&baseline.features),
        "MS2-on MS1 features == clean MS1-only baseline"
    );
    assert!(
        !on.ms2_hills.is_empty(),
        "MS2 hills present when emit_ms2 is on"
    );

    // MS2 *off* must never populate the MS2 side, regardless of input (the hook
    // is simply not called). MS1 output of the off path is unchanged from before
    // — it is pinned by the existing MS1-only parity tests above.
    let off = run_pipeline_from_spectra(
        ms1_only_iter(&spectra),
        &config,
        &PipelineOptions::default(),
    )
    .expect("ms2 off");
    assert!(
        off.ms2_hills.is_empty(),
        "no MS2 hills when emit_ms2 is off"
    );
    assert_eq!(
        hills_tsv(&off.hills),
        hills_tsv(&baseline.hills),
        "MS2-off over MS1-only input == baseline"
    );
}

/// MS1-only spectra iterator over a mixed DIA run.
fn ms1_only_iter(spectra: &[Spectrum]) -> impl Iterator<Item = Spectrum> + '_ {
    spectra.iter().filter(|s| s.ms_level != 2).cloned()
}

#[test]
fn on_ms2_hills_fires_once_after_features_and_groups_by_window() {
    let config = KothConfig::default();
    let spectra = synthetic_dia_run(11);

    let mut sink = Ms2Sink::default();
    run_pipeline_streaming_from_spectra(
        spectra.iter().cloned(),
        &config,
        &PipelineOptions {
            emit_ms2: true,
            ..PipelineOptions::default()
        },
        &mut sink,
    )
    .expect("stream emit_ms2");

    assert_eq!(sink.on_ms2_calls, 1, "on_ms2_hills fires exactly once");
    assert!(
        sink.ms2_after_features,
        "MS2 must be delivered after MS1 features"
    );
    assert!(!sink.ms2_hills.is_empty());

    // The window mapping koth_tracer relies on: two DIA channels, each hill
    // carrying its window, groupable deterministically.
    let groups = crate::pipeline::group_ms2_hills_by_window(sink.ms2_hills.clone());
    assert_eq!(groups.len(), 2, "two isolation windows");
    // Sorted by window key → targets ascending (500 then 700).
    assert!(groups[0].0.target < groups[1].0.target);
    let regrouped: usize = groups.iter().map(|(_, hs)| hs.len()).sum();
    assert_eq!(regrouped, sink.ms2_hills.len(), "grouping loses no hills");

    // Precursor→window mapping: an MS1 precursor at 500 m/z maps to the first
    // window (495–505), not the second (695–705).
    let precursor_mz = 500.0;
    let hit: Vec<_> = groups
        .iter()
        .filter(|(w, _)| w.lower <= precursor_mz && precursor_mz <= w.upper)
        .collect();
    assert_eq!(hit.len(), 1, "precursor maps to exactly one window");
    assert_eq!(hit[0].0.target, 500.0);
}
