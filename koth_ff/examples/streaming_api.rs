//! In-process feature-finding API — usage for the two real downstream tools.
//!
//! Run against any mzML (or Bruker `.d` with `--features tdf`):
//!
//! ```sh
//! cargo run --release -p koth_ff --example streaming_api -- path/to/run.mzML
//! ```
//!
//! None of the consumers below writes or reads any intermediate parquet/tsv
//! file: hills and features are obtained directly as owned Rust structs in
//! memory. The third consumer (`tracer_with_ms2`) also pulls the DIA MS2
//! fragment hills, partitioned by isolation window, in the same in-process call
//! — supported on mzML, Bruker **diaPASEF** `.d` (`--features tdf`, windows read
//! straight from the raw quadrupole settings), and Thermo **DIA** `.raw`
//! (`--features thermo`, one MS2 spectrum per scan stamped with its precursor
//! isolation window). On DDA acquisitions (ddaPASEF `.d` / DDA `.raw`) the MS2
//! set is empty (with a warning) and the MS1 side is unaffected.

use std::path::PathBuf;

use koth_ff::{
    config::KothConfig, group_ms2_hills_by_window, run_pipeline, run_pipeline_streaming,
    run_pipeline_with_ms2, Feature, Hill, IsolationWindow, PipelineOptions, PipelineSink,
    ScoredFeature,
};

// ---------------------------------------------------------------------------
// uno-style consumer: collect the final feature list into memory and search it.
// ---------------------------------------------------------------------------
fn uno_style(path: &std::path::Path, config: &KothConfig) {
    // One call, everything owned in memory. No files touched.
    let out = run_pipeline(path, config, &PipelineOptions::default()).expect("pipeline");

    // uno searches MS1 features by neutral mass + charge + RT.
    let searchable: Vec<&ScoredFeature> = out
        .features
        .iter()
        .filter(|f| f.feature.charge > 0) // same set the features file would hold
        .collect();

    println!(
        "[uno] {} hills, {} charged features available for MS1 search",
        out.hills.len(),
        searchable.len()
    );
    if let Some(top) = searchable.iter().max_by(|a, b| {
        a.feature
            .total_intensity()
            .partial_cmp(&b.feature.total_intensity())
            .unwrap()
    }) {
        println!(
            "[uno]   most intense: mass {:?} z{} score {:.3}",
            top.monoisotopic_neutral_mass(),
            top.feature.charge,
            top.combined_score
        );
    }
}

// ---------------------------------------------------------------------------
// koth_tracer-style consumer: take ownership of the hill set to build a trace
// index, then fold features incrementally without holding the whole vector.
// ---------------------------------------------------------------------------
#[derive(Default)]
struct TracerSink {
    hills: Vec<Hill>, // owned once, up front — used to trace precursor<->fragment
    n_features: usize,
    total_intensity: f64,
}

impl PipelineSink for TracerSink {
    fn on_hills(&mut self, hills: Vec<Hill>) {
        // koth_tracer indexes the finalized hill set here (m/z, RT, IM) so it
        // can attribute fragment hills to precursor features as they arrive.
        println!("[tracer] indexing {} finalized hills", hills.len());
        self.hills = hills;
    }

    fn on_feature(&mut self, feature: ScoredFeature) {
        // Consume-and-drop: fold each feature into whatever the tracer keeps.
        // The full feature vector is never materialized on the consumer side.
        let f: &Feature = &feature.feature;
        self.total_intensity += f.total_intensity();
        self.n_features += 1;
        // `feature` is dropped here, freeing its memory before the next arrives.
    }
}

fn tracer_style(path: &std::path::Path, config: &KothConfig) {
    let mut sink = TracerSink::default();
    run_pipeline_streaming(path, config, &PipelineOptions::default(), &mut sink)
        .expect("streaming pipeline");
    println!(
        "[tracer] streamed {} features (sum intensity {:.3e}) over {} indexed hills",
        sink.n_features,
        sink.total_intensity,
        sink.hills.len()
    );
}

// ---------------------------------------------------------------------------
// koth_tracer with MS2: get MS1 precursor features AND the DIA fragment hills,
// partitioned by isolation window, in one in-process call — no files. Then map
// each precursor to the window(s) that isolated it. MS2 side: mzML + Bruker
// diaPASEF `.d` + Thermo DIA `.raw`.
// ---------------------------------------------------------------------------
fn tracer_with_ms2(path: &std::path::Path, config: &KothConfig) {
    // `run_pipeline_with_ms2` forces MS2 emission on; equivalently, set
    // `PipelineOptions { emit_ms2: true, ..Default::default() }` on any entry.
    let out = run_pipeline_with_ms2(path, config, &PipelineOptions::default())
        .expect("pipeline with ms2");

    // Group the flat MS2 hill set into DIA channels (each hill already carries
    // its window in `hill.isolation_window`; this is just the convenience view).
    let windows: Vec<(IsolationWindow, Vec<Hill>)> = group_ms2_hills_by_window(out.ms2_hills);
    println!(
        "[tracer+ms2] {} MS1 features, {} DIA isolation windows",
        out.features.len(),
        windows.len()
    );

    // Precursor -> window mapping: a precursor feature is isolated by every
    // window whose [lower, upper] contains its (charge-1) precursor m/z.
    for f in out.features.iter().filter(|f| f.feature.charge > 0).take(3) {
        let mz = f.feature.monoisotopic_mz(); // precursor m/z of the feature
        let isolating: Vec<&(IsolationWindow, Vec<Hill>)> = windows
            .iter()
            .filter(|(w, _)| w.lower <= mz && mz <= w.upper)
            .collect();
        let n_frag: usize = isolating.iter().map(|(_, hs)| hs.len()).sum();
        println!(
            "[tracer+ms2]   precursor m/z {:.4} z{} -> {} window(s), {} candidate fragment hills",
            mz,
            f.feature.charge,
            isolating.len(),
            n_frag
        );
    }
}

fn main() {
    let path: PathBuf = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .expect("usage: streaming_api <input.mzML | input.d>");
    let config = KothConfig::default();

    tracer_style(&path, &config);
    uno_style(&path, &config);
    tracer_with_ms2(&path, &config);
}
