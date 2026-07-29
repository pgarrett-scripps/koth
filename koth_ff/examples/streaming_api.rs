//! In-process feature-finding API — usage for the two real downstream tools.
//!
//! Run against any mzML (or Bruker `.d` with `--features tdf`):
//!
//! ```sh
//! cargo run --release -p koth_ff --example streaming_api -- path/to/run.mzML
//! ```
//!
//! Neither consumer below writes or reads any intermediate parquet/tsv file:
//! hills and features are obtained directly as owned Rust structs in memory.

use std::path::PathBuf;

use koth_ff::{
    config::KothConfig, run_pipeline, run_pipeline_streaming, Feature, Hill, PipelineOptions,
    PipelineSink, ScoredFeature,
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
        sink.n_features, sink.total_intensity, sink.hills.len()
    );
}

fn main() {
    let path: PathBuf = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .expect("usage: streaming_api <input.mzML | input.d>");
    let config = KothConfig::default();

    tracer_style(&path, &config);
    uno_style(&path, &config);
}
