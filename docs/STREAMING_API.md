# In-process / streaming feature-finding API

`koth_ff` can hand its hills and isotope-chain features to another Rust tool
**directly in memory**, with no intermediate `hills.parquet` / `features.parquet`
written and re-read. This is an *additional* library seam over the same staged
functions the `koth_ff` binary uses — the CLI, its file output, and its byte
format are unchanged.

Everything lives in `koth_ff/src/pipeline.rs` and is re-exported at the crate
root (`koth_ff::run_pipeline`, `koth_ff::PipelineSink`, …). The in-memory result
types (`Hill`, `Feature`, `ScoredFeature`, `Spectrum`) are re-exported there too.

## API surface

```rust
// Options mirroring the binary's optional stages.
pub struct PipelineOptions { pub scoring: bool }   // default: scoring = true
                                                   // false == the binary's --no-scoring

// Owned collect-all result (the uno shape).
pub struct FeatureFindingOutput {
    pub hills:    Vec<Hill>,
    pub features: Vec<ScoredFeature>,
}

// Streaming sink (the koth_tracer shape). Both methods have no-op defaults.
pub trait PipelineSink {
    fn on_hills(&mut self, hills: Vec<Hill>) { … }        // once, whole finalized set, by value
    fn on_feature(&mut self, feature: ScoredFeature) { … } // once per feature, by value, in order
}

// --- streaming entry points (no files written) ---
pub fn run_pipeline_streaming<S: PipelineSink>(
    path: &Path, config: &KothConfig, opts: &PipelineOptions, sink: &mut S,
) -> Result<(), KothError>;

pub fn run_pipeline_streaming_from_spectra<I, S>(
    spectra: I, config: &KothConfig, opts: &PipelineOptions, sink: &mut S,
) -> Result<(), KothError>
where I: Iterator<Item = Spectrum>, S: PipelineSink;

// --- collect-all convenience (no files written) ---
pub fn run_pipeline(
    path: &Path, config: &KothConfig, opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError>;

pub fn run_pipeline_from_spectra<I: Iterator<Item = Spectrum>>(
    spectra: I, config: &KothConfig, opts: &PipelineOptions,
) -> Result<FeatureFindingOutput, KothError>;
```

Configuration is always passed in explicitly as a `&KothConfig` — the same
struct the binary loads from TOML.

## Streaming granularity, and why it is what it is

The pipeline has a hard data dependency that fixes the honest granularity:

* **Isotope-chain assembly needs the complete hill set.** `detect_features`
  builds an m/z-sorted index over *every* hill and grows chains across it, so no
  feature can be finalized before hill detection ends.
* **A hill is not finalized until end-of-read.** Co-elution *splitting*,
  baseline *filtering*, and stable *`hill_id` assignment* are global post-passes
  run after the last spectrum. There is no point mid-run at which a hill's final
  identity (its ID, its split boundaries) is known.

So the API does **not** pretend to emit finalized hills per RT-window. Instead:

1. Spectra stream in one at a time during detection — peak memory is
   `O(active hills)`, not `O(total peaks)`. (The mzML reader already does this,
   with a prefetch reader thread; Bruker `.d` / Thermo `.raw` batch-load because
   their readers expose no streaming API.)
2. `on_hills(Vec<Hill>)` fires **once**, with the whole finalized hill set, by
   value.
3. `on_feature(ScoredFeature)` fires **once per feature**, by value, in
   detection order, and each feature is freed as the consumer returns — so the
   consumer never holds the whole feature vector.

The feature order and set handed to the sink are exactly what the file writer
receives *before* its presentation step (which additionally drops `charge == 0`
and sorts by descending intensity). That is what makes parity provable.

## Parity with the file path

`pipeline.rs` calls the identical crate-level stages the binary calls
(`run_hills_streaming` / `detect_hills_from_iter` → `run_features` →
`run_scoring`), so the output is the same by construction. It is pinned by tests:

* `koth_ff/src/pipeline_tests.rs` (runs under default `cargo test -p koth_ff`):
  a synthetic multi-envelope run is put through both the manual staged path and
  the new API (collect-all *and* streaming sink), and the results are asserted
  byte-identical by serializing each through the **real** `write_hills_tsv` /
  `write_features_tsv` writers and comparing. Also covers `scoring: false`.
* `koth_ff/tests/pipeline_parity.rs` (`#[ignore]`, needs `--features tdf`):
  the same byte-identity check of `run_pipeline` vs the staged path on the real
  `tests/data/example_dda.d` Bruker fixture.

## How the two consumers plug in

### uno (MS1 search) — collect-all

```rust
use koth_ff::{config::KothConfig, run_pipeline, PipelineOptions};

let out = run_pipeline(path, &config, &PipelineOptions::default())?;
// out.features: Vec<ScoredFeature> with charge + averagine score + intensities.
let searchable: Vec<_> = out.features.iter()
    .filter(|f| f.feature.charge > 0)     // == the features-file row set
    .collect();
// index `searchable` by neutral mass / charge / RT and search.
```

### koth_tracer (DIA -> pseudo-DDA) — streaming sink

```rust
use koth_ff::{run_pipeline_streaming, Hill, PipelineOptions, PipelineSink, ScoredFeature};

#[derive(Default)]
struct TracerSink { hills: Vec<Hill> /* + trace index */ }

impl PipelineSink for TracerSink {
    fn on_hills(&mut self, hills: Vec<Hill>) {
        self.hills = hills;             // build the m/z x RT x IM trace index once
    }
    fn on_feature(&mut self, feature: ScoredFeature) {
        // attribute fragment hills to this precursor feature, then drop it —
        // the whole feature vector is never held on the consumer side.
    }
}

let mut sink = TracerSink::default();
run_pipeline_streaming(path, &config, &PipelineOptions::default(), &mut sink)?;
```

A runnable version of both is in `koth_ff/examples/streaming_api.rs`.
