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
pub struct PipelineOptions {
    pub scoring:  bool,   // default true;  false == the binary's --no-scoring
    pub emit_ms2: bool,   // default false; true also detects DIA MS2 hills (mzML + Bruker diaPASEF .d + Thermo DIA .raw)
}

// Owned collect-all result (the uno shape).
pub struct FeatureFindingOutput {
    pub hills:     Vec<Hill>,
    pub features:  Vec<ScoredFeature>,
    pub ms2_hills: Vec<Hill>,   // populated only when emit_ms2 = true, else empty
}

// Streaming sink (the koth_tracer shape). All methods have no-op defaults.
pub trait PipelineSink {
    fn on_hills(&mut self, hills: Vec<Hill>) { … }         // once, whole finalized MS1 set, by value
    fn on_feature(&mut self, feature: ScoredFeature) { … } // once per feature, by value, in order
    fn on_ms2_hills(&mut self, hills: Vec<Hill>) { … }     // at most once, only if emit_ms2, after the MS1 side
}

// Group MS2 hills into DIA channels (sorted by window key; MS1 hills skipped).
pub fn group_ms2_hills_by_window(hills: Vec<Hill>) -> Vec<(IsolationWindow, Vec<Hill>)>;

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

// Combined MS1 + MS2 collect-all convenience (forces emit_ms2 = true).
pub fn run_pipeline_with_ms2(
    path: &Path, config: &KothConfig, opts: &PipelineOptions,
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

## MS2 (DIA fragment hills) — opt-in (mzML + Bruker diaPASEF `.d`)

By default the API is MS1-only and byte-for-byte the historic behavior. Set
`PipelineOptions.emit_ms2 = true` (or call `run_pipeline_with_ms2`) to *also*
detect MS2 fragment hills over the same input and deliver them:

* streaming: `PipelineSink::on_ms2_hills(Vec<Hill>)` fires **at most once**,
  after `on_hills` and every `on_feature`, with the whole finalized MS2 hill set
  by value;
* collect-all: `FeatureFindingOutput.ms2_hills` holds the same set.

MS2 detection is [`run_ms2_hills_streaming`] — the DIA hill detector that
partitions MS2 spectra by their **precursor isolation window** and runs an
independent per-window detector, so a hill never spans two windows. It is
independent of MS1 feature finding; the two share only the input file.

**MS2 hill settings are independently configurable.** The pipeline detects MS2
hills with `KothConfig::ms2_hills()`, which resolves to `[hills]` plus any
`[hills_ms2]` overrides (see `docs/CONFIGURATION.md` §3.2). So one `KothConfig`
drives MS1 and MS2 hill detection with *different* settings — e.g. a smaller
`min_scans` for short fragment traces — without touching the MS1 path. With no
`[hills_ms2]` present, `ms2_hills()` returns `[hills]` unchanged, so the MS2 side
is byte-identical to using `[hills]` directly. An in-process caller that builds a
`KothConfig` programmatically can set `config.hills_ms2 = Some(HillsMs2Overrides
{ min_scans: Some(2), ..Default::default() })` to the same effect; the low-level
`run_ms2_hills_streaming(path, &HillsConfig, …)` still takes an explicit
`HillsConfig` for callers that want to pass one directly.

### Window carriage and the precursor → window mapping

Every MS2 `Hill` carries its isolation window in `hill.isolation_window`
(`Some(IsolationWindow { target, lower, upper })` — MS1 hills are `None`). The
window boundaries are absolute m/z (already expanded from any target offset).
`koth_tracer` maps an MS1 precursor feature to the DIA window(s) that isolated it
by testing containment of the precursor m/z:

```rust
let out = run_pipeline_with_ms2(path, &config, &PipelineOptions::default())?;
let windows = group_ms2_hills_by_window(out.ms2_hills);   // Vec<(IsolationWindow, Vec<Hill>)>

for f in out.features.iter().filter(|f| f.feature.charge > 0) {
    let mz = f.feature.monoisotopic_mz();
    let fragments = windows.iter()
        .filter(|(w, _)| w.lower <= mz && mz <= w.upper)   // isolating window(s)
        .flat_map(|(_, hs)| hs.iter());
    // reconstruct the pseudo-DDA spectrum for this precursor from `fragments`.
}
```

`group_ms2_hills_by_window` is a convenience over the same information already on
each hill — grouping is deterministic (sorted by `IsolationWindow::key`) and
loses no hills. Wide/overlapping DIA schemes can legitimately map a precursor to
more than one window; the filter above handles that naturally.

### Supported inputs and graceful degradation

MS2 detection supports **mzML** (`.mzML`, `.mzML.gz`), **Bruker diaPASEF `.d`**
(`--features tdf`), and **Thermo DIA `.raw`** (`--features thermo`):

* **mzML** — one MS2 spectrum per precursor isolation window, as recorded.
* **Bruker diaPASEF `.d`** — the fixed isolation windows are read directly from
  each MS2 frame's quadrupole settings (`isolation_mz` ± `isolation_width/2`,
  segmented by scan range). Each window segment of a frame is collapsed down the
  ion-mobility axis with the **same** dnoise vertical-IM filter + watershed
  centroiding as the MS1 Bruker reader, and emitted as one MS2 spectrum stamped
  with that window. No conversion to mzML is needed. This runs only when the
  caller opts into MS2; the MS1 Bruker output is byte-for-byte unchanged.
* **Thermo DIA `.raw`** — each MS2 scan's precursor isolation window is read from
  the Thermo `RawFileReader` (`precursor.isolation_window()` → absolute
  `lower`/`target`/`upper` m/z, falling back to the precursor m/z as center) and
  one centroided MS2 spectrum is emitted per scan, reusing the **same** scan→peaks
  conversion as the MS1 `.raw` reader. Orbitrap scans are 1-D centroids (no ion
  mobility), so there is no per-scan segmentation — one MS2 scan → one spectrum.
  DIA is detected first (below); only a DIA verdict emits MS2. Runs only under
  `emit_ms2`; the MS1 `.raw` output is byte-for-byte unchanged.

Everything else degrades gracefully — `emit_ms2` yields an **empty** MS2 set plus
a warning, leaving the MS1 result completely unaffected rather than erroring:

* **Bruker ddaPASEF `.d`** — MS2 frames describe per-precursor selection, not
  fixed windows; precursor reconstruction is out of scope, so no MS2 is emitted.
  The mode is detected via `FrameReader::get_acquisition()` (DIAPASEF vs
  DDAPASEF/Unknown).
* **Thermo DDA `.raw`** — detected heuristically: a first, signal-free pass
  collects the MS2 isolation windows and computes their **mean recurrence**
  (`n_ms2 / distinct_windows`). A fixed DIA schedule revisits a small window set
  every cycle (recurrence >> 1); DDA selects data-dependent precursors that each
  appear ~once (recurrence ~= 1). DIA requires a small distinct-window count and
  recurrence >= 3, with a floor of 8 MS2 scans; otherwise the file is treated as
  DDA and no MS2 is emitted (no DDA precursor reconstruction, MS1 untouched).
  *Known failure modes:* targeted **PRM/tMS2** (fixed repeating inclusion list)
  is classified as DIA; **all-ion / MSe / bbCID** shows up as a 1-2 channel DIA;
  highly-multiplexed / scanning-quad DIA with near-unique per-scan centers would
  be missed. See `io::thermo::is_dia_schedule` for the full rationale.

(This mirrors the binary's `[file] ms2_hills_enabled` behavior.) The
`run_pipeline_streaming_from_spectra` MS2 path splits a mixed-level in-memory
stream internally (MS2 scans route to the MS2 detector and never pollute MS1).

## Parity with the file path

`pipeline.rs` calls the identical crate-level stages the binary calls
(`run_hills_streaming` / `detect_hills_from_iter` → `run_features` →
`run_scoring`), so the output is the same by construction. It is pinned by tests:

* `koth_ff/src/pipeline_tests.rs` (runs under default `cargo test -p koth_ff`):
  a synthetic multi-envelope run is put through both the manual staged path and
  the new API (collect-all *and* streaming sink), and the results are asserted
  byte-identical by serializing each through the **real** `write_hills_tsv` /
  `write_features_tsv` writers and comparing. Also covers `scoring: false`.
  The MS2 tests use a synthetic mixed-level DIA run: `emit_ms2` yields exactly
  the hills `detect_ms2_hills_from_iter` produces (byte-identical through
  `write_ms2_hills_tsv`); the MS1 side with `emit_ms2` on is byte-identical to a
  clean MS1-only run (fragment scans never pollute MS1); `on_ms2_hills` fires
  exactly once, after the features; and the window grouping / precursor→window
  mapping is exercised.
* `koth_ff/tests/pipeline_parity.rs` (`#[ignore]`, needs `--features tdf`):
  the same byte-identity check of `run_pipeline` vs the staged path on the real
  `tests/data/example_dda.d` Bruker fixture; `emit_ms2_degrades_gracefully_on_bruker_d`
  asserts a **ddaPASEF** `.d` yields no MS2 with the MS1 side unchanged; and
  `bruker_diapasef_ms2` (points at a diaPASEF `.d`, override with
  `KOTH_DIAPASEF_D`; skips if absent) asserts the fixed windows are recovered,
  every MS2 hill carries its window, and the MS1 side stays byte-identical. The
  `QuadrupoleSettings` → window-segment mapping itself
  (`io::bruker::inner::window_segments`) is unit-tested under `--features tdf`
  (half-open scan ranges, degenerate-row rejection) with no `.d` file needed.
* `koth_ff/tests/thermo_dia_ms2.rs` (`#[ignore]`, needs `--features thermo`, a
  .NET 8 runtime, and a DIA `.raw` via `KOTH_DIA_RAW`; skips if absent):
  `thermo_dia_raw_ms2` asserts a DIA `.raw` recovers multiple isolation windows,
  every MS2 hill carries its window, and the MS1 side is unchanged. The pure DIA
  logic — isolation-window derivation and the DIA-vs-DDA recurrence heuristic
  (`io::thermo::{derive_isolation_window, is_dia_schedule}`) — is unit-tested in
  `io/thermo_tests.rs` under `--features thermo` with **no `.raw` file** (so it
  runs in CI without a .NET runtime): fixed-schedule vs data-dependent streams,
  the recurrence threshold boundary, all-ion single-window, and too-few-scans.

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
