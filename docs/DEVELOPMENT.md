# Development

## Recipes

```bash
just run-mzml data.mzML          # release build + run, output to ./out
just run-bruker data.d           # same for Bruker input
just run-config data.mzML cfg.toml  # with explicit config
just run-debug data.mzML         # debug log level
just run-no-score data.mzML      # skip scoring stage
just peek-hills                  # head -3 on hills.tsv
just peek-features               # head -3 on features.tsv
just count                       # row counts for both outputs
just test                        # cargo test
just lint                        # clippy -D warnings
```

## Library API

### Which crate

| Crate | Holds | Depend on it when |
|---|---|---|
| `koth-core` | hill detection, isotope features, averagine scoring, `KothConfig` (the same TOML as `koth_ff -c`), the in-memory pipeline | you have your own spectra and want no file readers, no SQLite and no native code (for example next to another `libsqlite3-sys`) |
| `koth-ms` | everything in `koth-core`, re-exported at the same `koth_ms::` paths, plus mzML/Bruker/Thermo readers, writers, alignment/LFQ and the executables | you want koth to read files or run LFQ |

`koth-core` from your own spectra:

```toml
[dependencies]
koth-core = "0.10"
```

```rust
use koth_core::{run_pipeline_from_spectra, KothConfig, PipelineOptions};

let config: KothConfig = KothConfig::from_toml(Path::new("koth_ff.toml"))?;
let out = run_pipeline_from_spectra(my_spectra.into_iter(), &config, &PipelineOptions::default())?;
```

`koth-ms` as a library without the executables:

```toml
[dependencies]
koth-ms = { version = "0.10", default-features = false, features = ["tdf", "thermo"] }
```

The `cli` feature (default on) builds `koth_ff` and `koth_align` and pulls in
`clap` and `env_logger`; library users can leave it off. `koth-core` functions
return `koth_core::Error`, which converts into `koth_ms::error::KothError` with
`?`.

### Single-run pipeline

`koth-ms` is also a library crate. The top-level functions mirror the CLI stages:

```rust
use std::path::Path;
use koth_ms::{run_hills_streaming, run_features, run_scoring, config::KothConfig};

let config = KothConfig::default();
let input = Path::new("data.mzML");

// Stage 1 — streams spectra under the default configuration
let hills = run_hills_streaming(input, &config.hills, &config.file)?;

// Stage 2
let features = run_features(&hills, &config.features, &config.file)?;

// Stage 3 — score using the configured isotope model and polarity
let scored = run_scoring(
    &features, &config.scoring, &config.features, config.file.polarity,
);
```

For non-streaming use (e.g. when you already have spectra in memory):

```rust
let spectra = koth_ms::read_spectra(input, &config.file)?;
let hills = koth_ms::run_hills(&spectra, &config.hills, &config.file);
```

### Alignment and LFQ

```rust
use koth_ms::{
    alignment::{align_runs, AlignmentConfig, RunInput},
    lfq::{quantify, LfqConfig},
};

// Build one RunInput per LC-MS run
let runs: Vec<RunInput> = vec![
    RunInput { name: "sample1".into(), features: scored1, hills: hills1, scan_times: vec![], rt_bounds: None },
    RunInput { name: "sample2".into(), features: scored2, hills: hills2, scan_times: vec![], rt_bounds: None },
];

// Stage 4: alignment
let alignment_config = AlignmentConfig::default();
let alignment = align_runs(&runs, &alignment_config);

// Stage 5: LFQ. Hills are fetched per run via the loader and dropped before
// the next run, so only one run's hills are resident at a time. Here they are
// already in memory; to stream from disk, read each run's hills in the closure.
let lfq_config = LfqConfig::default();
let matrix = quantify(&runs, &alignment, &lfq_config, |i| runs[i].hills.clone());

// Access results
println!("Reference run: {}", matrix.reference_run);
for feat in 0..matrix.n_features {
    for run in 0..matrix.n_runs {
        let intensity = matrix.intensity(feat, run);
        let qvalue   = matrix.q_value(feat, run);
    }
}
```

`scan_times` is a `Vec<f64>` mapping absolute scan index to retention time in minutes. Pass an
empty Vec to fall back to linear interpolation between each hill's `rt_start`/`rt_end` —
sufficient for most data. `rt_bounds` is only needed when hills are streamed and features
are absent.

Loading runs from `koth_ff` output files:

```rust
use koth_ms::input::{discover_runs, read_hills, read_features};

let run_paths = discover_runs(Path::new("batch/"))?;
for rp in &run_paths {
    let hills    = read_hills(&rp.hills_path)?;
    let features = read_features(&rp.features_path)?;
    // build RunInput ...
}
```

See also the [streaming API guide](STREAMING_API.md),
[contribution guide](../CONTRIBUTING.md), and [release process](../RELEASE.md).
