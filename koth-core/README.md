# koth-core

The in-memory LC–MS feature detection stages of
[koth](https://github.com/pgarrett-scripps/koth): chromatographic hill
detection, isotope-envelope feature assembly with charge assignment, and
averagine isotope-pattern scoring.

`koth-core` reads and writes no files and has no native dependencies (no
SQLite, no vendor readers, no Arrow/Parquet). Bring your own spectra from your
own mzML, Thermo or Bruker reader. Its normal dependencies are rayon, serde,
toml, log, thiserror, rand and rustc-hash. The `koth-ms` crate (the `koth_ff`
and `koth_align` executables) adds the readers, writers and cross-run LFQ, and
runs exactly these stages.

```toml
[dependencies]
koth-core = "0.11"
```

```rust
use koth_core::{run_pipeline_from_spectra, KothConfig, Peak, PipelineOptions, Spectrum};

let spectra: Vec<Spectrum> = (0..20)
    .map(|i| Spectrum {
        scan_index: i,
        retention_time: i as f64 * 0.05,
        peaks: vec![Peak { mz: 500.0, intensity: 1.0e5, ion_mobility: 0.0 }],
        ms_level: 1,
        isolation_window: None,
        faims_cv: None,
    })
    .collect();

// The configuration is the same struct, and the same TOML, as `koth_ff -c`.
let config = KothConfig::default();
let out = run_pipeline_from_spectra(spectra.into_iter(), &config, &PipelineOptions::default())?;
println!("{} hills, {} features", out.hills.len(), out.features.len());
# Ok::<(), koth_core::Error>(())
```

For large runs, `run_pipeline_streaming_from_spectra` with a `PipelineSink`
hands over features one at a time instead of collecting them. The stages can
also be called one by one: `run_hills`, `run_features`, `run_scoring`.

[Documentation](https://docs.rs/koth-core) ·
[Algorithm](https://github.com/pgarrett-scripps/koth/blob/master/docs/ALGORITHM.md) ·
[Source and issues](https://github.com/pgarrett-scripps/koth)

Licensed under MIT.
