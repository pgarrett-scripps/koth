# koth

[![crates.io](https://img.shields.io/crates/v/koth-ms.svg)](https://crates.io/crates/koth-ms) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22931619.svg)](https://doi.org/10.5281/zenodo.22931619)

**Fast, memory-efficient peptide feature detection for LC–MS data.**
koth reads mzML (including `.mzML.gz`), Bruker timsTOF `.d`, and Thermo
`.raw` files and reports isotope envelopes, charge states, retention times,
intensities, and, for timsTOF, Bruker-calibrated ion mobility (1/K0). Thermo
`.raw` is read natively by a pure-Rust reader, with no .NET runtime.

## Why koth?

On the paper's Orbitrap and timsTOF benchmarks, koth used less time and peak
memory than the compared feature finders at the tested settings.

![Per-run runtime and peak memory for koth and comparator tools on Orbitrap and timsTOF benchmarks](docs/assets/runtime-memory.png)

Lower is better on both axes. Points show individual runs and outlined points
show tool means: 20 Orbitrap runs and 18 timsTOF runs on the same workstation.
See [benchmark details and limitations](docs/BENCHMARKS.md)
for settings, provenance, and identified-peptide recall.

## Install

Download both executables from [GitHub Releases](https://github.com/pgarrett-scripps/koth/releases),
or build from source with Rust 1.88 or newer:

```bash
git clone https://github.com/pgarrett-scripps/koth.git
cd koth
cargo install --locked --path koth_ff
```

See [installation](docs/INSTALLATION.md)
for platform builds, crates.io installation after the first crate release,
and building without the Bruker or Thermo readers.

## Run

```bash
koth_ff sample.mzML --output results/
# Or: koth_ff sample.d --output results/     (Bruker timsTOF)
# Or: koth_ff sample.raw --output results/   (Thermo)
```

Results are written to `results/sample/`: chromatographic traces in
`hills.tsv`, isotope features in `features.tsv`, a `report.json` summary,
and the resolved `config.toml`. Parquet output is also available.

The included `koth_align` tool aligns multiple runs and produces an intensity
matrix (apex intensity per cell by default), optionally guided by Sage or other
peptide identifications. Its group and cell q-values are exploratory confidence measures,
not validated peptide-identification FDR; see the
[LFQ guide](docs/LFQ.md).

## Use as a library

Two crates, one version. `koth-core` is the feature-detection algorithm over
spectra you supply: no file readers, no SQLite, no native code. `koth-ms` adds
the mzML, Bruker and Thermo readers, output writers, alignment and LFQ, and
re-exports `koth-core`.

```toml
koth-core = "0.10"                                   # bring your own spectra
koth-ms = { version = "0.10", default-features = false, features = ["tdf", "thermo"] }  # readers, no CLI
```

```rust
let out = koth_core::run_pipeline_from_spectra(spectra, &config, &Default::default())?;
```

See [development](docs/DEVELOPMENT.md#library-api) and the
[streaming API](docs/STREAMING_API.md).

## Documentation

[Usage](docs/USAGE.md) ·
[Configuration](docs/CONFIGURATION.md) ·
[Outputs](docs/OUTPUTS.md) ·
[Algorithm](docs/ALGORITHM.md) ·
[All documentation](docs/README.md)

## Cite and contribute

Use [CITATION.cff](CITATION.cff)
to cite the software. Release archives use the accompanying
[Zenodo metadata](.zenodo.json).
The [manuscript and benchmark source](https://github.com/pgarrett-scripps/koth-paper)
are maintained separately and may remain private during preparation.

[Report an issue](https://github.com/pgarrett-scripps/koth/issues) ·
[Contributing](CONTRIBUTING.md) ·
[MIT license](LICENSE)
