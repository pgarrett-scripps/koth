# koth-ms

The `koth-ms` crate (binaries `koth_ff` and `koth_align`): fast, memory-efficient peptide feature detection for LC–MS data, written in Rust.
Reads mzML (including streamed `.mzML.gz`), Bruker timsTOF `.d` with
Bruker-calibrated 1/K0, and Thermo `.raw` through a pure-Rust reader (no .NET),
all in the default build.

```bash
cargo install --locked koth-ms
koth_ff sample.mzML --output results/
```

Installs both `koth_ff` and `koth_align`. Feature finding writes chromatographic
traces, isotope features, a summary, and the resolved configuration. The
alignment tool produces multi-run intensity matrices; its q-values are
exploratory confidence measures, not calibrated peptide-identification FDR.

[Documentation](https://github.com/pgarrett-scripps/koth/tree/master/docs) ·
[Benchmarks](https://github.com/pgarrett-scripps/koth/blob/master/docs/BENCHMARKS.md) ·
[Source and issues](https://github.com/pgarrett-scripps/koth) ·
[Citation](https://github.com/pgarrett-scripps/koth/blob/master/CITATION.cff)

Licensed under MIT.
