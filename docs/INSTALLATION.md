# Installation

The project requires Rust 1.88 or newer. To install both command-line tools from
a source checkout (clone it first):

```bash
git clone https://github.com/pgarrett-scripps/koth.git
cd koth
cargo install --locked --path koth_ff
```

Published GitHub releases produce archives containing both `koth_ff` and
`koth_align`, plus SHA-256 checksum files, for Linux x86_64, macOS x86_64 and
arm64, and Windows x86_64.

For a development build:

```bash
cargo build --release
# binaries are at target/release/koth_ff and target/release/koth_align
```

Both binaries are built together with the command above. Two input readers are
Cargo features, and both are on by default:

| Feature | Input | Dependencies |
|---|---|---|
| `tdf` | Bruker timsTOF `.d` | `timsrust`, `dnoise`, bundled SQLite |
| `thermo` | Thermo Fisher `.raw` | `opentfraw` (pure Rust) |

mzML and gzip-compressed `.mzML.gz` are always supported; `.mzML.gz` is
decompressed as it is parsed, not into memory first.

To build without one or both readers:

```bash
cargo build --release -p koth-ms --no-default-features --features tdf   # no Thermo reader
cargo build --release -p koth-ms --no-default-features                  # mzML only
```

### Native Thermo `.raw` input

Thermo `.raw` files are read directly, with no prior mzML conversion, by the
`thermo` feature. It uses [`opentfraw`](https://crates.io/crates/opentfraw), a
pure-Rust parser of the `.raw` format: no .NET runtime, vendor assemblies, or
system libraries are needed at build or run time. MS1 peaks are the
instrument's stored centroids; profile-only scans are centroided with mzdata's
peak picker, as for profile mzML.

Known limitation: on some Exploris 480 and Fusion Lumos DIA files the
isolation-window center cannot be recovered from the `.raw` (opentfraw issue
#44), so DIA MS2 detection on such a file emits no MS2 hills. MS1 is
unaffected. Convert to mzML for MS2 on those files.

## Install from crates.io

Once the first crate release has been published:

```bash
cargo install --locked koth-ms
```

This installs both `koth_ff` and `koth_align` with the default Bruker and
Thermo readers. Prebuilt release archives are built with the same default
features.

## Verify installation

```bash
koth_ff --version
koth_align --version
koth_ff --help
```

The version output includes the source commit when built from Git. Registry
and source-archive builds may report `unknown` for the commit.
