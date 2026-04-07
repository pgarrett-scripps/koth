# Changelog

All notable changes to `koth-ff` (the Python package) and `koth_ff` (the
Rust crate / CLI) are documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/),
and this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] — 2026-04-07

Initial release.

### Added
- `koth_ff` CLI for streaming hill detection, isotope-feature detection, and
  averagine scoring on mzML files and Bruker timsTOF `.d` directories.
- `koth-ff` Python package with PyO3 bindings (`detect_hills`,
  `detect_features`, `run_pipeline`) returning Polars DataFrames.
- abi3 wheel (`cp310-abi3`) covering Python 3.10+ on Linux x86_64/aarch64,
  macOS x86_64/aarch64, and Windows x86_64.
