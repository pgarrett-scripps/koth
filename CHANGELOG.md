# Changelog

All notable changes to `koth_ff` (the Rust crate / CLI) are documented in
this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/),
and this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- `koth_align`: `ppm_error` (and `im_delta`) in `lfq_details.tsv` are now
  computed against the per-row grid centre rather than the target consensus
  m/z. Previously, decoy rows reported the constant decoy m/z shift in ppm
  (~8000 ppm on Orbitrap data) instead of a real measurement error, which
  broke target/decoy comparability for downstream rescorers.

### Added
- `[lfq] decoy_mz_shift_da` and `[lfq] decoy_rt_shift_pct` config knobs that
  control the decoy grid offsets (previously hard-coded to 11 Da/charge and
  1 % of the RT span). Defaults preserve previous behaviour.

### Removed
- `koth-ff` PyO3 Python bindings (`koth_ff_py` crate) and the wheel-publish
  workflow. Use the `koth_ff` CLI instead.

## [0.1.0] — 2026-04-07

Initial release.

### Added
- `koth_ff` CLI for streaming hill detection, isotope-feature detection, and
  averagine scoring on mzML files and Bruker timsTOF `.d` directories.
