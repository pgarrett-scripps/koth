# Changelog

All notable changes to `koth_ff` (the Rust crate / CLI) are documented in
this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/),
and this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Release builds now use the published `dnoise` 0.1.0 crate instead of requiring
  a sibling source checkout. The minimum supported Rust version is 1.88 so the
  optional Thermo dependency graph and the default Bruker build share one MSRV.
- CI now checks formatting, Clippy, rustdoc, default and mzML-only builds across
  Linux, macOS, and Windows. Tagged releases verify the tag/version match, test
  and package the crate, then bundle both `koth_ff` and `koth_align` with SHA-256
  checksums for all supported targets.
- **The benchmark and the manuscript moved to their own repository**,
  [tacular-omics/koth-paper](https://github.com/tacular-omics/koth-paper), with
  their history. This repository is now the tool alone. Two consequences for
  anyone reading older entries below: paths of the form `benchmark/config/*.toml`
  and `benchmark/scripts/*` refer to that repository, and the build-time config
  parse test now covers only this repo's `example_config*.toml` rather than the
  tuned per-platform configs, which are validated where they are run.
- **Feature detection now defaults to the exhaustive non-destructive assembler
  with a seed-anchored isotope cosine** (previously a greedy resolver with an
  adjacent-anchored cosine). The exhaustive pool claims contested hills
  longest-envelope-first and truncates a partly-claimed candidate to its free
  prefix rather than dropping it; anchoring each isotope's chromatographic
  cosine to the monoisotope seed (as biosaur2/AlphaPept/Dinosaur do) recovers
  M+2/M+3 isotopes. On the 20-run PXD003881 Orbitrap cohort this lifts PSM
  recall 79.3 → 79.9 % with no quant regression (median CV/MV/FPR
  flat-to-better); on the 18-run timsTOF cohort it is a small regression (recall
  79.4 → 79.0 %, median CV 10.7 → 11.1 %) — accepted to keep a single shipped
  code path.
- **Exhaustive-assembler candidate generation is parallelized** over seeds (like
  the former greedy path), removing the ~2.5× slowdown of the initial serial
  implementation: **37.0 → 14.9 s per PXD003881 run** at
  `RAYON_NUM_THREADS = 4`, matching the former greedy speed. Output byte-identical.
- **Config surface slimmed (breaking for existing TOMLs).** Many dead or
  superseded experimental flags were removed. Because configs are parsed with
  `deny_unknown_fields`, a TOML that still sets a removed key will now error —
  see *Removed / renamed*. Feature output for the shipped configs is unchanged.
- **`koth_align` now streams hills one run at a time during LFQ** instead of
  loading every run's hills into memory up front. Peak memory is O(one run's
  hills) rather than O(all runs' hills), removing the main memory bottleneck on
  large cohorts. Only per-run *features* (small) are held for the whole run;
  each run's *hills* (the large intensity profiles) are loaded, quantified, and
  dropped before the next run. The `lfq::quantify` library API gains a
  `load_hills: impl Fn(usize) -> Vec<Hill>` parameter that fetches a run's hills
  on demand; `build_grid`/`SortedHills` now take hills directly rather than a
  `RunInput`.

- **`koth_align` MBR FDR uses a semi-supervised QDA rescorer.** Instead of ranking every cell by `hybrid_score`,
  it learns a quadratic discriminant over five symmetric per-cell features —
  `|ppm error|`, `|RT diff|`, spectral Bhattacharyya, isotope co-elution, and
  `|IM delta|` (inert on Orbitrap, active on timsTOF) — via a Percolator-style
  loop (iterative confident-positive selection, decoys as negatives, 3-fold
  cross-validation by feature). Every cell — detected and match-between-runs
  alike — competes against the decoys, so a background-contaminated cell in a
  depleted well is gatable regardless of how it was populated (required for
  fold-change rescue on large-dynamic-range designs). Improves target/decoy
  separation over the raw hybrid (MBR-vs-decoy AUROC 0.84 → 0.86). Fully
  deterministic; it is now the only MBR rescorer.
  NOTE: absolute q-value calibration is not yet independently validated.
- **Isotope scoring is no longer a config choice.** A cosine and a Bhattacharyya
  were historically conflated under one flag; the two *orthogonal* isotope
  signals are now always computed and named for what they are:
  - `spectral_bhattacharyya` — observed vs theoretical isotope *pattern* match
    (penalises expected-but-missing peaks); the term the hybrid uses.
  - `coelution` — cosine between the matched isotopes' XIC traces across RT
    (do they co-elute?).

  `hybrid_score = (rt × intensity × bhattacharyya × coelution)^¼`.
- `koth_align`: `ppm_error` (and `im_delta`) in `lfq_details.tsv` are now
  computed against the per-row grid centre rather than the target consensus
  m/z. Previously, decoy rows reported the constant decoy m/z shift in ppm
  (~8000 ppm on Orbitrap data) instead of a real measurement error, which
  broke target/decoy comparability for downstream rescorers.

### Added
- **Native Thermo Fisher `.raw` input** behind the opt-in `thermo` cargo
  feature (`cargo build -p koth_ff --features thermo`). Wraps Thermo's
  `RawFileReader` via a self-hosted .NET 8 runtime (required at build and run
  time; auto-detects `DOTNET_ROOT` under `~/.dotnet`, `/usr/share/dotnet`, …),
  so a `.raw` file can be fed to `koth_ff` directly with no mzML conversion.
  Off by default; `.raw` scans are read as MS1 centroids, matching the mzML
  path. A default (non-`thermo`) build gives a clear "rebuild with
  `--features thermo`" error for `.raw` input.
- **ID-free isotope-consistency m/z recalibration** (`[file] mz_recalibration`,
  default off; CLI `--recalibrate`). A pass-1 feature detection collects the
  signed ppm deviation of every adjacent isotope-hill spacing from its
  theoretical `neutron_mass / z` step and bins the residuals over (m/z, RT).
  Pass 2 shifts the *expected* isotope position during chain extension by the
  learned per-region median offset (hierarchical fallback: cell → m/z-marginal
  → global), so isotope hills are searched at their recalibrated location.
  Corrects the proportional, m/z-/RT-dependent mass-error component (dominant
  Orbitrap mode) without IDs or a lock mass; a constant Da offset stays
  unobservable by design. Deterministic (robust medians, fixed binning),
  skipped in decoy mode, and costs one extra cheap feature-assembly pass over
  the same hills. Tunable via `mz_recalibration_mz_bins` (20),
  `mz_recalibration_rt_bins` (8), `mz_recalibration_min_samples` (50).
  Inspired by Biosaur's per-isotope "smart" calibration and AlphaPept's
  multi-dimensional recalibration, adapted to koth_ff's ID-free feature stage.
  Benchmarked as a no-op on its own — well-calibrated instruments have little
  proportional error left to correct — but it is the scaffolding for the
  region-adaptive tolerance below, which *is* a win and is now applied together
  with it.
- **Region-adaptive isotope-match tolerance** — enabled automatically whenever
  `[file] mz_recalibration` is on (the former separate
  `mz_recalibration_adaptive_tol` flag and the `--adaptive-tol` / `--kish` CLI
  flags have been removed; see *Removed / renamed*). Replaces the fixed
  isotope-match ppm window with
  `clamp(tol_sigma_mult × σ(m/z, RT), tol_floor_ppm, mz_tolerance)`, where σ
  comes from the recalibration surface's per-region residual spread —
  tightening the search where the instrument is precise (rejecting false
  isotope matches) and relaxing it, up to the configured ceiling, where it
  isn't. **Validated on both benchmark platforms and shipped in
  `benchmark/config/koth_ff.toml` and `koth_ff_bruker.toml`:** on PXD003881
  (20-run Orbitrap cohort) it drops features −5.8% (spurious matches at the
  over-wide fixed 8 ppm window) at flat PSM recall, while median CV improves
  23.76→23.25%, MV 0.615→0.568%, CV@q≤0.05 20.64→20.47%, +70 complete-quant
  features @q≤0.05; on the 18-run timsTOF 15-min cohort (fixed 15 ppm vs a
  learned real spread of σ≈2.8 ppm) it improves PSM recall +1.06 pp
  (79.34→80.40%) *and* quant (median CV 13.89→13.76%, MV 3.12→3.08%,
  CV@q≤0.05 11.74→11.69%, +285 complete-quant features), with −0.2% features.
  No regression found on either platform or metric. Tunable via
  `mz_recalibration_tol_sigma_mult` (3.0) and `mz_recalibration_tol_floor_ppm`
  (1.0).
- `[lfq] rt_spread_scoring` (default `false`, experimental) — replaces the raw
  RT-closeness term with a σ-normalised Gaussian likelihood using the per-run
  post-warp RT-residual spread (region-aware RT scoring).
- New per-cell columns in `lfq_details.tsv`: `spectral_bhattacharyya`,
  `coelution`, `rt_score`, `int_score` — the individual hybrid components,
  exposed for downstream rescoring.
- `[lfq] decoy_mz_shift_da` and `[lfq] decoy_rt_shift_pct` config knobs that
  control the decoy grid offsets (previously hard-coded to 11 Da/charge and
  1 % of the RT span). Defaults preserve previous behaviour.

### Fixed
- Non-deterministic MS2 hill IDs: `detect_ms2_hills_from_iter` collected hills
  in `HashMap` iteration order before numbering them, so identical input
  produced different `hill_id`s across runs. Detectors are now sorted by
  isolation-window key first.
- NaN-safe sorts: several `partial_cmp(..).unwrap()` sorts (mzML/Bruker RT and
  peak m/z, output intensity, LFQ grid and feature m/z) could panic on a NaN
  value; all now fall back to `Ordering::Equal`.
- LFQ guards against a degenerate grid (`grid_cols`/`n_isotopes` = 0) that could
  underflow and panic; clamps to ≥ 1 with a warning.
- Scoring diagnostics that were printed unconditionally to stderr are now gated
  behind `--log-level debug`.

### Removed / renamed
- **The greedy feature assembler** and the `[features] exhaustive_assembly`
  flag. The exhaustive non-destructive assembler (now the default, see
  *Changed*) is the only feature-detection path.
- **Two superseded m/z-calibration mechanisms**, both fully replaced by
  `mz_recalibration`: the Kish mass-uncertainty model (`[file]
  mz_uncertainty_mode` + `mz_uncertainty_sigma_mult`) and the empirical
  adaptive-tolerance pass (`[file] adaptive_mz_tolerance` +
  `adaptive_mz_tolerance_pass1_multiplier` + `_sigma_mult`, CLI `--adaptive`).
  `mz_recalibration_adaptive_tol` is folded into `mz_recalibration` — a single
  `mz_recalibration = true` now enables recalibration *with* region-adaptive
  isotope-match tolerance (CLI `--adaptive-tol` / `--kish` removed).
- **The prominence hill splitter**: `[hills] split_algo`, `min_peak_distance`,
  `min_peak_height`, `min_prominence`. The persistence splitter is now the only
  splitter. The `benchmark/config/koth_ff_persist.toml` preset (which only
  selected persistence) is removed.
- **The legacy piecewise / linear RT-warp models** and their sigma-clip knobs
  (`[alignment] rt_warp_sigma_clip`, `rt_warp_clip_iters`). RANSAC is the only
  warp.
- `[features] cosine_intersection` (a documented no-op), `[features]
  cosine_hybrid_depth` + the `hybrid` cosine-anchor mode, `[file]
  n_most_abundant`, and the `[lfq] tdc_method` toggle (the QDA MBR rescorer is
  the only path; there is no `hybrid` fallback).
- `[scoring] isotope_offset_min` / `isotope_offset_max` collapsed to a single
  `[scoring] isotope_offset_enabled` bool (default `false`; `true` searches
  ±1 neutron offsets).
- Internal: the unreachable low-level adaptive-tolerance calibration primitives
  (`hills::calibration::MzDeltaHistogram` and the `HillDetector` calibration
  recording hooks), left over from the removed `adaptive_mz_tolerance` pass.
- `[lfq] spectral_cosine` output column removed (it held a Bhattacharyya value —
  the source of much confusion); replaced by the correctly named
  `spectral_bhattacharyya` and `coelution` columns.
- `[lfq] spectral_bhattacharyya` and `[lfq] spectral_coelution` *flags* removed
  — both metrics are always computed now. Leftover values in an existing TOML
  are silently ignored.
- `[lfq] spectral_cosine_min` renamed to `min_spectral_bhattacharyya` (it gates
  peak expansion on the Bhattacharyya score, never a cosine). The old name is
  ignored; the default (0.1) is unchanged.
- `koth-ff` PyO3 Python bindings (`koth_ff_py` crate) and the wheel-publish
  workflow. Use the `koth_ff` CLI instead.

## [0.2.0] — 2026-09-03

### Added
- Multi-CV FAIMS support for MS1 feature detection. mzML and native Thermo
  readers carry compensation voltage into independent, channel-local hill and
  feature pipelines; the feature TSV and Parquet schemas expose a nullable
  `FAIMS` column. Non-FAIMS inputs retain the existing single-channel detection
  path.

## [0.1.0] — 2026-04-07

Initial release.

### Added
- `koth_ff` CLI for streaming hill detection, isotope-feature detection, and
  averagine scoring on mzML files and Bruker timsTOF `.d` directories.
