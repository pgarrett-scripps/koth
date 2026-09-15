# Changelog

All notable changes to `koth_ff` (the Rust crate / CLI) are documented in
this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/),
and this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

_Nothing yet._

## [0.7.0] — 2026-09-14

### Feature finding and input streaming
- Rank isotope-chain candidates with tuned additive evidence and retain the
  best-scoring prefix. Per-run features must be regenerated for this release.
- Stream native Bruker and optional Thermo MS1 inputs through bounded buffers.
  Reader parity and error propagation have regression coverage.

### Release compatibility
- Preserve the 0.6.0 extraction defaults: RT half-window 0.005 of the observed
  run span and absolute IM half-window 0.015, independent of group limits.
- Regenerate both detection and LFQ outputs. Old consensus admission settings
  are rejected; migrate to `max_group_qvalue`. Long bundles use schema version 3.
- Group and cell q-values remain separate exploratory confidence measures;
  neither establishes calibrated peptide identification FDR. Known ownership
  and isotope-choice limitations remain documented below.

### Isotope-choice validation
- Add a calculated-peptide-mass audit of isotope competitors across both full
  cohorts. A shared-envelope selector worsened Bruker isotope choices despite
  nearly unchanged CV and was removed; the prior extraction behavior remains.
  See [the rejected trial](docs/lfq-isotope-choice.md) for results and provenance.

### Exclusive LFQ extraction
- Track original hill-profile samples and reserve native peak segments so
  competing groups cannot quantify the same signal or adjacent slices of it.
  Rank alternatives across runs, re-extract residuals, and mark unresolved
  residuals ambiguous. Targets and decoys use independent masks with identical
  rules; cell confidence is recomputed afterward.
- Add ownership provenance to LFQ details and long bundle schema version 3,
  synthetic regression tests, and full-cohort audits. Original group membership
  and permutation q-values remain unchanged. See
  [signal ownership](docs/lfq-signal-ownership.md) for results and limits.
  Full-cohort validation finds remaining isotope-choice errors and an Orbitrap
  precision regression. These remain documented limitations of the adopted
  signal-ownership approach, not claims of automatic isotope correction.

### LFQ group confidence
- Replace member/seed score floors, `min_group_size`, and
  `allow_replicated_weak_seeds` with `max_group_qvalue`. Continuous quality,
  cross-run coordinate agreement and ambiguity rank all valid candidates against
  ten independent RT-permutation control cohorts; at least two original runs
  supply support. Controls preserve RT density within run/charge/quality strata.
- Export separate group evidence and q-values in the consensus table and long
  bundle schema version 2. Per-cell extraction confidence remains independent.
- Add candidate auditing, CLI/export checks and an 18-run exploratory comparison.
  The group statistic is not a validated FDR guarantee, and the initial results
  show a coverage/quality tradeoff. The [permutation follow-up](docs/lfq-group-permutation.md)
  retains the original scorer and replaces unstable whole-run shifts; includes
  2-, 6-, 18- and 20-run checks, quantitative comparisons and member-link audits.
- The [integrity audit](docs/lfq-group-integrity.md) identifies unresolved cross-group
  splits, isotope aliases and shared-signal reporting; adds diagnostic tests and
  reproducible audits. Exclusive extraction addresses demonstrated signal reuse;
  cross-group membership and isotope identity remain unresolved in some cases.

## [0.6.0] — 2026-09-13

### Changed
- The default LFQ extraction RT half-window is 0.005 of each run's observed
  RT span (±0.5%), centred on the alignment-predicted native retention time.
- The default LFQ extraction ion-mobility half-window is 0.015 absolute 1/K0
  units. This setting is inert on inputs without ion mobility.
- Configuration documentation distinguishes extraction windows from the
  independent alignment and consensus-grouping tolerances.

### Compatibility and validation
- Explicit extraction settings continue to override the defaults. Per-run
  feature detection, alignment, consensus grouping, decoy offsets, and output
  schemas are unchanged. Regenerate LFQ matrices when adopting these defaults.
- The default extraction windows are covered by configuration-loading tests,
  including independence from consensus limits and explicit-setting overrides.

## [0.5.0] — 2026-09-13

### Added
- Optional long LFQ bundle (`[output] export_long = true`) containing original
  feature measurements, aligned coordinates, every original group member and
  alternative, recomputed LFQ cells, and a versioned provenance manifest.
  Original feature quality remains separate from extraction quality and q-values.
- Independent full-span mass, RT and IM limits in `[lfq.consensus]`, plus an
  optional fixed `[alignment] reference_run` for reproducible comparisons.
- Default-off `allow_replicated_weak_seeds` retains bounded groups supported by
  multiple original runs below the seed-quality floor. Other quality limits apply.

### Fixed
- Consensus grouping now uses deterministic quality order and bounded full-group
  spans. Minimum size and seed-quality filters run after membership is complete.
- Each run's primary observation is chosen by its own quality; alternatives are
  preserved. Missing ion mobility cannot bridge incompatible measured values.
- Long output reads original isotope counts and other feature descriptors from
  source tables, preserving measurements lost in reconstructed alignment hills.
- Synthetic extraction decoys have no original observation link or inherited
  target q-value. Q-values are unavailable when target/decoy scoring is disabled.

### Compatibility and validation
- Per-run feature-finding defaults are unchanged. Consensus grouping behavior
  changes, so LFQ matrices must be regenerated for version-matched comparisons.
- Wide matrices remain supported; long export and weak-group retention are opt-in.
- Extraction q-values assess extraction confidence, not peptide or protein
  identity. Absolute QDA q-value calibration remains an explicitly accepted
  limitation; no new absolute-calibration claim is made in this release.

## [0.4.0] — 2026-09-11

### Added
- **Non-peptide isotope models.** `[features] isotope_model` selects which
  analyte class's average composition the theoretical isotope pattern is built
  from: `"peptide"` (the default, Senko's averagine, and the only model the
  published benchmark exercises), `"rna"`, `"dna"`, or an explicit table
  (`{ residue_mass = 321.2916, c = 9.5, h = 11.75, n = 3.75, o = 7.0 }`). The
  nucleic-acid models are unweighted means of the four chain residues
  (a nucleoside monophosphate less one water): RNA C₉.₅H₁₁.₇₅N₃.₇₅O₇ per
  321.2916 Da, DNA C₉.₇₅H₁₂.₂₅N₃.₇₅O₆ per 308.8006 Da. Phosphorus is carried in
  the residue mass and nowhere else — ³¹P is the only stable phosphorus isotope,
  so it cannot shift a pattern, and excluding it from the convolution is exact
  rather than an approximation. Hill detection and chain assembly never assumed
  an analyte class; scoring was the only place one was hard-coded.
  `koth_align` takes the same key in `[lfq]`, and it must match the model the
  per-run features were detected with.

  A model without sulfur ignores `sulfur_offsets` entirely: varying an atom the
  analyte does not contain would hand every candidate, decoys included, a free
  maximum over templates.

  Measured on a PXD075396 RNase digest of *E. coli* rRNA (negative mode, 9601
  MS1 scans), RNA model vs peptide model: 15,578 vs 14,510 features, mean
  isotope score 0.848 vs 0.821, share scoring ≥ 0.90 up from 33.8 % to 42.9 %
  (on features with ≥ 4 isotopes, 0.926 vs 0.892).
- **`[file] polarity`** — `"positive"` (default, `M + zH`) or `"negative"`
  (`M − zH`). Nucleic acids are acquired in negative mode, where the previous
  fixed positive-mode arithmetic mis-massed every feature by 2·z·1.00728 Da, or
  8 Da at charge 4 — enough to defeat any downstream identification. koth does
  not read the polarity out of the file; set it alongside `isotope_model`.

### Fixed
- **The element cache clamped oxygen and nitrogen at 100 atoms**, sized for
  tryptic peptides. Oligonucleotides are oxygen-rich (7 O per 321 Da residue
  against a peptide's 1.48 per 111 Da), so a 5 kDa RNA needs 109 O and a 9 kDa
  one 196, and every count above the cap was silently scored against a truncated
  composition — a wrong isotope pattern, not a slow one. `MAX_O` is now 320 and
  `MAX_N` 220; the tables cost 80 bytes per count.

  This also reaches heavy peptides. On PXD003881 run B03_02, 3 of 135,085
  features change `isotope_score` (all at neutral mass ≥ 7.9 kDa, needing
  105–110 O): 0.7073 → 0.6998, 0.5729 → 0.5684, 0.5721 → 0.5684. No feature is
  gained or lost, hills are byte-identical, and every other feature row is
  unchanged. The published 0.3.0 numbers are unaffected — they were produced by
  the 0.3.0 binaries — and a re-run at 0.4.0 would move them by nothing
  measurable.

### Changed
- **Library API (breaking).** `run_scoring`, `scoring::score_features`,
  `write_features_tsv` and `write_features_parquet` take the ion polarity, and
  `Feature::monoisotopic_neutral_mass` / `ScoredFeature::monoisotopic_neutral_mass`
  take a `Polarity` rather than assuming protonation. There is deliberately no
  polarity-free overload: a silent positive-mode default is the bug this release
  fixes. The isotope model travels in the config structs, so the CLI is
  unaffected except for the two new keys.
- `scoring::averagine`'s mass-to-composition helpers moved onto the new
  `scoring::model::IsotopeModel` (`counts`, `distribution`,
  `distribution_with_sulfur`, `sulfur_ceil`); `lookup_template` is gone.
  `resolve_sulfur_counts` and `bhattacharyya_score_best_sulfur` take the model.

## [0.3.0] — 2026-09-08

### Changed
- **Benchmark operating point is now the default.** Every struct default (and
  every serde `default = ...` fallback) now equals the value the published
  benchmark configs set, so a run with no `--config` reproduces the paper's
  operating point and `example_config*.toml` really do show the defaults.
  Old → new:
  - `[hills]`: `max_gap` 0 → 1; `lfc_weight` 0.5 → 0.3; `split_valley_ratio`
    0.70 → 0.60; `split_sigma_mult` 4.0 → 5.0 (`split_height_frac` stays 0.10).
  - `[features]`: `min_charge` 1 → 2; `max_charge` 7 → 6; `min_chain_cosine`
    0.5 → 0.4; `min_isotope_score` 0.0 → 0.5; `min_isotope_step_ratio`
    (née `right_max_decrease`) 0.05 → 0.01.
  - `[file]`: `mz_recalibration` false → true (the `--recalibrate` flag is now
    a no-op unless a config turned it off).
  - `[lfq]`: `rt_window_pct` 0.02 → 0.01; `lone_coelution` 1.0 → 0.5;
    `decoy_own_template` false → true; `detected_use_grid` false → true (see
    below).
  - `[lfq.consensus]`: `min_group_size` 1 → 2; `min_member_combined_score`
    0.0 → 0.5 (now has a serde default, so it is no longer a required key);
    `min_seed_combined_score` 0.0 → 0.75.
  - `[alignment].rt_warp_kind` already defaulted to `"ransac"`; unchanged.
  A TOML that pinned the previous values still parses and behaves as before;
  only unset keys move.
- **`[features].right_max_decrease` is renamed `min_isotope_step_ratio`.** The
  semantics are unchanged: a heavier isotope must be at least this fraction of
  its predecessor's intensity or the chain stops. The old key is accepted as a
  serde alias, so existing configs keep working.
- **`koth_align` now defaults to `detected_use_grid = true`**: every consensus
  cell — detected and match-between-runs alike — is quantified by the same XIC
  grid re-integration, so the whole matrix sits on one intensity scale. The old
  default trusted the per-run feature's own intensity for detected cells, which
  mixes scales (on timsTOF the feature integrates ion mobility while the 2-D
  grid does not; replicate CV inflated ~38 → 13.7 %). All-grid quantification
  is also a measured win on Orbitrap (PXD003881: gated matrix CV 14.27 →
  12.31 %, ECOLI bias −0.067 → −0.020, HUMAN IQR 0.219 → 0.194; recall and
  missingness unchanged). Under all-grid, `quant_estimator = "sum"` beats
  `"apex"` on Orbitrap and ties it on timsTOF, so the shipped configs use
  `"sum"` on both platforms and the per-platform `[lfq]` split is gone. Set
  `detected_use_grid = false` only to reproduce older mixed-scale matrices.

## [0.2.0] — 2026-09-03

### Added
- Multi-CV FAIMS support for MS1 feature detection. mzML and native Thermo
  readers carry compensation voltage into independent, channel-local hill and
  feature pipelines; the feature TSV and Parquet schemas expose a nullable
  `FAIMS` column. Non-FAIMS inputs retain the existing single-channel detection
  path.

## [0.1.0] — 2026-08-18

Initial release, tagged from `992d2f4`. The sections below cover the whole
pre-release development history from 2026-04-07 to that tag, so *Changed*,
*Fixed* and *Removed / renamed* are relative to earlier development builds
rather than to any published version.

### Added
- `koth_ff` CLI for streaming hill detection, isotope-feature detection, and
  averagine scoring on mzML files and Bruker timsTOF `.d` directories.
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
