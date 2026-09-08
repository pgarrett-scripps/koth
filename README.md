# koth_ff

High-performance LC-MS feature finder and label-free quantifier for mzML, Bruker timsTOF (.d),
and (optionally) native Thermo Fisher (.raw) data, written in Rust.

Takes centroided MS1 data and produces two outputs: `hills.tsv` (chromatographic traces) and
`features.tsv` (isotope envelopes with charge states and averagine scores). The pipeline is
designed for large files — peak RSS is around 150 MB on a 1 GB mzML file.

The project also includes `koth_align`, a multi-run alignment and label-free quantification
pipeline. `koth_align` reads batch output from `koth_ff`, corrects systematic RT, mass, and
ion-mobility offsets between runs, and writes a feature × sample intensity matrix ready for
downstream statistical analysis.

## Install

The project requires Rust 1.88 or newer. To install both command-line tools from
a source checkout:

```bash
cargo install --locked --path koth_ff
```

Tagged releases also produce archives containing both `koth_ff` and
`koth_align`, plus SHA-256 checksum files, for Linux x86_64, macOS x86_64 and
arm64, and Windows x86_64.

For a development build:

```bash
cargo build --release
# binaries are at target/release/koth_ff and target/release/koth_align
```

Both binaries are built together with the command above. Bruker timsTOF support is compiled in by
default (requires `timsrust`). To build without it:

```bash
cargo build --release --no-default-features
```

### Native Thermo `.raw` input (optional)

Reading Thermo Fisher `.raw` files directly — no prior mzML conversion — is available behind the
`thermo` feature. It wraps Thermo's `RawFileReader` assemblies via a self-hosted **.NET 8 runtime**,
which must be installed at build and run time, so it is **off by default**:

```bash
cargo build --release -p koth_ff --features thermo
```

koth_ff auto-detects a .NET runtime in the usual locations (`~/.dotnet`, `/usr/share/dotnet`, …);
set `DOTNET_ROOT` explicitly if yours lives elsewhere. `.raw` reading is local-file only.

## Quick start

### Single-run feature finding (koth_ff)

```bash
# mzML input, results written to ./out/<stem>/
koth_ff data.mzML --output ./out

# Bruker .d directory
koth_ff data.d --output ./out

# Thermo .raw (requires a build with --features thermo)
koth_ff data.raw --output ./out

# With a custom config
koth_ff data.mzML --config my_config.toml --output ./out

# Skip the scoring stage (faster)
koth_ff data.mzML --output ./out --no-scoring
```

The output directory is `<output>/<input_stem>/` and always contains:

```
hills.tsv
features.tsv
report.json     # summary statistics (hill and feature counts, score distribution)
config.toml     # the config that was used, for reproducibility
```

### Multi-run alignment and LFQ (koth_align)

```bash
# Step 1 — process each sample with koth_ff
koth_ff sample1.mzML --output batch/
koth_ff sample2.mzML --output batch/
koth_ff sample3.mzML --output batch/
# produces batch/sample1/, batch/sample2/, batch/sample3/

# Step 2 — align and quantify
koth_align batch/                                          # output to batch/align_output/
koth_align batch/ --output results/                        # explicit output directory
koth_align batch/ --config example_config_align.toml --output results/
```

Output in `<output>/` (default: `<batch_dir>/align_output/`):

- `consensus_features.tsv` — one row per reference feature with mass, mz, charge, rtApex, im, score, seed_run, n_contributing_runs, n_runs_detected
- `intensity_matrix.tsv` — features × samples with integrated intensities (0 = not detected)
- `qvalue_matrix.tsv` — TDC q-values for each (feature, sample) cell; only written when `run_tdc = true`; use to filter intensity_matrix by FDR
- `align_config.toml` — copy of config used

## justfile recipes

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

## Algorithm

The single-run pipeline runs in three sequential stages. Multi-run alignment and LFQ add two
further stages that operate on the collected output from Stage 1–3.

### Stage 1: Hill detection

Hills are chromatographic traces — a single m/z signal tracked across consecutive MS1 scans.

The detector processes the file one spectrum at a time using a streaming iterator (`stream_mzml`).
No `Vec<Spectrum>` is ever held in memory; peak RSS is O(active hills) rather than
O(total peaks in the file).

For each scan the detector builds a scratch buffer of `(mz_mean, hill_id)` pairs sorted by m/z,
then binary-searches it for each incoming peak. Because peaks arrive sorted by m/z within a scan,
each binary search lands in a narrow window and the inner loop is short. Compared to a HashMap bin
approach there is no hash overhead and no bin-size tuning.

Matching rules:

- Closest unmatched hill within `mz_tolerance` (ppm or Da) is selected, with a combined m/z +
  intensity log-fold-change distance score weighted by `lfc_weight`.
- If ion mobility data is present, an IM tolerance check is applied and the combined
  m/z + IM distance is used to break ties.
- A matched peak extends the hill's running m/z mean (Welford online update) and appends
  the intensity to the profile.
- An unmatched peak starts a new hill.

After each peak is processed the detector periodically evicts stale hills — any hill whose
`last_scan_seen` is more than `max_gap` scans behind the current scan index. Evicted hills
shorter than `min_scans` are discarded; the rest are finalized and appended to the output list.

Intensity profiles are stored as `f32` vectors (4 bytes per position). Gap positions — scans
where no peak matched — are represented as `0.0` intensity rather than `Option<f64>`, cutting
per-element memory from 16 bytes to 4 bytes.

### Co-elution splitting

After all hills are finalized, any hill that contains multiple local intensity maxima is split
into sub-hills at the valley between them.

Split criteria (all must hold):

1. Two local maxima are at least `min_peak_distance` scans apart.
2. Each maximum reaches at least `min_peak_height` fraction of the hill's global maximum.
3. Each maximum has a prominence (peak height minus the highest valley between it and any taller
   neighbour) of at least `min_prominence` fraction of the global maximum. This prevents
   noise wiggles on a flank from triggering false splits.
4. Both resulting segments are at least `min_scans` scans long.

Splitting is controlled by `split_hills = true` in the config.

### Stage 2: Feature detection

Features are isotope envelopes — groups of hills whose m/z values are spaced by
`neutron_mass / charge` (where `neutron_mass` = 1.003354835 Da, the C13 offset).

Every hill is tried as a monoisotopic seed, at every charge state from `min_charge` to
`max_charge`. Chains extend **upward only** — seed → M+1 → M+2 → … — because the seed *is*
the monoisotopic hypothesis. There is no downward walk: a hill that one could reach is
itself a seed that builds the same envelope upward, with the averagine template indexed
from its own position.

A candidate partner must satisfy:

- m/z within `mz_tolerance` of the expected isotope position.
- Scan range overlaps with the reference hill.
- If IM data is present, IM within `im_tolerance` of the reference hill.
- `intensity_max >= ref_hill.intensity_max * min_isotope_step_ratio` (prevents linking to an
  implausibly weak signal).
- Cosine of the elution profile against the `cosine_anchor` reference hill (the seed by
  default) is at least `min_chain_cosine`.
- The apex intensity ratio against the chain predecessor matches the averagine ratio within
  ±`max_isotope_log2_ratio`.

This produces an over-complete `(seed, charge)` candidate pool, which is then resolved
**non-destructively**: contested hills are claimed longest-envelope-first, and a candidate
whose hills are partly claimed is truncated to its free monoisotope-anchored prefix,
re-scored, and re-queued rather than dropped. Feature detection also records
per-adjacent-pair cosine similarities and ppm errors for downstream quality filtering.

### Stage 3: Scoring

Each feature is compared against the theoretical isotope distribution predicted by the averagine
model (average amino acid composition: C₄.₉₃₈₄ H₇.₇₅₈₃ N₁.₃₅₇₇ O₁.₄₇₇₃ S₀.₀₄₁₇ per 111.13 Da).

The averagine distributions are precomputed at startup in 50 Da steps from 50 to 5050 Da and
stored in a static lookup table. Each distribution is built by convolving Poisson approximations
for each element, then normalized.

Scoring steps:

1. Look up the theoretical pattern for the feature's neutral mass.
2. Try neutron offsets in `[isotope_offset_min, isotope_offset_max]`. A nonzero offset shifts
   the observed pattern left or right, which corrects for cases where the detector assigned the
   wrong peak as the monoisotopic ion.
3. Score each offset with the Bhattacharyya coefficient between the (normalized) observed and
   theoretical patterns, penalized by the fraction of the theoretical distribution not covered
   by the observed peaks. Offset 0 receives a small `offset_zero_bonus` to prefer no reassignment
   when scores are close.
4. If the best score is below `min_score_threshold`, fall back to offset 0.

The final score is clamped to `[0, 1]`. A score near 1.0 means the observed isotope pattern
closely matches the averagine expectation for that mass. Features whose score falls below
`[features].min_score` are dropped before writing output.

### Stage 4: Multi-run alignment

The alignment module (`koth_ff/src/alignment/`) corrects systematic RT, mass, and ion-mobility
offsets between runs before quantification.

**Reference selection**: The run with the most high-confidence features (score ≥
`min_anchor_score`) is chosen automatically as the reference. All other runs are aligned to it.

**Anchor matching**: For each non-reference run, high-confidence features are matched to
reference features by (charge, monoisotopic m/z, RT). Matching uses binary search on (charge,
mz) sorted features and filters by:

- m/z within `anchor_mass_ppm` ppm
- Normalised RT within `rt_anchor_window` ([0, 1] space)
- Ion mobility within `im_tolerance` (when available)

Matching is 1:1 greedy (each reference feature claims at most one run feature, picked by lowest
PPM error).

**RT warp**: A sliding-window median of RT deltas is computed across the normalised RT axis
using a window of `rt_warp_bandwidth` width. Anchor pairs whose residuals exceed
`rt_warp_sigma_clip` × σ are clipped and the medians are refit (up to `rt_warp_clip_iters`
times). The resulting knots are interpolated piecewise-linearly to produce a continuous warp
function mapping run RT → reference RT. If fewer than `min_anchor_count` anchors are found the
run falls back to an identity warp.

**Mass and IM drift**: PPM errors and IM deltas across anchors are fit as linear functions of
RT (ordinary least squares with sigma-clipping). The fitted lines are used to correct m/z and
ion mobility values before LFQ extraction.

### Stage 5: Label-free quantification

The LFQ module (`koth_ff/src/lfq/`) re-extracts integrated intensities from the pre-computed
hill data for every (consensus feature, run) pair. This is done instead of matching scored
features, because a feature may not pass the scoring threshold in every run even when real
signal is present, and different runs may have different isotope coverage or noise.

**Consensus feature list**: All features from all runs are projected into reference-run coordinate
space and grouped by (charge, neutral mass ± ppm, aligned RT ± window, IM ± tolerance). The
highest-scoring feature in each group seeds the LFQ extraction. Groups can be filtered by
`[lfq.consensus]` settings.

**XIC grid construction**: For each consensus feature in each run, alignment corrections are
applied to obtain the expected RT, m/z, and IM in that run's coordinate space. A `grid_cols`-bin ×
`n_isotopes`-row grid is allocated covering ±`rt_window_pct` of the run's total gradient. Hills
from the run are matched to each isotopologue slot (M, M+1, M+2) within `mz_ppm` ppm and the
RT/IM window. Each matched hill's per-scan intensity profile is binned into the grid columns.

**Column scoring**: Each RT column is scored on four quality terms. Two of these are the
*orthogonal* isotope signals — a pattern match and a co-elution measure — that used to be
conflated under one "spectral" name:

- **RT score**: `1 − ∛(|col − centre| / half_cols)` — penalises columns far from the centre
- **Intensity score**: `√(col_total / max_col_total)` — relative intensity across the window
- **Bhattacharyya**: agreement between the observed (M, M+1, M+2) intensities and the theoretical
  averagine pattern, *penalising expected-but-missing peaks* (a lone monoisotope scores low, not
  a free 1.0). This is the "does the isotope pattern match theory?" signal.
- **Co-elution**: cosine similarity between the monoisotope's XIC trace and each higher isotope's
  trace across RT — "do the matched isotopes rise and fall together?". Orthogonal to the pattern
  match; 1.0 when fewer than two isotope rows carry signal.
- **Hybrid score**: `(RT × intensity × bhattacharyya × coelution)^¼` (default mode)

Both the Bhattacharyya and co-elution terms are always computed — there is no flag to swap one
metric for another.

**Peak integration**: The apex column (maximum hybrid score) is found. The peak is expanded
left and right until: 20 bins have been added, the spectral Bhattacharyya drops below
`min_spectral_bhattacharyya`, or the score falls below half the apex score. All intensities in
the expanded window across all isotopologue rows are summed.

**Target-decoy competition and MBR FDR**: A decoy is generated for each consensus feature by
shifting m/z by `decoy_mz_shift_da`/charge (default +11 Da) and RT back by `decoy_rt_shift_pct`
of the gradient, then extracted with identical logic. Q-values are then computed by
`tdc_method`:

- **`qda`** (default) — a semi-supervised quadratic-discriminant rescorer. Over all target cells
  (detected + match-between-runs) and the decoys, it learns a QDA over five symmetric per-cell
  features (`|ppm error|`, `|RT diff|`, Bhattacharyya, co-elution, `|IM delta|` — the last inert
  on Orbitrap, active on timsTOF) using a Percolator-style loop (iterative confident-positive
  selection, decoys as negatives, 3-fold cross-validation by feature index), then computes
  q-values from the held-out scores. Every cell is scored — detected cells get no free pass — so
  a background-contaminated cell in a depleted well is gatable regardless of how it was populated.
  Deterministic.
- **`hybrid`** — the legacy method: rank all target and decoy cells by `hybrid_score` and take
  the monotonised running FDR (n_decoy / n_target).

## Memory design

The streaming hill detector never accumulates spectra. It reads one spectrum, updates the active
hill map, evicts stale hills, and discards the spectrum. For a 1 GB mzML file this keeps peak
RSS around 150 MB.

Once hill detection is complete the hill list is passed to feature detection. Hills hold their
intensity profiles as `Arc<[f32]>`, so features can reference the same profile data without
copying. After feature detection the hill `Vec` is explicitly dropped; the `Arc` reference
counts keep the profiles alive inside each `Feature`.

`koth_align` applies the same discipline across runs: it holds only the per-run *features* (small)
in memory for the whole run, and **streams each run's hills one at a time** during LFQ — loading a
run's hills, quantifying every consensus feature against them, then dropping them before the next
run loads. Peak memory is therefore O(one run's hills) rather than O(all runs' hills), so a
20-run cohort costs about the same as a single run rather than 20×.

## Outputs

### hills.tsv

One row per chromatographic hill, sorted by `intensity_sum` descending.

| Column | Description |
|---|---|
| `mz` | Intensity-weighted mean m/z across all scans |
| `mz_std` | Intensity-weighted standard deviation of m/z |
| `rt` | Retention time at the apex scan (minutes) |
| `rt_start` / `rt_end` | Retention time at first and last scan |
| `rt_width` | `rt_end - rt_start` |
| `im` | Intensity-weighted mean ion mobility (0 if not available) |
| `im_std` | Intensity-weighted standard deviation of ion mobility |
| `scan_start` / `scan_apex` / `scan_end` | Absolute scan indices |
| `n_scans` | Total scans in profile (including gap positions) |
| `skipped_scans` | Number of gap positions (zero-intensity) |
| `intensity_sum` | Sum of all intensities in the profile |
| `intensity_max` | Peak intensity |
| `hill_score` | Shape score: fraction of scans monotone toward apex (0–1) |
| `intensity_profile` | JSON array of per-scan intensities (f32) |

### features.tsv

One row per isotope feature, sorted by `intensitySum` descending. Features with unknown charge
(charge = 0) are excluded.

| Column | Description |
|---|---|
| `massCalib` | Monoisotopic neutral mass (corrected for neutron offset) |
| `mz` | Monoisotopic m/z |
| `rtApex` / `rtStart` / `rtEnd` | Retention time extent |
| `intensityApex` | Sum of per-isotope intensities at the apex scan |
| `intensitySum` | Total intensity summed across all isotopes and scans |
| `charge` | Assigned charge state |
| `nIsotopes` | Number of isotope peaks in the envelope |
| `nScans` | Total scans covered across all hills in the feature |
| `im` | Ion mobility at apex (empty if not available) |
| `cosine_similarity` | Mean adjacent-pair cosine similarity of elution profiles |
| `ppm_error` | Mean m/z spacing error relative to the theoretical neutron step |
| `neutron_offset` | Offset applied to find the true monoisotopic peak (0, ±1) |
| `score` | Averagine Bhattacharyya score in [0, 1] |
| `theoretical_pattern` | JSON array: normalized averagine distribution |
| `isotope_profile` | JSON array: per-isotope apex intensities |
| `elution_profile` | JSON array: total intensity per scan across the feature |

### koth_align outputs

#### consensus_features.tsv

One row per consensus (reference) feature.

| Column | Description |
|---|---|
| `massCalib` | Monoisotopic neutral mass |
| `mz` | Monoisotopic m/z |
| `charge` | Charge state |
| `rtApex` | Retention time at apex in the reference run (minutes) |
| `im` | Ion mobility at apex (empty if not available) |
| `score` | Averagine isotope pattern score from the seed feature [0, 1] |
| `seed_run` | Name of the run that provided the seed feature for this row |
| `n_contributing_runs` | Runs that contributed a detection to this consensus group |
| `n_runs_detected` | Runs with intensity > 0 and q-value ≤ `max_qvalue` |

#### intensity_matrix.tsv

Feature metadata columns (massCalib, mz, charge, rtApex, im, score, seed_run,
n_contributing_runs) followed by one intensity column per run. Values are integrated intensities
from the XIC grid; 0 means no peak was found. No FDR filtering is applied — use
`qvalue_matrix.tsv` to filter downstream.

#### qvalue_matrix.tsv

Same layout as `intensity_matrix.tsv` but cells contain q-values (0–1) from the MBR rescorer
(`tdc_method`, default `qda`). Written only when `run_tdc = true`. A value of 1.0 means no signal
was found or TDC was not run; a detected (feature-supported) cell reports 0. Filter
intensity_matrix at q ≤ 0.01 for 1% FDR, for example.

## Configuration

Both binaries take a TOML config via `--config`. Every key is optional (defaults
apply); the fully-resolved config is written back out (`config.toml` /
`align_config.toml`) after each run. Parsing uses `deny_unknown_fields`, so a
stale or misspelled key is a hard error rather than a silent fallback.

**The complete, field-by-field reference — every setting, its default, what it
does, and what to set it to (including Orbitrap vs Bruker/timsTOF) — lives in
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md).** Read that first.

This repository ships two templates, [`example_config.toml`](example_config.toml)
(koth_ff) and [`example_config_align.toml`](example_config_align.toml)
(koth_align). Both list every knob at its struct default.

The tuned per-platform configurations behind the published benchmark
(`koth_ff.toml` / `koth_ff_bruker.toml` and `koth_align.toml` /
`koth_align_bruker.toml`, Orbitrap and timsTOF) live with the benchmark that
produced them, in [`benchmark/config/`](https://github.com/tacular-omics/koth-paper/tree/master/benchmark/config)
of the [koth-paper](https://github.com/tacular-omics/koth-paper) repository. They differ from the defaults in
24 settings between the two platforms, so start from the one matching your
instrument rather than from the templates if you are reproducing the paper.

koth_ff sections: `[file]` (reading + tolerances + Bruker front-end + m/z
recalibration), `[hills]` (trace detection + splitter), `[features]` (isotope
chains + charge + retention gates), `[scoring]`, `[output]`. koth_align sections:
`[alignment]` (RANSAC RT warp + mass/IM drift), `[lfq]` (XIC extraction +
integration + decoy/TDC), `[lfq.consensus]` (cross-run grouping), `[output]`. The
only per-platform differences are `[file].mz_tolerance` (8 ppm Orbitrap / 15 ppm
timsTOF) and, for koth_align, `[lfq].quant_estimator` and `[lfq].detected_use_grid`.

## Library API

### Single-run pipeline

`koth_ff` is also a library crate. The top-level functions mirror the CLI stages:

```rust
use std::path::Path;
use koth_ff::{run_hills_streaming, run_features, run_scoring, config::KothConfig};

let config = KothConfig::default();
let input = Path::new("data.mzML");

// Stage 1 — streaming; never holds Vec<Spectrum>
let hills = run_hills_streaming(input, &config.hills, &config.file)?;

// Stage 2
let features = run_features(&hills, &config.features, &config.file)?;

// Stage 3 — pass config.features.min_score to filter by isotope quality
let scored = run_scoring(&features, &config.scoring, config.features.min_score);
```

For non-streaming use (e.g. when you already have spectra in memory):

```rust
let spectra = koth_ff::read_spectra(input, &config.file)?;
let hills = koth_ff::run_hills(&spectra, &config.hills, &config.file);
```

### Alignment and LFQ

```rust
use koth_ff::{
    alignment::{align_runs, AlignmentConfig, RunInput},
    lfq::{quantify, LfqConfig},
};

// Build one RunInput per LC-MS run
let runs: Vec<RunInput> = vec![
    RunInput { name: "sample1".into(), features: scored1, hills: hills1, scan_times: vec![] },
    RunInput { name: "sample2".into(), features: scored2, hills: hills2, scan_times: vec![] },
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
sufficient for most data.

Loading runs from `koth_ff` output files:

```rust
use koth_ff::input::{discover_runs, read_hills, read_features};

let run_paths = discover_runs(Path::new("batch/"))?;
for rp in &run_paths {
    let hills    = read_hills(&rp.hills_path)?;
    let features = read_features(&rp.features_path)?;
    // build RunInput ...
}
```

## License

See [LICENSE](LICENSE).
