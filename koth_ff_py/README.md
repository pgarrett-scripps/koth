# koth_ff

High-performance LC-MS feature finder for mzML and Bruker timsTOF (.d) data, written in Rust.

Takes centroided MS1 data and produces two outputs: `hills.tsv` (chromatographic traces) and
`features.tsv` (isotope envelopes with charge states and averagine scores). The pipeline is
designed for large files — peak RSS is around 150 MB on a 1 GB mzML file.

## Install

```bash
cargo build --release
# binary is at target/release/koth_ff
```

Bruker timsTOF support is compiled in by default (requires `timsrust`). To build without it:

```bash
cargo build --release --no-default-features
```

## Quick start

```bash
# mzML input, results written to ./out/<stem>/
koth_ff data.mzML --output ./out

# Bruker .d directory
koth_ff data.d --output ./out

# With a custom config
koth_ff data.mzML --config my_config.toml --output ./out

# Skip the scoring stage (faster)
koth_ff data.mzML --output ./out --no-scoring

# Diagnostic: count scans and peaks without running the pipeline
koth_ff data.mzML --count-scans
```

The output directory is `<output>/<input_stem>/` and always contains:

```
hills.tsv
features.tsv
config.toml     # the config that was used, for reproducibility
```

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

The pipeline runs in three sequential stages.

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

- Closest unmatched hill within `mz_tolerance` (ppm or Da) is selected.
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
3. The valley between them is no deeper than `min_valley_ratio` times the shorter of the two
   flanking peaks. Shallower valleys are noise and are not split.
4. Both resulting segments are at least `min_scans` scans long.

Splitting is controlled by `split_hills = true` in the config.

### Stage 2: Feature detection

Features are isotope envelopes — groups of hills whose m/z values are spaced by
`neutron_mass / charge` (where `neutron_mass` = 1.003354835 Da, the C13 offset).

The algorithm seeds from the highest-intensity hill downward. For each unassigned seed it tries
every charge state from `max_charge` down to `min_charge` and searches for isotope partners
both to the right (M+1, M+2, ...) and to the left (M-1, M-2, ...) of the seed.

A candidate partner must satisfy:

- m/z within `mz_tolerance` of the expected isotope position.
- Scan range overlaps with the reference hill.
- If IM data is present, IM within `im_tolerance` of the reference hill.
- `intensity_max >= ref_hill.intensity_max * max_decrease` (prevents linking to an implausibly
  weak signal).
- Cosine similarity of the elution profile with the seed hill is at least `min_cosine_similarity`.

The charge state that produces the longest isotope chain is kept. All hills in the chain are
marked as assigned so they cannot be reused by a later seed. Feature detection also records
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
closely matches the averagine expectation for that mass.

## Memory design

The streaming hill detector never accumulates spectra. It reads one spectrum, updates the active
hill map, evicts stale hills, and discards the spectrum. For a 1 GB mzML file this keeps peak
RSS around 150 MB.

Once hill detection is complete the hill list is passed to feature detection. Hills hold their
intensity profiles as `Arc<[f32]>`, so features can reference the same profile data without
copying. After feature detection the hill `Vec` is explicitly dropped; the `Arc` reference
counts keep the profiles alive inside each `Feature`.

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
| `mono_hills_scan_lists` | JSON array of scan index arrays, one per isotope hill |
| `mono_hills_intensity_list` | JSON array of intensity arrays, one per isotope hill |

## Configuration

Configuration is a TOML file passed with `--config`. All values have defaults; you only need to
include the keys you want to change. Running the pipeline writes the resolved config to
`config.toml` in the output directory.

```toml
[hills]
mz_tolerance = 8.0          # m/z window for linking peaks to hills
mz_tolerance_type = "ppm"   # "ppm" or "da"
min_scans = 3               # discard hills shorter than this
max_gap = 1                 # consecutive missed scans before a hill is closed
split_hills = true          # split co-eluting hills at valleys
min_peak_distance = 10      # minimum scan separation between peaks when splitting
min_peak_height = 0.2       # minimum peak height relative to the hill maximum
min_valley_ratio = 0.6      # valley/peak ratio below which a split is applied
im_tolerance = 0.05         # ion mobility tolerance
im_tolerance_type = "relative"  # "relative" (fraction of IM value) or "absolute"
global_min_mz = 0.0         # ignore peaks below this m/z
global_max_mz = 10000.0     # ignore peaks above this m/z
bruker_mz_ppm = 5.0         # m/z tolerance for Bruker .d centroiding
bruker_im_pct = 3.0         # IM tolerance (%) for Bruker .d centroiding
# n_threads = 8             # parallelism; omit for all CPUs

[features]
mz_tolerance = 5.0          # m/z tolerance for isotope partner search
mz_tolerance_type = "ppm"
min_charge = 1
max_charge = 7
min_cosine_similarity = 0.5 # minimum elution profile cosine to extend an envelope
left_max_decrease = 0.9     # max allowed intensity drop on the low-m/z side
right_max_decrease = 0.9    # max allowed intensity drop on the high-m/z side
im_tolerance = 0.05
im_tolerance_type = "absolute"
max_isotopes = 6            # maximum isotope peaks per feature
neutron_mass = 1.003354835  # C13 mass offset in Da

[scoring]
isotope_offset_min = -1     # lower bound of neutron offset search
isotope_offset_max = 1      # upper bound of neutron offset search
offset_zero_bonus = 0.15    # score bonus for keeping offset = 0
min_score_threshold = 0.5   # features below this score keep offset = 0
```

A full annotated template is available at `example_config.toml`.

## Library API

`koth_ff` is also a library crate. The top-level functions mirror the CLI stages:

```rust
use std::path::Path;
use koth_ff::{run_hills_streaming, run_features, run_scoring, config::KothConfig};

let config = KothConfig::default();
let input = Path::new("data.mzML");

// Stage 1 — streaming; never holds Vec<Spectrum>
let hills = run_hills_streaming(input, &config.hills)?;

// Stage 2
let features = run_features(&hills, &config.features)?;

// Stage 3
let scored = run_scoring(&features, &config.scoring);
```

For non-streaming use (e.g. when you already have spectra in memory):

```rust
let spectra = koth_ff::read_spectra(input, &config.hills)?;
let hills = koth_ff::run_hills(&spectra, &config.hills);
```

## License

See [LICENSE](LICENSE).
