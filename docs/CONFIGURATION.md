# koth_ff / koth_align Configuration Reference

Authoritative, field-by-field reference for every configuration setting in the
two binaries. If you are an AI agent or a new user trying to understand *what a
setting does and what to set it to*, this is the file to read. It is kept in
sync with the config structs in `koth_ff/src/config/`, `koth_ff/src/lfq/`, and
`koth_ff/src/alignment/`.

> **Ground truth:** every TOML is parsed with `#[serde(deny_unknown_fields)]` at
> every level — a misspelled or stale key is a **hard load error**, not a silent
> fallback. A build-time test parses the shipped configs against the structs, so
> the field names here cannot drift from the code without breaking the build.
> Running either binary writes the fully-resolved config back out
> (`config.toml` / `align_config.toml`) for reproducibility.

---

## 1. Mental model — what the pipeline does

Two binaries, built together (`cargo build --release` → `target/release/koth_ff`
and `target/release/koth_align`).

### koth_ff — per-run feature finding

```
raw MS1  →  [file] read + tolerances + (Bruker front-end) + (m/z recalibration)
         →  [hills]    Stage 1: streaming chromatographic-trace ("hill") detection
         →  [features] Stage 2: isotope-chain assembly → charge + neutral mass + averagine score
         →  [scoring]  Stage 3 (optional): ±1-neutron monoisotope reassignment
         →  [output]   write hills.{tsv|parquet} + features.{tsv|parquet}
```

One input file → `hills` (chromatographic traces) + `features` (isotope
envelopes with charge and averagine/Bhattacharyya scores). Bruker `.d` input
adds a front-end (`dnoise` vertical-IM filter + watershed centroider) before
hill detection; on mzML those knobs are inert.

### koth_align — cross-run alignment + LFQ + FDR

```
batch of koth_ff run dirs
   →  [alignment]     pick a reference run; RANSAC-warp every other run's RT + fit mass/IM drift
   →  [lfq.consensus] project features into reference space; group the same peptide across runs
   →  [lfq]           per (consensus feature × run): build an XIC grid, score every RT column,
                      integrate the best peak; if run_tdc, repeat with a decoy and rank
   →  [output]        write consensus_features + intensity_matrix + qvalue_matrix (+ decoys, details)
```

Reads `koth_ff` output for N runs, corrects systematic RT/mass/IM offsets, and
emits a feature × run intensity matrix with per-cell target-decoy q-values.
Filter the matrix at e.g. `q ≤ 0.01` downstream for 1% FDR.

### The two things you actually tune per platform

Almost everything is platform-agnostic. The differences between the shipped
Orbitrap and Bruker/timsTOF configs are small and listed in
[§6 Platform presets](#6-platform-presets-orbitrap-vs-brukertimstof). The single
most important platform value is **`[file].mz_tolerance`** (8 ppm Orbitrap, 15
ppm timsTOF).

---

## 2. How configuration works

- Pass a TOML with `--config`. Every key is optional **except** where noted
  "required key"; omitted keys take the struct default. Because of
  `deny_unknown_fields`, you cannot leave a stale key lying around.
- A handful of CLI flags on `koth_ff` override the config for that run:
  `--no-scoring`, `--ms2`, `--recalibrate` (→ `[file].mz_recalibration`),
  `--filter-baseline-hills` (→ `[hills].filter_large_baseline_hills`), and
  `--threads` (→ `[file].n_threads`).
- This repository ships `example_config.toml` and `example_config_align.toml`,
  which list every knob at its struct default.
- The canonical worked examples are the tuned configs `koth_ff.toml` /
  `koth_ff_bruker.toml` (feature finding, Orbitrap / timsTOF) and
  `koth_align.toml` / `koth_align_bruker.toml` (alignment+LFQ). They live in
  `benchmark/config/` of the separate [koth-paper](https://github.com/tacular-omics/koth-paper) repository,
  alongside the benchmark that produced them. Prefer them over the templates
  when reproducing published results.

**Notation below:** `key` (type, `default`). Since 0.3.0 the struct defaults
*are* the benchmark's operating point, so the production Orbitrap config
differs from them only where a row says `default` → **`shipped`**. Platform
differences are called out inline.

---

## 3. `koth_ff` — feature-finding config (`KothConfig`)

Sections: `[file]`, `[hills]`, `[hills_ms2]` (optional MS2 overrides), `[features]`, `[scoring]`, `[output]`.

### 3.1 `[file]` — reading, shared tolerances, Bruker front-end, recalibration

Tolerances here are shared by **both** the hills and features stages.

#### Core tolerances
| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `mz_tolerance` | f64 | `8.0` | m/z match window for linking peaks into hills and isotopes. **The key platform value: 8.0 Orbitrap (Fusion OT is 2–5 ppm), 15.0 Bruker/timsTOF.** |
| `mz_tolerance_type` | enum | `"ppm"` | `"ppm"` or `"da"`. `"ppm"` on both platforms. Region-adaptive isotope tolerance is only honored for ppm. |
| `polarity` | enum | `"positive"` | `"positive"` (M + zH) or `"negative"` (M − zH) — how a neutral mass is recovered from an observed m/z. Peptides are positive; **nucleic acids are negative**, and the wrong setting shifts every reported `massCalib` by 2·z·1.00728 Da (8 Da at charge 4), enough to defeat any downstream identification. koth does not read polarity out of the file. Set it alongside `[features] isotope_model`. |
| `im_tolerance` | f64 | `0.05` | Ion-mobility tolerance. `0.05` everywhere; inert on Orbitrap but the key is schema-required. |
| `im_tolerance_type` | enum | `"relative"` | `"relative"` (fraction of IM value) or `"absolute"` (1/K0). `"relative"` on both. |
| `global_min_mz` | f64 | `0.0` | Ignore peaks below this m/z. Left at default. |
| `global_max_mz` | f64 | `inf` | Ignore peaks above this m/z. Left at default. |
| `intensity_coverage` | f64 | `1.0` | Fraction of hill intensity to retain. **Keep `1.0`** — validated strictly dominant (0.95→1.0 gave +2.9 pp recall at +0.5% features; 0.95 clips low-abundance apices). Only the AlphaPept-comparator config uses 0.95. |
| `n_threads` | Option\<usize\> | `None` (all cores) | Worker threads. **Unset in both shipped benchmark configs** (since 2026-08-27): v0.1.0 parsed this key but ignored it (always all cores), and every published timing ran uncapped. Current builds DO honor it, so a leftover value would silently cap a re-run. Set it only for deliberate throughput limiting. |

#### General noise / mode toggles
| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `noise_filter_sigma` | Option\<f64\> | `None` | Per-scan sigma-clip noise filter during hill detection (any format). Omit to disable (all shipped configs). Typical if used: `3.0`. |
| `decoy_mode` | bool | `false` | Shuffle MS1 scans before detection → null/decoy feature set. `false` normally; `true` only for FDR null-building. (Skips `mz_recalibration`.) |
| `ms2_hills_enabled` | bool | `false` | Also detect MS2 hills per isolation window, written to `hills_ms2.*`. Supports DIA mzML and diaPASEF; DDA inputs yield no MS2 hills. Overridable via `--ms2`. |

#### ID-free m/z recalibration + region-adaptive isotope tolerance
Enabled in both production configs (**code default is off; the tuned configs turn
it on**). Pass-1 detection learns a per-(m/z, RT)-region signed-ppm offset surface
from isotope-spacing residuals; pass-2 shifts the *expected* isotope position by
the learned offset and replaces the fixed isotope-match tolerance with
`clamp(sigma_mult × σ(m/z,RT), floor_ppm, mz_tolerance)`.

| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `mz_recalibration` | bool | `true` | Master switch. The **region-adaptive tolerance is the real win** (position-only recal is a documented no-op on well-calibrated instruments). It is a **precision** mechanism, not a recall one — see the measured ablation below. Overridable via `--recalibrate`. |
| `mz_recalibration_mz_bins` | usize | `20` | m/z bins in the recalibration surface. Default. |
| `mz_recalibration_rt_bins` | usize | `8` | RT bins in the recalibration surface. Default. |
| `mz_recalibration_min_samples` | usize | `50` | Min residuals a cell needs before its own median is trusted (else falls back to marginal → global median). Default. |
| `mz_recalibration_tol_sigma_mult` | f64 | `4.0` | The `N` in `N × σ` adaptive tolerance. **`4.0`**, finalized 2026-07-09 (was 3.0): a balanced single default within ~0.2 pp of each platform's optimum (5.0 Orbitrap / 3.0 Bruker). |
| `mz_recalibration_tol_floor_ppm` | f64 | `1.0` | Lower ppm bound so a tiny σ in a sparse region can't collapse the window. Default. |

**Measured on/off ablation (2026-07-30).** Both cohorts, all other settings at
the shipped defaults, recall scored under the shared 10 ppm + RT-interval join.
This supersedes an earlier entry here that claimed "Orbitrap recall flat" and
"Bruker recall +1.06 pp"; neither reproduced, and its CV figures (23.76→23.25 %
Orbitrap, 13.89→13.76 % Bruker) predate the config alignment and match no
current run.

| Metric | Orbitrap off | Orbitrap on | Bruker off | Bruker on |
|---|---|---|---|---|
| PSM recall | 79.85 % | **79.43 %** | 74.97 % | **75.13 %** |
| Median CV | 13.97 % | **13.94 %** | 9.95 % | **9.91 %** |
| HUMAN IQR | 0.287 | **0.285** | 0.399 | **0.398** |
| FFCR | 12.70 % | **12.50 %** | 16.42 % | **16.41 %** |
| Features | 2.78 M | 2.67 M | 2.03 M | 2.02 M |

Read this as a **precision/recall trade**, not a free win. Every quantitative
metric improves on both platforms, but by little (largest: Orbitrap FFCR,
0.20 pp). Recall moves in *opposite* directions — −0.42 pp Orbitrap, +0.16 pp
Bruker — both significant on an exact McNemar test over paired per-PSM outcomes
(p = 1.4e-62 and 1.2e-15). The Orbitrap recall loss is the price of the tighter
isotope-match window: it removes 3.9 % of features, some of which were matching
PSMs. Paired per-peptide ΔCV on Bruker is −0.0018 pp (Wilcoxon p = 0.007,
rank-biserial r = −0.009): significant and negligible at once, because n ≈ 82 k
paired cells. Keep it on for the precision and the lower spurious-feature count;
do not cite it as a recall win. Reproducers, all in the
[koth-paper](https://github.com/tacular-omics/koth-paper) repository: `benchmark/scripts/16_peptide_lfq.py`
(Orbitrap quant), `benchmark/scripts/bruker_validation.py` (Bruker), paired
tests in `paper/si/si-body.typ` @tab:si-recal-ablation.

#### Bruker vertical-IM filter (Stage 1) — inert on Orbitrap mzML
The `.d` front-end runs the `dnoise` vertical-IM feature filter over the raw
frames. All defaults are the tuned ddaPASEF values; leave them unless you know
the acquisition differs.

| Key | Type | Default | What it does |
|---|---|---|---|
| `bruker_filter_mz_half_width` | u32 | `2` | TOF-index half-width summed into each column profile. |
| `bruker_filter_max_internal_gap` | usize | `1` | Max empty scans tolerated inside a kept vertical run (morph-close radius). |
| `bruker_filter_min_feature_length` | usize | `5` | Min run span (scans) for a column feature to survive. |
| `bruker_filter_min_window_intensity` | u64 | `0` | Per-scan intensity floor for "occupied". |
| `bruker_filter_min_feature_intensity` | u64 | `0` | Total-intensity floor for a kept run. |
| `bruker_filter_num_iterations` | usize | `1` | Re-apply the filter to its own survivors N times (each pass stricter). |

#### Bruker watershed centroider (Stage 3) — inert on Orbitrap
| Key | Type | Default | What it does |
|---|---|---|---|
| `bruker_watershed_box_scan` | u32 | `10` | Neighbour reach on the scan (IM) axis. |
| `bruker_watershed_box_mz_idx` | u32 | `3` | Neighbour reach on the TOF (m/z index) axis. |
| `bruker_watershed_min_seed_intensity` | u64 | `0` | Intensity floor to seed a new group (else the orphan point is dropped). |
| `bruker_watershed_min_centroid_total` | u64 | `0` | Drop centroids whose summed group intensity is below this. |
| `bruker_watershed_max_tof_offset` | u32 | `10` | Cap on member distance from the group seed (TOF units); stops follower creep past the peak edge. |
| `bruker_noise_sigma` | Option\<f64\> | `None` | Per-frame MAD noise filter on centroided peaks before emitting a `Spectrum`. Omit (shipped). Typical if used: `3.0`. |

#### Bruker ion-mobility scale
| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `bruker_mobility_scale` | enum | `"calibrated"` | Scale of every reported 1/K0 (`im` columns) for `.d` input. `"calibrated"` converts each centroid's scan number with the run's acquisition calibration (`TimsCalibration` in `analysis.tdf`, ModelType 2), matching Bruker's timsdata SDK, DataAnalysis and SDK-based exports to within 1e-15 1/K0. `"linear"` restores the straight line between `OneOverK0AcqRange{Upper,Lower}` used up to 0.9.0, which differs by up to ~0.03 1/K0 (about 2%). A run whose calibration cannot be read, or uses another ModelType, is an error; set `"linear"` to process it. The chosen scale and the calibration rows are recorded in `report.json` under `mobility`, and `koth_align` refuses a batch that mixes scales (a timsTOF run with no recorded scale counts as linear). The streaming path's `bruker_ms1_polygon` gate still tests points on dnoise's linear scale. |

#### Bruker streaming path (experimental — see also [dnoise](#7-the-dnoise-integration-bruker-only))
| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `bruker_streaming` | bool | `false` | **Experimental / opt-in.** `true` = run dnoise's configured stages in-process (vertical filter → optional halo → optional MS1 polygon → watershed), with no denoised `.d` on disk. `false` = the historical local path (vertical + watershed, **no halo/polygon**) — **the paper's validated pipeline**. Keep `false` until re-validated on the Bruker cohort. |
| `bruker_halo` | bool | `true` | Streaming-only: apply the horizontal-halo filter after the vertical filter. No effect unless `bruker_streaming = true`. |
| `bruker_halo_peak_fraction` | f64 | `0.15` | Streaming-only: drop a peak below this fraction of the off-column box-max. |
| `bruker_halo_mz_idx_half_width` | u32 | `80` | Streaming-only: halo reference-box half-width along TOF index. |
| `bruker_halo_scan_half_width` | usize | `2` | Streaming-only: halo reference-box half-width along the IM scan axis. |
| `bruker_ms1_polygon` | bool | `false` | Streaming-only: apply a run's ddaPASEF IMS selection polygon to MS1 points. No-op for diaPASEF or runs without a polygon. Keep off until cohort-validated. |
| `bruker_ms1_polygon_mz_pad` | f64 | `0.0` | Expand each m/z edge of the selection polygon by this many Da, preserving isotope envelopes near an edge. |
| `bruker_ms1_polygon_im_pad` | f64 | `0.0` | Expand each ion-mobility edge of the selection polygon by this many 1/K0 units. |

### 3.2 `[hills]` — chromatographic-trace detection

| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_scans` | usize | `3` | Discard hills shorter than this. `3` on both. Interacts with `features.min_scan_overlap` (a 2-scan hill can't reach a 3-scan overlap). |
| `max_gap` | usize | `1` | Max internal scan gap within a hill. Every tuned config uses `1`. |
| `split_hills` | bool | `true` | Split merged traces at valleys (persistence splitter). `true` everywhere. |
| `split_valley_ratio` | f64 | `0.60` | A split survives only if the valley drops to ≤ this fraction of the *smaller* peak. `0.60` = precision-optimal; `0.70` leans recall. |
| `split_sigma_mult` | f64 | `5.0` | Notch-depth floor as a multiple of noise σ; rejects shallow noise notches. `5.0` = precision-optimal. |
| `split_height_frac` | f64 | `0.10` | Min peak height as fraction of the robust (95th-pct) max. `0.10` shipped; `0.05` leans recall (catches faint 10:1 co-eluters). |
| `lfc_weight` | f64 | `0.3` | Weight of the intensity log-fold-change term when matching a peak to a hill. `0.3` (helps noisy timsTOF). |
| `gap_fill_enabled` | bool | `false` | Interpolate intensity through internal zero gaps. **Keep `false`** on real configs (on+smoothing adds ~20% features for ~0.7 pp recall, worse redundancy). `relaxed` sets `true`. Warning: can create artificial maxima the splitter treats as new peaks. |
| `smoothing_enabled` | bool | `false` | Running-average the intensity profile. `false` on Orbitrap/Bruker; `true` in `relaxed`. |
| `smoothing_window` | usize | `1` | Half-width of the smoothing window (total = 2·w+1). `1`. Ignored when smoothing off. |
| `filter_large_baseline_hills` | bool | `false` | Drop baseline-like hills (≥`large_hill_min_scans` and no clear apex) — column bleed / plasticizers. **Comparator-only** (AlphaPept-ported); `true` only in `alphapept_like`. Overridable via `--filter-baseline`. |
| `large_hill_min_scans` | usize | `40` | Min scan span considered by the baseline-hill filter. Only meaningful when the filter is on. |
| `large_hill_peak_factor` | f64 | `2.0` | Required max/endpoint intensity ratio for a hill to count as "has an apex". Only meaningful when the filter is on. |
| `tic_norm_window` | usize | `0` (off) | **Experimental.** Pre-detection per-scan anomaly normalisation over a centred window (fixes brief ESI dropouts). `0` disables. Too small (≲ peak FWHM) flattens real apices; ~100 scans if used. |
| `tic_norm_min_scale` | f64 | `0.5` | Lower clamp on the per-scan multiplier. Active only when `tic_norm_window > 0`. |
| `tic_norm_max_scale` | f64 | `3.0` | Upper clamp on the per-scan multiplier. Active only when `tic_norm_window > 0`. |
| `tic_norm_mode` | String | `"median"` | Reference quantity: `"median"` (robust, recommended — fires only on whole-scan suppression) or `"tic"` (heavy-tailed, scales real apices down). |

#### `[hills_ms2]` — optional per-field overrides for MS2 (DIA fragment) hills

Optional table. Only consulted when `[file].ms2_hills_enabled = true` (DIA MS2
hill detection, written to `hills_ms2.*`). **Override semantics, not a fresh
config:** MS2 hill detection starts from your `[hills]` values and overwrites
**only** the keys you set under `[hills_ms2]`; every key you omit is inherited
from `[hills]` — it does **not** reset to the built-in default. If the whole
`[hills_ms2]` table is absent, MS2 hills use `[hills]` verbatim (byte-identical
to prior behavior — the default path is the regression guard).

Every key here is the `Option`-wrapped twin of the same-named `[hills]` key
(`min_scans`, `max_gap`, `split_hills`, `split_valley_ratio`, `split_sigma_mult`,
`split_height_frac`, `lfc_weight`, `gap_fill_enabled`, `smoothing_enabled`,
`smoothing_window`, `filter_large_baseline_hills`, `large_hill_min_scans`,
`large_hill_peak_factor`, `tic_norm_window`, `tic_norm_min_scale`,
`tic_norm_max_scale`, `tic_norm_mode`); see the `[hills]` table above for what
each does. `deny_unknown_fields` applies here too — a typo is a parse error. The
resolved MS2 config is computed by `KothConfig::ms2_hills()` and used by the
binary and every streaming pipeline entry point. Not set in any shipped config.

Typical use — fragment traces are often shorter than precursor traces, so relax
`min_scans` for MS2 only while keeping all other `[hills]` settings:

```toml
[hills_ms2]
min_scans = 2
```

### 3.3 `[features]` — isotope-chain assembly + charge + retention

| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_charge` | u8 | `2` | Min precursor charge. `2` (tryptic peptides are 2–5). |
| `max_charge` | u8 | `6` | Max precursor charge. `6`. |
| `min_chain_cosine` | f64 | `0.40` | Per-extension chromatographic-cosine gate while building a chain (vs the `cosine_anchor` reference). `0.40` recovers +1.24 pp recall at +19% features, neutral quant. **Not the contaminant filter** (co-eluters score ~0.74) — that's `min_isotope_score`. `relaxed` = 0.0; `alphapept_like` = 0.6. |
| `min_isotope_step_ratio` | f64 | `0.01` | Heavier isotope must be ≥ this fraction of its predecessor (chain-termination evidence). Chains extend **upward only**, so there is no downward counterpart (removed 2026-07-28). Named `right_max_decrease` before 0.3.0; the old key is still accepted as an alias. |
| `max_isotopes` | usize | `6` | Max isotope peaks added above the monoisotopic seed, so the envelope is at most `max_isotopes + 1` hills. `6`. (Before the downward walk was removed this bounded each direction separately, allowing envelopes up to `2·max_isotopes + 1`.) |
| `max_isotope_log2_ratio` | f64 | `1.5` | Intensity-ratio gate: after passing cosine, apex ratio vs predecessor must match averagine within ±this many log2 (≈2.83×). **This catches co-eluting contaminants** cosine misses. Default (not overridden). |
| `chain_predicted_intensity_gate` | bool | **`false`** | **Behavior-changing.** `true`: stop the chain when the averagine-*predicted* next-isotope intensity falls below the noise floor. `false` (default since the downward walk was removed, 2026-07-28): purely evidence-based termination. The gate only ever guarded the upward direction, so while the downward walk existed a dim monoisotope killed by it was still recovered by seeding its M+1 and stepping down. With one direction there is no second route, and leaving it on silently drops those features. |
| `min_scan_overlap` | usize | `3` | Min mutually-overlapping scans before two hills get a cosine (else 0, no extension). Keep `3` for PXD003881; lower to `2` on fast gradients (3–5-scan hills) to recover dim pairs. |
| `cosine_anchor` | String | `"seed"` | Which hill each isotope's cosine is measured against. **`"seed"`** (new 2026-07 default; anchor every isotope to the monoisotope — beat `"adjacent"` by +0.31 pp / +1565 PSMs on the 20-run cohort). `"adjacent"` reproduces pre-2026-07 paper output. Validated at load. |
| `sulfur_offsets` | list[int] | `[-1, 0, 1]` | Sulfur-count offsets, relative to the **ceiling** of the averagine-expected count, to score each chain against — one averagine template per offset, best Bhattacharyya kept. Corrects the systematic penalty on Cys/Met-rich peptides (³⁴S lifts M+2). Default resolves to `{0,1,2}` sulfurs up to 2665 Da, `{1,2,3}` to 5330, `{2,3,4}` to 7995. Negatives saturate at 0 and duplicates collapse. Widen to `[-2,-1,0,1]` to keep the no-sulfur template on large peptides (~32 % of 3 kDa peptides have none). **Empty list disables sulfur awareness.** Scoring is a max over templates, so a longer list can only raise scores incl. decoys — judge changes on recall, not the score distribution. `koth_ff_sulfur_{on,off}.toml` exist for A/B. |
| `isotope_model` | enum or table | `"peptide"` | Which analyte class's average composition the theoretical isotope pattern comes from. `"peptide"` is Senko's averagine (C₄.₉₃₈₄H₇.₇₅₈₃N₁.₃₅₇₇O₁.₄₇₇₃S₀.₀₄₁₇ / 111.1254 Da) and **the only model the published benchmark exercises**. `"rna"` (C₉.₅H₁₁.₇₅N₃.₇₅O₇ / 321.2916 Da) and `"dna"` (C₉.₇₅H₁₂.₂₅N₃.₇₅O₆ / 308.8006 Da) are unweighted means of the four chain residues; phosphorus is carried in the residue mass only, since ³¹P is monoisotopic and cannot shift a pattern. An explicit table overrides both: `isotope_model = { residue_mass = 321.2916, c = 9.5, h = 11.75, n = 3.75, o = 7.0 }` (`s` defaults to 0). A model without sulfur ignores `sulfur_offsets`. An unknown name is a parse error, never a silent fallback. On a PXD075396 RNase digest the RNA model raised the mean isotope score 0.821 → 0.848 and the share ≥0.90 from 33.8 % to 42.9 % against the peptide model. |
| `neutron_mass` | f64 | `1.003354835` | C13 mass offset for isotope-spacing targets. Default. |
| `exhaustive_min_isotope_score` | f64 | `0.0` | Minimum Bhattacharyya score a prefix needs to *claim* its hills. The resolver searches for the highest-evidence prefix that passes. `0.0` disables this claim gate; downstream retention filters still apply. |
| `isotope_evidence_ratio_sigma` | f64 | `0.75` | Signal standard deviation of seed-relative log2 apex-intensity-ratio errors, against a broad Normal(0, 2²) noise null. Must be finite and strictly between 0 and 2. Used for additive claim ranking and best-prefix selection; does not change the per-step ratio gate. |
| `isotope_evidence_cosine_shape` | f64 | `2.0` | Shape of the Beta(shape, 1) seed-anchored co-elution model against a uniform null. Must be finite and >1. Higher values favor tighter co-elution. These working likelihood models are not calibrated feature FDRs. |
| `exhaustive_isotope_priority` | bool | `false` | Legacy compatibility field; ignored. Contested-hill claims always rank by additive isotope log evidence. |

#### Final retention filters (AND-ed; drop the whole feature)
| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_isotope_score` | f64 | `0.5` | Drop features whose isotope (Bhattacharyya-vs-averagine) score is below this. **`0.5` is the real contaminant filter.** `0.0` in `relaxed`. |
| `min_cosine_score` | f64 | `0.0` | Drop features whose mean chromatographic-cosine score is below this. `0.0` (keep all). |
| `min_combined_score` | f64 | `0.0` | Drop features whose `isotope × cosine` combined score is below this. `0.0` (keep all). |

### 3.4 `[scoring]` — ±1-neutron monoisotope reassignment (optional)
Does **not** affect retention (use the `[features].min_*_score` gates for that).
Skippable entirely with `--no-scoring`.

| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `isotope_offset_enabled` | bool | `false` | Search neutron offsets [−1,+1] to reassign the monoisotope. `false` (test only offset 0) on all shipped configs. |
| `offset_zero_bonus` | f64 | `0.15` | Bhattacharyya bonus for keeping offset 0 when scores are close (internal only). Meaningful only when enabled. |
| `min_isotope_score_for_offset` | f64 | `0.5` | Features below this isotope score keep offset 0 (no reassignment). `0.5` (Orbitrap/Bruker/alphapept_like); `0.0` (relaxed). |

### 3.5 `[output]`
| Key | Type | Default → shipped | What it does |
|---|---|---|---|
| `format` | enum | `"tsv"` → **`"parquet"`** | Output format for hills+features. `"parquet"` (~10× smaller/faster for 20 runs); `"tsv"` for human-readable single-run debugging. |

---

## 4. `koth_align` — alignment + LFQ config (`AlignConfig`)

Sections: `[alignment]`, `[lfq]`, `[lfq.consensus]` (nested), `[output]`.

> **The align config is now IDENTICAL on Orbitrap and Bruker** (since
> 2026-08-27): `detected_use_grid = true` and `quant_estimator = "sum"` on both
> platforms. The former split (Bruker `apex` + all-grid, Orbitrap `sum` +
> detected-feature intensities) is gone — all-grid quantification was validated
> as a win on Orbitrap too, and sum-vs-apex is a wash on Bruker under all-grid.

### 4.1 `[alignment]` — RT/mass/IM alignment to an auto-selected reference
Anchor matching is coordinate-only (charge + ppm + RT window, no peptide ID), so
~25% of anchors are wrong-peptide matches — which is why the warp is RANSAC.

| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `reference_run` | optional string | unset | Exact run-directory name for a fixed reference. Unknown names are rejected; omission keeps automatic selection. |
| `anchor_mass_ppm` | f64 | `10.0` | ppm tolerance for forming an anchor pair (run feature ↔ reference feature). `10.0` both. Widen only for poorer mass accuracy. |
| `rt_anchor_window` | f64 | `0.05` | Normalised RT half-window `[0,1]` for candidate anchors (±5% ≈ ±6 min on a 2 hr gradient). `0.05` both. |
| `im_tolerance` | f64 | `0.05` | IM tolerance (1/K0) for anchors. Inert on Orbitrap; `0.05` both. |
| `min_anchor_combined_score` | f64 | `0.5` | Min feature `combined_score` (isotope × cosine) to be an anchor; also picks the reference run (most features clearing this). `0.5`. (Renamed from `min_anchor_score`.) |
| `min_anchor_count` | usize | `10` | Min anchors to fit a warp; below → identity warp + warning. `10` (safe floor). |
| `rt_warp_bandwidth` | f64 | `0.1` | Sliding-window width (fraction of RT range) for the piecewise-linear warp knots fit on RANSAC inliers. `0.1`. |
| `rt_warp_kind` | enum | `"ransac"` | **The only supported value** — `WarpKind` is single-variant. RANSAC line consensus (4000 deterministic iterations) → ≤4 refinement passes of median piecewise-linear + inlier re-collection → PAVA isotonic (monotone) warp. Matches the old piecewise fit on clean data but degrades far less under contamination (3× window: piecewise p95 RT residual ~8 min vs RANSAC ~0.7 min). |
| `rt_warp_ransac_thresh` | f64 | `0.01` | RANSAC inlier half-band (normalised RT; ≈1.4 min on a 142-min gradient). Omitted from shipped configs (uses default). Tighten to reject contamination, loosen if anchors are sparse. |

*Compile-time constants (not in TOML):* `RANSAC_ITERS=4000`, `RANSAC_REFINE_PASSES=4`;
and the RT-residual σ model (`RtSigmaModel`, consumed by `[lfq].rt_spread_scoring`):
`RT_SIGMA_BINS=12`, `RT_SIGMA_SHRINK_N0=20.0`, `RT_SIGMA_FLOOR_FRAC=0.25`.

### 4.2 `[lfq]` — XIC extraction, integration, decoy, TDC
For each consensus feature × run, build a `grid_cols × n_isotopes` XIC grid at the
alignment-predicted (m/z, RT, IM), score every RT column, integrate the best peak,
and (if `run_tdc`) repeat with a decoy.

| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `mz_ppm` | f64 | `10.0` | ppm tolerance for hill lookup per isotopologue. Consensus grouping has its own mass-span limit in `[lfq.consensus]`. |
| `rt_window_pct` | f64 | `0.005` | XIC extraction half-window, ±0.5% of each run's observed RT span, centred on the alignment-predicted native RT. A 142-minute span gives ±42.6 seconds. This is independent of the consensus RT-span limit. |
| `im_tolerance` | f64 | `0.015` | Absolute IM half-window (1/K0) for hill lookup, not a percentage. Inert on Orbitrap. Alignment and consensus use their own independent tolerances. |
| `n_isotopes` | usize | `3` | Grid isotope rows (1=M, 2=M+M1, 3=M+M1+M2). `3`. Clamped to ≥1. |
| `grid_cols` | usize | `100` | RT bins per grid. `100`. Clamped to ≥1. More = finer RT at higher cost. |
| `min_spectral_bhattacharyya` | f64 | `0.1` | Min spectral Bhattacharyya (observed vs averagine pattern) for a column to keep **extending** the integration peak (a gate, not a cosine). `0.1`. |
| `score_mode` | enum | `"hybrid"` | Which per-column score picks the integration window + feeds TDC ranking. **`"hybrid"`** = geo-mean `(rt·intensity·bhattacharyya·coelution)^¼`. `"rt"`/`"intensity"`/`"spectral"` use one signal (diagnostic). Note: the *apex column* is chosen by raw intensity regardless. |
| `run_tdc` | bool | `true` | Run target-decoy competition + compute per-cell q-values. `true`. `false` → q-values all 1.0, decoys zero. |
| `decoy_mz_shift_da` | f64 | `11.0` | Decoy m/z = target + `shift/charge`; must clear any real isotopologue/adduct. `11.0`. Tunable null-model knob (not shipped explicitly). |
| `decoy_rt_shift_pct` | f64 | `0.01` | Decoy RT = target RT − `shift × run_rt_span`. With the default extraction half-window of 0.005, the two RT windows meet at one boundary. A shift greater than twice the half-window fully separates the intervals. |
| `decoy_own_template` | bool | `true` | Score the +`decoy_mz_shift_da` decoy against **its own** averagine template (from the shifted mass), not the target's. Correctness fix (audit A2); measured neutral on its own but paired with `lone_coelution`. Target scoring + reported intensities unchanged. `false` reproduces pre-0.3.0 q-values. |
| `lone_coelution` | f64 | `0.5` | Co-elution value for a cell with <2 isotope rows carrying signal (a lone monoisotope — nothing to co-elute), applied identically to target and decoy. `1.0` (the pre-0.3.0 default) hands noise-grabbing lone-hill decoys a free target-like coordinate on the QDA co-elution feature; **`0.5` neutralises that freebie** — validated q-calibration win (audit A4: q-AUROC 0.934→0.938, +141 PSMs at q≤0.05, no quant cost). |
| `isotope_model` | enum or table | `"peptide"` | As `[features] isotope_model`, for the LFQ consensus templates. **Must match** the model the per-run features were detected with, or every cell is scored against a pattern the detector never used. |
| `normalize` | String | `"none"` | Cross-run matrix normalisation. **`"none"` for the paper** (the benchmark normalises every tool identically downstream, so an in-binary median-of-ratios would double-normalise unfairly). `"median_ratios"` = DESeq/edgeR size factors — a **validated option for standalone use** where you consume the matrix directly. |
| `quant_estimator` | String | `"sum"` (both platforms) | Per-cell estimator over the grid. Under all-grid quantification `"sum"` beats `"apex"` on Orbitrap (CV 12.31 vs 13.19 %, HUMAN IQR 0.194 vs 0.206) and the two are a wash on Bruker (CV 9.06 vs 9.19 %). The old "apex wrecks IQR 0.227→0.413" result was measured in the mixed detected-feature/grid regime and no longer applies. |
| `detected_use_grid` | bool | `true` (both platforms; **code default flipped 2026-08-27**) | Quantify **every** cell (detected + MBR) by the same grid re-integration — one estimator, one scale. Essential on timsTOF (feature integrates IM, 2-D grid doesn't; mixed scales inflated CV 38→13.7%) and a validated win on Orbitrap too (gated CV 14.27→12.31 %, ECOLI bias −0.067→−0.020, HUMAN IQR 0.219→0.194). `false` now warns and uses grid intensities because exclusive native-signal ownership requires grid quantification. |
| `rt_spread_scoring` | bool | `false` | Replace the raw RT term with a σ-normalised Gaussian likelihood using the per-run post-warp RT-residual spread (region-aware; strict where alignment is confident). Applied to target+decoy so TDC stays calibrated. **Validated but default-off**; omit unless experimenting. |
| `averagine_projection` | bool | `false` | Report the averagine matched-filter projection per cell instead of the raw box-sum (keeps on-pattern signal, rejects orthogonal contamination). **Tested negative on Orbitrap** (CV +3.1 pp, IQR +0.036, FFCR +1.9 pp, recall flat — see [§5](#5-experimental-knob-status-do-not-re-litigate)); untested on Bruker (its background-floor regime is where it might help). Byte-identical when off. **Keep `false`.** |

*Removed keys that now error:* `spectral_cosine_min` (→ `min_spectral_bhattacharyya`),
`spectral_bhattacharyya`, `spectral_coelution` (both signals always on),
`min_feature_score` (removed; see the group confidence gate below).

**q-value mechanics.** With `run_tdc=true` the active path is a semi-supervised,
cross-validated **QDA rescorer** (`lfq/rescore.rs`, Percolator-style protocol)
over five symmetric per-cell features: |ppm error|, |RT diff|, spectral
Bhattacharyya, co-elution cosine, |IM delta|. Every target cell (detected **or**
MBR) competes against decoys — a detected cell in a depleted well is still
gatable. Deterministic (no RNG). The simpler ranker in `lfq/tdc.rs` is retained
as reference but not used by the pipeline.

### 4.3 `[lfq.consensus]` — experimental cross-run confidence

The [independent-permutation audit](lfq-group-permutation.md) replaces the initial
whole-run shifts with stable controls on the tested cohorts. The group gate is
still experimental: `0.05` is not a validated 5% identification or member-link FDR.

All finite scored features with known charge and valid projected coordinates can
enter bounded candidate groups, including scores below 0.5. The highest-quality
member still provides reference coordinates, but its score is not an admission
threshold. Each run supplies at most one primary observation; alternatives stay
available in the long bundle and dilute ambiguous group support.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `mz_ppm` | f64 | `20.0` | Maximum full group mass span in ppm. |
| `rt_window_pct` | f64 | `0.02` | Maximum full group RT span as a fraction of reference run span. |
| `im_tolerance` | f64 | `0.05` | Maximum full IM span; missing IM cannot bridge incompatible values. |
| `max_group_qvalue` | f64 | `0.05` | Experimental permuted-RT group gate. `1.0` keeps every replicated candidate for auditing. |

`min_member_combined_score`, `min_seed_combined_score`, `min_group_size` and
`allow_replicated_weak_seeds` are removed and rejected by the parser. Delete them
from old configs and set `max_group_qvalue`. Two distinct original runs are
required for cross-run support; a single-run dataset yields no consensus groups.

Group evidence combines continuous detector quality, pairwise mass/RT/IM
agreement and ambiguity from within-run alternatives. Ten deterministic,
independent RT permutations within each run, charge and detector-score bin
(width 0.1) preserve that stratum's exact RT distribution. Mass, IM, quality and
feature provenance remain unchanged. The controls undergo identical grouping
and scoring; fixed points are allowed. This replaces the initial whole-run
circular translations, whose control counts were unstable with small run counts.
At each score cutoff, estimated false groups equal one plus the average null
count; division by the target count and a reverse cumulative minimum give the
group q-value. The log reports retained counts from the two five-control halves
and the pooled ten-control estimate. These internal controls add no user settings.

This is not a learned probability model or a validated FDR guarantee. Permutations
do not preserve within-stratum mass–RT or IM–RT dependence, or within-run group
structure. Concentrated RT distributions can provide no useful discrimination;
recurring artifacts and false attachments to real groups require separate
validation. Group confidence does not identify peptides or establish every link.

The existing per-cell extraction and gate still apply. `group_score` and
`group_qvalue` are exported in the consensus table and long bundle, separately
from `lfq_q_value`. No group-confidence value is reused as cell confidence.
See [lfq-group-permutation.md](lfq-group-permutation.md) for the current experiment.

### 4.4 `[output]`
| Key | Type | Default → shipped | What it does |
|---|---|---|---|
| `export_long` | bool | `false` | Write schema-versioned original feature observations, recomputed LFQ cells and a manifest with source/output hashes. See [long-matrix.md](long-matrix.md). |
| `format` | enum | `"tsv"` → **`"parquet"`** | Matrix/consensus output format. `"parquet"` shipped. |
| `max_qvalue` | f64 | `1.0` | Only write matrix entries with q ≤ this. **Required key** (no serde default). `1.0` = emit all, filter downstream. |
| `export_decoys` | bool | `true` | Write `decoy_intensity_matrix`. No effect when `run_tdc=false`. `true`. |
| `export_details` | bool | `true` | Write `lfq_details` (long-format: one row per feature×run×is_decoy, with scores/RT-diff/observed m/z+IM). `true`. |

---

## 5. Experimental knob status — do not re-litigate

These have been tested; the verdicts are settled. **Keep them at their defaults**
unless you are deliberately re-opening the experiment (and re-benchmarking).

| Knob | Section | Verdict | Default |
|---|---|---|---|
| `mz_recalibration` (adaptive tol) | `[file]` | ✅ **Validated win**, both platforms (the adaptive tolerance, not position-only) | on; code default `true` since 0.3.0 |
| `detected_use_grid` | `[lfq]` | ✅ **Validated win on both platforms** (timsTOF CV 38→13.7%; Orbitrap gated CV 14.27→12.31%) | On everywhere; code default `true` |
| `normalize = "median_ratios"` | `[lfq]` | ✅ Valid for **standalone** use; off in paper for fairness | `"none"` |
| `rt_spread_scoring` | `[lfq]` | ✅ Validated, but kept **default-off** | off |
| `cosine_anchor = "seed"` | `[features]` | ✅ **Won**, now the default (+0.31 pp) | `"seed"` |
| `sulfur_offsets` | `[features]` | ✅ Sulfur-aware scoring is a small win (+0.06–0.2 pp) and stays on. Offsets are ceiling-relative as of 2026-07-28 — the old fixed set `{0, avg, avg+2, avg+4}` resolved to `{0,1,3,5}` over 1332–3997 Da and **skipped n_S = 2** (~7.6 % of peptides) while spending a slot on n_S = 5 (~0.02 %). Requires the pending data re-run (paper `TODO.md` item 6) to requantify. | `[-1, 0, 1]` |
| `averagine_projection` | `[lfq]` | ❌ **Dud on Orbitrap** (CV +3.1 pp, IQR +0.036, FFCR +1.9 pp, recall flat); Bruker untested | off |
| `chain_predicted_intensity_gate = false` | `[features]` | ✅ **Now the default** — required once the downward chain walk was removed | `false` |
| `bruker_streaming` | `[file]` | ⚪ Experimental; parity-verified vs local path but not yet cohort-validated | off |
| `tic_norm_*` | `[hills]` | ⚪ Experimental (ESI-dropout salvage) | off |
| `gap_fill_enabled` / `smoothing_enabled` | `[hills]` | ❌ Net-negative for quant (adds features, hurts redundancy); only in `relaxed` | off |
| position-only m/z recalibration | `[file]` | ❌ No-op on well-calibrated instruments | (subsumed by adaptive) |

---

## 6. Platform presets: Orbitrap vs Bruker/timsTOF

Start from the shipped configs; the deltas are minimal.

| Setting | Orbitrap (`koth_ff.toml` / `koth_align.toml`) | Bruker/timsTOF (`*_bruker.toml`) | Why |
|---|---|---|---|
| `[file].mz_tolerance` | `8.0` | `15.0` | timsTOF MS1 mass accuracy is looser |
| `[file].n_threads` | unset | unset | all cores on both platforms; v0.1.0 ignored the key, current builds honor it |
| `[file].bruker_*` front-end | inert (mzML) | active (`.d` filter + watershed) | Bruker raw-frame denoising |
| `[lfq].quant_estimator` | `"sum"` | `"sum"` | identical since 2026-08-27 (sum-vs-apex is a wash under all-grid) |
| `[lfq].detected_use_grid` | `true` | `true` | identical since 2026-08-27; all-grid everywhere |

Everything else — splitter params, `min_chain_cosine=0.40`, recalibration,
`sulfur_offsets`, `min_isotope_score=0.5`, all alignment/consensus/TDC
settings — is **identical across platforms**.

---

## 7. The dnoise integration (Bruker only)

koth_ff depends on the published `dnoise` 0.1 crate (behind the `tdf` feature)
for Bruker denoising. Two paths:

- **Default (paper-validated):** the `.d` is denoised **externally** by the
  `dnoise` CLI first, then koth_ff reads the denoised `.d` and applies its
  built-in vertical filter + watershed (`bruker_streaming = false`). koth_ff only
  exposes the vertical-filter (`bruker_filter_*`) and watershed
  (`bruker_watershed_*`) knobs.
- **Streaming (opt-in):** `bruker_streaming = true` drives dnoise's `RunContext`
  in-process (vertical → optional halo → optional MS1 selection polygon →
  watershed in one pass, no denoised `.d` on disk), adding the `bruker_halo*`
  and `bruker_ms1_polygon*` knobs. Keep off until cohort-re-validated.

---

## 8. Quick recipes

- **Maximise depth (more features / recall):** lower `[features].min_chain_cosine`
  (→ 0.0), zero the retention gates (`min_isotope_score=0`), enable
  `[hills].gap_fill_enabled` + `smoothing_enabled`, lower `split_height_frac` to
  0.05. This is essentially `koth_ff_relaxed.toml`. Expect worse quant/redundancy.
- **Tighter quant (lower CV):** keep the shipped precision-optimal splitter
  (`split_valley_ratio=0.60`, `split_sigma_mult=5.0`, `split_height_frac=0.13`),
  `intensity_coverage=1.0`, `min_isotope_score=0.5`. For koth_align keep
  `quant_estimator="sum"` (Orbitrap). The default extraction half-window is
  `rt_window_pct=0.005`; evaluate precision and completeness for the intended use.
- **Fast/short gradients (3–5-scan hills):** lower `[hills].min_scans` to 2 and
  `[features].min_scan_overlap` to 2; re-benchmark quant.
- **MS1-search feature depth (dim 2+/3+):** `[features].chain_predicted_intensity_gate`
  is already `false` by default; setting it `true` re-enables the predicted-intensity
  early-stop and will drop dim monoisotopes.
- **Match AlphaPept for comparison:** use `koth_ff_alphapept_like.toml`
  (`intensity_coverage=0.95`, `min_chain_cosine=0.6`, `filter_large_baseline_hills=true`).

---

*Generated from the config structs and shipped TOMLs. If you add or rename a
config field, update this file, `example_config*.toml` here, and the tuned
`benchmark/config/*.toml` in [koth-paper](https://github.com/tacular-omics/koth-paper). The
`deny_unknown_fields` parse test fails the build if the templates in THIS repo
disagree with the structs; it can no longer see the tuned configs, and it cannot
check this doc — keep both current by hand.*

### Experimental exclusive LFQ extraction

The group-confidence worktree reserves native hill-profile peak segments across
competing LFQ groups and rebuilds residual extractions before cell scoring. This
requires grid quantification: `detected_use_grid=false` now warns and uses the
grid estimator. There are no additional user settings. See
[signal ownership](lfq-signal-ownership.md) for behavior, provenance and limits.

### Search-guided CLI options (0.8.0)

Experimental `[lfq]` settings in the LFQ-performance worktree:

- `search_scoring = "legacy"` (default) preserves release scoring exactly.
  `"peptide_grouped"` assigns all runs and charges of each modified peptide to
  the same deterministic fold and fits imputation/scaling on training data only.
  `"charge_stratified"` additionally fits 2+, 3+, and other charges separately
  within direct and transferred evidence classes. If any present charge bin
  lacks 100 positive targets, 100 positive decoys, or 20 peptide groups in any
  training fold, the entire evidence class falls back to peptide-grouped pooled
  scoring. Candidate q-values count score ties together and use a +1 decoy
  correction. Compare grouped versus stratified to isolate charge splitting;
  the grouped-versus-legacy contrast includes the validation/calibration changes.
- `search_scoring = "charge_signed"` uses signed mass/RT residuals in the same
  charge-stratified, peptide-grouped QDA. `"charge_signed_quality"` also adds
  log1p isotope count, log1p nonnegative peak width, and preceding-isotope signal
  fraction. Features use each target or decoy's own extraction coordinates;
  preprocessing remains training-fold-only. Sparse strata use the same pooled
  fallback with the selected feature definition. Both are experimental,
  default-off, development candidates. Rounded replay motivates native testing;
  it does not validate transfer FDR or establish FlashLFQ parity.
- `search_rt_rescue = false` (default) preserves target rejection. When true,
  non-IM targets can retain supported native RT clusters: an ambiguous run needs
  a unique bounded cluster containing a strict majority of its PSMs and at least
  two observations. Unsupported runs are withheld. Any excluded PSM or remaining
  cross-run RT conflict blocks transfers for that target; unambiguous same-run
  IDs can still be extracted. `search_manifest.json` records excluded source rows
  and cross-run ambiguity. Mass-conflicting or IM-bearing targets retain the
  legacy rejection policy. This conservative first rescue does not resolve
  chromatographic isomers or enable transfers from ambiguous targets.

- `search_expand_charges = false` (default). When enabled, unambiguous non-IM
  peptides are queried at charges 2–4 within the run's observed feature/hill m/z
  range. This conservative range may be narrower than the instrument's acquired
  range. Native PSM RTs at other charges anchor additional same-run queries even
  with `--no-mbr`; inferred charges are never direct MS2 identifications. Native
  exports label `inferred_charge` or `mbr_inferred_charge` and retain real seed
  PSMs at their original charges in the manifest. Require at least two isotopes
  and co-elution >=0.5 for inferred targets and paired decoys before scoring.
  Reuse signal ownership to avoid duplicate signal. Other observed charges remain
  available. This bounded ablation is not an exact FlashLFQ charge-search replica.

These are development ablations, not validated improvements or calibrated FDR
claims. Identification-free scoring is unchanged. Apex-versus-sum comparisons
for search-guided LFQ reopen that estimator question only in this new workflow.

`koth_align --sage-psms results.sage.tsv` (or `.parquet`) and
`--targets targets.tsv` select optional identification-guided quantification.
`--max-id-qvalue` and `--max-extraction-qvalue` both default to `0.01`;
`--no-mbr` disables transfers and `--ignore-target-im` ignores imported IM.
These are CLI options, not TOML fields. Existing LFQ extraction settings apply;
search targets bypass identification-free group support and group-q gating.
The consensus coordinate-span limits still reject ambiguous peptide targets.
See [search-guided LFQ](search-guided-lfq.md) for schemas, output confidence
semantics, alignment requirements, and validation limits.
