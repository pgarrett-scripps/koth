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
  `--no-scoring`, `--ms2-hills`, `--recalibrate` (→ `[file].mz_recalibration`),
  `--filter-baseline` (→ `[hills].filter_large_baseline_hills`).
- The canonical worked examples are the shipped configs in `benchmark/config/`:
  `koth_ff.toml` / `koth_ff_bruker.toml` (feature finding, Orbitrap / timsTOF)
  and `koth_align.toml` / `koth_align_bruker.toml` (alignment+LFQ). Prefer these
  over the older `example_config*.toml` templates.

**Notation below:** `key` (type, `default` → `shipped`) where "shipped" is the
value in the production Orbitrap config when it differs from the struct default.
Platform differences are called out inline.

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
| `im_tolerance` | f64 | `0.05` | Ion-mobility tolerance. `0.05` everywhere; inert on Orbitrap but the key is schema-required. |
| `im_tolerance_type` | enum | `"relative"` | `"relative"` (fraction of IM value) or `"absolute"` (1/K0). `"relative"` on both. |
| `global_min_mz` | f64 | `0.0` | Ignore peaks below this m/z. Left at default. |
| `global_max_mz` | f64 | `inf` | Ignore peaks above this m/z. Left at default. |
| `intensity_coverage` | f64 | `1.0` | Fraction of hill intensity to retain. **Keep `1.0`** — validated strictly dominant (0.95→1.0 gave +2.9 pp recall at +0.5% features; 0.95 clips low-abundance apices). Only the AlphaPept-comparator config uses 0.95. |
| `n_threads` | Option\<usize\> | `None` (all cores) | Worker threads. Shipped 4 (Orbitrap) / 8 (Bruker) for benchmark reproducibility; omit for max throughput. |

#### General noise / mode toggles
| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `noise_filter_sigma` | Option\<f64\> | `None` | Per-scan sigma-clip noise filter during hill detection (any format). Omit to disable (all shipped configs). Typical if used: `3.0`. |
| `decoy_mode` | bool | `false` | Shuffle MS1 scans before detection → null/decoy feature set. `false` normally; `true` only for FDR null-building. (Skips `mz_recalibration`.) |
| `ms2_hills_enabled` | bool | `false` | Also detect MS2 hills per isolation window (DIA), written to `hills_ms2.*`. `false` — shipped data is DDA. Bruker `.d` MS2 not yet supported. Overridable via `--ms2-hills`. |

#### ID-free m/z recalibration + region-adaptive isotope tolerance
Enabled in both production configs (**code default is off; the tuned configs turn
it on**). Pass-1 detection learns a per-(m/z, RT)-region signed-ppm offset surface
from isotope-spacing residuals; pass-2 shifts the *expected* isotope position by
the learned offset and replaces the fixed isotope-match tolerance with
`clamp(sigma_mult × σ(m/z,RT), floor_ppm, mz_tolerance)`.

| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `mz_recalibration` | bool | `false` → **`true`** | Master switch. The **region-adaptive tolerance is the real win** (position-only recal is a documented no-op on well-calibrated instruments). Orbitrap: −5.8% spurious features, recall flat, CV 23.76→23.25%. Bruker (bigger win, since fixed 15 ppm ≫ real ~2.8 ppm spread): recall +1.06 pp, CV 13.89→13.76%. Overridable via `--recalibrate`. |
| `mz_recalibration_mz_bins` | usize | `20` | m/z bins in the recalibration surface. Default. |
| `mz_recalibration_rt_bins` | usize | `8` | RT bins in the recalibration surface. Default. |
| `mz_recalibration_min_samples` | usize | `50` | Min residuals a cell needs before its own median is trusted (else falls back to marginal → global median). Default. |
| `mz_recalibration_tol_sigma_mult` | f64 | `4.0` | The `N` in `N × σ` adaptive tolerance. **`4.0`**, finalized 2026-07-09 (was 3.0): a balanced single default within ~0.2 pp of each platform's optimum (5.0 Orbitrap / 3.0 Bruker). |
| `mz_recalibration_tol_floor_ppm` | f64 | `1.0` | Lower ppm bound so a tiny σ in a sparse region can't collapse the window. Default. |

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

#### Bruker streaming path (experimental — see also [dnoise](#7-the-dnoise-integration-bruker-only))
| Key | Type | Default | What it does / what to set |
|---|---|---|---|
| `bruker_streaming` | bool | `false` | **Experimental / opt-in.** `true` = run dnoise's full pipeline (vertical filter → halo → watershed) in-process, no denoised `.d` on disk. `false` = the historical local path (vertical + watershed, **no halo**) — **the paper's validated pipeline**. Keep `false` until re-validated on the Bruker cohort. |
| `bruker_halo` | bool | `true` | Streaming-only: apply the horizontal-halo filter after the vertical filter. No effect unless `bruker_streaming = true`. |
| `bruker_halo_peak_fraction` | f64 | `0.15` | Streaming-only: drop a peak below this fraction of the off-column box-max. |
| `bruker_halo_mz_idx_half_width` | u32 | `80` | Streaming-only: halo reference-box half-width along TOF index. |
| `bruker_halo_scan_half_width` | usize | `2` | Streaming-only: halo reference-box half-width along the IM scan axis. |

### 3.2 `[hills]` — chromatographic-trace detection

| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_scans` | usize | `3` | Discard hills shorter than this. `3` on both. Interacts with `features.min_scan_overlap` (a 2-scan hill can't reach a 3-scan overlap). |
| `max_gap` | usize | `0` → **`1`** | Max internal scan gap within a hill. Every tuned config uses `1`. |
| `split_hills` | bool | `true` | Split merged traces at valleys (persistence splitter). `true` everywhere. |
| `split_valley_ratio` | f64 | `0.70` → **`0.60`** | A split survives only if the valley drops to ≤ this fraction of the *smaller* peak. `0.60` = precision-optimal; `0.70` leans recall. |
| `split_sigma_mult` | f64 | `4.0` → **`5.0`** | Notch-depth floor as a multiple of noise σ; rejects shallow noise notches. `5.0` = precision-optimal. |
| `split_height_frac` | f64 | `0.10` → **`0.13`** | Min peak height as fraction of the robust (95th-pct) max. `0.13` shipped; `0.05` leans recall (catches faint 10:1 co-eluters). |
| `lfc_weight` | f64 | `0.5` → **`0.3`** | Weight of the intensity log-fold-change term when matching a peak to a hill. `0.3` (helps noisy timsTOF). |
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
| `min_charge` | u8 | `1` → **`2`** | Min precursor charge. `2` (tryptic peptides are 2–5). |
| `max_charge` | u8 | `7` → **`6`** | Max precursor charge. `6`. |
| `min_chain_cosine` | f64 | `0.5` → **`0.40`** | Per-extension chromatographic-cosine gate while building a chain (vs the `cosine_anchor` reference). `0.40` recovers +1.24 pp recall at +19% features, neutral quant. **Not the contaminant filter** (co-eluters score ~0.74) — that's `min_isotope_score`. `relaxed` = 0.0; `alphapept_like` = 0.6. |
| `right_max_decrease` | f64 | `0.05` → **`0.01`** | Heavier isotope must be ≥ this fraction of its predecessor (chain-termination evidence). `0.01` (production). Chains extend **upward only**, so there is no `left_max_decrease` counterpart (removed 2026-07-28). |
| `max_isotopes` | usize | `6` | Max isotope peaks added above the monoisotopic seed, so the envelope is at most `max_isotopes + 1` hills. `6`. (Before the downward walk was removed this bounded each direction separately, allowing envelopes up to `2·max_isotopes + 1`.) |
| `max_isotope_log2_ratio` | f64 | `1.5` | Intensity-ratio gate: after passing cosine, apex ratio vs predecessor must match averagine within ±this many log2 (≈2.83×). **This catches co-eluting contaminants** cosine misses. Default (not overridden). |
| `chain_predicted_intensity_gate` | bool | **`false`** | **Behavior-changing.** `true`: stop the chain when the averagine-*predicted* next-isotope intensity falls below the noise floor. `false` (default since the downward walk was removed, 2026-07-28): purely evidence-based termination. The gate only ever guarded the upward direction, so while the downward walk existed a dim monoisotope killed by it was still recovered by seeding its M+1 and stepping down. With one direction there is no second route, and leaving it on silently drops those features. |
| `min_scan_overlap` | usize | `3` | Min mutually-overlapping scans before two hills get a cosine (else 0, no extension). Keep `3` for PXD003881; lower to `2` on fast gradients (3–5-scan hills) to recover dim pairs. |
| `cosine_anchor` | String | `"seed"` | Which hill each isotope's cosine is measured against. **`"seed"`** (new 2026-07 default; anchor every isotope to the monoisotope — beat `"adjacent"` by +0.31 pp / +1565 PSMs on the 20-run cohort). `"adjacent"` reproduces pre-2026-07 paper output. Validated at load. |
| `sulfur_offsets` | list[int] | `[-1, 0, 1]` | Sulfur-count offsets, relative to the **ceiling** of the averagine-expected count, to score each chain against — one averagine template per offset, best Bhattacharyya kept. Corrects the systematic penalty on Cys/Met-rich peptides (³⁴S lifts M+2). Default resolves to `{0,1,2}` sulfurs up to 2665 Da, `{1,2,3}` to 5330, `{2,3,4}` to 7995. Negatives saturate at 0 and duplicates collapse. Widen to `[-2,-1,0,1]` to keep the no-sulfur template on large peptides (~32 % of 3 kDa peptides have none). **Empty list disables sulfur awareness.** Scoring is a max over templates, so a longer list can only raise scores incl. decoys — judge changes on recall, not the score distribution. `koth_ff_sulfur_{on,off}.toml` exist for A/B. |
| `neutron_mass` | f64 | `1.003354835` | C13 mass offset for isotope-spacing targets. Default. |
| `exhaustive_min_isotope_score` | f64 | `0.0` | Exhaustive-resolver knob (unmerged experiment): min Bhattacharyya a candidate needs to *claim* its hills. `0.0` = no gate. Not set in shipped configs. |
| `exhaustive_isotope_priority` | bool | `false` | Exhaustive-resolver knob: order contested-hill claims by envelope length → isotope score → composite. `false`. Not set in shipped configs. |

#### Final retention filters (AND-ed; drop the whole feature)
| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_isotope_score` | f64 | `0.0` → **`0.5`** | Drop features whose isotope (Bhattacharyya-vs-averagine) score is below this. **`0.5` is the real contaminant filter.** `0.0` in `relaxed`. |
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

> **The only Orbitrap-vs-Bruker differences in the whole align config are two
> `[lfq]` keys:** Bruker sets `quant_estimator = "apex"` and
> `detected_use_grid = true`. Everything else is identical.

### 4.1 `[alignment]` — RT/mass/IM alignment to an auto-selected reference
Anchor matching is coordinate-only (charge + ppm + RT window, no peptide ID), so
~25% of anchors are wrong-peptide matches — which is why the warp is RANSAC.

| Key | Type | Default | What it does / what to set |
|---|---|---|---|
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
| `mz_ppm` | f64 | `10.0` | ppm tolerance for hill lookup per isotopologue. `10.0` both. **Also reused by consensus grouping at 2× (20 ppm)** to absorb alignment drift. |
| `rt_window_pct` | f64 | `0.02` → **`0.01`** | XIC grid half-window as a fraction of RT span. **Shipped `0.01` (±~1.2 min), halved from the code default:** the generous ±2.4 min window summed contaminant hills (+11.5% MBR bias, 2× scatter); `0.01` cuts bias to −2.1% and matrix CV 16.17→15.33 with no MV loss; `0.005` over-tightens (−24%). Also reused by consensus grouping at 2×. |
| `im_tolerance` | f64 | `0.05` | IM tolerance (1/K0) for hill lookup. Inert on Orbitrap. Reused by consensus grouping at 1×. |
| `n_isotopes` | usize | `3` | Grid isotope rows (1=M, 2=M+M1, 3=M+M1+M2). `3`. Clamped to ≥1. |
| `grid_cols` | usize | `100` | RT bins per grid. `100`. Clamped to ≥1. More = finer RT at higher cost. |
| `min_spectral_bhattacharyya` | f64 | `0.1` | Min spectral Bhattacharyya (observed vs averagine pattern) for a column to keep **extending** the integration peak (a gate, not a cosine). `0.1`. |
| `score_mode` | enum | `"hybrid"` | Which per-column score picks the integration window + feeds TDC ranking. **`"hybrid"`** = geo-mean `(rt·intensity·bhattacharyya·coelution)^¼`. `"rt"`/`"intensity"`/`"spectral"` use one signal (diagnostic). Note: the *apex column* is chosen by raw intensity regardless. |
| `run_tdc` | bool | `true` | Run target-decoy competition + compute per-cell q-values. `true`. `false` → q-values all 1.0, decoys zero. |
| `decoy_mz_shift_da` | f64 | `11.0` | Decoy m/z = target + `shift/charge`; must clear any real isotopologue/adduct. `11.0`. Tunable null-model knob (not shipped explicitly). |
| `decoy_rt_shift_pct` | f64 | `0.01` | Decoy RT = target − `shift × rt_span`. At `0.01` the decoy window still overlaps the target; increase to separate. Tunable null-model knob. |
| `decoy_own_template` | bool | `false` (code) / **`true` (shipped)** | Score the +`decoy_mz_shift_da` decoy against **its own** averagine template (from the shifted mass), not the target's. Correctness fix (audit A2); measured neutral on its own but paired with `lone_coelution`. Code default `false` for byte-safety; enabled in the shipped `koth_align.toml`. Target scoring + reported intensities unchanged. |
| `lone_coelution` | f64 | `1.0` (code) / **`0.5` (shipped)** | Co-elution value for a cell with <2 isotope rows carrying signal (a lone monoisotope — nothing to co-elute), applied identically to target and decoy. `1.0` hands noise-grabbing lone-hill decoys a free target-like coordinate on the QDA co-elution feature; **`0.5` neutralises that freebie** — validated q-calibration win (audit A4: q-AUROC 0.934→0.938, +141 PSMs at q≤0.05, no quant cost). Code default `1.0` (byte-safe); shipped `0.5`. |
| `normalize` | String | `"none"` | Cross-run matrix normalisation. **`"none"` for the paper** (the benchmark normalises every tool identically downstream, so an in-binary median-of-ratios would double-normalise unfairly). `"median_ratios"` = DESeq/edgeR size factors — a **validated option for standalone use** where you consume the matrix directly. |
| `quant_estimator` | String | `"sum"` → **Bruker `"apex"`** | Per-cell estimator. **Orbitrap `"sum"`** (integrated area — best fold-change accuracy; `"apex"` cuts CV but wrecks IQR 0.227→0.413 because MBR apex columns sit at the jittery consensus RT). **Bruker `"apex"`** — *for comparator fairness* (AlphaPept's `.d` output exposes only apex), not quality. |
| `detected_use_grid` | bool | `false` → **Bruker `true`** | Quantify **every** cell (detected + MBR) by the same grid re-integration. **Bruker-specific validated win** (CV 38→13.7%): on timsTOF the per-run feature integrates IM but the 2-D grid doesn't, so mixing scales inflates CV. Orbitrap `false` (trust the detected feature's own intensity). |
| `rt_spread_scoring` | bool | `false` | Replace the raw RT term with a σ-normalised Gaussian likelihood using the per-run post-warp RT-residual spread (region-aware; strict where alignment is confident). Applied to target+decoy so TDC stays calibrated. **Validated but default-off**; omit unless experimenting. |
| `averagine_projection` | bool | `false` | Report the averagine matched-filter projection per cell instead of the raw box-sum (keeps on-pattern signal, rejects orthogonal contamination). **Tested negative on Orbitrap** (CV +3.1 pp, IQR +0.036, FFCR +1.9 pp, recall flat — see [§5](#5-experimental-knob-status-do-not-re-litigate)); untested on Bruker (its background-floor regime is where it might help). Byte-identical when off. **Keep `false`.** |

*Removed keys that now error:* `spectral_cosine_min` (→ `min_spectral_bhattacharyya`),
`spectral_bhattacharyya`, `spectral_coelution` (both signals always on),
`min_feature_score` (→ `[lfq.consensus].min_member_combined_score`).

**q-value mechanics.** With `run_tdc=true` the active path is a semi-supervised,
cross-validated **QDA rescorer** (`lfq/rescore.rs`, Percolator-style protocol)
over five symmetric per-cell features: |ppm error|, |RT diff|, spectral
Bhattacharyya, co-elution cosine, |IM delta|. Every target cell (detected **or**
MBR) competes against decoys — a detected cell in a depleted well is still
gatable. Deterministic (no RNG). The simpler ranker in `lfq/tdc.rs` is retained
as reference but not used by the pipeline.

### 4.3 `[lfq.consensus]` — cross-run grouping quality filters
Tolerances are **not** here — grouping reuses `[lfq]`'s `mz_ppm`/`rt_window_pct`
(at 2×) and `im_tolerance` (at 1×). Features are projected into reference space
and single-linkage grouped by (charge, neutral mass, aligned RT, IM).

| Key | Type | Default → shipped | What it does / what to set |
|---|---|---|---|
| `min_member_combined_score` | f64 | `0.0` → **`0.5`** | Pre-grouping filter: a feature's `combined_score` must clear this to join/seed a group (else excluded, doesn't count toward `n_contributing_runs`). **Required key when the `[lfq.consensus]` table is present** (no serde default). Shipped `0.5`. |
| `min_group_size` | usize | `1` → **`2`** | Min distinct runs that must detect a feature to keep the group. **`2`** (of 20) shipped. `1` = full MBR; raise for stricter reproducibility. |
| `min_seed_combined_score` | f64 | `0.0` → **`0.75`** | Post-grouping filter: drop a group whose best member (seed) is below this. Shipped `0.75`. |

### 4.4 `[output]`
| Key | Type | Default → shipped | What it does |
|---|---|---|---|
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
| `mz_recalibration` (adaptive tol) | `[file]` | ✅ **Validated win**, both platforms (the adaptive tolerance, not position-only) | on in shipped configs |
| `detected_use_grid` | `[lfq]` | ✅ **Validated win on timsTOF** (CV 38→13.7%) | Bruker on, Orbitrap off |
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
| `[file].n_threads` | `4` | `8` | (benchmark reproducibility only) |
| `[file].bruker_*` front-end | inert (mzML) | active (`.d` filter + watershed) | Bruker raw-frame denoising |
| `[lfq].quant_estimator` | `"sum"` | `"apex"` | match what comparators expose per platform |
| `[lfq].detected_use_grid` | `false` | `true` | IM-vs-2D-grid scale commensurability |

Everything else — splitter params, `min_chain_cosine=0.40`, recalibration,
`sulfur_offsets`, `min_isotope_score=0.5`, all alignment/consensus/TDC
settings — is **identical across platforms**.

---

## 7. The dnoise integration (Bruker only)

koth_ff depends on the sibling `dnoise` crate (`../../d_noise`, behind the `tdf`
feature) for Bruker denoising. Two paths:

- **Default (paper-validated):** the `.d` is denoised **externally** by the
  `dnoise` CLI first, then koth_ff reads the denoised `.d` and applies its
  built-in vertical filter + watershed (`bruker_streaming = false`). koth_ff only
  exposes the vertical-filter (`bruker_filter_*`) and watershed
  (`bruker_watershed_*`) knobs — **not** dnoise's full setting surface (halo,
  smoothing, polygon gate, etc.).
- **Streaming (opt-in):** `bruker_streaming = true` drives dnoise's `RunContext`
  in-process (vertical → halo → watershed in one pass, no denoised `.d` on disk),
  adding the `bruker_halo*` knobs. Keep off until cohort-re-validated.

---

## 8. Quick recipes

- **Maximise depth (more features / recall):** lower `[features].min_chain_cosine`
  (→ 0.0), zero the retention gates (`min_isotope_score=0`), enable
  `[hills].gap_fill_enabled` + `smoothing_enabled`, lower `split_height_frac` to
  0.05. This is essentially `koth_ff_relaxed.toml`. Expect worse quant/redundancy.
- **Tighter quant (lower CV):** keep the shipped precision-optimal splitter
  (`split_valley_ratio=0.60`, `split_sigma_mult=5.0`, `split_height_frac=0.13`),
  `intensity_coverage=1.0`, `min_isotope_score=0.5`. For koth_align keep
  `quant_estimator="sum"` (Orbitrap) and `rt_window_pct=0.01`.
- **Fast/short gradients (3–5-scan hills):** lower `[hills].min_scans` to 2 and
  `[features].min_scan_overlap` to 2; re-benchmark quant.
- **MS1-search feature depth (dim 2+/3+):** `[features].chain_predicted_intensity_gate`
  is already `false` by default; setting it `true` re-enables the predicted-intensity
  early-stop and will drop dim monoisotopes.
- **Match AlphaPept for comparison:** use `koth_ff_alphapept_like.toml`
  (`intensity_coverage=0.95`, `min_chain_cosine=0.6`, `filter_large_baseline_hills=true`).

---

*Generated from the config structs and shipped TOMLs. If you add or rename a
config field, update this file and the shipped `benchmark/config/*.toml`; the
`deny_unknown_fields` parse test will fail the build if the TOMLs and structs
disagree, but it cannot check this doc — keep it current by hand.*
