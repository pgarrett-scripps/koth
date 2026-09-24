# Output reference

## Single-run files

The default TSV output is written to `<output>/<input_stem>/`. Set
`[output].format = "parquet"` to use Parquet; column names are shared between
both writers. Each run also writes `report.json` and the resolved `config.toml`.
`report.json` records the 1/K0 scale (`mobility.scale`) and any Bruker
calibration coefficients applied.
Optional MS2 hill output adds isolation-window columns.

### Hills

One row represents one chromatographic trace, sorted by `intensity_sum`
descending. Times are in minutes and profiles include zero-filled scan gaps.

| Columns | Meaning |
|---|---|
| `hill_id` | Identifier used by feature-to-hill references |
| `mz`, `mz_std`, `mz_se` | Weighted m/z, spread, and standard error |
| `rt`, `rt_start`, `rt_end`, `rt_width` | Apex time, boundaries, and width |
| `im`, `im_std` | Ion mobility (1/K0) and spread; zero when unavailable. Bruker `.d` values use the acquisition-calibrated scale by default |
| `scan_start`, `scan_apex`, `scan_end` | Absolute scan indices |
| `n_scans`, `skipped_scans` | Profile length and zero-filled gaps |
| `intensity_sum`, `intensity_max` | Integrated and maximum intensity |
| `hill_score` | Chromatographic shape score |
| `intensity_profile` | Per-scan intensity array |

MS2 hills append `iso_target_mz`, `iso_lower_mz`, and `iso_upper_mz`.

### Features

One row represents a scored isotope envelope, sorted by `intensitySum`
descending. Unknown-charge features are excluded. Array-valued fields are JSON-encoded strings in both TSV and Parquet;
missing-value representation follows the writer.

| Columns | Meaning |
|---|---|
| `massCalib`, `mz` | Corrected monoisotopic neutral mass and m/z |
| `rtApex`, `rtStart`, `rtEnd` | Apex and retention-time boundaries, in minutes |
| `FAIMS` | Compensation voltage, when present |
| `intensityApex`, `intensitySum` | Apex and total intensity across isotopes |
| `charge`, `nIsotopes`, `nScans` | Charge, isotope count, and scan count |
| `im` | Apex ion mobility; empty in TSV when unavailable |
| `cosine_score`, `ppm_error` | Chromatographic co-elution score and isotope-spacing error |
| `neutron_offset` | Selected monoisotope reassignment offset |
| `isotope_score`, `combined_score` | Averagine score and combined isotope/co-elution score |
| `theoretical_pattern`, `isotope_profile` | Predicted pattern and observed isotope apex intensities |
| `elution_profile`, `hill_ids` | Summed elution profile and constituent hill identifiers |
| `intensityApexParab` | Apex intensity refined by parabolic interpolation |
| `intensityScattered5` | Sum of the five largest profile values |
| `intensityConsec5` | Largest sum over five consecutive profile values |

The column constants and writers in
[`koth_ff/src/output/mod.rs`](https://github.com/pgarrett-scripps/koth/blob/master/koth_ff/src/output/mod.rs) define the exact
schema and serialization behavior.

## Multi-run files

`koth_align` writes to `<batch_dir>/align_output/` unless `--output` is set.
Its matrices are TSV; `[output].format = "parquet"` is accepted but falls
back to TSV with a warning.

| File | Contents |
|---|---|
| `consensus_features.tsv` | Coordinates, seed score, run support, and group confidence |
| `intensity_matrix.tsv` | Feature metadata followed by unfiltered intensity per run |
| `qvalue_matrix.tsv` | Per-cell target-decoy confidence, when `run_tdc = true` |
| `decoy_intensity_matrix.tsv` | Decoy-pass intensities, when `run_tdc` and `[output].export_decoys` are set |
| `lfq_details.tsv` | Long-format per-cell diagnostics, when `[output].export_details = true` |
| `align_config.toml` | Resolved configuration |
| `align_report.json` | Alignment and LFQ diagnostics and timing |

Each intensity cell uses the `[lfq].quant_estimator`, apex intensity by default
since 0.10.0 (`"sum"` gives integrated area). Search-guided runs
(`--sage-psms`, `--targets`) write additional peptide-level outputs described in
[search-guided LFQ](search-guided-lfq.md).

Consensus metadata columns are `massCalib`, `mz`, `charge`, `rtApex`, `im`,
`combined_score`, `seed_run`, `n_contributing_runs`, `n_runs_detected`,
`group_score`, and `group_qvalue`. `n_runs_detected` counts positive-intensity
cells with q-value at or below `max_qvalue`.

The intensity and q-value matrices share the first eight metadata columns
(through `n_contributing_runs`), followed by run-name columns. Zero intensity
means no extracted signal; q-value 1 is the default when a cell has no scored
signal. Detected cells are rescored rather than automatically assigned zero.

**Group and cell q-values are exploratory measures, not validated
peptide-identification FDR.** Filtering at 0.01 does not establish a 1% error
rate. See [LFQ limitations](LFQ.md#confidence-and-interpretation). Optional
decoy/detail exports and the [long-format bundle](long-matrix.md) retain
additional diagnostic evidence. The exact matrix writers are in
[`align_main.rs`](https://github.com/pgarrett-scripps/koth/blob/master/koth_ff/src/align_main.rs).
