# Search-guided LFQ and match between runs

Version 0.8.0 adds optional peptide targets to `koth_align`. A search engine
identifies a modified peptide and charge in a donor run; Koth integrates its
MS1 isotope signal in that run and, with MBR enabled, in aligned recipient runs.
A complete feature detection for that peptide is not required. Extraction uses
the existing hill grids, isotope scoring, exclusive signal ownership, and
intensity estimator. Alignment still uses Koth feature anchors.

## Run it

First run `koth_ff` on every sample to create a batch directory containing
`hills.{tsv,parquet}` and `features.{tsv,parquet}` per run. Then:

```bash
# Sage TSV or Parquet; extract direct IDs and attempt cross-run transfers.
koth_align batch/ --sage-psms results.sage.tsv --output quant/
koth_align batch/ --sage-psms results.sage.parquet --output quant/

# Direct identification-guided extraction only.
koth_align batch/ --sage-psms results.sage.tsv --no-mbr --output direct/

# Identifications from another search engine, converted to the schema below.
koth_align batch/ --targets targets.tsv --output quant/
```

`--targets` and `--sage-psms` are mutually exclusive. Omit both for the existing
identification-free workflow. Search mode quantifies only imported targets;
it does not append them to the identification-free consensus list.

| Option | Default | Meaning |
| --- | --- | --- |
| `--max-id-qvalue` | `0.01` | Maximum imported ID q-value; Sage must pass both spectrum and peptide q-values. |
| `--max-extraction-qvalue` | `0.01` | Exploratory cell-score cutoff for the peptide intensity matrix. This is not a validated 1% transfer FDR. |
| `--no-mbr` | off | Extract only in runs with an accepted same-run identification. |
| `--ignore-target-im` | off | Ignore imported ion mobility if it is not expressed in native Koth coordinates. Recorded in the search manifest. |

All thresholds must be finite and between zero and one. Existing TOML LFQ
extraction settings apply, including normalization. Search options are CLI
arguments, not additional TOML keys. The peptide matrix gate is independent of
`[output].max_qvalue`, which only affects legacy consensus summary counts.

## Generic identification TSV

One row is an accepted MS2 identification in a named run. The caller must remove
search decoys and non-primary/ambiguous PSMs and apply its search engine's
identification validation before export. Supply an ID q-value that satisfies
both the intended PSM and peptide confidence requirements. A predicted peptide
library with no observed donor identification is not this input format.

| Column | Required | Meaning |
| --- | --- | --- |
| `modified_peptide` | yes | Stable peptide string including all modifications and their positions. Exact strings define identity; Koth does not translate modification notation. |
| `run` | yes | Run directory name or raw filename/path. |
| `charge` | yes | Positive precursor charge, integer 1–255. |
| `neutral_mass` | yes | Calculated neutral monoisotopic peptide mass in Da, including modifications; not m/z, MH+, or an isotope-shifted measured precursor mass. |
| `rt_minutes` | yes | Observed native retention time in minutes; not aligned/normalized/predicted RT. |
| `id_qvalue` | yes | Validated input identification q-value in [0, 1]. |
| `im` | no | Native inverse reduced mobility, 1/K0; empty or 0 means unavailable. |

```text
modified_peptide	run	charge	neutral_mass	rt_minutes	id_qvalue
PEPTIDE	A.mzML	2	799.359964	12.3	0.001
```

The example uses literal tab-separated fields. Mass is supplied by the caller;
Koth does not calculate or verify it from the sequence. Sequence-derived exact
isotope distributions and protein aggregation are outside this release.

Paths can use Unix or Windows separators. Known `.gz`, `.mzML`, `.raw`, and `.d`
suffixes are removed for matching; remaining sample names are case-sensitive.
This also matches Koth's `sample.mzML` output names for compressed mzML input.
Unknown retained donor runs or colliding normalized run names fail explicitly;
filter a larger search cohort to the batch being quantified first.

## Sage import

The adapter reads TSV and Parquet. Required columns are `peptide`, `filename`,
`charge`, `calcmass`, `rt`, `rank`, `spectrum_q`, `peptide_q`, and either
`is_decoy` (boolean or 0/1) or `label` (1 target, -1 decoy). Only rank-1 targets
passing both q-value gates are retained. The reported `id_qvalue` is the maximum
of `spectrum_q` and `peptide_q`, used as a conservative admission summary; it
is not a newly estimated identification q-value. `expmass`, `aligned_rt`, and
predicted coordinates are not used.

Optional `ion_mobility` must already be native 1/K0. Since 0.10.0, koth's
native timsTOF 1/K0 is Bruker's acquisition-calibrated scale (the scale the
timsdata SDK and SDK-based mzML exports report); features written with
`bruker_mobility_scale = "linear"` or by koth 0.9.0 and earlier use a linear
scale instead. Some converted timsTOF
search inputs have different mobility coordinates; correct them before import
or use `--ignore-target-im`. Koth does not infer a coordinate conversion.

## Target selection and transfer eligibility

Each exact modified-peptide/charge pair produces at most one target. Repeated
PSMs do not create duplicate matrix rows. Koth retains the lowest-q accepted
PSM per run, breaking ties deterministically by native coordinates and source
row. The seed is the best accepted ID in a run connected to the reference by a
usable alignment, or the best donor if none is connected.

PSMs with inconsistent calculated masses or incompatible aligned RT/IM are
rejected as a target and listed in `search_rejections.tsv`. Compatibility uses
`[lfq.consensus]` mass, RT-span, and IM-span limits; repeated native donor
observations are checked as well. This conservative release does not split one
peptide/charge identity into chromatographic isomers or mobility conformers.
The identification-free group support requirement and group-q gate do not
apply: a single donor ID is sufficient to define a search target.

Same-run extraction uses the selected PSM's native RT and mobility. Recipient
RT/mass/IM come from the seed projected into reference coordinates and the
existing alignment model. Missing target mobility stays missing. Mass centres
use the calculated monoisotopic mass with positive protonation and the existing
relative run mass correction; absolute reference-run mass bias is not fitted
from IDs. Averagine templates remain the configured LFQ isotope model.

Transfers require an aligned seed donor and recipient, with the seed RT inside
the measured alignment range. Out-of-range donor IDs remain direct-only to
avoid transferring a clipped boundary coordinate. Non-reference runs need
at least `[alignment].min_anchor_count` retained RT inliers (at least two even
if configured lower). Runs falling back to identity alignment receive direct
extraction only. Featureless runs can still quantify direct IDs using their
hill RT bounds. No signal remains zero; MBR does not impute an intensity.

## Outputs and confidence

| File | Contents |
| --- | --- |
| `peptide_quant.tsv` | Every target/run opportunity: identity, `direct_ms2` or `mbr` evidence, raw/accepted intensities, status, extraction q-value, same-run ID q-value, seed run/q-value, and original input row references. |
| `peptide_intensity_matrix.tsv` | Target × run matrix after the explicit extraction gate. Includes target ID, modified peptide, and charge. |
| `search_rejections.tsv` | Peptide/charge pairs excluded for inconsistent mass or ambiguous RT/IM. |
| `search_manifest.json` | Input hash, import/filter counts, options, run order, retained donor records, rejected targets, and peptide-output hashes. |

`target_id` is a zero-based consensus row index and links to existing extraction
outputs. Input source rows are one-based data-row numbers, excluding the TSV
header or spanning all Parquet batches. A seed is a target-coordinate source;
it is not proof that its peptide exists in every recipient.

Statuses are `accepted`, `rejected`, `no_signal`, `not_attempted`, and `unscored`.
`direct_ms2` records accepted identification evidence even when no quantitative
signal was extracted. `mbr` denotes an opportunity without a same-run accepted
ID; only `accepted` with positive intensity is a retained transfer. Disabled or
unalignable opportunities are `not_attempted`, never successful transfers.

Koth ranks direct and transferred extraction cells separately, each with its
own paired synthetic decoys and feature-group cross-validation. Search decoys
are excluded from donor targets; synthetic extraction decoys are generated
by the existing coordinate-shift method. ID q-values are never features in
the extraction rescorer and never substitute for recipient confidence.
These extraction q-values remain exploratory, including in small target panels;
independent absent-species/entrapment validation is needed before interpreting
a numerical cutoff as calibrated peptide-transfer FDR.

Legacy `intensity_matrix.tsv` and `lfq_details.tsv` retain raw quantities, even
when the peptide matrix rejects a cell. Existing normalization applies to the
wide intensity matrices and peptide output; detail entries remain pre-normalization.
With `run_tdc = false`, peptide cells with signal are `unscored`, extraction
q-values are blank, and the accepted peptide matrix is zero; raw quantities
remain available. No confidence is invented when rescoring is disabled.

With `export_long = true`, search mode emits long-bundle schema 4: original
feature membership and seed-observation links are empty for peptide targets,
and `lfq_is_mbr` means no accepted same-run ID. The peptide target IDs link to
`peptide_quant.tsv` and `search_manifest.json`. Detector/group scores are not
applicable and are blank in the long bundle (`NaN` in legacy wide metadata).
Identification-free mode keeps schema 3 and its existing feature-detection
meaning of MBR. Neither mode silently promotes an MS1 extraction to an MS2 ID.

## Validation scope

Automated tests cover Sage TSV/Parquet filtering, malformed inputs, run mapping,
modification/charge identity, duplicate PSMs, ambiguous coordinates, extraction
without complete features, aligned transfers, missing signal, MBR disabling,
weak alignment, separate confidence provenance, and long-format export.
These establish software behavior; quantitative accuracy and transfer-error
calibration require the planned independent LFQ-paper experiments.

Local validation on September 16, 2026 passed formatting, Clippy with warnings
as errors, and rustdoc with warnings as errors. The default test suite passed
231 tests (six pre-existing tests remain ignored); the no-default-features
suite passed 225 (one pre-existing test remains ignored). Nineteen new tests
exercise the importer, extraction, provenance, and CLI. The identification-free
CLI fixture reproduced all seven compared TSV outputs byte-for-byte against
the 0.7.0 release binary at `88f9e7c54c9a`. Native vendor-fixture tests and
biological LFQ/transfer-error benchmarks were not rerun for this change.
