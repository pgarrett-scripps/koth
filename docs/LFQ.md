# Alignment and label-free quantification

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

- `consensus_features.tsv` — one row per reference feature with coordinates, seed score, run support, and group confidence
- `intensity_matrix.tsv` — features × samples; each cell is the `[lfq].quant_estimator` value, apex intensity by default since 0.10.0 (`"sum"` gives integrated area); 0 = no extracted signal
- `qvalue_matrix.tsv` — TDC q-values for each (feature, sample) cell; only written when `run_tdc = true`; an exploratory confidence measure for filtering intensities
- `decoy_intensity_matrix.tsv`, `lfq_details.tsv` — decoy and per-cell diagnostics (see [outputs](OUTPUTS.md))
- `align_config.toml` — copy of config used
- `align_report.json` — alignment and LFQ diagnostics

All runs in a batch must report ion mobility on the same 1/K0 scale. Each run's
scale is read from its `report.json`; a timsTOF run with none recorded (koth
0.9.0 or earlier) counts as linear, and a mixed batch is refused with an error
naming each run. Re-run detection, or process every input with
`[file] bruker_mobility_scale = "linear"`.

## Search-guided LFQ and match-between-runs

Accepted peptide identifications can target the same extraction engine:

```bash
koth_align batch/ --sage-psms results.sage.tsv --output search_quant/
koth_align batch/ --targets targets.tsv --output search_quant/
koth_align batch/ --sage-psms results.sage.parquet --no-mbr --output direct_quant/
```

Sage TSV/Parquet and a generic identification TSV are supported.
`--max-id-qvalue` and `--max-extraction-qvalue` both default to `0.01`;
`--no-mbr` disables transfers between runs. Peptide outputs separate same-run
MS2 evidence from transferred extraction and report identification and
exploratory extraction confidence separately. Omitting these options keeps
identification-free LFQ. See [search-guided LFQ](search-guided-lfq.md) for
schemas, outputs, and validation limits.

## Confidence and interpretation

Group and per-cell q-values are separate exploratory confidence measures.
They have not been independently calibrated as peptide-identification FDR;
`q <= 0.01` must not be described as a demonstrated 1% error rate. Detected
cells also undergo rescoring and do not automatically receive a zero q-value.
See the [output reference](OUTPUTS.md) for how to interpret exported matrices.

## Alignment

`koth_align` selects a reference run using high-confidence feature support,
then matches anchors by charge, mass, normalized retention time, and mobility
where available. A RANSAC inlier selection followed by a median piecewise warp
and isotonic monotonicity corrects retention time. Mass and mobility drift are
fit across anchors; insufficient anchors fall back to an identity RT warp.
See [configuration](CONFIGURATION.md) for the current matching and warp settings.

## Quantification

Features are projected into reference coordinates and grouped across runs.
Group confidence uses RT-permutation controls; per-cell confidence uses the
configured target-decoy rescorer. Each run's hills are loaded for extraction,
then released, so all runs' hill profiles need not be resident together.

Extraction scores chromatographic signal and isotope-pattern agreement in
an RT/mass/mobility window. Exclusive signal ownership prevents competing
groups from simply reusing the same native peak samples. This does not resolve
all isotope-choice errors, and the current validation records an Orbitrap
precision regression. Review these limitations when selecting this pipeline.

The detailed development evidence is retained in [group confidence](lfq-group-confidence.md),
[permutation controls](lfq-group-permutation.md), [group integrity](lfq-group-integrity.md),
[signal ownership](lfq-signal-ownership.md), and [isotope choice](lfq-isotope-choice.md).
The optional [long-format bundle](long-matrix.md) exposes evidence for auditing.
