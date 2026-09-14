# Original additive isotope evidence pilot

This records the initial implementation before best-prefix selection and parameter
tuning. See [the tuning report](tuning/README.md) for the final implementation and
held-out validation.

The first fixed model improved matched peptide-charge recall from **83.9527% to
84.2803% (+0.3276 percentage points, +333 IDs)** on four IonStar files. Every
file improved. Quantification error on common complete peptide-charge pairs was
essentially unchanged; E. coli within-level CV increased slightly. This is a
promising small experiment, not a calibrated noise model or a full-cohort result.

## Scope and implementation

- Branch: `codex/additive-isotope-evidence`, based on `0908b50` (current master at
  worktree creation), not the older manuscript-pinned build.
- `rescore_chain` now returns summed isotope log evidence for claim priority.
  `HeapItem::cmp` uses that score alone, followed only by deterministic mono-index
  and charge ties. There is no envelope-length bonus or isotope-score override.
- Existing candidate-generation gates, Bhattacharyya claim gate and downstream
  scoring/filtering are retained. The exported `combined_score` is still the
  downstream score, not this internal claim-priority score.
- Sulfur offsets still select the best **whole-chain** template. The legacy
  `exhaustive_isotope_priority` option is accepted for compatibility but ignored.
- Only four raw files were processed once per binary. Paper data, manuscript
  results, alignment, LFQ and Sage searches were not rerun or modified.

For each isotope after the mono seed, let `r` be the observed log2 apex-intensity
ratio to the seed minus its theoretical log2 ratio, and let `c` be its
seed-anchored chromatographic cosine. The working models are:

| Evidence | Isotope model | Noise null |
| --- | --- | --- |
| Ratio residual `r` | Normal(0, 0.5²) | Normal(0, 2²) |
| Cosine `c` | Beta(4, 1) | Uniform(0, 1) |

The natural-log contribution is `ln(4) - 1.875*r² + ln(4) + 3*ln(c)`.
The chain score is the sum over M+1 onward, maximized over the allowed sulfur
templates. The seed contributes zero. Exact-fit, perfectly co-eluting isotopes
add about 2.773 each; a 1.5-log2 ratio error contributes about -1.446 even with
perfect co-elution. Thus the tested six-hill envelope with two such contaminants
scores about 5.425 and loses to a clean three-hill envelope scoring about 5.545.

Retained per-isotope observations do not depend on chain length or its combined
apex. Truncating positive-evidence tails lowers the score; removing negative
evidence can raise it. A prefix is always rescored and requeued at its own score.
This change only reranks the candidate pool: it does not generate every possible
prefix proactively or reject all negative-score uncontested chains.

The distribution shapes and widths were fixed before reading the A/B results;
there was no parameter sweep. They are explicit starting assumptions, not learned
from decoys. Ratio and shape observations share the same seed and are only
approximately conditionally independent. Candidate-generation gates also select
the observations before ranking. The score therefore is not a calibrated feature
probability, Bayes factor for the full detection process, or FDR estimate.

## Data and matching

Inputs are the existing uncompressed mzML files from PXD003881 under
`/home/patrick-garrett/Repos/koth-paper/analysis/data/ionstar_plain`. They were
chosen as A replicates 1–2 (B03_10, B03_11) and E replicates 1–2 (B03_05, B03_06)
before evaluation. Both binaries used the same copied paper Orbitrap config,
release builds and eight worker threads. The saved hill Parquet files are
byte-identical between baseline and additive for all four inputs.

Sage targets at peptide-q ≤ 0.01 are deduplicated per run/peptide/charge, selecting
the lowest peptide-q PSM with stable tie handling. A match requires identical
charge, 10 ppm m/z agreement and PSM retention time inside the feature interval;
features are restricted to charges 2–6. An apex ±0.5-minute sensitivity check is
also recorded. Missing MS/MS identification is not treated as a false feature.

| File | Baseline matches | Additive matches | Recall change (percentage points) |
| --- | ---: | ---: | ---: |
| B03_10, A1 | 21,181 | 21,237 | +0.223 |
| B03_11, A2 | 21,170 | 21,209 | +0.153 |
| B03_05, E1 | 21,405 | 21,534 | +0.506 |
| B03_06, E2 | 21,571 | 21,680 | +0.425 |
| Total / pooled | 85,327 | 85,660 | +0.328 |

| Metric | Baseline | Additive |
| --- | ---: | ---: |
| Matched IDs / 101,637 anchors | 85,327 | 85,660 |
| Pooled interval recall | 83.9527% | 84.2803% |
| Pooled apex-window recall | 81.2873% | 81.5077% |
| Features, charge 2–6 | 538,024 | 541,568 (+0.66%) |
| Features per matched anchor | 1.001535 | 1.001599 |
| Human median absolute log2 ratio error | 0.321535 | 0.321495 |
| E. coli median absolute log2 ratio error | 0.173683 | 0.173683 |
| Human median within-level CV | 8.5237% | 8.5223% |
| E. coli median within-level CV | 7.7095% | 7.7946% |

The net +333 comprises **1,905 gained and 1,572 lost** anchor matches. The changed
assignments are therefore more numerous than the net gain. Median isotope count
is not used to assess accuracy; mean chain length fell from about 2.71–2.72 to
2.65–2.66 hills per feature.

Quantification uses the exact same 12,640 peptide-charge pairs matched in all four
files by both scorers, including 11,745 human and 891 E. coli pairs (four other
pairs are excluded from species summaries). Intensities are the strongest
interval-matching feature's `intensitySum`; A/E ratios use the two within-level
medians and expected ratios 1 and 1/3. These are raw, PSM-anchored feature
intensities, without alignment, normalization or missing-value recovery. The
human log2 bias is about +0.309 in both versions. Two replicates per level and
selection for complete observations limit the CV comparison.

Observed total wall time was 94.5 s baseline and 93.9 s additive. These sequential
runs are not a controlled performance benchmark; the difference should not be
interpreted as a speed improvement.

## Validation and reproducibility

Passed: workspace tests with default and no-default features, Clippy on all
workspace targets with warnings denied, rustdoc with warnings denied, optional
Thermo feature compilation, format and diff whitespace checks. Seven new
regressions cover evidence signs, additive clean extensions, the contaminated-six
comparison, prefix heap ordering, deterministic score-only ties, invalid values,
and whole-chain sulfur selection. Existing assembly/determinism tests also pass. The default suite retains four
pre-existing ignored data-dependent tests; these were not enabled.

`compare.py` reproduces only this four-file experiment. It refuses to overwrite
an existing per-file output, and `score` requires all four outputs for both
variants. The exact commands and binary/config SHA-256 values are saved in
`artifacts/{baseline,additive}/manifest.json`; detailed matched anchors, logs,
binaries and raw feature outputs remain in the ignored `artifacts/` directory.
The small reviewable results are `per_file.csv` and `summary.json` beside this
README.

From the worktree root, after building and preserving the corresponding binaries:

```bash
/home/patrick-garrett/Repos/koth-paper/analysis/.venv/bin/python experiments/additive-isotope-evidence/compare.py run baseline experiments/additive-isotope-evidence/artifacts/bin/baseline
/home/patrick-garrett/Repos/koth-paper/analysis/.venv/bin/python experiments/additive-isotope-evidence/compare.py run additive experiments/additive-isotope-evidence/artifacts/bin/additive
/home/patrick-garrett/Repos/koth-paper/analysis/.venv/bin/python experiments/additive-isotope-evidence/compare.py score
```

The saved baseline was built from unmodified `0908b50` with
`cargo build --release -p koth_ff --bin koth_ff --locked` and copied before source
edits. The additive binary was built with the same command after the scoring
change. Both builds reused `CARGO_TARGET_DIR=/home/patrick-garrett/Repos/koth_rust/target`.
