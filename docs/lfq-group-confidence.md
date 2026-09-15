# LFQ group-confidence experiment

**Historical report:** the prototype now uses independent RT permutations after
the original whole-run controls proved unstable. See the
[current estimator and validation](lfq-group-permutation.md); the report below
describes the original implementation and its first cohort.

The prototype uses weaker evidence with one quality gate, but the first real-data
trial shows a coverage/quality tradeoff rather than a clear overall improvement.

Branch: `codex/lfq-group-confidence`, based on `0908b50e5f3f`. This is an isolated
prototype; the main checkout and paper results are not updated. Experimental
artifacts are under `out/group-confidence/` in this worktree.

## Configuration and behavior

The prototype replaces `min_member_combined_score`, `min_seed_combined_score`,
`min_group_size`, and `allow_replicated_weak_seeds` with one setting:

```toml
[lfq.consensus]
mz_ppm = 20.0
rt_window_pct = 0.02
im_tolerance = 0.05
max_group_qvalue = 0.05
```

The removed keys cause a configuration error rather than being silently ignored.
Mass, RT and IM span limits remain independent of extraction windows. All valid
finite-scored features can enter candidate groups, including scores below 0.5.
Two distinct original runs are required for cross-run evidence, so singletons
cannot pass even at `max_group_qvalue = 1.0`. That value is useful for auditing
all replicated candidates. A single-run dataset has no cross-run consensus.

The best-scoring member still supplies reference coordinates and the isotope
template; its score is metadata, not an admission threshold. Group assembly
remains deterministic and bounded, so low-quality bridges cannot create groups
that exceed the full coordinate-span limits. Each run has one primary observation;
alternatives retain their original provenance.

## Evidence and controls

This first prototype uses a transparent fixed score, not a learned classifier.
For run `r`, let `s_r` be its best member's detector score clipped to [0, 1]. Its
weight is `w_r = s_r² / sum(s)` over the group's alternatives from that run.
This uses detector quality continuously and penalizes ambiguity without counting
duplicate observations as extra runs. The detector score is not interpreted as
a probability.

For a group containing `n` distinct runs, the evidence is
`sum(sqrt(w_r*w_s) * exp(-2*d_rs²)) / (n-1)` over distinct run pairs. Here `d_rs²`
is squared mass/RT/IM distance normalized by the full grouping tolerances; absent
IM contributes no term. The normalization makes evidence scale approximately
linearly with run support. Singleton and zero-evidence groups receive q = 1.

Five deterministic null cohorts independently translate each run's RT positions
on a circular reference coordinate interval. Mass, IM, detector quality and
within-run feature structure remain intact. Each cohort undergoes the same
bounded grouping, primary selection, ambiguity penalty and evidence calculation.
At a threshold, the estimated false-group count is one plus the mean null count,
divided by the target count. Equal-score blocks and a reverse cumulative minimum
produce `group_qvalue`. Zero-evidence candidates cannot dilute that denominator.
The +1 correction means fewer than 20 replicated hypotheses cannot pass the
default 0.05 gate, even if no null groups score as highly.

The group q-value only selects hypotheses for extraction. The existing LFQ grid
integration, QDA cell rescorer and cell gate remain separate. Strong cross-run
support cannot substitute for measured signal in a particular run. Output adds
`group_score` and `group_qvalue` to `consensus_features.tsv` and to long-bundle
schema version 2; the wide intensity and cell-q matrix schemas are unchanged.

## Validation limits

These group q-values are experimental estimates, not validated FDR guarantees.
Circular translations assume enough RT stationarity to represent accidental
grouping. Local density changes, shared alignment fitting, recurring artifacts,
correlated runs and incorrect attachments to otherwise genuine groups can defeat
that assumption. An accepted group is not a calibrated set of correct member
links, a peptide identification, or proof that every extraction is correct.
Global group confidence also does not guarantee error control specifically among
the newly rescued weak groups. Existing cell q-values retain their prior
calibration limitations.

The real-data experiment is exploratory on an existing benchmark, not a held-out
validation. Lower missingness can result from dropping sparse rows, so the
comparison script additionally evaluates CV on exactly the same positive gated
peptide/run cells. Any broader adoption needs independent entrapment/absence
controls and a second cohort, including a check of the weak-only subset.

## Reproduction

From this worktree, build and check the implementation:

```bash
cargo build --offline --release -j 2 --bin koth_align --example lfq_group_audit
cargo test --offline --workspace
cargo clippy --offline --workspace --all-targets -- -D warnings
python scripts/lfq_confidence_smoke.py target/release/koth_align
```

Run the experiment using the saved config (the original baseline config with
only grouping filters replaced, details export disabled, and explicit TSV output):

```bash
RAYON_NUM_THREADS=4 target/release/koth_align \
  /home/patrick-garrett/Repos/koth-paper/benchmark/data/bruker_15min/koth \
  --config out/group-confidence/config.toml \
  --output out/group-confidence/align
RUST_LOG=info RAYON_NUM_THREADS=4 target/release/examples/lfq_group_audit \
  /home/patrick-garrett/Repos/koth-paper/benchmark/data/bruker_15min/koth \
  out/group-confidence/config.toml out/group-confidence/audit
```

The audit table includes rejected candidates and supports a group-gate count
sweep without re-extraction. The `q020` trial repeats extraction/rescoring with
`max_group_qvalue = 0.2` as an exploratory sensitivity check, not a recommended
FDR operating point. Both trials retain the existing per-cell gate of 0.05.

Evaluate each trial with `analysis/scripts/bruker_koth_align.py` from the paper
checkout, passing its isolated `--align` directory, a new `--json` output path and
`--sage-psms /home/patrick-garrett/Data/dnoise/results/dda_15min/original/results.sage.tsv`.
Use the paper's `analysis/.venv/bin/python`, which contains its existing analysis
dependencies. The final comparison command is:

```bash
/home/patrick-garrett/Repos/koth-paper/analysis/.venv/bin/python \
  scripts/lfq_confidence_compare.py \
  --paper-scripts /home/patrick-garrett/Repos/koth-paper/analysis/scripts \
  --baseline /home/patrick-garrett/Repos/koth-paper/analysis/results/bruker_15min/lfq_im_sensitivity_20260913/im015/align \
  --baseline-metrics /home/patrick-garrett/Repos/koth-paper/analysis/results/bruker_15min/lfq_im_sensitivity_20260913/im015/metrics.json \
  --sage-psms /home/patrick-garrett/Data/dnoise/results/dda_15min/original/results.sage.tsv \
  --experiment out/group-confidence
```

`comparison.json`, `paired_cv.parquet`, per-trial `metrics.json`, logs, configs
and `provenance.json` hold the measurements. Results are summarized below.


## Results on the 18-run timsTOF HYE cohort

All three trials use identical fitted alignments and a per-cell gate of 0.05.
The baseline is the cached 0.5 member / 0.75 seed / two-run configuration used
for the existing IM=0.015 trial. Neither sensitivity setting was selected on
held-out data; the prototype default remains 0.05.

| Metric | Baseline | Group gate 0.05 | Group gate 0.20 |
|---|---:|---:|---:|
| Consensus groups | 107,215 | 80,666 | 92,927 |
| PSM recall (%) | 66.283 | 65.995 | 66.799 |
| Median CV (%) | 10.332 | 10.329 | 10.490 |
| Human ratio IQR (log2) | 0.3972 | 0.3911 | 0.4008 |
| Human false fold-change rate (%) | 15.465 | 15.435 | 15.915 |
| Missingness within selected anchors (%) | 29.294 | 25.992 | 27.846 |
| Mean absolute nonhuman species bias (log2) | 0.2891 | 0.2714 | 0.2820 |

At the default group gate, 27,170 retained groups include 46,989 observations
below the old 0.5 member cutoff (including alternative observations). There are
4,176 retained groups whose best member is below 0.75, including six whose best
member is below 0.5. Of the MS2 peptide/charge anchors, 721 map to groups with a
best score below 0.75; none map to the six all-below-0.5 groups. This demonstrates
access to weaker evidence, but does not establish that all rescued observations
are correct or yield new peptide identifications.

Relative to baseline, the 0.05 trial gains 692 anchors with any gated signal and
loses 2,585. At the peptide/run-cell level it gains 13,197 cells and loses 18,415.
The 0.20 sensitivity trial gains 1,071 anchors and loses 1,365; it gains 14,351
cells and loses 9,318. These counts describe assignment and extraction changes,
not independently validated discoveries. The looser gate increases PSM recall
by 0.516 percentage points but also increases aggregate CV and false fold-change
rate. Its nominal group q-value must not be presented as calibrated 20% FDR.

The identical-cell control contains 504,725 peptide/run cells and 97,461 paired
peptide/condition CV observations. Median CV is 10.1650%, 10.1636% and 10.1637%
for baseline, gate 0.05 and gate 0.20. The median paired changes are effectively
zero. Thus the apparent differences in aggregate CV/missingness are chiefly
selection effects; they do not demonstrate improved intensity estimation.

The first full 0.05 trial took 119.92 seconds with four Rayon threads and peak
RSS of 1,708,364 KiB (about 1.63 GiB). These are local warm/cold mixed-cache
observations, not a controlled performance comparison with the baseline.

## Verification and assessment

The 168 library tests and one doc test pass. Four pre-existing large-fixture
integration tests remain ignored under their normal configuration. Clippy passes
with warnings denied, and Rust formatting and diff whitespace checks pass.
Tests cover bounded grouping, missing IM, duplicate support, ambiguity, input/run
reordering, tied q-values, zero-evidence denominator protection, removed config
keys and invalid gate values.

In a controlled synthetic mixture across 20 deterministic trials, the prototype
recovers all 2,000 planted weak groups (each supported by six runs at score 0.35)
and accepts 96 false groups, a pooled false discovery proportion of 4.58% at the
0.05 gate. An independent uniform-noise fixture retains no groups. These toy
results test the implementation under the modeled null; they do not validate
the real-data null assumptions.

The CLI/export smoke test retains 200 replicated weak groups but, with empty
hill inputs, produces zero signal in all 800 target/decoy cells. It verifies
schema version 2, output hashes and separation of group confidence from cell
confidence. A final-build repeat of the 0.05 cohort is checked against the
measured trial by hashes of all four consensus/intensity/q-value/decoy tables;
`final-verification.json` records the result.

This branch is suitable for further experimentation, not replacement of the
paper's defaults. The mechanics achieve the requested simplification and weak
feature inclusion. The current fixed score has not established a better overall
operating point, and calibrated member-link/weak-subset confidence remains open.
