# LFQ group confidence with independent RT permutations

Final adoption decision: [signal ownership](lfq-signal-ownership.md) records the
retained improvements and known isotope-choice and precision limitations. The
recommendations below describe the earlier experimental stages.

**Subsequent integrity audit:** [split/isotope/shared-signal checks](lfq-group-integrity.md)
found unresolved duplicate-reporting cases. The estimator remains the preferred
prototype, but group reconciliation and signal ownership need resolution before
default adoption. The results below concern the estimator comparison.

**Decision: replace whole-run RT shifts with independent feature-level RT
permutations in the experimental worktree.** The replacement resolves the large
control-set swings observed previously and passes the preset quantitative
comparison limits. It keeps the original evidence score and ambiguity penalty,
removes the same four legacy settings, and introduces no new production settings.
The main checkout and paper results remain unchanged.

## What changed

For each original run, features are partitioned by charge and detector-score bin
(width 0.1). Their existing retention times are independently shuffled within
each partition using deterministic Fisher–Yates permutations. The partition's
exact RT distribution is preserved, including crowded regions and ties. Each
feature retains its original mass, IM, quality, charge and provenance; only its
RT assignment changes in the artificial control cohort.

This differs from both whole-run translations and the earlier circular rank
shifts: a single displacement no longer moves thousands of features coherently.
Fixed points are allowed rather than forcing control features away from possible
matches, which could underestimate accidental agreement. Sparse or constant-RT
partitions can therefore have little or no discriminating power.

Ten control cohorts undergo the same bounded grouping and evidence scoring as
the real features. At each cutoff, the estimated false-group fraction remains
`(1 + mean control count) / real count`, with tied-score blocks and a reverse
cumulative minimum. The pooled ten-control estimate is now used by the normal
`koth_align` path. The log also reports the accepted counts from each disjoint
five-control half, making sensitivity visible without another setting.

## Stability and weak-feature recovery

All methods below use the same score and group gate of 0.05. The two control
halves differ only in their deterministic shuffle seeds; set agreement is
intersection divided by union, not an accuracy measurement.

| Cohort | First five controls | Second five controls | Set agreement | Pooled ten |
|---|---:|---:|---:|---:|
| Bruker, 2 runs | 31,606 | 31,682 | 99.760% | 31,642 |
| Bruker, 6 runs | 55,942 | 55,955 | 99.977% | 55,947 |
| Bruker, 18 runs | 77,239 | 77,293 | 99.930% | 77,262 |
| Orbitrap, 20 runs | 78,962 | 79,054 | 99.884% | 79,004 |

The previous whole-run controls accepted 1,647 versus 43,306 groups in the same
two-run comparison. Independent permutation passes the predeclared 95% set
agreement requirement on every tested cohort. Its two-run accepted set includes
2,230 groups with a best feature below the old 0.75 seed floor. The 18-run set
contains 3,432 such groups and 44,088 original observations below 0.5; the latter
includes alternatives retained only for provenance.

## Quantification compared with the preceding prototype

The scorer and both 0.05 gates stayed fixed. Full extraction and cell rescoring
were repeated for the newly accepted groups. The preset limits were no recall
loss greater than 0.5 percentage points, no CV increase greater than 0.2 points,
and no human false-fold-change-rate increase greater than 0.5 points.

| Cohort / estimator | Groups | PSM recall | Median CV | Human false-fold-change rate |
|---|---:|---:|---:|---:|
| Bruker / whole-run shifts | 80,666 | 65.995% | 10.329% | 15.435% |
| Bruker / independent permutations | 77,262 | 65.658% | 10.274% | 15.273% |
| Orbitrap / whole-run shifts | 84,684 | 49.777% | 16.136% | 9.283% |
| Orbitrap / independent permutations | 79,004 | 49.385% | 16.039% | 8.989% |

Both cohorts pass those limits. Recall declines by 0.34 and 0.39 percentage
points, respectively, with slightly better aggregate CV and false-fold-change
rate. On shared cells, Bruker CV is essentially unchanged (10.2731% versus
10.2743%, 521,325 cells); Orbitrap is identical at the reported precision
(16.0101%, 470,448 cells). The aggregate improvements do not establish better
measurement of the same cells. Relative to the original legacy score-floor
workflow, the broader coverage/precision tradeoff reported earlier still exists.

The cached independent MS2 labels were reused only after verifying identical
candidate coordinates, indices and run ordering. Among assessed primary weak
observations with positive extraction and cell q <= 0.05, Bruker has 18 discordant
links among 1,766 (1.02%; Wilson 95% interval 0.65–1.61%), compared with 18/1,779
previously. Orbitrap has 14/586 (2.39%; interval 1.43–3.97%), compared with 15/587.
There is no supported increase in this observed discordance proxy. It excludes
unlabelled features and groups without a strong labelled donor, so it does not
measure all rescued-group or transfer FDR.

## Synthetic checks and remaining limits

On the existing mixed fixture with independent RT background correlated with
feature quality, permutations recover all 500 planted groups across five trials
and accept 17 false groups: 3.29% false discoveries, versus 11.5% for whole-run
shifts. The pure independent-noise regression admits no groups. These limited
fixtures support the replacement but do not prove general FDR control.

The concentrated-background fixture still rejects every group, including genuine
signals. Recurrent MS1 artifacts are still accepted, and a strong genuine group
can still contain an erroneous member: the absence fixture retains 62 spurious
absent-run links among 1,562 original links in 500 genuine groups. Confidence in
group existence cannot establish peptide identity or every member's correctness.
The separate cell-extraction gate remains essential.

Permutations assume RT assignments are exchangeable within each partition. They
do not preserve within-partition mass–RT or IM–RT dependence, or within-run
alternative-feature structure. A biological entrapment/absence validation is
still unavailable locally, and the already inspected Orbitrap cohort is a
secondary check rather than a new held-out dataset. The configured 0.05 must
continue to be described as an experimental group statistic, not a validated
5% identification or member-link FDR. Published MBR work such as
[IonQuant](https://pmc.ncbi.nlm.nih.gov/articles/PMC8131922/) models identification
transfer using different decoys and calibration; that work does not validate
this group-level estimator.

## Implementation and reproduction

`group_confidence.rs` contains the production control generator and diagnostics.
`lfq_permutation_audit` shares production projection, grouping, scoring and
extraction. The older `lfq_refinement` example remains a historical comparison
of the original shift/rank methods. Both research examples share the same TSV
writer. No changes were made to the cell rescorer or extraction algorithm.

The no-signal CLI fixture now spreads planted weak groups over RT so permutation
controls have information to discriminate them; it still supplies no measured
hills and checks that group acceptance cannot invent intensity. A separate
constant-RT regression requires q = 1 when permutations cannot discriminate.
Tests also verify exact within-partition RT preservation, unchanged original
feature attributes, reproducibility under reordered observations/runs, and
rejection of pure independent noise.

The final checks passed: 174 library tests, one documentation test, Clippy with
warnings denied, formatting, and the no-signal CLI/export smoke check. The
synthetic audit was run explicitly; four existing external-fixture tests remain
ignored. All four standard TSV outputs are byte-identical between the two-run
research extraction and the final production CLI.

Artifacts, protocol, control counts, selected candidates, metric tables and
source/input hashes are under `out/group-permutation/`. Build with:

```bash
cargo build --offline --release -j 2 --bin koth_align --example lfq_permutation_audit
cargo test --offline --workspace
cargo clippy --offline --workspace --all-targets -- -D warnings
python scripts/lfq_confidence_smoke.py target/release/koth_align
```

Run `lfq_permutation_audit BATCH CONFIG OUTPUT`, optionally with `--take 2` or
`--take 6`; `--quantify` extracts the pooled accepted set. Candidate q columns
are `average_first`, `average_repeat` and `average_permuted`. Evaluate quantification
with `scripts/lfq_refinement_metrics.py --control permuted --variants average`
and the paper's Python environment and analysis paths, as in the previous report.
Its `legacy` result key denotes the provided `--baseline` path: in this experiment
that is the preceding prototype, not the original score-floor workflow.
