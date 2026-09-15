# LFQ group confidence: bounded refinement and landing decision

**Historical decision:** the subsequent [independent-permutation experiment](lfq-group-permutation.md)
resolved the demonstrated control instability and replaced the original null in
the experimental worktree. Results below describe the preceding refinement pass.

**Decision: retain the original evidence formula, but do not promote this
prototype to the default LFQ workflow yet.** The proposed simplification still
removes the four legacy switches, including `allow_replicated_weak_seeds`, and
lets detector quality contribute continuously. However, this validation found
that the current shifted-run group threshold is not reliable across background
patterns and small run counts. Further formula tuning is not the priority;
null-model calibration is the remaining scientific blocker.

All changes remain in `codex/lfq-group-confidence`. The main checkout and the
paper's canonical results were not updated. The original report remains in
[lfq-group-confidence.md](lfq-group-confidence.md); this report supersedes its
landing assessment. Reproduction artifacts are in `out/group-refinement/`.

## Fixed comparison and implementation

The protocol was saved before the new experiments. Four scorers were compared
on the existing 18-run Bruker HYE cohort: the original pair average, that average
without the within-run ambiguity penalty, a maximum spanning-tree sum, and that
tree without the penalty. The tree can add a weak run without diluting an
existing strong core. Both group and cell gates stayed at 0.05; no threshold
sweep or learned classifier was introduced.

The original scorer was retained before evaluating the independent 20-run
PXD003881 Orbitrap cohort. Research code now shares exactly the production
projection and extraction pipeline, exports candidate/member provenance, and
supports the fixed scorer/control comparisons. Production does not calculate
the alternative tree scores or allocate their pair matrices. No new production
configuration switches were added during refinement.

## Development results

Percentages below use the existing paper analysis functions. Human FFCR is the
fraction of assessed human fold changes exceeding 1.5-fold despite an expected
ratio of one; it is not an identification FDR.

| Bruker variant | Groups | PSM recall | Median CV | Human FFCR |
|---|---:|---:|---:|---:|
| Legacy score floors | 107,215 | 66.283% | 10.332% | 15.465% |
| Original pair average | 80,666 | 65.995% | 10.329% | 15.435% |
| Average, no ambiguity penalty | 79,679 | 65.977% | 10.322% | 15.378% |
| Tree | 75,074 | 65.300% | 10.202% | 15.060% |
| Tree, no ambiguity penalty | 73,965 | 65.181% | 10.179% | 14.989% |

None of the challengers improved coverage, which was a prespecified requirement.
The apparent precision gain from the tree accompanies fewer retained rows.
Across 503,856 cells shared by all five methods, median CV is effectively
unchanged: 10.0172% for legacy, 10.0168% for the original scorer and 10.0156%
for the tree. There is no demonstrated benefit sufficient to justify changing
the evidence formula or removing the ambiguity penalty.

## Confidence controls: the blocker

The 18-run original scorer retains 80,666 groups with the first five global RT
shifts and 87,812 with a disjoint five shifts at the same gate. Pooling all ten
retains 84,011. More seriously, the two-run subset retains **1,647 versus 43,306**
groups with those disjoint control sets; pooling ten retains only 3,378. This is
control-model sensitivity, not nondeterminism: repeating the same controls
reproduces the same result.

The prespecified synthetic audit deliberately tested departures from the global
shift assumptions. With feature quality correlated with RT in otherwise uniform
background, the original scorer retained 500 planted groups and 65 false groups
across five trials: 11.5% false discoveries at a nominal 5% gate. With independent
background features concentrated in a narrow RT interval, it retained 500 planted
groups and 2,435 false groups. These are constructed counterexamples, not an
estimate of either real cohort's FDR, but they invalidate a general 5% guarantee.

RT-rank controls preserving each run's exact RT distribution did not solve the
problem. Following that failure, one documented diagnostic amendment added rank
shifts within fixed detector-score bins of width 0.1, pooled over ten controls.
That reduced the uniform quality-correlated case to 20 false groups alongside
500 planted groups (3.85%). Both rank-control approaches rejected every planted
group in the dense-background case, and both retained recurrent artifacts that
look like genuine repeated MS1 signals. MS1 recurrence does not establish
peptide identity.

On Bruker, the quality-conditioned diagnostic retained 74,156 groups, with recall
65.197%, CV 10.227%, and human FFCR 15.068%. Its shared-cell CV was again unchanged
from the original scorer (10.2250% versus 10.2248%). It rejected every group in
the two-run subset. This repair was therefore not adopted. Increasing the number
of global shifts alone does not address the demonstrated mismatch in the null.

The six-run subset retained 58,840 groups, including 3,631 weak-seed groups,
and 268,760 positive cells passing the existing cell gate. Both subset runs
completed with finite intensities and q-values. These are operating checks;
they do not establish acceptable small-cohort calibration.

## Independent cohort and weak-member audit

| Orbitrap variant | Groups | PSM recall | Median CV | Human FFCR |
|---|---:|---:|---:|---:|
| Legacy score floors | 47,491 | 42.251% | 13.926% | 5.284% |
| Original pair average | 84,684 | 49.777% | 16.136% | 9.283% |

The new approach gains 7.53 percentage points of recall on Orbitrap, but its
expanded output has worse aggregate CV and false-fold-change rate. On the
367,952 shared cells, median CV is nearly identical (13.7525% legacy versus
13.7526% prototype). This supports a coverage/quality tradeoff rather than a
claim that quantification improved. No scorer or threshold was retuned using
these confirmation results.

The member audit assigned unambiguous independent MS2 peptide/charge labels to
original features, treating I/L as indistinguishable. A tested observation was
compared only with unanimous strong-primary donor labels from other runs.
After selecting primary weak observations and requiring a positive extraction
passing cell q <= 0.05, the original scorer had:

| Cohort | Assessed weak links | Discordant | Observed fraction | Wilson 95% interval |
|---|---:|---:|---:|---:|
| Bruker | 1,779 | 18 | 1.01% | 0.64–1.59% |
| Orbitrap | 587 | 15 | 2.56% | 1.55–4.17% |

These are MS2-observed discordance proxies, not calibrated FDR estimates. They
exclude unlabelled/ambiguous features and cannot validate all-weak groups lacking
a strong labelled donor. Alternatives retained only for provenance are excluded
from this reported-primary analysis. Conditioning on the existing cell gate
substantially improves the observed subset; it does not certify every member.

In a strengthened synthetic absence fixture, 500 genuine groups survived with
62 spurious absent-run attachments among 1,562 original member links. This
shows why group existence cannot establish every member link. The separate
CLI smoke fixture supplies no measured hills: it retains 200 weak groups but
emits zero invented intensity cells, with cell q = 1. A genuine biological
absence/entrapment dataset was not available in the inspected local benchmark
inputs. Depleted samples and missing MS2 were not relabelled as known absences.

## Verification and reproduction

The final workspace checks passed: 171 library tests, one documentation test,
Clippy with warnings denied, formatting, and the CLI/export smoke check. The
ignored synthetic audit was run explicitly and its measurements are reported
above; a successful test-process exit does not mean its calibration passed.
Four existing external-fixture tests remain ignored. All four standard TSVs
from the two-run research extraction are byte-identical to the production CLI.
The original Bruker metrics also exactly reproduce the preceding prototype.

Build with `cargo build --offline --release -j 2 --bin koth_align --example
lfq_refinement`. The example accepts `BATCH CONFIG OUTPUT`, `--quantify 0,1,2,3`,
`--take 2` or `--take 6`, and `--quick` for just the first five global controls.
Without `--quick`, it evaluates 25 controls: two global sets, marginal ranks,
and two quality-conditioned rank sets. `--conditioned` chooses the diagnostic
pooled conditioned gate for extraction; this is a research executable option,
not a production setting. Unrun sensitivity columns are NaN in the final writer.

`scripts/lfq_member_audit.py` produces the label/provenance audit;
`scripts/lfq_refinement_metrics.py` evaluates output using the paper's existing
analysis modules. They require the paper's Python environment and
`--paper-scripts` path. For Orbitrap, the metric script's `--psms` argument is the
Sage directory, while the member script takes its parquet file; the metric
script also needs `--files analysis/config/files.tsv`. The conditioned run uses
`--control conditioned` in the metric script. Protocol, selection, metrics,
source snapshots, input hashes and equivalence evidence are recorded under
`out/group-refinement/`.

The first audit build wrote six rather than four decimal places for cell q-values;
metrics round them to production precision before applying the cell gate. Early
quick audits copied global q-values into unused sensitivity columns; those
columns were excluded from analysis and the final writer emits NaN instead.
The conditioned run's initially hardcoded control label was corrected in its
metadata and explicitly annotated. These reporting corrections do not alter
selected candidates or measured intensities.

The remaining landing requirement is a null model that behaves acceptably under
nonstationary background and small run counts, checked against independent
known-negative/absence data. This bounded pass stops without adding another
scorer or relaxing a threshold to conceal the failures.
