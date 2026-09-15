# Exclusive LFQ signal ownership

The subsequent [shared-envelope isotope trial](lfq-isotope-choice.md) was
rejected after worsening the calculated-mass audit. The implementation described
below remains active; the rejected selector is archived separately.

**Final adoption decision:** retain group confidence and exclusive native peak
ownership, with the failed shared-envelope selector excluded. This prevents the
demonstrated signal reuse; it does not guarantee the correct mono hypothesis.
The isotope-choice errors and Orbitrap precision regression below remain known
limitations. Earlier recommendations to defer adoption are preserved as the
experimental record and superseded by this scoped decision.

This worktree now evaluates competing extraction hypotheses before assigning
signal. Raw hill profiles remain immutable. The original consensus membership,
full-span limits, evidence scores, and ten permutation controls are unchanged;
this layer reconciles competing **quantifications**, without unioning groups or
claiming that an isotope-shifted mass is automatically the correct identity.

## Implementation

Quantification streams the runs twice. The first pass records native isotope
pattern and co-elution evidence for each target and decoy hypothesis. Each fit is reduced by the fraction of co-eluting signal at the preceding
isotope position that this mono hypothesis leaves unexplained. The mean adjusted
fit across signal-bearing runs determines a fixed preference across runs.
Ties use detection count and measured coordinates. This prevents the relative
preference between two aliases from switching with run processing order. It is
an evidence ranking, not a probability or a confirmed monoisotope assignment.

The second pass quantifies in that preference order. Every grid records original
`(hill index, profile index)` identities, and integration identifies which
samples support the selected peak. Ownership reserves the native chromatographic
segment containing each used sample. A hill with one peak is one segment; a long
hill is split at a valley below half both neighboring peak heights, after a
three-point triangular smoothing used only to determine ownership boundaries.
The stored intensities are never smoothed or subtracted. Reserving a segment is
necessary because a sparse grid may integrate only one scan; claiming that scan
alone allowed another group to extract the adjacent slice of the same peak.

When a candidate's grid includes reserved signal, its grid, isotope scores,
co-elution and integration are rebuilt with those original samples excluded.
A residual must contain at least two isotope rows. If most of the originally
selected samples were reserved, a surviving second peak must have a disjoint
integration interval and lie closer to the candidate's predicted RT than the
original shared peak. Otherwise the residual is marked ambiguous and its
intensity is zero. Separate m/z, charge, IM and resolved chromatographic signals
remain available when they do not use the same native segment.

Targets and decoys have separate ownership masks and independently computed
preferences, using the same code and rules. A target never depletes a decoy's
input. Cell QDA and q-values are recomputed after conflict resolution. This is a
necessary symmetry safeguard, not a new validation of FDR calibration or of
identification correctness. Group q-values remain those of their original
candidate groups; no confidence is inherited from a pooled alias group.

## Outputs and compatibility

`lfq_details.tsv` and long-bundle schema version 3 expose ownership status,
number of integrated owned samples, number of grid samples excluded during
reassessment, the highest-preference competing consensus ID, and the preceding co-eluting
signal fraction used to compare mono hypotheses. Statuses are
`exclusive`, `residual`, `shared_signal_only`, and `ambiguous_residual`. The last
two have zero intensity and target q = 1. Counts of integrated samples differ
from the larger native segments reserved to prevent adjacent-scan reuse.

Original group rows and observation links remain available, including fully
suppressed rows, rather than silently merging their provenance or changing their
group confidence. The research audit additionally exports `signal_aliases.tsv`
for entirely suppressed groups consistently explained by one competitor in at
least two runs. Mixed-owner cases are not labeled a single canonical identity.
A retained residual row represents independently extracted signal, not proof of
a distinct peptide. Detector scores are not used as an ownership admission gate.

Exclusive extraction requires grid intensities in every cell. The historical
`detected_use_grid=false` setting produces a warning and uses grid intensities;
feature-supplied intensity overrides cannot establish native sample ownership.
No new user configuration setting was added.

## Limits

A reserved native segment is a conservative conflict unit. Very close unresolved
co-eluting species, small shoulders, and isotope/charge interference can remain
ambiguous; this implementation does not deconvolve them or correct their masses
from an isotope-step difference alone. A fixed global preference can favor a
wrong explanation, and preserving distinct quantities is not the same as proving
distinct identities. The original geometric group-splitting and mixed-member
identity issues are not repaired by rewriting membership in this pass.

An intermediate native-segment version still preferred the +1-isotope candidate
in the inspected Bruker isotope case. The independently identified peptide has
calculated neutral mass 1798.9827 Da; the lower candidate at 1798.9903 Da is the
compatible mono hypothesis. This motivated the general preceding-isotope check,
which uses only native MS1 signal at runtime, not the peptide identity or mass.
Tests also require the higher-mass candidate to win when the lower hypothesis
has no M peak, and impose no penalty for a preceding peak that does not co-elute.

The initial exact-sample-only prototype is retained in `out/group-ownership/`.
It removed identical intensity duplication but left many pairs reporting adjacent
slices of the same peak. The intermediate native-segment revision is in `out/group-ownership-peaks/`.
The revision with preceding-isotope evidence and deduplicated segment reservation
is evaluated in `out/group-ownership-final/`.

## Reproduction

Run `lfq_permutation_audit BATCH CONFIG OUTPUT --quantify` for each cohort, then
`scripts/lfq_ownership_audit.py --root OUTPUT --baseline out/group-permutation/COHORT
--integrity out/group-integrity/COHORT`. The audit verifies all original candidate
coordinates and group q-values are unchanged, compares the previously identified
pairs, and reports ownership status counts and the two raw-signal witnesses.
`scripts/lfq_refinement_metrics.py` performs the existing fixed-anchor coverage,
CV and fold-change comparisons against the permutation baseline.

## Final benchmark results

The final version uses native-segment ownership and preceding-isotope evidence.
Both cohorts retain exactly the baseline candidate set, coordinates, group scores
and group q-values. Every final positive cell passes a runtime assertion that its
integrated original samples have not already been reserved by another cell in
the same run and target/decoy world. All 79 input-file hashes remain unchanged.

| Metric | Bruker baseline | Bruker exclusive | Orbitrap baseline | Orbitrap exclusive |
| --- | ---: | ---: | ---: | ---: |
| Groups (audit rows retained) | 77,262 | 77,262 | 79,004 | 79,004 |
| PSM recall | 65.658% | 65.801% | 49.385% | 49.957% |
| Median CV | 10.274% | 10.331% | 16.039% | 16.422% |
| Human false fold-change rate | 15.273% | 15.099% | 8.989% | 9.544% |
| Reported cells at cell q <= .05 | 951,589 | 946,874 | 1,251,652 | 1,264,796 |

The final shared-cell CV comparison is 10.221% → 10.247% for Bruker and
15.995% → 16.184% for Orbitrap. Thus the Orbitrap regression is not solely a
change in which cells are present. Recall gains after recomputing cell q-values
do not establish improved FDR calibration. The peak-segment ownership stage can
remove decoy as well as target signals; procedural symmetry alone does not prove
that their selected distributions remain a calibrated error model.

Among the previously screened Bruker same-mass pairs with supported matching MS2
identities, pairs reported together in at least three runs drop from 88 to 2;
the 87 pairs with near-identical intensities in at least three runs drop to zero.
The analogous Orbitrap co-reported count drops from two to zero. Bruker one-step
isotope pairs reported together in at least three runs drop from 21 to 5 and
two-step pairs from four to three. These counts describe selected warning cases,
not a measured duplicate-error rate or evidence that every surviving pair is wrong.

The reference-run duplicate witness changes from 130,585 / 130,585 to
130,585 / 0. The isotope witness changes from 142,555 / 112,940 to 142,555 / 0,
retaining the candidate compatible with the identified peptide's neutral mass.
The audit links 938 fully suppressed Bruker groups and 1,218 Orbitrap groups to
a consistent sole competitor across at least two runs; these are reporting
aliases, without merged group membership or pooled confidence.

## Why this is not a landing recommendation

Fixing the inspected isotope witness was insufficient. A broader Bruker audit
screens 27 one/two-isotope pairs sharing supported MS2 labels. Of these, 25 pairs
have one candidate compatible with the identified peptide's calculated neutral
mass within 10 ppm, giving 450 pair/run comparisons. At the cell-q gate:

| Outcome relative to peptide mass | Baseline | Native segments only | Segments + preceding evidence |
| --- | ---: | ---: | ---: |
| Only compatible candidate reported | 69 | 222 | 199 |
| Only incompatible candidate reported | 36 | 77 | 83 |
| Both reported | 310 | 97 | 81 |
| Neither reported | 35 | 54 | 87 |

These correlated, selectively labelled comparisons are not an FDR estimate.
They nevertheless show that the current ranking can choose the incompatible
mass when signal becomes exclusive. The preceding-isotope term fixes the chosen
witness but does not improve this broader comparison over segment ownership
alone. It must not be described as a validated automatic isotope correction.
A reliable comparison of competing isotope envelopes, with explicit ambiguity
when the mass cannot be determined, remains a release blocker. Original group
membership/outlier reassignment is also still unresolved; this pass reconciles
quantification, not the underlying membership lists.

Keep exclusive signal accounting as the demonstrated mechanism, but do not
promote the current global winner ranking to the default on these results.
The native-segment boundary convention also needs evaluation on closer peaks and
co-eluting mixtures beyond the synthetic preservation tests used here.

## Verification and artifacts

All 188 library tests and the documentation test pass; one existing research
stress test remains ignored. Eleven new ownership tests cover repeated queries,
sparse-scan slices, isotope shifts in both directions, non-co-eluting preceding
peaks, distinct charge/mass signals, a recovered second peak on the same hills,
run/candidate ordering, target/decoy symmetry, input preservation and empty data.
Clippy with warnings denied, formatting and CLI long-export smoke checks pass.
The smoke test retains 200 weak groups and emits 800 cells without invented signal.

`out/group-ownership-final/` contains full-cohort matrices, ownership records,
metrics, pair audits, calculated-mass audit, test logs and provenance. The broader
mass check is reproducible with its saved `mono_identity_audit.py`; the peptide
labels are audit evidence only and are never used by the extraction algorithm.
The preceding term and segment reservation optimization were evaluated together
in the final full-cohort run. Subsequent cleanup removes an unused intensity
optional override and fixes the long-export diagnostic column; those changes
were covered by the final unit, Clippy and CLI export checks.
