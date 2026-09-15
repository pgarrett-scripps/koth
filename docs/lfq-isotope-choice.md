# Shared-envelope isotope comparison: rejected experiment

The proposed isotope selector made the labelled Bruker comparison worse and was
removed from the active pipeline. The previous exclusive-signal implementation
is restored; no settings or export schema changes from this experiment remain.
The experiment and its evidence are saved under `out/isotope-competition/`.

## What was tested

Competing groups were compared only when they selected the same native hill peak
segment, had the same charge, and differed by one or two isotope steps within
the existing mass tolerance. Both hypotheses were scored against one common XIC
observation. Isotope rows were projected onto a chromatographic trace made from
the native segments shared by both hypotheses, then compared with shifted
theoretical envelopes using squared cosine similarity. This replaced the prior
preceding-isotope penalty; independent envelope/co-elution ranking was the fallback.

A preference required at least two runs and had to survive removal of any one
run. Stable preferences constrained the global ownership order. Conflicting or
unstable preferences were flagged; fallback quantities were retained rather than
silently excluded from evaluation. Targets and decoys used separate comparisons
with identical rules. No peptide labels or new user settings entered extraction.

## Results

Both complete cohorts were rerun with fixed configurations and reference runs:
18 Bruker runs and 20 Orbitrap runs. All 77,262 and 79,004 selected groups,
respectively, retained their previous membership and group q-values. The table
compares the preceding-isotope version with the rejected shared-envelope trial.

| Metric | Bruker previous | Bruker trial | Orbitrap previous | Orbitrap trial |
| --- | ---: | ---: | ---: | ---: |
| PSM recall | 65.801% | 66.094% | 49.957% | 50.049% |
| Median CV | 10.331% | 10.349% | 16.422% | 16.427% |
| Human false fold-change rate | 15.099% | 15.122% | 9.544% | 9.508% |

The original Bruker isotope screen contains 25 pairs compatible with the
calculated peptide mass within 10 ppm, across 18 runs (450 comparisons). Each
group must have the same independently assigned MS2 identity on high-quality
primary members in at least two runs. Cells below include only positive
quantities passing the existing cell q-value cutoff of 0.05, without removing
ambiguity-flagged cases.

| Outcome | Previous | Trial |
| --- | ---: | ---: |
| Correct candidate only | 199 | 119 |
| Wrong candidate only | 83 | 134 |
| Both candidates | 81 | 104 |
| Neither candidate | 87 | 93 |

The new native-signal screen found 32 labelled Bruker pairs, including nine
additional pairs outside the previous screen. On all 32, correct-only counts
fell from 228 to 149 and wrong-only counts rose from 85 to 136. The nine
additional pairs showed no output change. Only 13 of 27 stable pair preferences
agreed with the calculated mass; 28 of 32 preferences agreed between disjoint
even/odd run sets. Repeatability therefore did not establish correctness, and
the ambiguity flag failed to capture many wrong preferences. These selective,
correlated pair/run observations are diagnostics, not an isotope error-rate or
FDR estimate; the additional pairs are not an independent biological holdout.

Orbitrap supplied only two labelled pairs in the original screen, so its isotope
check is too small to establish generality. Their outcomes were unchanged: 20
correct-only, zero wrong-only, 20 both, zero neither. One pair had actual shared
native-segment comparisons, with the correct stable preference.

Exclusive signal accounting remained intact in the rejected trial. The 87
Bruker same-identity pairs previously sharing identical quantities in at least
three runs still fell to zero; both previously inspected witnesses remained
resolved. This does not rescue the selector: choosing a consistent wrong mass
can leave aggregate CV almost unchanged.

## Decision and verification

Reject this selector and retain the prior quantification behavior. The useful
addition is a reproducible isotope-choice audit, which prevents unchanged CV
from being mistaken for correct monoisotope assignment. The cause of the
real-data envelope-ranking errors is not established by this experiment;
interference, incomplete isotope observations, and incorrect group membership
remain hypotheses to investigate, rather than reasons to add another cutoff.

The proposal passed 194 library tests and one documentation test, with one
existing research stress test ignored. Clippy and the schema-4 trial export smoke
test passed, demonstrating that synthetic checks alone missed the scientific
regression. After removing the proposal, the original production sources were
checked against the previous recorded SHA-256 hashes, and the original tests,
Clippy, release build, and schema-3 export smoke check were repeated.

Artifacts include cohort matrices, `isotope_choice_summary.json`,
`isotope_truth_pairs.tsv`, `isotope_decisions.tsv`, ownership audits, quantitative
metrics, and test logs. `proposed-source/` and `rejected-proposal.patch` preserve
the rejected code separately from the active pipeline. The new audit is
`scripts/lfq_isotope_choice_audit.py`; it uses peptide labels only for evaluation.
Nothing was merged or pushed, and canonical paper outputs were not modified.

## Reproduce the audit

Use the paper analysis environment, which supplies pandas and the existing
validation helpers. For the saved Bruker trial, run from this worktree:

```bash
/home/patrick-garrett/Repos/koth-paper/analysis/.venv/bin/python scripts/lfq_isotope_choice_audit.py \
  --root out/isotope-competition/bruker \
  --baseline out/group-ownership-final/bruker \
  --labels-root out/group-refinement/bruker \
  --integrity out/group-integrity/bruker \
  --psms /home/patrick-garrett/Data/dnoise/results/dda_15min/original/results.sage.tsv \
  --platform bruker \
  --paper-scripts /home/patrick-garrett/Repos/koth-paper/analysis/scripts
```

For Orbitrap, replace cohort paths with `orbitrap`, use `--platform orbitrap`,
and supply `benchmark/data/sage/results.sage.parquet` in the paper repository as
the PSM file. The saved trial matrices require the archived proposal, not the
restored production binary; source and result hashes are recorded in provenance.
