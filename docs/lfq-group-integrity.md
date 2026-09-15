# Group integrity: splits, isotope aliases and shared signal

Follow-up: [exclusive signal ownership](lfq-signal-ownership.md) implements
conflict resolution and documents the remaining membership/identity limits.
The findings below describe the pre-ownership baseline; its recommendations are
superseded by the scoped adoption decision in that follow-up.

The previous validation checked within-group duplicate handling, bounded spans,
confidence stability and quantification metrics. It did not establish one
consensus group per underlying chromatographic feature. This dedicated audit
finds concrete unresolved split/alias cases, so the recommendation to land the
whole implementation was premature. Keep the independent-permutation estimator,
but resolve cross-group identity and shared-signal conflicts before default adoption.

## Existing protections and gaps

Every projected original observation is assigned to one candidate group. Within
each group, only one primary observation per run counts as independent support;
alternatives retain provenance and reduce the evidence score. Full mass, RT and
IM span limits prevent a chain of nearby observations from bridging arbitrarily
far apart. Input projection uses the feature finder's corrected monoisotopic m/z.
These protect record ownership and local geometry, not analyte uniqueness.

Group assembly is greedy and quality-ordered, and groups are never reconsidered
or reconciled after formation. A feature at the edge of a group's allowed span
can prevent otherwise coherent observations from joining it. A residual isotope
assignment error changes the reported mass by approximately 1.003355 Da per
isotope in neutral mass, or that value divided by charge in m/z; the grouping
code does not compare alternative monoisotope hypotheses across runs. The
assembler itself documents that an M+1 peak can be emitted as a mono hypothesis
when the true mono's chain fails, and downstream correction cannot recover a
missing chain member merely from its recorded offset.

Each consensus group subsequently queries the hill database independently.
`build_grid` accumulates every matching hill in the group/run/isotope box; there
is no ownership or competition between neighboring group queries. Thus a unique
original feature assignment does not prevent two output rows from integrating
shared input signal. Neither the group statistic nor independent per-cell q-values
establish that two accepted rows are distinct analytes.

## Retained-group screen

The fixed screen uses the current pooled-permutation accepted groups at q <= 0.05,
same charge, a 20 ppm neutral-mass residual, RT difference <= 0.1 minutes and IM
difference <= 0.015 where available. It tests mass differences of 0, one and two
C13 isotope steps. These are candidate-pair screens, not measured error rates;
the isotope screen intentionally uses the broad existing grouping tolerance.

| Cohort | Nearby same-mass pairs | One-isotope-step pairs | Two-isotope-step pairs | Reused original observations |
|---|---:|---:|---:|---:|
| Bruker, 77,262 groups | 713 | 3,504 | 4,385 | 0 |
| Orbitrap, 79,004 groups | 72 | 1,287 | 1,148 | 0 |

For an independent supporting-identity check, each member group must have one
unanimous MS2 peptide/charge label among strong primary observations in at least
two original runs, using the previously cached label audit. Among Bruker's
same-mass candidate pairs, 88 share that supported identity; all 88 are reported
in at least three common runs, and 87 have near-identical reported intensities
(relative tolerance 1e-5) in at least three runs. Orbitrap has three such identity
pairs, with two reported in at least three common runs and one near-identical
in at least three runs. The identity evidence is selective and does not label
all remaining candidates as errors or validate every unidentified feature.

For one-isotope-step pairs, 23 Bruker pairs and two Orbitrap pairs share the
supported identity; 21 and two respectively are reported in at least three
common runs. Four Bruker two-step pairs share a supported identity. These are
stronger warning signs than proximity alone, but still need chromatographic and
isotope-envelope adjudication.

The audit also finds multiple strong primary MS2 labels within 4,557 of 30,025
Bruker groups with at least two strongly labelled runs, and 556 of 11,803 such
Orbitrap groups. These are possible mixed-identity groups or label/matching
errors, not a calibrated false-merge rate. They reinforce that avoiding splits
by indiscriminately merging neighbors would be unsafe.

## Inspected raw-signal examples

In the unwarped Bruker reference run, candidates **99522 and 99523** have masses
1758.9189668 and 1758.9199844 Da and RTs 10.9504962 and 10.9498787 minutes. Both
report intensity 130,585 with cell q = 0.0001. Their extraction boxes admit the
same three raw hills and 41 positive hill-profile samples, and the groups share
a supported MS2 identity. This is a concrete duplicate-reporting concern rather
than a claim inferred only from the number of groups.

Reconstructing all original members explains why the greedy grouping split them.
The first group's mass span is 19.720 ppm, the second's is 0.483 ppm, and their
union is 20.574 ppm, exceeding the 20 ppm full-span bound. RT and IM union spans
remain within bounds. A high-scoring observation (quality 0.981246) at mass
1758.884552 Da stretches the first group toward the low-mass edge. Most other
members cluster near 1758.918–1758.921 Da. A detector score floor would not remove
this high-scoring outlier; robust membership/reassignment is needed. Simply
merging the two current membership lists would violate the protective span bound.

Candidates **102260 and 102329** have masses 1798.9902573 and 1799.9924815 Da,
approximately one isotope step apart. Both groups have a supported matching MS2
identity and all inspected input features have recorded neutron offset zero.
Their reference-run intensities are 142,555 and 112,940 with cell q = 0.0006 and
0.0192. Their boxes admit three shared hills and 12 positive profile samples.
This is a supported isotope-alias concern that the current per-cell gate allows.
The raw-sample audit checks input eligibility for both grids; it does not rebuild
the final integrated peak boundaries or quantify exactly how much shared signal
is included in both reported values.

## Proposed handling before landing

1. **Reconcile neighboring candidate groups using coherent membership.** Inspect
   close same-charge mass/RT/IM groups together. Compare robust agreement of their
   original members and native elution traces, identify unsupported outliers,
   and allow reassignment while preserving the full-span protection. Do not
   widen all tolerances or merge by transitive proximity.
2. **Compare explicit isotope hypotheses.** Test zero, plus/minus one and, where
   justified, two isotope-step alternatives using accurate mass, charge spacing,
   native isotope envelopes, co-elution and evidence from other runs. Prefer a
   canonical mono only when the measurements support it; retain correction and
   ambiguity provenance. The lowest-mass candidate is not automatically correct.
3. **Resolve shared-signal conflicts during extraction.** Track original hill
   and peak/sample provenance across nearby group queries. Duplicate or shifted
   hypotheses explaining the same chromatographic envelope should compete for
   one explanation. Unresolved overlap should be marked ambiguous rather than
   confidently reported as two independent measurements. Distinct charge states
   and genuinely resolved RT/IM peaks should remain distinct measurements.
4. **Recalibrate after structural changes.** Apply candidate reconciliation to
   controls as well as targets, recalculate group confidence, and repeat cell
   extraction/scoring. Add known split, wrong-mono, wrong-charge, neighboring
   analyte and shared-signal cases, and rerun the existing cohort comparisons.

A mass difference near one Da is not sufficient evidence for an isotope error.
For example, deamidation adds 0.984016 Da ([Unimod](https://www.unimod.org/modifications_view.php?editid1=7)),
whereas a C13 step is about 1.003355 Da. At typical peptide masses the broad
20 ppm screen can include both; [Matrix Science's discussion](https://www.matrixscience.com/help/the_plus_one_dilemma.html)
explains this ambiguity. The screen therefore must not become an automatic
merge rule. Similar caution applies to isomers, conformers and partially
co-eluting different peptides.

## Scope and verification

This pass adds an integrity screen, raw shared-signal inspection, aligned-member
geometry export, and three diagnostic regression tests. The tests reproduce RT
fragmentation, separate grouping after a residual mono error, and reuse of one
hill by two isotope-offset queries. They demonstrate current limitations rather
than implement their resolution. No automatic merging, isotope correction,
production grouping change, or main-branch merge was performed in this pass.
All 43 LFQ-targeted tests pass; the existing ignored research stress test remains
separate. Clippy and formatting checks pass.

Artifacts are in `out/group-integrity/`. Run `scripts/lfq_group_integrity_audit.py`
with `--root out/group-permutation/COHORT --labels-root out/group-refinement/COHORT
--output out/group-integrity/COHORT`. Run `scripts/lfq_shared_signal_audit.py` with
that root and output for the inspected reference-run examples. Both use the
paper's existing Python environment. The `lfq_member_geometry` Rust example
reconstructs original aligned member coordinates for selected cached candidate
IDs using the original batch, config and members table.
