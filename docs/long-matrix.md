# Long-format LFQ evidence

Enable the versioned long bundle in the `koth_align` configuration:

```toml
[output]
export_long = true

[lfq.consensus]
mz_ppm = 20.0
rt_window_pct = 0.02
im_tolerance = 0.05
max_group_qvalue = 0.05
```

Grouping tolerances are independent of extraction tolerances. These defaults
retain the old effective pairwise limits, but apply them to the complete group
span. Missing IM cannot bridge incompatible known IM values. Groups are formed
in deterministic quality order; cross-run evidence and the experimental group
q-value gate replace member and seed floors. The primary observation in each run is selected using
its own score. Alternatives remain in the observation table.

Set `[alignment] reference_run = "exact run directory name"` to fix the reference
for a comparison. An unknown name is rejected. If omitted, the existing automatic
reference selection is used.

The additional outputs are:

| File | Row identity | Purpose |
| --- | --- | --- |
| `feature_observations.tsv` | run ID + original row ID | Original feature measurements, aligned coordinates, group membership, seed/primary flags and exclusion reasons |
| `lfq_matrix.long.tsv` | consensus ID + run ID + extraction kind | Quantitative extraction and grid quality, with original/seed observation links |
| `matrix_manifest.json` | bundle | Schema version 3, units, runs and source hashes, reference, configuration, and output hashes |

The observation table streams the source TSV or Parquet measurements directly:
alignment's reconstructed single-hill placeholders do not contain the original
isotope counts. Each original source row is retained, including alternatives,
failed group members and charge-zero observations. Source row IDs are scoped to
the source hashes in this bundle. Consensus IDs can change after regrouping.

Original `feature_*` measurements are distinct from `lfq_*` grid measurements.
A target cell's optional `observation_id` links to its actual original primary
feature. An MBR-only target and a synthetic decoy have no such original feature.
The seed link describes the consensus reference and must not be confused with a
direct observation in that cell's run.

Quantitative intensity is written once per cell, using the same final normalized
or unnormalized value as the wide export. Empty numeric fields mean unavailable;
they do not mean zero. Target extraction q-values are omitted when TDC is off;
synthetic decoys do not inherit their paired target's q-value. Signal/no-signal
status describes an attempted extraction. This is extraction confidence, not
peptide or protein identification confidence.

The manifest is installed only after both tables have been flushed and hashed.
Wide exports and the existing `lfq_details.tsv` remain available for compatibility.
# Experimental group confidence

The old member/seed/minimum-size settings and `allow_replicated_weak_seeds` are
removed. Every valid candidate can contribute evidence, with two original runs
required for cross-run support. `max_group_qvalue` is the sole group-quality gate.

Schema version 2 adds `group_score` and `group_qvalue` to `lfq_matrix.long.tsv`.
These describe the parent target consensus group, including on synthetic decoy
extraction rows; they are not confidences for the decoy signal. `lfq_q_value`
remains separate and is blank for decoys. Ungrouped observations use
`group_confidence_or_singleton` or `invalid_feature_score` instead of the removed
member-quality exclusion. See [lfq-group-confidence.md](lfq-group-confidence.md).

# Exclusive signal provenance

Schema version 3 adds `lfq_ownership_status`, `lfq_owned_samples`,
`lfq_excluded_samples`, `lfq_competing_consensus_id`, and
`lfq_preceding_signal_fraction`. Suppressed shared or
ambiguous residual signal has intensity zero; original observations and group
statistics are retained. See [lfq-signal-ownership.md](lfq-signal-ownership.md)
for native peak ownership, control symmetry and limitations.
