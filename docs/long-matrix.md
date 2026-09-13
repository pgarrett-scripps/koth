# Long-format LFQ evidence

Enable the versioned long bundle in the `koth_align` configuration:

```toml
[output]
export_long = true

[lfq.consensus]
mz_ppm = 20.0
rt_window_pct = 0.02
im_tolerance = 0.05
min_member_combined_score = 0.5
min_seed_combined_score = 0.75
min_group_size = 2
```

Grouping tolerances are independent of extraction tolerances. These defaults
retain the old effective pairwise limits, but apply them to the complete group
span. Missing IM cannot bridge incompatible known IM values. Groups are formed
in deterministic quality order; minimum size and seed score are applied after
membership is complete. The primary observation in each run is selected using
its own score. Alternatives remain in the observation table.

Set `[alignment] reference_run = "exact run directory name"` to fix the reference
for a comparison. An unknown name is rejected. If omitted, the existing automatic
reference selection is used.

The additional outputs are:

| File | Row identity | Purpose |
| --- | --- | --- |
| `feature_observations.tsv` | run ID + original row ID | Original feature measurements, aligned coordinates, group membership, seed/primary flags and exclusion reasons |
| `lfq_matrix.long.tsv` | consensus ID + run ID + extraction kind | Quantitative extraction and grid quality, with original/seed observation links |
| `matrix_manifest.json` | bundle | Schema version 1, units, runs and source hashes, reference, configuration, and output hashes |

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
# Experimental retention of replicated moderate-quality features

Set `[lfq.consensus] allow_replicated_weak_seeds = true` to retain a bounded group
when at least two original runs support it, even if its best member falls below
`min_seed_combined_score`. The existing member-quality floor, full mass/RT/IM
span limits and `min_group_size` still apply. Multiple alternatives from one run
do not satisfy the two-run requirement, and weak singleton groups remain excluded.

The option defaults to false. It changes the final retention decision after
grouping, so existing retained groups keep their members and seed. Added groups
receive LFQ extraction and confidence estimates through the same path, with all
original measurements retained in the long bundle for joint identification.
