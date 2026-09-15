"""Compare exclusive extraction with its fixed, unchanged candidate-set baseline."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--baseline', type=Path, required=True)
    ap.add_argument('--integrity', type=Path, required=True)
    args = ap.parse_args()
    old = pd.read_csv(args.baseline / 'candidates.tsv', sep='\t')
    new = pd.read_csv(args.root / 'candidates.tsv', sep='\t')
    pd.testing.assert_frame_equal(old, new)
    selected = new[(new.n_contributing_runs >= 2) & (new.average_permuted <= .05)]
    rows = dict(zip(selected.candidate_id.astype(int), range(len(selected))))
    folder = args.root / 'average'
    table = pd.read_csv(folder / 'ownership.tsv', sep='\t')
    meta = ['massCalib','mz','charge','rtApex','im','combined_score','seed_run','n_contributing_runs']
    intensity = pd.read_csv(folder / 'intensity_matrix.tsv', sep='\t').drop(columns=meta)
    q = pd.read_csv(folder / 'qvalue_matrix.tsv', sep='\t').drop(columns=meta)
    old_i = pd.read_csv(args.baseline / 'average/intensity_matrix.tsv', sep='\t').drop(columns=meta)
    old_q = pd.read_csv(args.baseline / 'average/qvalue_matrix.tsv', sep='\t').drop(columns=meta)
    reported = (intensity.to_numpy() > 0) & (q.to_numpy() <= .05)
    old_reported = (old_i.to_numpy() > 0) & (old_q.to_numpy() <= .05)
    summary = {'unchanged_group_candidates_and_qvalues': True, 'n_groups': len(selected),
               'reported_cells_before': int(old_reported.sum()), 'reported_cells_after': int(reported.sum()),
               'all_zero_target_rows': int((intensity == 0).all(axis=1).sum())}
    for side, data in table.groupby('is_decoy'):
        name = 'decoy' if side else 'target'
        summary[name] = {'statuses': data.ownership_status.value_counts().to_dict(),
                         'sample_exclusion_events': int(data.excluded_samples.sum()),
                         'owned_samples': int(data.owned_samples.sum())}
        assert not ((data.intensity == 0) & (data.owned_samples != 0)).any()
        assert not ((data.intensity > 0) & (data.owned_samples == 0)).any()
    pairs = pd.read_csv(args.integrity / 'candidate_pairs.tsv', sep='\t')
    pairs['common_reported_runs_after'] = [int((reported[int(a)] & reported[int(b)]).sum()) for a,b in zip(pairs.row_a,pairs.row_b)]
    pairs['identical_intensity_runs_after'] = [int((reported[int(a)] & reported[int(b)] & np.isclose(intensity.iloc[int(a)].to_numpy(),intensity.iloc[int(b)].to_numpy(),rtol=1e-5,atol=0)).sum()) for a,b in zip(pairs.row_a,pairs.row_b)]
    pairs.to_csv(args.root / 'pairs_after.tsv',sep='\t',index=False)
    for step in range(3):
        part = pairs[(pairs.isotope_step == step) & pairs.same_supported_identity]
        summary[f'same_identity_step_{step}'] = {'n_pairs': len(part),
            'reported_3plus_before': int((part.common_reported_runs >= 3).sum()),
            'reported_3plus_after': int((part.common_reported_runs_after >= 3).sum()),
            'identical_3plus_before': int((part.identical_intensity_runs >= 3).sum()),
            'identical_3plus_after': int((part.identical_intensity_runs_after >= 3).sum())}
    witnesses = []
    reference = json.loads((args.root/'audit.json').read_text())['reference']
    run_index = list(intensity.columns).index(reference)
    for a,b in [(99522,99523),(102260,102329)]:
        if a not in rows or b not in rows: continue
        if 'bruker' not in str(args.root): continue
        indices = [rows[a], rows[b]]
        records = table[(~table.is_decoy) & table.feature_idx.isin(indices) & (table.run_idx == run_index)]
        witnesses.append({'candidate_ids':[a,b], 'reference_run':reference,
            'before':[float(old_i.iloc[r,run_index]) for r in indices],
            'after':[float(intensity.iloc[r,run_index]) for r in indices],
            'q_after':[float(q.iloc[r,run_index]) for r in indices],
            'cells':json.loads(records.to_json(orient='records'))})
    summary['witnesses'] = witnesses
    # An entirely suppressed group can be linked only when the same measured
    # competitor explains it in >=2 runs. Do not merge members or pool group q.
    targets = table[~table.is_decoy]
    active = targets.groupby('feature_idx').intensity.max() > 0
    aliases = []
    for g, part in targets[targets.excluded_samples > 0].groupby('feature_idx'):
        owners = part.competing_feature.dropna().unique()
        if not active.loc[g] and len(part) >= 2 and len(owners) == 1:
            aliases.append({'suppressed_group':int(g),'preferred_group':int(owners[0]),'runs':len(part)})
    pd.DataFrame(aliases, columns=['suppressed_group','preferred_group','runs']).to_csv(args.root/'signal_aliases.tsv',sep='\t',index=False)
    summary['fully_suppressed_single_owner_groups'] = len(aliases)
    (args.root/'ownership_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
