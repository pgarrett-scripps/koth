"""Screen retained groups for close duplicates and isotope aliases; not error labels."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

C13 = 1.003354835


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--labels-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    all_c = pd.read_csv(args.root/'candidates.tsv', sep='\t')
    selected = (all_c.n_contributing_runs >= 2) & (all_c.average_permuted <= .05)
    c = all_c.loc[selected].reset_index(drop=True)
    m = pd.read_parquet(args.labels_root/'member_labels.parquet')
    ownership_conflicts = int(m.duplicated(['run_id', 'feature_idx']).sum())
    m = m[m.candidate_id.isin(c.candidate_id)].copy()
    primary = m[m.is_primary]
    strong = primary[(primary.label >= 0) & (primary.quality >= .75)]
    labels = strong.groupby('candidate_id').label.agg(lambda x: frozenset(x))
    n_label_runs = strong.groupby('candidate_id').run_id.nunique()
    supported_labels = {int(i): next(iter(ls)) for i, ls in labels.items()
                       if len(ls) == 1 and n_label_runs.loc[i] >= 2}
    counts = strong.groupby('candidate_id').agg(labelled_runs=('run_id','nunique'), identities=('label','nunique'))
    masks = primary.groupby('candidate_id').run_id.agg(lambda x: sum(1 << int(i) for i in set(x)))
    intensity = pd.read_csv(args.root/'average/intensity_matrix.tsv',sep='\t')
    q = pd.read_csv(args.root/'average/qvalue_matrix.tsv',sep='\t')
    assert len(intensity) == len(c)
    assert np.allclose(intensity.massCalib, c.massCalib, rtol=0, atol=0.00000051)
    cols = list(intensity.columns[8:])
    values, qualities = intensity[cols].to_numpy(), q[cols].to_numpy()
    passed = (values > 0) & np.isfinite(qualities) & (qualities <= .05)
    rows = []
    # A deliberately tight co-elution screen, fixed before reading its results.
    # 20 ppm neutral-mass residual; same charge; <=0.1 min RT; <=0.015 IM if known.
    for charge, sub in c.groupby('charge', sort=True):
        ids = sub.sort_values('massCalib').index.to_numpy()
        masses = c.massCalib.to_numpy()[ids]
        rts = c.rtApex.to_numpy()[ids]
        ims = c.im.to_numpy()[ids]
        for k in [0,1,2]:
            expected = masses + k*C13
            width = masses * 20e-6
            lo = np.searchsorted(masses, expected-width)
            hi = np.searchsorted(masses, expected+width, side='right')
            for x in range(len(ids)):
                possible = np.arange(max(int(lo[x]), x+1), hi[x])
                drt = np.abs(rts[possible]-rts[x])
                dim = np.abs(ims[possible]-ims[x])
                valid_im = (ims[possible] == 0) | (ims[x] == 0) | (dim <= .015)
                for y in possible[(drt <= .1) & valid_im]:
                    a,b = int(ids[x]),int(ids[y])
                    ca,cb = int(c.candidate_id[a]),int(c.candidate_id[b])
                    common = passed[a] & passed[b]
                    n = int(common.sum())
                    v,w = values[a,common],values[b,common]
                    corr = float(np.corrcoef(np.log2(v),np.log2(w))[0,1]) if n>=3 and np.std(np.log2(v))>0 and np.std(np.log2(w))>0 else None
                    la,lb = supported_labels.get(ca),supported_labels.get(cb)
                    mask_a,mask_b = int(masks.get(ca,0)),int(masks.get(cb,0))
                    rows.append(dict(row_a=a,row_b=b,candidate_a=ca,candidate_b=cb,charge=int(charge),isotope_step=k,
                        mass_a=float(masses[x]),mass_b=float(masses[y]),residual_ppm=float((masses[y]-expected[x])/masses[x]*1e6),
                        rt_delta=float(abs(rts[y]-rts[x])),im_delta=float(abs(ims[y]-ims[x])),
                        primary_runs_a=mask_a.bit_count(),primary_runs_b=mask_b.bit_count(),shared_primary_runs=(mask_a&mask_b).bit_count(),
                        common_reported_runs=n,identical_intensity_runs=int(np.isclose(v,w,rtol=1e-5,atol=0).sum()),log_intensity_correlation=corr,
                        identity_a=la,identity_b=lb,same_supported_identity=la is not None and lb is not None and la==lb))
    pairs = pd.DataFrame(rows)
    pairs.to_csv(args.output/'candidate_pairs.tsv',sep='\t',index=False)
    result = dict(retained_groups=len(c),original_observation_ownership_conflicts=ownership_conflicts,
        groups_with_two_or_more_strong_labelled_runs=int((counts.labelled_runs>=2).sum()),
        groups_with_multiple_strong_primary_identities=int(((counts.labelled_runs>=2)&(counts.identities>=2)).sum()),
        screening=dict(mass_ppm=20,rt_minutes=.1,im_absolute=.015,same_charge=True),pairs={})
    for k, subset in pairs.groupby('isotope_step'):
        result['pairs'][str(k)] = dict(candidate_pairs=len(subset),unique_groups=len(set(subset.candidate_a)|set(subset.candidate_b)),
            both_reported_in_three_or_more_runs=int((subset.common_reported_runs>=3).sum()),
            near_identical_intensity_in_three_or_more_runs=int((subset.identical_intensity_runs>=3).sum()),
            same_supported_identity=int(subset.same_supported_identity.sum()),
            disjoint_primary_runs=int((subset.shared_primary_runs==0).sum()))
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
