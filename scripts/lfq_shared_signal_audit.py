"""Inspect two flagged pairs in the unwarped reference run against raw hill samples."""
import argparse
import json
import tomllib
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    root,out = args.root,args.output
    pairs = pd.read_csv(out/'candidate_pairs.tsv',sep='\t')
    c = pd.read_csv(root/'candidates.tsv',sep='\t')
    c = c[(c.n_contributing_runs>=2)&(c.average_permuted<=.05)].reset_index(drop=True)
    i = pd.read_csv(root/'average/intensity_matrix.tsv',sep='\t')
    q = pd.read_csv(root/'average/qvalue_matrix.tsv',sep='\t')
    meta = json.loads((root/'audit.json').read_text())
    runs = json.loads((root/'runs.json').read_text())
    ref = next(r for r in runs if r['name']==meta['reference'])
    config = tomllib.loads((root/'average/align_config.toml').read_text())['lfq']
    f = pd.read_parquet(ref['features'],columns=['charge','rtStart','rtEnd'])
    f = f[f.charge!=0]
    half = config['rt_window_pct']*(f.rtEnd.max()-f.rtStart.min())
    results=[]
    for step in [0,1]:
        subset = pairs[(pairs.isotope_step==step)&pairs.same_supported_identity&(pairs.common_reported_runs>=3)]
        choices=[]
        for p in subset.itertuples():
            a,b=int(p.row_a),int(p.row_b)
            if i.loc[a,ref['name']]>0 and i.loc[b,ref['name']]>0 and q.loc[a,ref['name']]<=.05 and q.loc[b,ref['name']]<=.05:
                choices.append(p)
        if not choices:continue
        chosen=max(choices,key=lambda p:(p.identical_intensity_runs if step==0 else p.log_intensity_correlation,-abs(p.residual_ppm)))
        a,b=int(chosen.row_a),int(chosen.row_b)
        queries=[]
        for row in [a,b]:
            g=c.iloc[row]
            for isotope in range(config['n_isotopes']):
                mass=g.mz+isotope*1.003354835/g.charge
                tol=mass*config['mz_ppm']/1e6
                queries.append([('mz','>=',mass-tol),('mz','<=',mass+tol)])
        hills=pq.read_table(Path(ref['features']).with_name('hills.parquet'),filters=queries,
            columns=['hill_id','mz','rt','rt_start','rt_end','im','intensity_profile']).to_pandas()
        signals=[]
        for row in [a,b]:
            g=c.iloc[row];samples={};hill_ids=set()
            for h in hills.itertuples():
                mass_match=any(abs(h.mz-(g.mz+k*1.003354835/g.charge)) <= (g.mz+k*1.003354835/g.charge)*config['mz_ppm']/1e6 for k in range(config['n_isotopes']))
                im_match=g.im==0 or h.im==0 or abs(g.im-h.im)<=config['im_tolerance']
                if not mass_match or not im_match:continue
                profile=json.loads(h.intensity_profile)
                # read_hills gives no scan_times: production uses this same interpolation.
                times=np.linspace(h.rt_start,h.rt_end,len(profile)) if len(profile)>1 else [h.rt]
                for index,(rt,value) in enumerate(zip(times,profile)):
                    if value>0 and g.rtApex-half<=rt<g.rtApex+half:
                        samples[(int(h.hill_id),index)]=value;hill_ids.add(int(h.hill_id))
            signals.append((samples,hill_ids))
        common=signals[0][0].keys()&signals[1][0].keys()
        results.append(dict(isotope_step=step,candidate_a=int(c.candidate_id[a]),candidate_b=int(c.candidate_id[b]),
            mass_a=float(c.massCalib[a]),mass_b=float(c.massCalib[b]),rt_a=float(c.rtApex[a]),rt_b=float(c.rtApex[b]),
            reference_run=ref['name'],reference_intensity_a=float(i.loc[a,ref['name']]),reference_intensity_b=float(i.loc[b,ref['name']]),
            reference_q_a=float(q.loc[a,ref['name']]),reference_q_b=float(q.loc[b,ref['name']]),
            same_supported_ms2_identity=True,shared_hills=len(signals[0][1]&signals[1][1]),
            shared_positive_samples=len(common),shared_raw_input_intensity=float(sum(signals[0][0][s] for s in common)),
            scope='Input samples eligible for both grids; not a reconstruction of final integrated peak boundaries.'))
    (out/'shared_signal_examples.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
