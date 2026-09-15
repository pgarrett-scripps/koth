"""Independent peptide-mass checks for shared-envelope isotope choices.

Labels are evaluation-only. Report both raw output choices and ambiguity;
never improve the headline by silently removing flagged or missing cells.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

META=['massCalib','mz','charge','rtApex','im','combined_score','seed_run','n_contributing_runs']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--labels-root',type=Path,required=True)
    ap.add_argument('--integrity',type=Path,required=True)
    ap.add_argument('--psms',type=Path,required=True)
    ap.add_argument('--platform',choices=['bruker','orbitrap'],required=True)
    ap.add_argument('--paper-scripts',type=Path,required=True)
    args=ap.parse_args()
    if args.platform=='bruker':
        sys.path.insert(0,str(args.paper_scripts))
        import bruker_validation as bv
        p=bv.load_sage_psms(args.psms)
        raw=pd.read_csv(args.psms,sep='\t',usecols=['peptide','charge','calcmass','peptide_q','label'])
        raw=raw[(raw.peptide_q<=.01)&(raw.label==1)]
    else:
        raw=pd.read_parquet(args.psms,columns=['filename','peptide','charge','calcmass','expmass','rt','peptide_q','is_decoy'])
        raw=raw[(raw.peptide_q<=.01)&~raw.is_decoy].copy()
        p=raw.copy();p['stem']=p.filename.str.removesuffix('.gz')
        p['precursor_mz']=p.expmass/p.charge+1.007276466621
        p=p.groupby(['stem','peptide','charge'],sort=False).agg(precursor_mz=('precursor_mz','median'),rt=('rt','median')).reset_index()
    p['identity']=pd.factorize(p.peptide.str.replace('I','L',regex=False)+'/'+p.charge.astype(str))[0]
    names=p[['peptide','charge','identity']].drop_duplicates()
    truth=raw.groupby(['peptide','charge']).calcmass.median().reset_index().merge(names,on=['peptide','charge']).groupby('identity').calcmass.median()
    c=pd.read_csv(args.root/'candidates.tsv',sep='\t')
    assert c.equals(pd.read_csv(args.baseline/'candidates.tsv',sep='\t')), 'Candidate IDs changed between comparisons'
    selected=c[(c.n_contributing_runs>=2)&(c.average_permuted<=.05)].reset_index(drop=True)
    labels=pd.read_parquet(args.labels_root/'member_labels.parquet')
    strong=labels[labels.is_primary&(labels.quality>=.75)&(labels.label>=0)]
    g=strong.groupby('candidate_id').agg(n_labels=('label','nunique'),n_runs=('run_id','nunique'),identity=('label','first'))
    identity=g[(g.n_labels==1)&(g.n_runs>=2)].identity
    label_by_row=selected.candidate_id.map(identity)
    old_pairs=pd.read_csv(args.integrity/'candidate_pairs.tsv',sep='\t')
    screened=old_pairs[old_pairs.same_supported_identity&(old_pairs.isotope_step>0)]
    seen=set(zip(screened.row_a.astype(int),screened.row_b.astype(int)))
    comparisons=pd.read_csv(args.root/'average/isotope_comparisons.tsv',sep='\t')
    comparisons=comparisons[~comparisons.is_decoy]
    actual=comparisons[['lower','upper']].drop_duplicates()
    records=[]
    for a,b in actual.itertuples(index=False,name=None):
        la,lb=label_by_row.iloc[a],label_by_row.iloc[b]
        if pd.isna(la) or pd.isna(lb) or la!=lb:continue
        mass=float(truth.loc[int(la)])
        da,db=abs(selected.iloc[a].massCalib-mass),abs(selected.iloc[b].massCalib-mass)
        if min(da,db)/mass*1e6>10:continue
        good,bad=(a,b) if da<db else (b,a)
        records.append(dict(lower=int(a),upper=int(b),good=int(good),bad=int(bad),identity=int(la),theoretical_mass=mass,previously_screened=(a,b) in seen))
    pairs=pd.DataFrame(records,columns=['lower','upper','good','bad','identity','theoretical_mass','previously_screened'])
    pairs.to_csv(args.root/'isotope_truth_pairs.tsv',sep='\t',index=False)
    ownership=pd.read_csv(args.root/'average/ownership.tsv',sep='\t',usecols=['feature_idx','run_idx','is_decoy','isotope_assignment_ambiguous'])
    ownership=ownership[~ownership.is_decoy]
    n_runs=len(json.loads((args.root/'runs.json').read_text()))
    ambiguity=np.zeros((len(selected),n_runs),dtype=bool)
    ambiguity[ownership.feature_idx,ownership.run_idx]=ownership.isotope_assignment_ambiguous
    matrices={}
    for name,root in [('previous',args.baseline),('shared_envelope',args.root)]:
        i=pd.read_csv(root/'average/intensity_matrix.tsv',sep='\t').drop(columns=META).to_numpy()
        q=pd.read_csv(root/'average/qvalue_matrix.tsv',sep='\t').drop(columns=META).to_numpy()
        matrices[name]=(i>0)&(q<=.05)
    summary={'method':'Same strong-primary MS2 identity in >=2 runs per group; calculated peptide mass compatible within 10 ppm. Pair/run counts are correlated diagnostic observations, not FDR.', 'actual_shared_signal_pairs':len(pairs),'previously_screened_pairs':int(pairs.previously_screened.sum())}
    for subset,frame in [('all',pairs),('additional_pairs',pairs[~pairs.previously_screened])]:
        part={}
        for name,keep in matrices.items():
            counts=dict(pairs=len(frame),correct_only=0,wrong_only=0,both=0,neither=0)
            for r in frame.itertuples():
                good,bad=keep[r.good],keep[r.bad]
                counts['correct_only']+=int((good&~bad).sum());counts['wrong_only']+=int((bad&~good).sum())
                counts['both']+=int((good&bad).sum());counts['neither']+=int((~good&~bad).sum())
            part[name]=counts
        summary[subset]=part
    # Preserve the exact previous screen as a separate denominator, including
    # pairs for which no common native segment was found by the new comparison.
    legacy=[]
    for r in screened.itertuples():
        mass=float(truth.loc[int(r.identity_a)])
        da,db=abs(r.mass_a-mass),abs(r.mass_b-mass)
        if min(da,db)/mass*1e6>10:continue
        legacy.append((int(r.row_a),int(r.row_b)) if da<db else (int(r.row_b),int(r.row_a)))
    summary['original_screen']={}
    for name,keep in matrices.items():
        counts=dict(pairs=len(legacy),correct_only=0,wrong_only=0,both=0,neither=0)
        for good,bad in legacy:
            a,b=keep[good],keep[bad]
            counts['correct_only']+=int((a&~b).sum());counts['wrong_only']+=int((b&~a).sum());counts['both']+=int((a&b).sum());counts['neither']+=int((~a&~b).sum())
        summary['original_screen'][name]=counts
    decisions=[]
    for r in pairs.itertuples():
        evidence=comparisons[(comparisons.lower==r.lower)&(comparisons.upper==r.upper)]
        delta=(evidence.lower_fit-evidence.upper_fit).to_numpy()
        stable=len(delta)>=2 and ((delta.sum()-delta.max()>64*np.finfo(float).eps*len(delta)) or (delta.sum()-delta.min() < -64*np.finfo(float).eps*len(delta)))
        winner=r.lower if delta.sum()>0 else r.upper
        # Split by run index solely for evaluation; no parameter is fitted.
        first=delta[evidence.run_idx.to_numpy()%2==0];second=delta[evidence.run_idx.to_numpy()%2==1]
        decisions.append(dict(lower=r.lower,upper=r.upper,good=r.good,n_runs=len(delta),stable=bool(stable),correct=bool(winner==r.good),
            ambiguous=bool(ambiguity[r.lower].any() or ambiguity[r.upper].any()),
            disjoint_halves_agree=bool(len(first) and len(second) and first.sum()*second.sum()>0),previously_screened=r.previously_screened))
    pd.DataFrame(decisions).to_csv(args.root/'isotope_decisions.tsv',sep='\t',index=False)
    summary['pair_preferences']={'total':len(decisions),'stable':sum(d['stable'] for d in decisions),'stable_correct':sum(d['stable'] and d['correct'] for d in decisions), 'flagged_ambiguous':sum(d['ambiguous'] for d in decisions), 'disjoint_halves_agree':sum(d['disjoint_halves_agree'] for d in decisions)}
    (args.root/'isotope_choice_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
