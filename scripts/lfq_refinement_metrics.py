"""Evaluate fixed LFQ variants with the existing paper metrics and shared-cell CV."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--platform",choices=["bruker","orbitrap"],required=True)
    ap.add_argument("--psms",type=Path,required=True)
    ap.add_argument("--paper-scripts",type=Path,required=True)
    ap.add_argument("--files",type=Path)
    ap.add_argument("--variants",default="average,average_no_penalty,tree,tree_no_penalty")
    ap.add_argument("--baseline",type=Path)
    ap.add_argument("--control",choices=["global","conditioned","permuted"],default="global")
    args=ap.parse_args()
    sys.path.insert(0,str(args.paper_scripts))
    import bruker_koth_align as ba
    import bruker_validation as bv
    orbit=importlib.import_module("18_koth_align_lfq")
    if args.platform=="bruker":
        psms=bv.load_sage_psms(args.psms)
        if psms.rt.max()>60:psms["rt"]/=60
        anchors=psms.groupby(["peptide","charge"],sort=False).agg(precursor_mz=("precursor_mz","median"),im=("im","median"),rt=("rt","median"),species=("species","first")).reset_index()
    else:
        anchors=orbit.load_anchors(args.psms)
        psms=orbit.load_psms_per_run(args.psms)
        files=orbit.parse_files_tsv(args.files)
    folders={v:args.root/v for v in args.variants.split(',')}
    if args.baseline:folders={"legacy":args.baseline,**folders}
    frames={};summary={};indices={}
    index=pd.MultiIndex.from_frame(anchors[["peptide","charge"]])
    species=pd.Series(anchors.species.to_numpy(),index=index)
    for variant,folder in folders.items():
        intensity,q,cols=ba.load_consensus(folder)
        q[cols]=q[cols].round(4)  # Match production TSV precision in the first audit build.
        cons=pd.read_csv(folder/'consensus_features.tsv',sep='\t')
        if args.platform=="bruker":
            row=ba.match_anchors(anchors,intensity,10.0,.5)
        else:
            row=orbit.match_anchors_to_consensus(anchors,cons,10.0,.5).consensus_row.to_numpy()
        found=row>=0
        values=np.full((len(anchors),len(cols)),np.nan)
        v=intensity[cols].to_numpy(float)[row[found]]
        quality=q[cols].to_numpy(float)[row[found]]
        values[found]=np.where((v>0)&np.isfinite(quality)&(quality<=.05),v,np.nan)
        names=[c.removesuffix('.d') for c in cols]
        frame=pd.DataFrame(values,index=index,columns=names)
        frames[variant]=frame;indices[variant]=row
        if args.platform=='bruker':
            m=bv.matrix_metrics(frame.loc[found & anchors.species.isin(bv.SPECIES).to_numpy()],species.loc[found & anchors.species.isin(bv.SPECIES).to_numpy()])
            pi=index.get_indexer(pd.MultiIndex.from_frame(psms[['peptide','charge']]))
            ri=pd.Index(names).get_indexer(psms.stem)
            valid=(pi>=0)&(ri>=0)
            covered=int(np.isfinite(values[pi[valid],ri[valid]]).sum())
            m['recall']={'psm_total':len(psms),'psm_covered':covered,'psm_recall':covered/len(psms)}
        else:
            mat=frame.loc[found].fillna(0).reset_index()
            sp=species.loc[found].reset_index(drop=True)
            normalized=orbit.median_normalize_wide(mat,names)
            ratios=orbit.per_pair_log2_wide(normalized,sp,files,names)
            cv=orbit.per_level_cv_wide(normalized,sp,files,names)
            human=ratios.loc[ratios.species=='HUMAN','observed_log2']
            ecoli=ratios.loc[ratios.species=='ECOLI']
            recall,_=orbit.psm_recall_metrics(psms,cons,intensity,q,cols,10.0,.5,.05)
            m=dict(n_anchors=int(found.sum()),median_cv=float(cv.cv.median()),human_iqr=float(human.quantile(.75)-human.quantile(.25)),
                ffcr_human=float((human.abs()>np.log2(1.5)).mean()),mv_rate=float(frame.loc[found].isna().to_numpy().mean()),
                species_bias=float((ecoli.observed_log2-ecoli.expected_log2).median()),recall=recall)
        m['consensus_groups']=len(cons)
        summary[variant]=m
        if variant!='legacy':(folder/'metrics.json').write_text(json.dumps({'koth_align':m},indent=2)+'\n')
        print(variant, 'groups',len(cons),'recall',m['recall']['psm_recall'],'CV',m['median_cv'],'FFCR',m['ffcr_human'],flush=True)
    shared=np.logical_and.reduce([f.notna().to_numpy() for f in frames.values()])
    cv_frames={}
    for variant,frame in frames.items():
        f=frame.where(shared)
        if args.platform=='bruker':
            cv=bv.per_peptide_cv_frame(f.loc[anchors.species.isin(bv.SPECIES).to_numpy()]).set_index(['peptide','charge','cond']).cv
        else:
            mat=orbit.median_normalize_wide(f.fillna(0).reset_index(),list(f.columns))
            cv=orbit.per_level_cv_wide(mat,anchors.species,files,list(f.columns)).set_index(['peptide','charge','level']).cv
        cv_frames[variant]=cv
    paired=pd.concat(cv_frames,axis=1).dropna()
    paired.to_parquet(args.root/'paired_cv.parquet')
    control={'shared_cells':int(shared.sum()),'paired_CVs':len(paired),'median_CV':{v:float(paired[v].median()) for v in frames}}
    result={'variants':summary,'common_cells':control}
    (args.root/'metrics_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('COMMON',json.dumps(control))
    # Post-cell-gate link audit: only primary weak observations linked to a
    # reported extraction. It remains an MS2-observed discordance proxy.
    member_path=args.root/'member_labels.parquet'
    if member_path.exists():
        members=pd.read_parquet(member_path)
        candidates=pd.read_csv(args.root/'candidates.tsv',sep='\t')
        from lfq_member_audit import wilson
        links={}
        for variant in folders:
            if variant=='legacy':continue
            keep=(candidates.n_contributing_runs>=2)&(candidates[variant+'_'+args.control]<=.05)
            mapping=np.full(len(candidates),-1,dtype=int);mapping[np.flatnonzero(keep)]=np.arange(keep.sum())
            intensity,q,cols=ba.load_consensus(folders[variant])
            q[cols]=q[cols].round(4)
            assert len(intensity)==keep.sum()
            subset=members[(members.quality<.5)&members.is_primary&members.assessed].copy()
            row=mapping[subset.candidate_id.to_numpy()];valid=row>=0
            row=row[valid];run=subset.run_id.to_numpy()[valid]
            accepted=(intensity[cols].to_numpy()[row,run]>0)&(q[cols].to_numpy()[row,run]<=.05)
            k=int(subset.loc[valid,'discordant'].to_numpy()[accepted].sum());n=int(accepted.sum())
            links[variant]={'assessed':n,'discordant':k,'fraction':k/n if n else None,'wilson95':wilson(k,n)}
        (args.root/'reported_weak_links.json').write_text(json.dumps(links,indent=2)+'\n')
        print('REPORTED_WEAK_LINKS',json.dumps(links))


if __name__=='__main__':main()
