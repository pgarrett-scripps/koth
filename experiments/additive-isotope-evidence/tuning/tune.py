"""Predefined tuning and untouched-file validation for additive isotope evidence."""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import compare as pilot
PROTOCOL = json.loads((HERE / 'protocol.json').read_text())
ART = pilot.ARTIFACTS / 'tuning'
REPO = HERE.parents[2]
TARGET = Path('/home/patrick-garrett/Repos/koth_rust/target/release')
COLS = ['mz', 'charge', 'rtStart', 'rtEnd', 'rtApex', 'intensitySum', 'nIsotopes']


def files(split):
    return [(level, f'B03_{slot:02d}_150304_human_ecoli_{level}_3ul_3um_column_95_HCD_OT_2hrs_30B_9B')
            for level, slot in PROTOCOL[split]]


def grid():
    return [(f's{s:g}_c{c:g}', s, c) for s in PROTOCOL['ratio_sigma'] for c in PROTOCOL['cosine_shape']]


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def prepare():
    (ART / 'configs').mkdir(parents=True, exist_ok=True)
    (ART / 'bin').mkdir(exist_ok=True)
    base = (pilot.PAPER / 'config/koth_ff.toml').read_text()
    for label, sigma, shape in grid():
        (ART / 'configs' / f'{label}.toml').write_text(base + f'\n[features]\nisotope_evidence_ratio_sigma = {sigma}\nisotope_evidence_cosine_shape = {shape}\n')
    for name, source in [('replay', TARGET / 'examples/replay_features'), ('candidate', TARGET / 'koth_ff'),
                         ('baseline', pilot.ARTIFACTS / 'bin/baseline')]:
        dest = ART / 'bin' / name
        if dest.exists():
            raise RuntimeError(f'Preserved binary already exists: {dest}')
        shutil.copy2(source, dest)
    provenance = {'created_utc': datetime.now(timezone.utc).isoformat(),
                  'protocol_sha256': sha(HERE / 'protocol.json'),
                  'binaries': {p.name: sha(p) for p in (ART/'bin').iterdir()},
                  'source_sha256': {str(p.relative_to(REPO)):sha(p) for p in [
                      REPO/'koth_ff/src/features/assemble.rs', REPO/'koth_ff/src/config/features.rs']}}
    (ART / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


def replay_training():
    for _, stem in files('training'):
        command = [str(ART/'bin/replay'), str(pilot.ARTIFACTS/'baseline'/stem/'hills.parquet'),
                   '--output', str(ART/'training'/stem), '--threads', str(PROTOCOL['threads'])]
        for label, _, _ in grid():
            command += ['--config', str(ART/'configs'/f'{label}.toml')]
        subprocess.run(command, check=True)
        print(f'Training {stem[:6]} finished', flush=True)


def feature_path(split, label, stem):
    if split == 'training' and label == 'baseline':
        return pilot.ARTIFACTS / 'baseline' / stem / 'features.parquet'
    if split == 'training' and label == 'original_additive':
        return pilot.ARTIFACTS / 'additive' / stem / 'features.parquet'
    if split == 'heldout':
        return ART / 'heldout' / label / stem / 'features.parquet'
    return ART / 'training' / stem / label / 'features.parquet'


def diagnostics(p, f, counts):
    """Suspicious assignments / chance-match proxies, never called feature FDR."""
    pmz, prt, pz = p.mz.to_numpy(), p.rt.to_numpy(), p.charge.to_numpy()
    fmz, fz, r0, r1 = (f[c].to_numpy() for c in ['mz','charge','rtStart','rtEnd'])
    wrong = np.zeros(len(p), bool)
    lo, hi = np.searchsorted(fmz, pmz*(1-10e-6)), np.searchsorted(fmz, pmz*(1+10e-6), side='right')
    for k,(a,b) in enumerate(zip(lo,hi)):
        wrong[k] = counts[k] == 0 and np.any((fz[a:b]!=pz[k]) & (r0[a:b]<=prt[k]) & (prt[k]<=r1[a:b]))
    shifted = p.copy()
    shifted['mz'] = pmz + 17.331281/pz  # fixed non-isotopic neutral-mass displacement
    shifted_counts = pilot.match(shifted, f)[0]
    isotope = np.zeros(len(p), bool)
    for k in [-2,-1,1,2]:
        shifted['mz'] = pmz + k*1.003354835/pz
        isotope |= pilot.match(shifted, f)[0] > 0
    return wrong, isotope & (counts==0), shifted_counts > 0


def evaluate(split, labels):
    pilot.FILES = files(split)
    anchors = pilot.anchors()  # only this split reaches the scoring code
    rows, details, matrices = [], {}, {}
    for label in labels:
        pieces = []
        for level, stem in files(split):
            f = pq.read_table(feature_path(split,label,stem), columns=COLS).to_pandas()
            f = f[f.charge.between(2,6)].sort_values('mz').reset_index(drop=True)
            p = anchors[anchors.stem == stem].copy()
            counts, apex, intensity = pilot.match(p,f)
            wrong, isotope, shifted = diagnostics(p,f,counts)
            p['count'], p['apex_count'], p['intensity'], p['level'] = counts, apex, intensity, level
            p['wrong_charge_only'], p['isotope_shift_only'], p['shifted_mass_match'] = wrong,isotope,shifted
            pieces.append(p)
            rows.append(dict(label=label, stem=stem, anchors=len(p), features=len(f),
                matched=int((counts>0).sum()), recall=float((counts>0).mean()),
                mean_isotopes=float(f.nIsotopes.mean()),
                wrong_charge_only_rate=float(wrong.mean()), isotope_shift_only_rate=float(isotope.mean()),
                shifted_mass_match_rate=float(shifted.mean())))
        d = pd.concat(pieces).sort_index()
        details[label] = d
        dest=ART/split/'matched'; dest.mkdir(parents=True,exist_ok=True)
        d.to_parquet(dest/f'{label}.parquet',index=False)
        matrices[label] = d.pivot(index=['peptide','charge','species'],columns='stem',values='intensity').replace(0,np.nan).reindex(columns=[s for _,s in files(split)])
        print(f'Evaluated {split} {label}: recall {(d["count"]>0).mean():.6f}',flush=True)
    # One fixed complete-pair cohort shared by all variants being compared.
    common = matrices[labels[0]].dropna().index
    for label in labels[1:]: common = common.intersection(matrices[label].dropna().index)
    per_file=pd.DataFrame(rows)
    summary={'split':split,'common_complete_pairs':len(common),'variants':{}}
    baseline=details['baseline']
    for label,d in details.items():
        m=matrices[label].loc[common]
        a=m[[s for l,s in files(split) if l=='A']]
        e=m[[s for l,s in files(split) if l=='E']]
        ratios=np.log2(a.median(axis=1)/e.median(axis=1))
        q={}
        for sp,expected in [('HUMAN',0),('ECOLI',np.log2(1/3))]:
            take=m.index.get_level_values('species')==sp
            errors=ratios[take]-expected
            cvs=pd.concat([a[take].std(axis=1)/a[take].mean(axis=1),e[take].std(axis=1)/e[take].mean(axis=1)])
            q[sp]=dict(n=int(take.sum()),median_absolute_log2_error=float(errors.abs().median()),
                       median_cv=float(cvs.median()),fraction_abs_log2_error_gt_0_5=float((errors.abs()>.5).mean()))
        summary['variants'][label]=dict(features=int(per_file[per_file.label==label].features.sum()),
            anchors=len(d),matched=int((d['count']>0).sum()),recall=float((d['count']>0).mean()),
            apex_recall=float((d.apex_count>0).mean()), features_per_matched_anchor=float(d.loc[d['count']>0,'count'].mean()),
            gained=int(((baseline['count']==0)&(d['count']>0)).sum()), lost=int(((baseline['count']>0)&(d['count']==0)).sum()),
            wrong_charge_only_rate=float(d.wrong_charge_only.mean()),isotope_shift_only_rate=float(d.isotope_shift_only.mean()),
            shifted_mass_match_rate=float(d.shifted_mass_match.mean()),quant=q)
    per_file.to_csv(HERE/f'{split}_per_file.csv',index=False)
    (HERE/f'{split}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary,per_file


def guardrails(summary,per_file,label):
    g=PROTOCOL['guardrails_vs_baseline']; b=summary['variants']['baseline']; v=summary['variants'][label]
    failures=[]
    bp=per_file[per_file.label=='baseline'].set_index('stem')
    vp=per_file[per_file.label==label].set_index('stem')
    if (vp.recall-bp.recall < g['minimum_per_file_recall_delta']).any(): failures.append('per-file recall')
    if v['features']/b['features']-1>g['maximum_feature_count_relative_increase']: failures.append('feature count')
    for field in ['shifted_mass_match_rate','wrong_charge_only_rate','isotope_shift_only_rate']:
        if v[field]-b[field]>g[f'maximum_{field}_increase']: failures.append(field)
    for sp in ['HUMAN','ECOLI']:
        if v['quant'][sp]['n']==0: failures.append(f'{sp} no quant pairs')
        for key in ['median_absolute_log2_error','median_cv','fraction_abs_log2_error_gt_0_5']:
            if v['quant'][sp][key]-b['quant'][sp][key]>g[f'maximum_species_{key}_increase']: failures.append(f'{sp} {key}')
    return failures


def select():
    if (HERE/'selection.json').exists(): raise RuntimeError('Selection already frozen')
    summary,pf=evaluate('training',['baseline']+[l for l,_,_ in grid()])
    eligible=[]; rows=[]
    for label,sigma,shape in grid():
        fails=guardrails(summary,pf,label)
        rows.append(dict(label=label,sigma=sigma,shape=shape,recall=summary['variants'][label]['recall'],failures=fails))
        if not fails: eligible.append(rows[-1])
    (HERE/'training_selection_audit.json').write_text(json.dumps(rows,indent=2)+'\n')
    if not eligible: raise RuntimeError('No candidate passed the predefined guardrails; do not select on held-out data')
    best=max(r['recall'] for r in eligible)
    near=[r for r in eligible if r['recall']>=best-.0005]
    rs,rc=PROTOCOL['reference']
    chosen=min(near,key=lambda r: ((np.log(r['sigma']/rs)**2+np.log(r['shape']/rc)**2),-r['recall'],r['label']))
    frozen=dict(chosen,created_utc=datetime.now(timezone.utc).isoformat(),
                protocol_sha256=sha(HERE/'protocol.json'),training_summary_sha256=sha(HERE/'training_summary.json'),
                config_sha256=sha(ART/'configs'/f'{chosen["label"]}.toml'))
    (HERE/'selection.json').write_text(json.dumps(frozen,indent=2)+'\n')
    print('FROZEN',json.dumps(frozen),flush=True)


def heldout():
    selection=json.loads((HERE/'selection.json').read_text())
    cfg=ART/'configs'/f'{selection["label"]}.toml'
    assert sha(cfg)==selection['config_sha256']
    manifest=[]
    for _,stem in files('heldout'):
        for label,binary,config in [('baseline',ART/'bin/baseline',pilot.PAPER/'config/koth_ff.toml'),
                                    ('selected',ART/'bin/candidate',cfg)]:
            out=ART/'heldout'/label
            out.mkdir(parents=True,exist_ok=True)
            if (out/stem).exists(): raise RuntimeError('Refusing to overwrite heldout result')
            command=[str(binary),str(pilot.PAPER/'data/ionstar_plain'/f'{stem}.mzML'),
                     '--config',str(config),'--output',str(out),'--threads',str(PROTOCOL['threads'])]
            import time
            start=time.monotonic()
            with (out/f'{stem}.log').open('w') as log: subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            manifest.append(dict(label=label,stem=stem,wall_s=time.monotonic()-start,command=command,binary_sha256=sha(binary)))
            (ART/'heldout'/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            print('Heldout',stem[:6],label,f'{manifest[-1]["wall_s"]:.1f}s',flush=True)
    summary,pf=evaluate('heldout',['baseline','selected'])
    failures=guardrails(summary,pf,'selected')
    result={'guardrail_failures':failures,'pass':not failures,'selection_sha256':sha(HERE/'selection.json')}
    (HERE/'heldout_decision.json').write_text(json.dumps(result,indent=2)+'\n')
    print('HELDOUT',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','training','select','heldout'])
    args=parser.parse_args()
    {'prepare':prepare,'training':replay_training,'select':select,'heldout':heldout}[args.action]()
