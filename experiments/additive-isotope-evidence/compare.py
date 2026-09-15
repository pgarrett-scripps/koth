"""Bounded four-file IonStar A/B experiment; never changes paper outputs.

Run each binary with `run LABEL BINARY`, then use `score`. Requires numpy,
pandas and pyarrow. See README.md for exact commands and model assumptions.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ARTIFACTS = HERE / 'artifacts'
PAPER = Path('/home/patrick-garrett/Repos/koth-paper/analysis')
FILES = [(level, f'B03_{slot:02d}_150304_human_ecoli_{level}_3ul_3um_column_95_HCD_OT_2hrs_30B_9B')
         for level, slot in [('A', 10), ('A', 11), ('E', 5), ('E', 6)]]
PROTON = 1.007276466621


def run(label: str, binary: Path) -> None:
    out = ARTIFACTS / label
    out.mkdir(parents=True, exist_ok=True)
    config = out / 'config.toml'
    config.write_text((PAPER / 'config/koth_ff.toml').read_text())
    manifest = dict(label=label, binary=str(binary.resolve()),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                    threads=8, runs=[])
    for level, stem in FILES:
        source = PAPER / 'data/ionstar_plain' / f'{stem}.mzML'
        if (out / stem / 'features.parquet').exists():
            raise RuntimeError(f'Refusing to overwrite previous output for {stem}')
        command = [str(binary.resolve()), str(source), '--config', str(config),
                   '--output', str(out), '--threads', '8']
        start = time.monotonic()
        with (out / f'{stem}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        elapsed = time.monotonic() - start
        manifest['runs'].append(dict(level=level, stem=stem, input=str(source),
                                     size=source.stat().st_size, wall_s=elapsed, command=command))
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        print(f'{label} {stem[:6]} complete: {elapsed:.1f}s', flush=True)


def anchors() -> pd.DataFrame:
    p = pq.read_table(PAPER / 'data/sage_perz/results.sage.parquet', columns=[
        'filename', 'peptide', 'charge', 'expmass', 'rt', 'peptide_q', 'is_decoy', 'proteins']).to_pandas()
    p['stem'] = p.filename.str.removesuffix('.mzML.gz')
    p = p[(p.peptide_q <= .01) & (~p.is_decoy) & p.stem.isin([s for _, s in FILES])].copy()
    p = p.sort_values('peptide_q', kind='stable').drop_duplicates(['stem', 'peptide', 'charge'])
    p['mz'] = (p.expmass + p.charge * PROTON) / p.charge
    human = p.proteins.str.contains('_HUMAN', regex=False)
    ecoli = p.proteins.str.contains('_ECOLI', regex=False)
    p['species'] = np.where(human & ~ecoli, 'HUMAN', np.where(ecoli & ~human, 'ECOLI', 'OTHER'))
    return p.reset_index(drop=True)


def match(p: pd.DataFrame, f: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mz = f.mz.to_numpy()
    z = f.charge.to_numpy()
    r0, r1, ra = (f[c].to_numpy() for c in ['rtStart', 'rtEnd', 'rtApex'])
    intensity = f.intensitySum.to_numpy()
    counts = np.zeros(len(p), dtype=int)
    apex_counts = counts.copy()
    intensities = np.zeros(len(p))
    lo = np.searchsorted(mz, p.mz.to_numpy() * (1 - 10e-6))
    hi = np.searchsorted(mz, p.mz.to_numpy() * (1 + 10e-6), side='right')
    for k, (a, b, charge, rt) in enumerate(zip(lo, hi, p.charge, p.rt)):
        ok = (z[a:b] == charge) & (r0[a:b] <= rt) & (rt <= r1[a:b])
        counts[k] = ok.sum()
        apex_counts[k] = ((z[a:b] == charge) & (np.abs(ra[a:b] - rt) <= .5)).sum()
        if ok.any():
            intensities[k] = intensity[a:b][ok].max()
    return counts, apex_counts, intensities


def score() -> None:
    p = anchors()
    details, per_file = {}, []
    for label in ['baseline', 'additive']:
        rows = []
        for level, stem in FILES:
            f = pq.read_table(ARTIFACTS / label / stem / 'features.parquet', columns=[
                'mz', 'charge', 'rtStart', 'rtEnd', 'rtApex', 'intensitySum', 'nIsotopes']).to_pandas()
            f = f[f.charge.between(2, 6)].sort_values('mz').reset_index(drop=True)
            sub = p[p.stem == stem].copy()
            counts, apex_counts, intensity = match(sub, f)
            sub['level'], sub['count'], sub['intensity'] = level, counts, intensity
            sub['apex_count'] = apex_counts
            rows.append(sub)
            per_file.append(dict(label=label, stem=stem, features=len(f), mean_isotopes=float(f.nIsotopes.mean()), anchors=len(sub),
                                 matched=int((counts > 0).sum()), recall=float((counts > 0).mean()),
                                 apex_recall=float((apex_counts > 0).mean()),
                                 features_per_matched_anchor=float(counts[counts > 0].mean())))
        details[label] = pd.concat(rows).sort_index()
        details[label].to_parquet(ARTIFACTS / f'{label}_matched.parquet', index=False)
    per_file = pd.DataFrame(per_file)
    per_file.to_csv(HERE / 'per_file.csv', index=False)
    # Quantification is compared on the exact same peptide/charge pairs,
    # observed in all four files by both scorers: avoids a changing cohort.
    matrices = {label: d.pivot(index=['peptide', 'charge', 'species'], columns='stem', values='intensity')
                .replace(0, np.nan).reindex(columns=[s for _, s in FILES]) for label, d in details.items()}
    common = matrices['baseline'].dropna().index.intersection(matrices['additive'].dropna().index)
    summary = {'anchor_definition': '1% peptide-q target Sage IDs, deduplicated per run/peptide/charge',
               'matching': 'same charge, 10 ppm, PSM RT inside feature interval; apex +/-0.5 min sensitivity',
               'files': [s for _, s in FILES], 'common_complete_quant_pairs': len(common), 'variants': {}}
    for label, d in details.items():
        stats = per_file[per_file.label == label]
        m = matrices[label].loc[common]
        a = m[[s for l, s in FILES if l == 'A']]
        e = m[[s for l, s in FILES if l == 'E']]
        ratio = np.log2(a.median(axis=1) / e.median(axis=1))
        quant = {}
        for sp, expected in [('HUMAN', 0.0), ('ECOLI', np.log2(1/3))]:
            take = m.index.get_level_values('species') == sp
            err = ratio[take] - expected
            # Two-replicate within-level CV is descriptive, not a precision estimate.
            cvs = pd.concat([a[take].std(axis=1)/a[take].mean(axis=1),
                             e[take].std(axis=1)/e[take].mean(axis=1)])
            quant[sp] = dict(n=int(take.sum()), median_abs_log2_error=float(err.abs().median()),
                             median_log2_bias=float(err.median()), median_within_level_cv=float(cvs.median()))
        summary['variants'][label] = dict(features=int(stats.features.sum()), anchors=len(d),
             matched=int((d['count'] > 0).sum()), recall=float((d['count'] > 0).mean()),
             apex_recall=float((d.apex_count > 0).mean()),
             features_per_matched_anchor=float(d.loc[d['count'] > 0, 'count'].mean()), quant=quant)
    b, a = details['baseline'], details['additive']
    summary['paired_recall'] = dict(gained=int(((b['count'] == 0) & (a['count'] > 0)).sum()),
                                    lost=int(((b['count'] > 0) & (a['count'] == 0)).sum()))
    (HERE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    run_args = sub.add_parser('run')
    run_args.add_argument('label', choices=['baseline', 'additive'])
    run_args.add_argument('binary', type=Path)
    sub.add_parser('score')
    args = parser.parse_args()
    if args.action == 'run':
        run(args.label, args.binary)
    else:
        score()
