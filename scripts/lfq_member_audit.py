"""PSM-observed weak-member discordance, without treating missing MS2 as absence."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

VARIANTS = ["average", "average_no_penalty", "tree", "tree_no_penalty"]


def wilson(k, n):
    if not n:
        return [None, None]
    z = 1.96
    p = k / n
    d = 1 + z*z/n
    mid = (p + z*z/(2*n))/d
    half = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return [float(mid-half), float(mid+half)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--psms", type=Path, required=True)
    ap.add_argument("--platform", choices=["bruker", "orbitrap"], required=True)
    ap.add_argument("--paper-scripts", type=Path, required=True)
    args = ap.parse_args()
    sys.path.insert(0, str(args.paper_scripts))
    if args.platform == "bruker":
        import bruker_validation as bv
        psms = bv.load_sage_psms(args.psms)
        if psms.rt.max() > 60:
            psms["rt"] /= 60
    else:
        psms = pd.read_parquet(args.psms, columns=["filename", "peptide", "charge", "expmass", "rt", "peptide_q", "is_decoy"])
        psms = psms[(psms.peptide_q <= .01) & ~psms.is_decoy].copy()
        psms["stem"] = psms.filename.str.removesuffix(".gz")
        psms["precursor_mz"] = psms.expmass / psms.charge + 1.007276466621
        psms = psms.groupby(["stem", "peptide", "charge"], sort=False).agg(precursor_mz=("precursor_mz", "median"), rt=("rt", "median")).reset_index()
    # I/L are indistinguishable by MS1. Other ambiguous candidate identities are excluded.
    psms["label"] = pd.factorize(psms.peptide.str.replace("I", "L", regex=False) + "/" + psms.charge.astype(str))[0]
    runs = json.loads((args.audit / "runs.json").read_text())
    labels = []
    for r in runs:
        f = pd.read_parquet(r["features"], columns=["mz", "charge", "rtApex", "rtStart", "rtEnd", "im", "combined_score"])
        f = f[f.charge != 0].reset_index(drop=True)
        ids = np.argsort(f.mz.to_numpy())
        mz = f.mz.to_numpy()[ids]
        charge = f.charge.to_numpy()
        rt, lo_rt, hi_rt = (f[c].to_numpy() for c in ["rtApex", "rtStart", "rtEnd"])
        im = f.im.to_numpy()
        lab = np.full(len(f), -1, dtype=np.int64)
        name = r["name"].removesuffix(".d").removesuffix(".gz")
        ps = psms[psms.stem == name]
        low = np.searchsorted(mz, ps.precursor_mz.to_numpy()*(1-10e-6))
        high = np.searchsorted(mz, ps.precursor_mz.to_numpy()*(1+10e-6), side="right")
        for j, p in enumerate(ps.itertuples()):
            candidates = ids[low[j]:high[j]]
            ok = (charge[candidates] == p.charge) & (np.abs(rt[candidates]-p.rt) <= .1)
            ok &= (lo_rt[candidates]-.02 <= p.rt) & (hi_rt[candidates]+.02 >= p.rt)
            if args.platform == "bruker":
                ok &= np.abs(im[candidates]-p.im) <= .015
            hit = candidates[ok]
            current = lab[hit]
            lab[hit] = np.where((current == -1) | (current == p.label), p.label, -2)
        labels.append(pd.DataFrame(dict(run_id=r["run_id"], feature_idx=np.arange(len(f)), label=lab, quality=f.combined_score)))
        print(f"{name}: {sum(lab>=0):,}/{len(f):,} unambiguous feature labels", flush=True)
    observations = pd.concat(labels, ignore_index=True)
    observations.to_parquet(args.audit / "feature_labels.parquet", index=False)
    members = pd.read_csv(args.audit / "members.tsv", sep="\t")
    m = members.merge(observations, on=["run_id", "feature_idx"], validate="many_to_one")
    candidates = pd.read_csv(args.audit / "candidates.tsv", sep="\t")
    donors = {}
    for p in m[(m.is_primary) & (m.label >= 0) & (m.quality >= .75)].itertuples():
        donors.setdefault(p.candidate_id, {}).setdefault(p.label, set()).add(p.run_id)
    assessed, wrong = np.zeros(len(m), bool), np.zeros(len(m), bool)
    for p in m[m.label >= 0].itertuples():
        other = [label for label, run_ids in donors.get(p.candidate_id, {}).items() if run_ids - {p.run_id}]
        if len(other) == 1:
            assessed[p.Index] = True
            wrong[p.Index] = p.label != other[0]
    m["assessed"], m["discordant"] = assessed, wrong
    m.to_parquet(args.audit / "member_labels.parquet", index=False)
    result = {"method": "Leave-one-run-out unanimous strong-primary donor labels; PSM-observed discordance only, not FDR", "variants": {}}
    n_controls = len(json.loads((args.audit / "audit.json").read_text())["controls"])
    for variant in VARIANTS:
        keep = (candidates.n_contributing_runs >= 2) & (candidates[f"{variant}_global"] <= .05)
        selected = m[m.candidate_id.isin(candidates.loc[keep, "candidate_id"])]
        groups = candidates[keep]
        stats = dict(groups=len(groups), weak_observations=int((selected.quality < .5).sum()), weak_seed_groups=int((groups.combined_score < .75).sum()))
        for name, rows in [("weak", selected[selected.quality < .5]), ("other", selected[selected.quality >= .5]),
                           ("weak_primary", selected[(selected.quality < .5) & selected.is_primary])]:
            n, k = int(rows.assessed.sum()), int(rows.discordant.sum())
            stats[name] = dict(observations=len(rows), labelled=int((rows.label >= 0).sum()), assessed=n, discordant=k,
                               fraction=k/n if n else None, wilson95=wilson(k, n))
        stability = {}
        control_names = (["repeat", "density", "pooled"] if n_controls >= 15 else []) + (["conditioned"] if n_controls >= 25 else [])
        for control in control_names:
            alternate = (candidates.n_contributing_runs >= 2) & (candidates[f"{variant}_{control}"] <= .05)
            stability[control] = dict(groups=int(alternate.sum()), retained_fraction=float((keep & alternate).sum()/max(keep.sum(), 1)),
                                      jaccard=float((keep & alternate).sum()/max((keep | alternate).sum(), 1)))
        stats["controls"] = stability
        result["variants"][variant] = stats
    (args.audit / "member_audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
