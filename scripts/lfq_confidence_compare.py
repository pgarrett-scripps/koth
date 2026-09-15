"""Compare the isolated group-confidence trials to a cached LFQ baseline.

Uses the paper's existing PSM matching and CV code. Paths are explicit; outputs
are written only inside the experiment directory. See docs/lfq-group-confidence.md.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper-scripts", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--baseline-metrics", type=Path, required=True)
    parser.add_argument("--sage-psms", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.paper_scripts.resolve()))
    import bruker_koth_align as ba
    import bruker_validation as bv

    root = args.experiment
    candidates = pd.read_csv(root / "audit/candidates.tsv", sep="\t")
    summary = {"group_gate_sweep": [], "trials": {}}
    for gate in [0.01, 0.025, 0.05, 0.1, 0.2, 1.0]:
        rows = candidates[(candidates.n_contributing_runs >= 2) & (candidates.group_qvalue <= gate)]
        summary["group_gate_sweep"].append(dict(gate=gate, groups=len(rows),
            seeds_below_075=int((rows.combined_score < 0.75).sum()),
            seeds_below_05=int((rows.combined_score < 0.5).sum()),
            groups_with_weak_members=int((rows.weak_members > 0).sum()),
            weak_observations=int(rows.weak_members.sum())))

    psms = bv.load_sage_psms(args.sage_psms)
    if len(psms) and psms.rt.max() > 60:
        psms = psms.assign(rt=psms.rt / 60.0)
    anchors = psms.groupby(["peptide", "charge"], sort=False).agg(
        precursor_mz=("precursor_mz", "median"), im=("im", "median"),
        rt=("rt", "median"), species=("species", "first")).reset_index()
    index = pd.MultiIndex.from_frame(anchors[["peptide", "charge"]])
    frames = {}
    baseline_report = json.loads((args.baseline / "align_report.json").read_text())
    paths = {"baseline": (args.baseline, args.baseline_metrics),
             "q005": (root / "align", root / "metrics.json"),
             "q020": (root / "q020/align", root / "q020/metrics.json")}
    for name, (folder, metric_path) in paths.items():
        intensity, qvalue, run_cols = ba.load_consensus(folder)
        assert qvalue is not None
        row = ba.match_anchors(anchors, intensity, 10.0, 0.5)
        found = row >= 0
        signal = np.full((len(anchors), len(run_cols)), np.nan)
        values = intensity[run_cols].to_numpy(float)[row[found]]
        q = qvalue[run_cols].to_numpy(float)[row[found]]
        signal[found] = np.where((values > 0) & np.isfinite(q) & (q <= 0.05), values, np.nan)
        frames[name] = pd.DataFrame(signal, index=index,
            columns=[c.removesuffix(".d") for c in run_cols])
        m = json.loads(metric_path.read_text())["koth_align"]
        trial = {k: m[k] for k in ["n_anchors", "median_cv", "human_iqr", "ffcr_human", "mv_rate", "recall"]}
        trial["mean_abs_species_bias_log2"] = float(np.mean([
            abs(v["bias"]) for k, v in m["ratios"].items() if not k.endswith("|HUMAN")]))
        trial["consensus_groups"] = len(intensity)
        trial["matched_anchors_before_cell_gate"] = int(found.sum())
        report = json.loads((folder / "align_report.json").read_text())
        trial["same_alignment_as_baseline"] = (
            report["reference_run"] == baseline_report["reference_run"]
            and all(a["name"] == b["name"] and a.get("alignment") == b.get("alignment")
                    for a, b in zip(report["runs"], baseline_report["runs"], strict=True)))
        if name != "baseline":
            mapped_scores = intensity.combined_score.to_numpy()[row[found]]
            trial["anchors_assigned_to_seed_below_075"] = int((mapped_scores < 0.75).sum())
            trial["anchors_assigned_to_seed_below_05"] = int((mapped_scores < 0.5).sum())
            base_cells = frames["baseline"].notna().to_numpy()
            cells = np.isfinite(signal)
            trial["gained_anchor_run_cells"] = int((cells & ~base_cells).sum())
            trial["lost_anchor_run_cells"] = int((base_cells & ~cells).sum())
            trial["gained_anchors_with_any_signal"] = int((cells.any(axis=1) & ~base_cells.any(axis=1)).sum())
            trial["lost_anchors_with_any_signal"] = int((base_cells.any(axis=1) & ~cells.any(axis=1)).sum())
        summary["trials"][name] = trial

    # Same peptide/run cells in all trials, before normalization: selection
    # alone cannot make this CV comparison look better by dropping difficult rows.
    shared = np.logical_and.reduce([frame.notna().to_numpy() for frame in frames.values()])
    species_ok = anchors.species.isin(bv.SPECIES).to_numpy()
    cv = {name: bv.per_peptide_cv_frame(frame.where(shared).loc[species_ok])
          .set_index(["peptide", "charge", "cond"])["cv"] for name, frame in frames.items()}
    paired = pd.concat(cv, axis=1).dropna()
    paired.to_parquet(root / "paired_cv.parquet")
    summary["common_cell_control"] = {
        "shared_cells": int(shared[species_ok].sum()), "paired_cv_observations": len(paired),
        "median_cv_pct": {name: float(paired[name].median() * 100) for name in frames},
        "median_paired_cv_change_pp": {name: float((paired[name] - paired.baseline).median() * 100) for name in frames}}
    summary["inputs"] = {"baseline": str(args.baseline.resolve()), "sage_psms": str(args.sage_psms.resolve()),
        "baseline_metrics_sha256": hashlib.sha256(args.baseline_metrics.read_bytes()).hexdigest(),
        "paper_code_sha256": {name: hashlib.sha256((args.paper_scripts / name).read_bytes()).hexdigest()
                              for name in ["bruker_koth_align.py", "bruker_validation.py"]}}
    (root / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
