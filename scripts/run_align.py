"""
Align feature tables from a batch run and produce a quantification matrix.

Usage:
    uv run python scripts/run_align.py <batch_output_dir> [aligned_output_dir]
                                       [--fdr] [--sweep]
                                       [--fdr-threshold FLOAT]
                                       [--mass-ppm FLOAT] [--rt-window FLOAT]
                                       [--im-tolerance FLOAT]
                                       [--rt-anchor-window FLOAT]
                                       [--min-feature-score FLOAT]

Expects <batch_output_dir> to contain one subdirectory per sample, each with
a features.parquet file (as written by run_batch.py).

Outputs (in <aligned_output_dir>/):
    consensus_features.tsv        – one row per consensus feature
    intensity_matrix.tsv          – consensus features × sample intensities
    alignment_report.html         – visual QC report
    <run>_calibration.json        – per-run alignment curves (one file per non-reference run):
                                    MZ drift (linear ppm coefficients),
                                    IM drift (linear coefficients),
                                    RT warp (dense lookup table + spline control points)

Target / decoy flags:
    --fdr     Estimate false discovery rate using a mass-shifted decoy reference.
              Adds per-run FDR estimates to the report and console output.
    --sweep   Grid search over (mass_ppm, rt_window, im_tolerance) to find the
              tightest tolerances that keep FDR below --fdr-threshold (default 1%).
              Implies --fdr. Adds heatmap plots to the report.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from koth_ff.align import align_runs
from alignment_report import generate_report


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("batch_dir", help="Directory containing per-sample subdirectories")
    p.add_argument("out_dir", nargs="?", help="Output directory (default: <batch_dir>/../aligned)")
    p.add_argument("--fdr", action="store_true",
                   help="Run target/decoy FDR estimation for each run")
    p.add_argument("--sweep", action="store_true",
                   help="Grid-search tolerance parameters (implies --fdr)")
    p.add_argument("--fdr-threshold", type=float, default=0.01, metavar="FLOAT",
                   help="Maximum acceptable competition FDR (default: 0.01 = 1%%)")
    p.add_argument("--mass-shift-da", type=float, default=100.0, metavar="FLOAT",
                   help="Mass shift for decoy reference in Da (default: 100)")
    # Matching tolerances
    p.add_argument("--mass-ppm", type=float, default=10.0, metavar="FLOAT",
                   help="Mass tolerance for final feature matching in ppm (default: 10)")
    rt_group = p.add_mutually_exclusive_group()
    rt_group.add_argument("--rt-window", type=float, default=None, metavar="FLOAT",
                          help="RT tolerance for final matching in minutes (default: 0.5)")
    rt_group.add_argument("--rt-window-pct", type=float, default=None, metavar="FLOAT",
                          help="RT tolerance as %% of reference gradient length "
                               "(e.g. 1.5 means 1.5%% of the gradient). "
                               "Scales automatically across different gradient lengths.")
    p.add_argument("--im-tolerance", type=float, default=0.05, metavar="FLOAT",
                   help="Ion-mobility tolerance in 1/K0 units (default: 0.05). "
                        "Skipped automatically for features without ion mobility.")
    p.add_argument("--rt-anchor-window", type=float, default=0.05, metavar="FLOAT",
                   help="RT tolerance for anchor-pair finding, as a fraction of the "
                        "gradient length (default: 0.05 = ±5%% of gradient). "
                        "This is the only anchor-specific tolerance; tighten it if "
                        "the warp fit looks noisy, loosen it if too few anchors are found.")
    p.add_argument("--min-anchor-score", type=float, default=0.75, metavar="FLOAT",
                   help="Minimum Bhattacharyya score for a feature to be used as a "
                        "warp anchor (default: 0.75)")
    p.add_argument("--min-anchor-count", type=int, default=20, metavar="INT",
                   help="Minimum anchor pairs required to fit a warp; fewer falls back "
                        "to identity (no RT correction) with a warning (default: 20)")
    warp_group = p.add_argument_group("RT warp tuning")
    warp_group.add_argument("--rt-warp-bandwidth", type=float, default=0.05, metavar="FLOAT",
                            help="Sliding-window width (in normalised RT [0,1]) used to compute "
                                 "per-window median deltas before fitting the spline. "
                                 "Smaller = finer resolution but noisier; larger = smoother "
                                 "(default: 0.05 = 5%% of gradient)")
    warp_group.add_argument("--rt-warp-sigma-clip", type=float, default=3.0, metavar="FLOAT",
                            help="Sigma threshold for iterative outlier rejection of anchor pairs "
                                 "during warp fitting. Lower values remove more outliers "
                                 "(default: 3.0)")
    warp_group.add_argument("--rt-warp-clip-iters", type=int, default=5, metavar="INT",
                            help="Maximum iterations of sigma-clipping during warp fitting "
                                 "(default: 5)")
    warp_group.add_argument("--rt-warp-lam", type=float, default=None, metavar="FLOAT",
                            help="Smoothing spline lambda (regularisation strength). "
                                 "None = auto-select via GCV (default). "
                                 "Larger values produce a smoother/stiffer curve.")
    p.add_argument("--min-feature-score", type=float, default=0.7, metavar="FLOAT",
                   help="Minimum Bhattacharyya score for a reference feature to appear "
                        "as a row in the output matrix. Higher-scoring features are "
                        "processed first and can still match lower-scoring counterparts "
                        "in other runs (default: 0.7). Pass 0 to keep all features.")
    return p.parse_args()


REQUIRED_COLUMNS = {"rt_apex", "charge", "mz", "intensity_sum"}

# Mapping from koth_ff CLI camelCase column names to snake_case expected by align_runs
_COLUMN_RENAMES = {
    "rtApex": "rt_apex",
    "rtStart": "rt_start",
    "rtEnd": "rt_end",
    "intensityApex": "intensity_apex",
    "intensitySum": "intensity_sum",
    "massCalib": "mass",
    "nIsotopes": "n_isotopes",
    "nScans": "n_scans",
}


def _write_calibration_json(
    run_name: str,
    run_diag: dict,
    reference: str,
    out_dir: Path,
    n_rt_samples: int = 200,
) -> None:
    """Write per-run alignment calibration summary as JSON."""
    ms = run_diag["match_stats"]
    n_ref = ms["n_ref"]
    n_matched = ms["n_matched"]

    summary = {
        "n_anchors": run_diag["n_anchors"],
        "identity_fallback": run_diag.get("identity_fallback", False),
        "n_ref_features": n_ref,
        "n_run_features": ms["n_run"],
        "n_matched": n_matched,
        "match_rate_pct": round(n_matched / n_ref * 100, 2) if n_ref else 0.0,
        "after_charge_mass_rt": ms["after_charge_mass_rt"],
        "after_im": ms["after_im"],
    }

    mz_calibration = {
        "type": "linear_ppm",
        "slope_ppm_per_min": run_diag.get("mass_drift_slope_ppm_per_min", 0.0),
        "intercept_ppm": run_diag.get("mass_drift_intercept_ppm", 0.0),
        "description": (
            "ppm_correction = slope_ppm_per_min * rt_apex + intercept_ppm; "
            "corrected_mass = raw_mass / (1 + ppm_correction * 1e-6)"
        ),
    }

    im_calibration = {
        "type": "linear",
        "slope_per_min": run_diag.get("im_drift_slope_per_min", 0.0),
        "intercept": run_diag.get("im_drift_intercept", 0.0),
        "description": (
            "im_correction = slope_per_min * rt_apex + intercept; "
            "corrected_im = raw_im - im_correction"
        ),
    }

    run_rt_min = run_diag["run_rt_min"]
    run_rt_max = run_diag["run_rt_max"]
    identity_fallback = run_diag.get("identity_fallback", False)

    if identity_fallback or "warp_fn" not in run_diag:
        rt_calibration = {
            "type": "identity",
            "run_rt_min": run_rt_min,
            "run_rt_max": run_rt_max,
            "ref_rt_min": run_diag["ref_rt_min"],
            "ref_rt_max": run_diag["ref_rt_max"],
            "description": "No RT correction applied (identity warp).",
        }
    else:
        warp_fn = run_diag["warp_fn"]
        input_rt = np.linspace(run_rt_min, run_rt_max, n_rt_samples)
        warped_rt = warp_fn(input_rt)
        rt_calibration = {
            "type": "spline_lookup",
            "run_rt_min": run_rt_min,
            "run_rt_max": run_rt_max,
            "ref_rt_min": run_diag["ref_rt_min"],
            "ref_rt_max": run_diag["ref_rt_max"],
            # Dense lookup table: interpolate linearly between adjacent points
            "input_rt": [round(float(v), 4) for v in input_rt],
            "warped_rt": [round(float(v), 4) for v in warped_rt],
            # Spline control points in normalised RT space (for reference/reconstruction)
            "window_centers_norm": [round(float(v), 6) for v in run_diag.get("window_centers", [])],
            "window_medians_delta_norm": [round(float(v), 6) for v in run_diag.get("window_medians", [])],
            "description": (
                "input_rt (run space, absolute minutes) -> warped_rt (ref space, absolute minutes). "
                "Interpolate linearly between adjacent sample points. "
                "window_centers_norm and window_medians_delta_norm are the spline control points "
                "in normalised RT space [0,1] used to fit the warp."
            ),
        }

    out = {
        "run_name": run_name,
        "reference": reference,
        "summary": summary,
        "mz_calibration": mz_calibration,
        "im_calibration": im_calibration,
        "rt_calibration": rt_calibration,
    }

    out_path = out_dir / f"{run_name}_calibration.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)


def load_runs(batch_dir: Path) -> dict[str, pl.DataFrame]:
    runs: dict[str, pl.DataFrame] = {}
    for sub in sorted(batch_dir.iterdir()):
        if not sub.is_dir():
            continue
        if (pq := sub / "features.parquet").exists():
            df = pl.read_parquet(pq)
        elif (tsv := sub / "features.tsv").exists():
            df = pl.read_csv(tsv, separator="\t")
        else:
            continue

        renames = {k: v for k, v in _COLUMN_RENAMES.items() if k in df.columns}
        if renames:
            df = df.rename(renames)

        missing = REQUIRED_COLUMNS - set(df.columns)
        if missing:
            print(f"Warning: skipping run '{sub.name}': missing columns {sorted(missing)}", file=sys.stderr)
            continue

        if len(df) == 0:
            print(f"Warning: skipping run '{sub.name}': features file is empty", file=sys.stderr)
            continue

        runs[sub.name] = df
    return runs


def main(args: argparse.Namespace) -> None:
    batch_dir = Path(args.batch_dir)
    out_dir = Path(args.out_dir) if args.out_dir else batch_dir.parent / "aligned"

    if not batch_dir.is_dir():
        print(f"Error: {batch_dir} is not a directory")
        sys.exit(1)

    print(f"Loading features from {batch_dir}")
    runs = load_runs(batch_dir)

    if len(runs) < 2:
        print(f"Found {len(runs)} run(s) — need at least 2 to align.")
        sys.exit(1)

    print(f"Found {len(runs)} runs:")
    for name, df in runs.items():
        print(f"  {name}: {len(df):,} features")

    # Resolve RT window: --rt-window-pct scales by the longest run's gradient span
    if args.rt_window_pct is not None:
        def _span(df: pl.DataFrame) -> float:
            rt = df["rt_apex"].cast(pl.Float64).to_numpy()
            return float(rt.max()) - float(rt.min())
        max_span = max(_span(df) for df in runs.values())
        rt_window = max_span * args.rt_window_pct / 100.0
        print(f"\nRT window: {args.rt_window_pct}% of {max_span:.1f} min gradient = {rt_window:.3f} min")
    else:
        rt_window = args.rt_window if args.rt_window is not None else 0.5

    run_fdr = args.fdr or args.sweep

    min_feature_score_val = args.min_feature_score if args.min_feature_score > 0 else None
    print("\nSettings:")
    print(f"  mass tolerance (matching):  {args.mass_ppm} ppm")
    print(f"  RT window (matching):       {rt_window:.3f} min", end="")
    if args.rt_window_pct is not None:
        print(f"  ({args.rt_window_pct}% of gradient)", end="")
    print()
    print(f"  IM tolerance:               {args.im_tolerance} 1/K0")
    print(f"  RT anchor window:           {args.rt_anchor_window*100:.1f}% of gradient")
    print(f"  min anchor score:           {args.min_anchor_score}")
    print(f"  min anchor count:           {args.min_anchor_count}")
    lam_str = str(args.rt_warp_lam) if args.rt_warp_lam is not None else "auto (GCV)"
    print(f"  RT warp bandwidth:          {args.rt_warp_bandwidth*100:.1f}% of gradient")
    print(f"  RT warp sigma clip:         {args.rt_warp_sigma_clip}")
    print(f"  RT warp clip iters:         {args.rt_warp_clip_iters}")
    print(f"  RT warp spline lam:         {lam_str}")
    print(f"  min feature score:          {min_feature_score_val if min_feature_score_val is not None else 'none (all features)'}")
    if run_fdr:
        print(f"  FDR mass shift:             {args.mass_shift_da} Da")
        print(f"  FDR threshold:              {args.fdr_threshold*100:.1f}%")
    if args.sweep:
        print("  parameter sweep:            enabled")

    t0 = time.perf_counter()
    consensus_df, matrix_df, diagnostics = align_runs(
        runs,
        mass_ppm=args.mass_ppm,
        anchor_mass_ppm=args.mass_ppm,
        rt_window=rt_window,
        rt_anchor_window=args.rt_anchor_window,
        im_tolerance=args.im_tolerance,
        min_anchor_score=args.min_anchor_score,
        min_feature_score=args.min_feature_score if args.min_feature_score > 0 else None,
        min_anchor_count=args.min_anchor_count,
        rt_warp_bandwidth=args.rt_warp_bandwidth,
        rt_warp_sigma_clip=args.rt_warp_sigma_clip,
        rt_warp_clip_iters=args.rt_warp_clip_iters,
        rt_warp_spline_lam=args.rt_warp_lam,
        intensity_col="intensity_sum",
        return_diagnostics=True,
        run_fdr_estimation=run_fdr,
        run_parameter_sweep=args.sweep,
        fdr_threshold=args.fdr_threshold,
        mass_shift_da=args.mass_shift_da,
    )
    elapsed = time.perf_counter() - t0

    out_dir.mkdir(parents=True, exist_ok=True)

    consensus_df.write_csv(out_dir / "consensus_features.tsv", separator="\t")
    matrix_df.write_csv(out_dir / "intensity_matrix.tsv", separator="\t")
    generate_report(consensus_df, matrix_df, runs, diagnostics, out_dir / "alignment_report.html")

    for run_name, run_diag in diagnostics["runs"].items():
        _write_calibration_json(run_name, run_diag, diagnostics["reference"], out_dir)

    print(f"\nReference run: {diagnostics['reference']}")
    for run_name, run_diag in diagnostics["runs"].items():
        ms = run_diag["match_stats"]
        n_ref   = ms["n_ref"]
        n_run   = ms["n_run"]
        n_cmr   = ms["after_charge_mass_rt"]
        n_im    = ms["after_im"]
        n_match = ms["n_matched"]
        fallback = " [identity RT warp]" if run_diag.get("identity_fallback") else f" [{run_diag['n_anchors']} anchors]"
        print(f"\n  {run_name}{fallback}")
        mass_slope = run_diag.get("mass_drift_slope_ppm_per_min", 0.0)
        mass_intercept = run_diag.get("mass_drift_intercept_ppm", 0.0)
        im_slope = run_diag.get("im_drift_slope_per_min", 0.0)
        im_intercept = run_diag.get("im_drift_intercept", 0.0)
        print(f"    mass drift:              {mass_slope:+.4f} ppm/min  (intercept {mass_intercept:+.3f} ppm)")
        if im_slope != 0.0 or im_intercept != 0.0:
            print(f"    IM drift:                {im_slope:+.5f} 1/K0/min  (intercept {im_intercept:+.4f} 1/K0)")
        print(f"    ref features:            {n_ref:>8,}")
        print(f"    run features:            {n_run:>8,}")
        print(f"    ref with candidate       {n_cmr:>8,}  ({n_cmr/n_ref*100:.1f}%)  [charge+mass+RT]")
        print(f"    ref with candidate       {n_im:>8,}  ({n_im/n_ref*100:.1f}%)  [after IM filter, {n_cmr - n_im:,} lost]")
        print(f"    matched (best per ref):  {n_match:>8,}  ({n_match/n_ref*100:.1f}%)  [{n_ref - n_match:,} ref features unmatched]")

        if "fdr" in run_diag:
            fdr = run_diag["fdr"]
            flag = "  *** ABOVE THRESHOLD" if fdr["fdr"] > args.fdr_threshold else ""
            print(f"    FDR estimate:")
            print(f"      target matches:         {fdr['n_target']:>8,}")
            print(f"      decoy matches:          {fdr['n_decoy']:>8,}")
            print(f"      FDR (decoy/target):     {float(fdr['fdr'])*100:>7.2f}%{flag}")

        if "sweep" in run_diag:
            sweep_df: pl.DataFrame = run_diag["sweep"]
            best = sweep_df.filter(pl.col("passes_fdr")).head(5)
            if len(best):
                print(f"\n    Top parameter combinations (FDR ≤ {args.fdr_threshold*100:.1f}%):")
                print(f"    {'mass_ppm':>8}  {'rt_window':>9}  {'im_tol':>6}  {'n_target':>8}  {'FDR%':>6}")
                for row in best.iter_rows(named=True):
                    print(f"    {row['mass_ppm']:>8.1f}  {row['rt_window']:>9.2f}  "
                          f"{row['im_tolerance']:>6.3f}  {row['n_target']:>8,}  "
                          f"{row['fdr']*100:>5.2f}%")
                r = best.row(0, named=True)
                print(f"\n    Recommended: --mass-ppm {r['mass_ppm']} --rt-window {r['rt_window']} "
                      f"--im-tolerance {r['im_tolerance']}")
            else:
                print(f"\n    No parameter combination achieved FDR ≤ {args.fdr_threshold*100:.1f}%")
                print(f"    Lowest FDR found: {sweep_df['fdr'].to_numpy().min()*100:.2f}%")

    n_detected = consensus_df["n_runs_detected"]
    print(f"\nAlignment complete in {elapsed:.1f}s")
    print(f"  Consensus features:   {len(consensus_df):,}")
    print(f"  Detected in all runs: {(n_detected == len(runs)).sum():,}")
    print(f"  Detected in ≥2 runs:  {(n_detected >= 2).sum():,}")
    print(f"\nOutputs written to {out_dir}/")
    print("  consensus_features.tsv")
    print("  intensity_matrix.tsv")
    print("  alignment_report.html")
    for run_name in diagnostics["runs"]:
        print(f"  {run_name}_calibration.json")


if __name__ == "__main__":
    main(_parse_args())
