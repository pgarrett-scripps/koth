"""
Align feature tables from a batch run and produce a quantification matrix.

Usage:
    uv run python scripts/run_align.py <batch_output_dir> [aligned_output_dir]
                                       [--fdr] [--sweep]
                                       [--fdr-threshold FLOAT]
                                       [--mass-ppm FLOAT] [--rt-window FLOAT]
                                       [--im-tolerance FLOAT]

Expects <batch_output_dir> to contain one subdirectory per sample, each with
a features.parquet file (as written by run_batch.py).

Outputs (in <aligned_output_dir>/):
    consensus_features.parquet   – one row per consensus feature
    intensity_matrix.parquet     – consensus features × sample intensities
    intensity_matrix.tsv         – same, as tab-separated text (no list cols)
    alignment_report.html        – visual QC report

Target / decoy flags:
    --fdr     Estimate false discovery rate using a mass-shifted decoy reference.
              Adds per-run FDR estimates to the report and console output.
    --sweep   Grid search over (mass_ppm, rt_window, im_tolerance) to find the
              tightest tolerances that keep FDR below --fdr-threshold (default 1%).
              Implies --fdr. Adds heatmap plots to the report.
"""

import argparse
import sys
import time
from pathlib import Path

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
    p.add_argument("--mass-ppm", type=float, default=10.0, metavar="FLOAT",
                   help="Mass tolerance for final feature matching in ppm (default: 20)")
    rt_group = p.add_mutually_exclusive_group()
    rt_group.add_argument("--rt-window", type=float, default=None, metavar="FLOAT",
                          help="RT tolerance for final matching in minutes")
    rt_group.add_argument("--rt-window-pct", type=float, default=None, metavar="FLOAT",
                          help="RT tolerance as %% of reference gradient length "
                               "(e.g. 1.5 means 1.5%% of the gradient). "
                               "Scales automatically across different gradient lengths.")
    p.add_argument("--im-tolerance", type=float, default=1.0, metavar="FLOAT",
                   help="Ion-mobility tolerance in 1/K0 units (default: 1.0)")
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


def load_runs(batch_dir: Path) -> dict[str, pl.DataFrame]:
    runs: dict[str, pl.DataFrame] = {}
    for sub in sorted(batch_dir.iterdir()):
        if not sub.is_dir():
            continue
        if (pq := sub / "features.parquet").exists():
            df = pl.read_parquet(pq)
        elif (tsv := sub / "features.tsv").exists():
            df = pl.read_csv(tsv, separator="\t")
            renames = {k: v for k, v in _COLUMN_RENAMES.items() if k in df.columns}
            if renames:
                df = df.rename(renames)
        else:
            continue

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
        rt_window = args.rt_window if args.rt_window is not None else 1.0

    run_fdr = args.fdr or args.sweep
    if run_fdr:
        print(f"\nTarget/decoy FDR enabled (mass shift = {args.mass_shift_da} Da, "
              f"threshold = {args.fdr_threshold*100:.1f}%)")
    if args.sweep:
        print("Parameter sweep enabled — this will take a moment...")

    t0 = time.perf_counter()
    consensus_df, matrix_df, diagnostics = align_runs(
        runs,
        anchor_mass_ppm=20.0,
        mass_ppm=args.mass_ppm,
        rt_window=rt_window,
        rt_anchor_window=0.05,
        im_tolerance=args.im_tolerance,
        min_anchor_score=0.7,
        min_anchor_count=10,
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


if __name__ == "__main__":
    main(_parse_args())
