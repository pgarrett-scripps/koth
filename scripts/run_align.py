"""
Align feature tables from a batch run and produce a quantification matrix.

Usage:
    uv run python scripts/run_align.py <batch_output_dir> [aligned_output_dir]

Expects <batch_output_dir> to contain one subdirectory per sample, each with
a features.parquet file (as written by run_batch.py).

Outputs (in <aligned_output_dir>/):
    consensus_features.parquet   – one row per consensus feature
    intensity_matrix.parquet     – consensus features × sample intensities
    intensity_matrix.tsv         – same, as tab-separated text (no list cols)
    alignment_report.html        – visual QC report
"""

import sys
import time
from pathlib import Path

import polars as pl

from koth_ff.align import align_runs
from alignment_report import generate_report


def load_runs(batch_dir: Path) -> dict[str, pl.DataFrame]:
    runs: dict[str, pl.DataFrame] = {}
    for sub in sorted(batch_dir.iterdir()):
        pq = sub / "features.parquet"
        if sub.is_dir() and pq.exists():
            df = pl.read_parquet(pq)
            runs[sub.name] = df
    return runs


def main(batch_dir: Path, out_dir: Path) -> None:
    print(f"Loading features from {batch_dir}")
    runs = load_runs(batch_dir)

    if len(runs) < 2:
        print(f"Found {len(runs)} run(s) — need at least 2 to align.")
        sys.exit(1)

    print(f"Found {len(runs)} runs:")
    for name, df in runs.items():
        print(f"  {name}: {len(df):,} features")

    t0 = time.perf_counter()
    consensus_df, matrix_df, diagnostics = align_runs(
        runs,
        anchor_mass_ppm=15.0,    # tight — reliable landmark pairs only for warp
        mass_ppm=15.0,           # looser — final matching after RT correction
        rt_window=2.0,
        rt_anchor_window=0.1,
        im_tolerance=1.0,
        min_anchor_score=0.5,
        min_anchor_count=10,
        intensity_col="intensity_sum",
        return_diagnostics=True,
    )
    elapsed = time.perf_counter() - t0

    out_dir.mkdir(parents=True, exist_ok=True)

    consensus_df.write_parquet(out_dir / "consensus_features.parquet")
    matrix_df.write_parquet(out_dir / "intensity_matrix.parquet")
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

    n_detected = consensus_df["n_runs_detected"]
    print(f"\nAlignment complete in {elapsed:.1f}s")
    print(f"  Consensus features:   {len(consensus_df):,}")
    print(f"  Detected in all runs: {(n_detected == len(runs)).sum():,}")
    print(f"  Detected in ≥2 runs:  {(n_detected >= 2).sum():,}")
    print(f"\nOutputs written to {out_dir}/")
    print("  consensus_features.parquet")
    print("  intensity_matrix.parquet")
    print("  intensity_matrix.tsv")
    print("  alignment_report.html")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    batch_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else batch_dir.parent / "aligned"

    if not batch_dir.is_dir():
        print(f"Error: {batch_dir} is not a directory")
        sys.exit(1)

    main(batch_dir, out_dir)
