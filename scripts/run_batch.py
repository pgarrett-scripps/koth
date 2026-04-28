"""
Batch hill/feature extraction for all mzML(.gz) files in a directory.

Usage:
    uv run python scripts/run_batch.py <input_dir> [output_dir]

Outputs per run (in <output_dir>/<stem>/):
    hills.parquet / features.parquet   – full data, all columns
    hills.tsv / features.tsv           – first 50 + last 50 rows, list columns as JSON
"""

import json
import sys
import time
from pathlib import Path

import polars as pl

import koth_ff


def _head_tail(df: pl.DataFrame, n: int = 50) -> pl.DataFrame:
    """Return first n + last n rows (or the full df if shorter)."""
    if len(df) <= 2 * n:
        return df
    return pl.concat([df.head(n), df.tail(n)])


def _flatten_lists(df: pl.DataFrame) -> pl.DataFrame:
    """Serialize any List columns to JSON strings so TSV writer can handle them."""
    list_cols = [c for c, d in zip(df.columns, df.dtypes) if isinstance(d, pl.List)]
    if not list_cols:
        return df
    return df.with_columns([
        pl.col(c).map_elements(lambda x: json.dumps(x.to_list()), return_dtype=pl.String)
        for c in list_cols
    ])


def run_batch(input_dir: Path, output_dir: Path) -> None:
    files = sorted(
        p for p in input_dir.iterdir()
        if p.name.lower().endswith(".mzml") or p.name.lower().endswith(".mzml.gz")
    )

    if not files:
        print(f"No mzML files found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(files)} file(s) in {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, path in enumerate(files, 1):
        stem = path.name.removesuffix(".gz").removesuffix(".mzML").removesuffix(".mzml")
        out = output_dir / stem
        out.mkdir(exist_ok=True)

        print(f"\n[{i}/{len(files)}] {path.name}")
        t0 = time.perf_counter()

        result = koth_ff.run_pipeline(
            str(path),
            mz_tolerance=5.0,
            mz_tolerance_type="ppm",
            min_scans=4,
            max_gap=0,
            split_hills=False,
            min_charge=1,
            max_charge=9,
            min_cosine_similarity=0.3,
            max_isotopes=6,
            min_score_threshold=0.3,
        )

        hills: pl.DataFrame    = result["hills"]
        features: pl.DataFrame = result["features"]

        hills.write_parquet(out / "hills.parquet")
        features.write_parquet(out / "features.parquet")

        _flatten_lists(_head_tail(hills)).write_csv(out / "hills.tsv", separator="\t")
        _flatten_lists(_head_tail(features)).write_csv(out / "features.tsv", separator="\t")

        elapsed = time.perf_counter() - t0
        print(f"  hills:    {len(hills):>8,}  →  {out / 'hills.parquet'}")
        print(f"  features: {len(features):>8,}  →  {out / 'features.parquet'}")
        print(f"  elapsed:  {elapsed:.1f}s")

    print(f"\nDone. Results in {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_dir  = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else input_dir.parent / "koth_ff_output"

    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory")
        sys.exit(1)

    run_batch(input_dir, output_dir)
