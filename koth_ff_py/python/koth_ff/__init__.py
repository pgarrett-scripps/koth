"""
koth_ff — high-performance LC-MS feature finder.

Functions
---------
detect_hills(path, **kwargs) -> polars.DataFrame
    Stage 1: detect chromatographic hills from an mzML or Bruker .d file.

detect_features(hills, **kwargs) -> polars.DataFrame
    Stages 2–3: group hills into isotope features and score them.
    Accepts either the DataFrame returned by detect_hills or a list[dict].

run_pipeline(path, **kwargs) -> dict[str, polars.DataFrame]
    Run all three stages. Returns {"hills": DataFrame, "features": DataFrame}.

All functions accept keyword arguments to override any config field:
    mz_tolerance, min_scans, max_gap, intensity_coverage, split_hills,
    im_tolerance, global_min_mz, global_max_mz,   # hills stage
    min_charge, max_charge, min_cosine_similarity,  # features stage
    min_score_threshold, offset_zero_bonus          # scoring stage
"""

from __future__ import annotations

import polars as pl

from .koth_ff import (
    detect_hills    as _detect_hills_raw,
    detect_features as _detect_features_raw,
    run_pipeline    as _run_pipeline_raw,
)


def detect_hills(path: str, **kwargs) -> pl.DataFrame:
    """Detect chromatographic hills from *path* (mzML or Bruker .d).

    Returns a Polars DataFrame with one row per hill.
    ``intensity_profile`` is a ``List(Float32)`` column.
    """
    return pl.DataFrame(_detect_hills_raw(path, **kwargs))


def detect_features(
    hills: pl.DataFrame | list[dict],
    **kwargs,
) -> pl.DataFrame:
    """Detect isotope features from *hills* and score them.

    *hills* may be the DataFrame returned by :func:`detect_hills` or the
    raw ``list[dict]``.  Returns a Polars DataFrame with one row per feature.
    List columns: ``elution_profile``, ``isotope_profile``, ``theoretical_pattern``.
    """
    if isinstance(hills, pl.DataFrame):
        hills = hills.to_dicts()
    return pl.DataFrame(_detect_features_raw(hills, **kwargs))


def run_pipeline(path: str, **kwargs) -> dict[str, pl.DataFrame]:
    """Run the full hill → feature → scoring pipeline.

    Returns a dict with two Polars DataFrames::

        result = koth_ff.run_pipeline("data.mzML", mz_tolerance=8.0)
        hills_df    = result["hills"]
        features_df = result["features"]
    """
    raw = _run_pipeline_raw(path, **kwargs)
    return {
        "hills":    pl.DataFrame(raw["hills"]),
        "features": pl.DataFrame(raw["features"]),
    }


__all__ = ["detect_hills", "detect_features", "run_pipeline"]
