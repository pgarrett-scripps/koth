"""
koth_ff — high-performance LC-MS feature finder.

Functions
---------
detect_hills(path, ...) -> polars.DataFrame
    Stage 1: detect chromatographic hills from an mzML or Bruker .d file.

detect_features(hills, ...) -> polars.DataFrame
    Stages 2–3: group hills into isotope features and score them.
    Accepts either the DataFrame returned by detect_hills or a list[dict].

run_pipeline(path, ...) -> dict[str, polars.DataFrame]
    Run all three stages. Returns {"hills": DataFrame, "features": DataFrame}.

align_runs(runs, ...) -> tuple[polars.DataFrame, polars.DataFrame]
    Align features from multiple runs (potentially different gradients) and
    return (consensus_features, intensity_matrix).
"""

from __future__ import annotations

from typing import Literal

import polars as pl

from .koth_ff import (
    detect_hills    as _detect_hills_raw,
    detect_features as _detect_features_raw,
    run_pipeline    as _run_pipeline_raw,
)

_KEEP_F64 = frozenset({"mz", "mass"})


def _downcast(df: pl.DataFrame) -> pl.DataFrame:
    """Cast Float64 columns to Float32, keeping mz/mass at full precision."""
    casts = [
        pl.col(c).cast(pl.Float32)
        for c in df.columns
        if df[c].dtype == pl.Float64 and c not in _KEEP_F64
    ]
    return df.with_columns(casts) if casts else df


def detect_hills(
    path: str,
    *,
    mz_tolerance: float = 8.0,
    mz_tolerance_type: Literal["ppm", "da"] = "ppm",
    min_scans: int = 3,
    max_gap: int = 0,
    split_hills: bool = True,
    im_tolerance: float = 0.05,
    im_tolerance_type: Literal["relative", "absolute"] = "relative",
    intensity_coverage: float = 1.0,
    global_min_mz: float = 0.0,
    global_max_mz: float = float("inf"),
    bruker_mz_ppm: float = 5.0,
    bruker_im_pct: float = 3.0,
) -> pl.DataFrame:
    """Detect chromatographic hills from *path* (mzML or Bruker .d).

    Returns a Polars DataFrame with one row per hill.
    ``intensity_profile`` is a ``List(Float32)`` column.

    Parameters
    ----------
    path:
        Path to an mzML file or a Bruker .d directory.
    mz_tolerance:
        m/z tolerance value.
    mz_tolerance_type:
        ``"ppm"`` or ``"da"``.
    min_scans:
        Minimum number of scans per hill.
    max_gap:
        Maximum allowed gap in scans within a hill.
    split_hills:
        Whether to split multi-apex hills.
    im_tolerance:
        Ion-mobility tolerance value.
    im_tolerance_type:
        ``"relative"`` or ``"absolute"``.
    intensity_coverage:
        Fraction of total hill intensity to retain, 0–1.
    global_min_mz:
        Discard peaks below this m/z.
    global_max_mz:
        Discard peaks above this m/z.
    bruker_mz_ppm:
        m/z tolerance in ppm used when centroiding Bruker .d data.
    bruker_im_pct:
        Ion-mobility tolerance in percent used when centroiding Bruker .d data.
    """
    return _downcast(pl.DataFrame(
        _detect_hills_raw(
            path,
            mz_tolerance=mz_tolerance,
            mz_tolerance_type=mz_tolerance_type,
            min_scans=min_scans,
            max_gap=max_gap,
            split_hills=split_hills,
            im_tolerance=im_tolerance,
            im_tolerance_type=im_tolerance_type,
            intensity_coverage=intensity_coverage,
            global_min_mz=global_min_mz,
            global_max_mz=global_max_mz,
            bruker_mz_ppm=bruker_mz_ppm,
            bruker_im_pct=bruker_im_pct,
        )
    ))


def detect_features(
    hills: pl.DataFrame | list[dict],
    *,
    mz_tolerance: float = 5.0,
    min_charge: int = 1,
    max_charge: int = 7,
    min_cosine_similarity: float = 0.5,
    im_tolerance: float = 0.05,
    max_isotopes: int = 6,
    isotope_offset_min: int = -1,
    isotope_offset_max: int = 1,
    offset_zero_bonus: float = 0.15,
    min_score_threshold: float = 0.5,
) -> pl.DataFrame:
    """Detect isotope features from *hills* and score them.

    *hills* may be the DataFrame returned by :func:`detect_hills` or the
    raw ``list[dict]``.  Returns a Polars DataFrame with one row per feature.
    List columns: ``elution_profile``, ``isotope_profile``, ``theoretical_pattern``.

    Parameters
    ----------
    hills:
        Hill data as a Polars DataFrame or list of dicts.
    mz_tolerance:
        m/z tolerance for isotope grouping in ppm.
    min_charge:
        Minimum charge state to consider.
    max_charge:
        Maximum charge state to consider.
    min_cosine_similarity:
        Minimum cosine similarity for isotope pattern matching.
    im_tolerance:
        Ion-mobility tolerance for grouping hills.
    max_isotopes:
        Maximum number of isotopes per feature.
    isotope_offset_min:
        Minimum neutron offset to test during scoring.
    isotope_offset_max:
        Maximum neutron offset to test during scoring.
    offset_zero_bonus:
        Score bonus applied when offset == 0.
    min_score_threshold:
        Features scoring below this keep offset=0.
    """
    if isinstance(hills, pl.DataFrame):
        hills = hills.to_dicts()
    return _downcast(pl.DataFrame(
        _detect_features_raw(
            hills,
            mz_tolerance=mz_tolerance,
            min_charge=min_charge,
            max_charge=max_charge,
            min_cosine_similarity=min_cosine_similarity,
            im_tolerance=im_tolerance,
            max_isotopes=max_isotopes,
            isotope_offset_min=isotope_offset_min,
            isotope_offset_max=isotope_offset_max,
            offset_zero_bonus=offset_zero_bonus,
            min_score_threshold=min_score_threshold,
        )
    ))


def run_pipeline(
    path: str,
    *,
    # Hills stage
    mz_tolerance: float = 8.0,
    mz_tolerance_type: Literal["ppm", "da"] = "ppm",
    min_scans: int = 3,
    max_gap: int = 0,
    split_hills: bool = True,
    im_tolerance: float = 0.05,
    im_tolerance_type: Literal["relative", "absolute"] = "relative",
    intensity_coverage: float = 1.0,
    global_min_mz: float = 0.0,
    global_max_mz: float = float("inf"),
    bruker_mz_ppm: float = 5.0,
    bruker_im_pct: float = 3.0,
    # Features stage
    features_mz_tolerance: float = 5.0,
    min_charge: int = 1,
    max_charge: int = 7,
    min_cosine_similarity: float = 0.5,
    max_isotopes: int = 6,
    # Scoring stage
    isotope_offset_min: int = -1,
    isotope_offset_max: int = 1,
    offset_zero_bonus: float = 0.15,
    min_score_threshold: float = 0.5,
) -> dict[str, pl.DataFrame]:
    """Run the full hill → feature → scoring pipeline.

    Returns a dict with two Polars DataFrames::

        result = koth_ff.run_pipeline("data.mzML", mz_tolerance=8.0)
        hills_df    = result["hills"]
        features_df = result["features"]

    Parameters
    ----------
    path:
        Path to an mzML file or a Bruker .d directory.
    mz_tolerance:
        m/z tolerance for the hills stage.
    mz_tolerance_type:
        ``"ppm"`` or ``"da"`` — hills stage only.
    min_scans:
        Minimum scans per hill.
    max_gap:
        Maximum scan gap within a hill.
    split_hills:
        Whether to split multi-apex hills.
    im_tolerance:
        Ion-mobility tolerance (shared by hills and features stages).
    im_tolerance_type:
        ``"relative"`` or ``"absolute"`` — hills stage only.
    intensity_coverage:
        Fraction of total hill intensity to retain, 0–1.
    global_min_mz:
        Discard peaks below this m/z.
    global_max_mz:
        Discard peaks above this m/z.
    bruker_mz_ppm:
        m/z tolerance in ppm for Bruker .d centroiding.
    bruker_im_pct:
        Ion-mobility tolerance in percent for Bruker .d centroiding.
    features_mz_tolerance:
        m/z tolerance in ppm for the features stage.
    min_charge:
        Minimum charge state.
    max_charge:
        Maximum charge state.
    min_cosine_similarity:
        Minimum cosine similarity for isotope matching.
    max_isotopes:
        Maximum isotopes per feature.
    isotope_offset_min:
        Minimum neutron offset tested during scoring.
    isotope_offset_max:
        Maximum neutron offset tested during scoring.
    offset_zero_bonus:
        Score bonus for zero offset.
    min_score_threshold:
        Features below this score keep offset=0.
    """
    raw = _run_pipeline_raw(
        path,
        mz_tolerance=mz_tolerance,
        mz_tolerance_type=mz_tolerance_type,
        min_scans=min_scans,
        max_gap=max_gap,
        split_hills=split_hills,
        im_tolerance=im_tolerance,
        im_tolerance_type=im_tolerance_type,
        intensity_coverage=intensity_coverage,
        global_min_mz=global_min_mz,
        global_max_mz=global_max_mz,
        bruker_mz_ppm=bruker_mz_ppm,
        bruker_im_pct=bruker_im_pct,
        min_charge=min_charge,
        max_charge=max_charge,
        min_cosine_similarity=min_cosine_similarity,
        max_isotopes=max_isotopes,
        isotope_offset_min=isotope_offset_min,
        isotope_offset_max=isotope_offset_max,
        offset_zero_bonus=offset_zero_bonus,
        min_score_threshold=min_score_threshold,
    )
    return {
        "hills":    _downcast(pl.DataFrame(raw["hills"])),
        "features": _downcast(pl.DataFrame(raw["features"])),
    }


from .align import align_runs

__all__ = ["detect_hills", "detect_features", "run_pipeline", "align_runs"]
