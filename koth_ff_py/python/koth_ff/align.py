"""
align — multi-run feature alignment and quantification matrix builder.

Aligns LC-MS feature sets from runs with potentially different gradients by
normalising RT to [0, 1] before finding anchor pairs, then estimating a smooth
RT correction curve via Gaussian kernel regression (running weighted average of
RT deltas), and finally doing a mass/charge/RT match to build a feature × sample
intensity matrix.
"""

from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def align_runs(
    runs: dict[str, pl.DataFrame],
    *,
    reference: str | None = None,
    mass_ppm: float = 10.0,
    anchor_mass_ppm: float = 5.0,
    rt_window: float = 0.5,
    rt_anchor_window: float = 0.05,
    im_tolerance: float = 0.05,
    min_anchor_score: float = 0.7,
    min_anchor_count: int = 20,
    rt_warp_bandwidth: float = 0.15,
    intensity_col: str = "intensity_sum",
    return_diagnostics: bool = False,
) -> tuple[pl.DataFrame, pl.DataFrame] | tuple[pl.DataFrame, pl.DataFrame, dict]:
    """Align features from multiple runs and return a quantification matrix.

    Parameters
    ----------
    runs:
        Mapping of sample name → features DataFrame (from ``run_pipeline()``).
    reference:
        Name of the reference run. ``None`` auto-selects the run with the most
        high-confidence features (score >= min_anchor_score).
    mass_ppm:
        Mass tolerance for final feature matching (parts-per-million).
    anchor_mass_ppm:
        Mass tolerance used only for finding RT anchor pairs. Should be tighter
        than ``mass_ppm`` to avoid false coincidental mass matches polluting the
        warp estimate. Default 5 ppm.
    rt_window:
        RT tolerance for final matching, in the reference run's absolute
        minutes (applied after RT warping).
    rt_anchor_window:
        RT window for anchor-pair finding in **normalised** RT space [0, 1].
        0.05 means ±5 % of the total gradient length.
    im_tolerance:
        Ion-mobility tolerance in 1/K₀ units. Applied only when both features
        have ``im > 0``; ignored otherwise (mixed-IM datasets).
    min_anchor_score:
        Minimum Bhattacharyya score for a feature to be used as an anchor.
    min_anchor_count:
        Minimum number of anchor pairs required to fit a warp. Runs with fewer
        anchors fall back to identity (no correction) with a warning.
    intensity_col:
        Column from the features DataFrame to use as the quantification value.
        Typically ``"intensity_sum"`` or ``"intensity_apex"``.

    Returns
    -------
    consensus_df:
        One row per consensus feature with columns:
        ``mass, mz, charge, rt_apex, im, n_runs_detected``.
    matrix_df:
        Same rows as ``consensus_df`` plus one column per sample containing
        the intensity value (``0.0`` when not detected in that run).
    """
    if len(runs) < 2:
        raise ValueError("align_runs requires at least 2 runs.")

    ref_name = reference if reference is not None else _pick_reference(runs, min_anchor_score)
    ref_df = runs[ref_name]

    matched: dict[str, pl.DataFrame] = {ref_name: ref_df}
    diag_runs: dict[str, dict] = {}

    ref_rt_min = float(ref_df["rt_apex"].min())
    ref_rt_max = float(ref_df["rt_apex"].max())

    for name, df in runs.items():
        if name == ref_name:
            continue
        warp_fn, warp_diag = _build_warp(
            ref_df, df,
            mass_ppm=anchor_mass_ppm,
            rt_anchor_window=rt_anchor_window,
            min_anchor_score=min_anchor_score,
            min_anchor_count=min_anchor_count,
            bandwidth=rt_warp_bandwidth,
        )
        matched_df, match_stats = _match_features(ref_df, df, warp_fn, mass_ppm, rt_window, im_tolerance, intensity_col)
        diag_runs[name] = {**warp_diag, "match_stats": match_stats}
        matched[name] = matched_df

    consensus_df, matrix_df = _build_matrix(ref_df, matched, intensity_col)

    if return_diagnostics:
        diagnostics = {
            "reference": ref_name,
            "ref_rt_min": ref_rt_min,
            "ref_rt_max": ref_rt_max,
            "runs": diag_runs,
        }
        return consensus_df, matrix_df, diagnostics
    return consensus_df, matrix_df


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pick_reference(runs: dict[str, pl.DataFrame], min_score: float) -> str:
    best_name = max(
        runs,
        key=lambda n: (runs[n]["score"] >= min_score).sum(),
    )
    return best_name


def _select_anchors(df: pl.DataFrame, min_score: float) -> pl.DataFrame:
    return df.filter(pl.col("score") >= min_score)


def _match_anchors(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    mass_ppm: float,
    rt_norm_window: float,
    ref_rt_min: float,
    ref_rt_max: float,
    run_rt_min: float,
    run_rt_max: float,
) -> list[tuple[float, float]]:
    """Return (ref_rt_norm, run_rt_norm) anchor pairs.

    Normalisation uses the FULL run RT range (not the anchor subset range) so
    that the normalised coordinates are consistent with the warp function.
    """
    ref_span = ref_rt_max - ref_rt_min or 1.0
    run_span = run_rt_max - run_rt_min or 1.0

    ref_anchors = ref_df.with_columns(
        ((pl.col("rt_apex").cast(pl.Float64) - ref_rt_min) / ref_span).alias("rt_norm")
    )
    run_anchors = run_df.with_columns(
        ((pl.col("rt_apex").cast(pl.Float64) - run_rt_min) / run_span).alias("rt_norm")
    )

    # Inner-join on charge first (O(n)), then filter mass and RT
    ref_s = ref_anchors.select(
        pl.col("mass").alias("ref_mass"),
        pl.col("charge"),
        pl.col("rt_norm").alias("ref_rt_norm"),
    )
    run_s = run_anchors.select(
        pl.col("mass").alias("run_mass"),
        pl.col("charge"),
        pl.col("rt_norm").alias("run_rt_norm"),
    )
    pairs = ref_s.join(run_s, on="charge", how="inner").filter(
        (((pl.col("ref_mass") - pl.col("run_mass")).abs() / pl.col("ref_mass") * 1e6) <= mass_ppm)
        & ((pl.col("ref_rt_norm") - pl.col("run_rt_norm")).abs() <= rt_norm_window)
    )

    if len(pairs) == 0:
        return []

    # Keep best mass match per reference feature
    pairs = pairs.with_columns(
        (((pl.col("ref_mass") - pl.col("run_mass")).abs() / pl.col("ref_mass")) * 1e6).alias("_ppm")
    ).sort("_ppm").unique(subset=["ref_rt_norm"], keep="first")

    return list(zip(
        pairs["ref_rt_norm"].to_list(),
        pairs["run_rt_norm"].to_list(),
    ))


def _fit_rt_warp(
    anchor_pairs: list[tuple[float, float]],
    ref_rt_min: float,
    ref_rt_max: float,
    run_rt_min: float,
    run_rt_max: float,
    bandwidth: float = 0.15,
) -> Callable[[np.ndarray], np.ndarray]:
    """Build a warp function using Gaussian kernel regression on RT deltas.

    For each query RT, the correction is a weighted average of anchor deltas
    (ref_norm - run_norm), with weights falling off as a Gaussian in normalised
    RT distance.  This is more robust than exact interpolation because noisy
    individual anchors are smoothed out rather than fitted exactly.
    """
    run_norms = np.array([p[1] for p in anchor_pairs])
    deltas = np.array([p[0] - p[1] for p in anchor_pairs])  # ref_norm - run_norm

    ref_span = ref_rt_max - ref_rt_min
    run_span = run_rt_max - run_rt_min or 1.0

    def warp(rt_abs: np.ndarray) -> np.ndarray:
        rt_norm = (np.asarray(rt_abs, dtype=float) - run_rt_min) / run_span
        # (n_query, n_anchors) Gaussian weights
        w = np.exp(-0.5 * ((rt_norm[:, np.newaxis] - run_norms[np.newaxis, :]) / bandwidth) ** 2)
        w_sum = w.sum(axis=1)
        # Where no anchors are nearby, fall back to nearest anchor delta
        no_support = w_sum < 1e-10
        w_sum = np.where(no_support, 1.0, w_sum)
        delta = (w * deltas[np.newaxis, :]).sum(axis=1) / w_sum
        if no_support.any():
            nearest = np.abs(rt_norm[no_support, np.newaxis] - run_norms[np.newaxis, :]).argmin(axis=1)
            delta[no_support] = deltas[nearest]
        return (rt_norm + delta) * ref_span + ref_rt_min

    return warp


def _build_warp(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    mass_ppm: float,
    rt_anchor_window: float,
    min_anchor_score: float,
    min_anchor_count: int,
    bandwidth: float = 0.15,
) -> tuple[Callable[[np.ndarray], np.ndarray], dict]:
    _ref_rt = ref_df["rt_apex"].cast(pl.Float64).to_numpy()
    _run_rt = run_df["rt_apex"].cast(pl.Float64).to_numpy()
    ref_rt_min, ref_rt_max = float(_ref_rt.min()), float(_ref_rt.max())
    run_rt_min, run_rt_max = float(_run_rt.min()), float(_run_rt.max())

    ref_anchors = _select_anchors(ref_df, min_anchor_score)
    run_anchors = _select_anchors(run_df, min_anchor_score)

    pairs = _match_anchors(
        ref_anchors, run_anchors, mass_ppm, rt_anchor_window,
        ref_rt_min, ref_rt_max, run_rt_min, run_rt_max,
    )

    diag = {
        "anchor_pairs": pairs,
        "n_anchors": len(pairs),
        "bandwidth": bandwidth,
        "run_rt_min": run_rt_min,
        "run_rt_max": run_rt_max,
        "ref_rt_min": ref_rt_min,
        "ref_rt_max": ref_rt_max,
        "identity_fallback": False,
    }

    if len(pairs) < min_anchor_count:
        warnings.warn(
            f"Only {len(pairs)} anchor pairs found (need {min_anchor_count}); "
            "falling back to identity RT warp (no correction).",
            RuntimeWarning,
            stacklevel=3,
        )
        diag["identity_fallback"] = True
        return lambda rt: rt, diag  # identity

    warp_fn = _fit_rt_warp(pairs, ref_rt_min, ref_rt_max, run_rt_min, run_rt_max, bandwidth)
    diag["warp_fn"] = warp_fn
    return warp_fn, diag


def _match_features(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    warp_fn: Callable[[np.ndarray], np.ndarray],
    mass_ppm: float,
    rt_window: float,
    im_tolerance: float,
    intensity_col: str,
) -> tuple[pl.DataFrame, dict]:
    """Return (DataFrame with ref feature rows + warped intensity, match stats dict)."""
    n_ref = len(ref_df)
    n_run = len(run_df)

    run_rt_warped = warp_fn(run_df["rt_apex"].cast(pl.Float64).to_numpy())

    run_aligned = run_df.with_columns(
        pl.Series("rt_warped", run_rt_warped, dtype=pl.Float64)
    ).select(["mass", "charge", "im", intensity_col, "rt_warped"])

    ref_keyed = ref_df.select(["mass", "charge", "im", "rt_apex"]).with_row_index("_ref_idx")
    run_keyed = run_aligned.with_row_index("_run_idx")

    cross = ref_keyed.join(run_keyed, how="cross")

    cross = cross.filter(
        (pl.col("charge") == pl.col("charge_right"))
        & (
            ((pl.col("mass") - pl.col("mass_right")).abs() / pl.col("mass") * 1e6) <= mass_ppm
        )
        & ((pl.col("rt_apex") - pl.col("rt_warped")).abs() <= rt_window)
    )
    # Count distinct ref features with ≥1 candidate (not raw pair count)
    n_after_charge_mass_rt = cross["_ref_idx"].n_unique()

    # Apply IM filter only when both features have im > 0
    cross = cross.filter(
        (pl.col("im") == 0.0)
        | (pl.col("im_right") == 0.0)
        | ((pl.col("im") - pl.col("im_right")).abs() <= im_tolerance)
    )
    n_after_im = cross["_ref_idx"].n_unique()

    empty_result = ref_df.select(["mass", "charge", "rt_apex", "im"]).with_columns(
        pl.lit(0.0).cast(pl.Float64).alias(intensity_col)
    )

    if len(cross) == 0:
        stats = {
            "n_ref": n_ref, "n_run": n_run,
            "after_charge_mass_rt": 0, "after_im": 0,
            "after_dedup": 0, "n_matched": 0,
        }
        return empty_result, stats

    # Keep best mass match per reference feature
    best = (
        cross
        .with_columns(
            (((pl.col("mass") - pl.col("mass_right")).abs() / pl.col("mass")) * 1e6).alias("_ppm")
        )
        .sort("_ppm")
        .unique(subset=["_ref_idx"], keep="first")
        .select(["_ref_idx", intensity_col])
    )
    n_after_dedup = len(best)  # already one row per ref feature

    result = (
        ref_df.select(["mass", "charge", "rt_apex", "im"])
        .with_row_index("_ref_idx")
        .join(best, on="_ref_idx", how="left")
        .with_columns(pl.col(intensity_col).fill_null(0.0))
        .drop("_ref_idx")
    )
    n_matched = int((result[intensity_col] > 0).sum())

    stats = {
        "n_ref": n_ref,
        "n_run": n_run,
        "after_charge_mass_rt": n_after_charge_mass_rt,
        "after_im": n_after_im,
        "after_dedup": n_after_dedup,
        "n_matched": n_matched,
    }
    return result, stats


def _build_matrix(
    ref_df: pl.DataFrame,
    matched: dict[str, pl.DataFrame],
    intensity_col: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    consensus = ref_df.select(["mass", "mz", "charge", "rt_apex", "im"]).with_row_index("_idx")

    intensity_frames = []
    for name, df in matched.items():
        col = df.select(intensity_col).rename({intensity_col: name})
        intensity_frames.append(col)

    n_detected = pl.Series(
        "n_runs_detected",
        sum(
            (matched[n].select(intensity_col).to_series() > 0).to_numpy()
            for n in matched
        ),
    )

    consensus_out = consensus.drop("_idx").with_columns(n_detected)

    matrix_df = pl.concat(
        [consensus.drop(["mz", "im"])] + intensity_frames,
        how="horizontal",
    ).drop("_idx")

    return consensus_out, matrix_df
