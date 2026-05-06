"""
align — multi-run feature alignment and quantification matrix builder.

Aligns LC-MS feature sets from runs with potentially different gradients by
normalising RT to [0, 1] before finding anchor pairs, then estimating a smooth
RT correction curve via sliding-window medians (iteratively sigma-clipped) fit
with a smoothing spline, and finally doing a mass/charge/RT match to build a
feature × sample intensity matrix.
"""

from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
import polars as pl
from scipy.interpolate import make_smoothing_spline


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
    min_anchor_score: float = 0.75,
    min_feature_score: float | None = 0.7,
    min_anchor_count: int = 20,
    rt_warp_bandwidth: float = 0.05,
    rt_warp_sigma_clip: float = 3.0,
    rt_warp_clip_iters: int = 5,
    rt_warp_spline_lam: float | None = None,
    intensity_col: str = "intensity_sum",
    return_diagnostics: bool = False,
    run_fdr_estimation: bool = False,
    run_parameter_sweep: bool = False,
    fdr_threshold: float = 0.01,
    mass_shift_da: float = 100.0,
    sweep_mass_ppm: list[float] | None = None,
    sweep_rt_window: list[float] | None = None,
    sweep_im_tolerance: list[float] | None = None,
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
    min_feature_score:
        Minimum Bhattacharyya score for a reference feature to appear as a row
        in the output matrix. Features below this threshold are dropped from the
        matrix entirely — they cannot initiate a match. When ``None`` all
        reference features are kept (default behaviour). Features in non-reference
        runs may still have any score; a qualifying reference feature can match a
        lower-scoring counterpart in another run. Matching is performed greedily
        in descending reference-score order so the highest-confidence features
        claim their best run-counterpart first.
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
    ref_df_full = runs[ref_name]

    # Filter which reference features become matrix rows. The full reference is
    # still used for warp building so the RT span is accurate.
    if min_feature_score is not None and "score" in ref_df_full.columns:
        ref_df = ref_df_full.filter(pl.col("score") >= min_feature_score)
    else:
        ref_df = ref_df_full

    matched: dict[str, pl.DataFrame] = {ref_name: ref_df}
    diag_runs: dict[str, dict] = {}

    ref_rt_min = float(ref_df_full["rt_apex"].min())
    ref_rt_max = float(ref_df_full["rt_apex"].max())

    for name, df in runs.items():
        if name == ref_name:
            continue
        warp_fn, mass_corr_fn, im_corr_fn, warp_diag = _build_warp(
            ref_df_full, df,
            mass_ppm=anchor_mass_ppm,
            rt_anchor_window=rt_anchor_window,
            min_anchor_score=min_anchor_score,
            min_anchor_count=min_anchor_count,
            bandwidth=rt_warp_bandwidth,
            sigma_clip=rt_warp_sigma_clip,
            clip_iters=rt_warp_clip_iters,
            spline_lam=rt_warp_spline_lam,
        )
        matched_df, match_stats = _match_features(
            ref_df, df, warp_fn, mass_ppm, rt_window, im_tolerance, intensity_col,
            mass_corr_fn=mass_corr_fn, im_corr_fn=im_corr_fn,
        )
        diag_runs[name] = {**warp_diag, "match_stats": match_stats}
        matched[name] = matched_df

        if run_fdr_estimation or run_parameter_sweep:
            diag_runs[name]["fdr"] = estimate_match_fdr(
                ref_df, df, warp_fn, mass_ppm, rt_window, im_tolerance,
                intensity_col, mass_shift_da=mass_shift_da,
            )
        if run_parameter_sweep:
            diag_runs[name]["sweep"] = sweep_match_parameters(
                ref_df, df, warp_fn, intensity_col,
                mass_ppm_values=sweep_mass_ppm or [5.0, 10.0, 15.0, 20.0, 30.0],
                rt_window_values=sweep_rt_window or [0.25, 0.5, 0.75, 1.0, 1.5],
                im_tolerance_values=sweep_im_tolerance or [0.03, 0.05, 0.1, 0.5, 1.0],
                fdr_threshold=fdr_threshold,
                mass_shift_da=mass_shift_da,
            )

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
) -> tuple[list[tuple[float, float]], np.ndarray, np.ndarray, np.ndarray]:
    """Return RT anchor pairs and per-anchor drift data.

    Normalisation uses the FULL run RT range (not the anchor subset range) so
    that the normalised coordinates are consistent with the warp function.

    Returns
    -------
    pairs : list of (ref_rt_norm, run_rt_norm)
    ref_rt_abs : absolute RT of each matched ref feature (minutes)
    ppm_errors : (run_mass - ref_mass) / ref_mass * 1e6 per pair
    im_deltas  : run_im - ref_im per pair; 0.0 when either side has no IM
    """
    ref_span = ref_rt_max - ref_rt_min or 1.0
    run_span = run_rt_max - run_rt_min or 1.0

    ref_mass     = ref_df["mass"].cast(pl.Float64).to_numpy()
    ref_charge   = ref_df["charge"].to_numpy()
    ref_rt_abs_arr = ref_df["rt_apex"].cast(pl.Float64).to_numpy()
    ref_rt_norm  = (ref_rt_abs_arr - ref_rt_min) / ref_span
    ref_im       = ref_df["im"].cast(pl.Float64).to_numpy() if "im" in ref_df.columns else np.zeros(len(ref_df))

    run_mass   = run_df["mass"].cast(pl.Float64).to_numpy()
    run_charge = run_df["charge"].to_numpy()
    run_rt_norm = (run_df["rt_apex"].cast(pl.Float64).to_numpy() - run_rt_min) / run_span
    run_im     = run_df["im"].cast(pl.Float64).to_numpy() if "im" in run_df.columns else np.zeros(len(run_df))

    # Build per-charge sorted index into run (by mass) for binary search
    run_index: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for z in np.unique(run_charge):
        mask = run_charge == z
        idx = np.where(mask)[0]
        order = np.argsort(run_mass[idx])
        sorted_idx = idx[order]
        run_index[int(z)] = (run_mass[sorted_idx], run_rt_norm[sorted_idx], sorted_idx)

    # ref_rt_norm -> (run_rt_norm, best_ppm, ref_idx, run_idx)
    best_per_ref: dict[float, tuple[float, float, int, int]] = {}

    for i in range(len(ref_df)):
        z = int(ref_charge[i])
        if z not in run_index:
            continue

        sorted_masses, sorted_rts, sorted_idx = run_index[z]

        m = ref_mass[i]
        delta = m * mass_ppm * 1e-6
        lo = int(np.searchsorted(sorted_masses, m - delta))
        hi = int(np.searchsorted(sorted_masses, m + delta, side="right"))
        if lo >= hi:
            continue

        rrt = ref_rt_norm[i]
        rt_mask = np.abs(rrt - sorted_rts[lo:hi]) <= rt_norm_window
        if not rt_mask.any():
            continue

        cands_idx = sorted_idx[lo:hi][rt_mask]
        ppms = np.abs(m - run_mass[cands_idx]) / m * 1e6
        best = int(np.argmin(ppms))
        best_ppm = float(ppms[best])
        best_run_idx = int(cands_idx[best])

        prev = best_per_ref.get(rrt)
        if prev is None or best_ppm < prev[1]:
            best_per_ref[rrt] = (float(run_rt_norm[best_run_idx]), best_ppm, i, best_run_idx)

    if not best_per_ref:
        empty = np.array([], dtype=np.float64)
        return [], empty, empty, empty

    pairs = [(rrt, v[0]) for rrt, v in best_per_ref.items()]
    ri = np.array([v[2] for v in best_per_ref.values()], dtype=int)
    rj = np.array([v[3] for v in best_per_ref.values()], dtype=int)

    ref_rt_abs_out = ref_rt_abs_arr[ri]
    ppm_errors = (run_mass[rj] - ref_mass[ri]) / ref_mass[ri] * 1e6
    im_deltas = np.where(
        (ref_im[ri] != 0.0) & (run_im[rj] != 0.0),
        run_im[rj] - ref_im[ri],
        0.0,
    )

    return pairs, ref_rt_abs_out, ppm_errors, im_deltas


def _sliding_window_medians(
    run_norms: np.ndarray,
    deltas: np.ndarray,
    window_width: float,
    step: float,
    min_points: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Sliding-window median of RT deltas in normalised RT space.

    Returns (centers, medians) for windows containing >= min_points anchors.
    Windows advance by `step` (50% overlap when step == window_width / 2).
    """
    centers: list[float] = []
    medians: list[float] = []
    lo = 0.0
    while lo < 1.0:
        mask = (run_norms >= lo) & (run_norms < lo + window_width)
        if mask.sum() >= min_points:
            centers.append(lo + window_width / 2.0)
            medians.append(float(np.median(deltas[mask])))
        lo += step
    return np.array(centers), np.array(medians)


def _sigma_clip_mask(residuals: np.ndarray, sigma_clip: float) -> np.ndarray:
    """Inlier mask using a lower-quartile sigma estimate.

    Estimating sigma from the 25th percentile of |residuals| is resistant up to
    75% contamination — as long as true anchors make up >25% of the set, the
    threshold correctly separates them from scatter.
    """
    abs_res = np.abs(residuals - np.median(residuals))
    p25 = float(np.percentile(abs_res, 25))
    # For a Gaussian, the 25th percentile of |X - median| ≈ 0.3186 * sigma
    robust_sigma = p25 / 0.3186 if p25 > 1e-10 else float(np.percentile(abs_res, 75)) / 1.4826
    if robust_sigma < 1e-10:
        return np.ones(len(residuals), dtype=bool)
    return abs_res <= sigma_clip * robust_sigma


def _fit_linear_drift(
    rt: np.ndarray,
    drift: np.ndarray,
    sigma_clip: float = 3.0,
    min_points: int = 4,
) -> tuple[Callable[[np.ndarray], np.ndarray], float, float]:
    """Fit a sigma-clipped linear drift model (drift vs absolute RT).

    Returns (fn, slope, intercept) where fn(rt_abs) -> drift estimate.
    When there are too few points the function returns zeros and slope/intercept
    are both 0.0.
    """
    valid = np.isfinite(rt) & np.isfinite(drift)
    rt_v, drift_v = rt[valid], drift[valid]
    if len(rt_v) < min_points:
        return (lambda x: np.zeros_like(np.asarray(x, dtype=float))), 0.0, 0.0
    coef = np.polyfit(rt_v, drift_v, 1)
    residuals = drift_v - np.polyval(coef, rt_v)
    mask = _sigma_clip_mask(residuals, sigma_clip)
    if mask.sum() >= min_points:
        coef = np.polyfit(rt_v[mask], drift_v[mask], 1)
    slope, intercept = float(coef[0]), float(coef[1])
    return (lambda x: np.polyval(coef, np.asarray(x, dtype=float))), slope, intercept


def _sorted_unique(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort by x and average y for any duplicate x values."""
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    unique_x, inverse = np.unique(xs, return_inverse=True)
    unique_y = np.bincount(inverse, weights=ys) / np.bincount(inverse)
    return unique_x, unique_y


def _fit_rt_warp(
    anchor_pairs: list[tuple[float, float]],
    ref_rt_min: float,
    ref_rt_max: float,
    run_rt_min: float,
    run_rt_max: float,
    window_width: float = 0.05,
    sigma_clip: float = 3.0,
    clip_iters: int = 5,
    spline_lam: float | None = None,
) -> tuple[Callable[[np.ndarray], np.ndarray], np.ndarray, np.ndarray]:
    """Build a warp function from sliding-window medians with iterative pruning.

    Steps:
    1. Compute per-window median RT delta (robust seed, no Gaussian weighting).
    2. Fit a smoothing spline through (window_center, median) points.
    3. Sigma-clip raw anchor pairs against the spline residuals.
    4. Recompute window medians on survivors and refit — repeat until stable.

    This combines the interpretability of plain sliding-window medians with the
    iterative outlier pruning that cleans up the final spline.
    """
    run_norms = np.array([p[1] for p in anchor_pairs])
    deltas    = np.array([p[0] - p[1] for p in anchor_pairs])

    ref_span = ref_rt_max - ref_rt_min
    run_span = run_rt_max - run_rt_min or 1.0

    def _fit_spline(norms: np.ndarray, ds: np.ndarray):
        centers, medians = _sliding_window_medians(norms, ds, window_width, step=window_width / 2)
        if len(centers) < 2:
            med = float(np.median(ds))
            return (lambda x: np.full_like(np.asarray(x, dtype=float), med)), centers, medians
        cx, cy = _sorted_unique(centers, medians)
        return make_smoothing_spline(cx, cy, lam=spline_lam), centers, medians

    spl, centers, medians = _fit_spline(run_norms, deltas)
    mask = _sigma_clip_mask(deltas - spl(run_norms), sigma_clip)

    for _ in range(clip_iters):
        if mask.sum() < 4:
            break
        spl, centers, medians = _fit_spline(run_norms[mask], deltas[mask])
        residuals = deltas[mask] - spl(run_norms[mask])
        inlier = _sigma_clip_mask(residuals, sigma_clip)
        new_mask = mask.copy()
        new_mask[np.where(mask)[0][~inlier]] = False
        if new_mask.sum() == mask.sum():
            break
        mask = new_mask

    def warp(rt_abs: np.ndarray) -> np.ndarray:
        rt_norm = (np.asarray(rt_abs, dtype=float) - run_rt_min) / run_span
        return (rt_norm + spl(rt_norm)) * ref_span + ref_rt_min

    return warp, centers, medians


def _build_warp(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    mass_ppm: float,
    rt_anchor_window: float,
    min_anchor_score: float,
    min_anchor_count: int,
    bandwidth: float = 0.05,
    sigma_clip: float = 3.0,
    clip_iters: int = 5,
    spline_lam: float | None = None,
) -> tuple[Callable, Callable, Callable, dict]:
    """Return (rt_warp_fn, mass_corr_fn, im_corr_fn, diag).

    mass_corr_fn(rt_abs) -> ppm shift to subtract from run features.
    im_corr_fn(rt_abs)   -> IM delta to subtract from run features.
    """
    _ref_rt = ref_df["rt_apex"].cast(pl.Float64).to_numpy()
    _run_rt = run_df["rt_apex"].cast(pl.Float64).to_numpy()
    ref_rt_min, ref_rt_max = float(_ref_rt.min()), float(_ref_rt.max())
    run_rt_min, run_rt_max = float(_run_rt.min()), float(_run_rt.max())

    ref_anchors = _select_anchors(ref_df, min_anchor_score)
    run_anchors = _select_anchors(run_df, min_anchor_score)

    pairs, ref_rt_abs, ppm_errors, im_deltas = _match_anchors(
        ref_anchors, run_anchors, mass_ppm, rt_anchor_window,
        ref_rt_min, ref_rt_max, run_rt_min, run_rt_max,
    )

    def _zero(rt: np.ndarray) -> np.ndarray:
        return np.zeros_like(np.asarray(rt, dtype=float))

    diag: dict = {
        "anchor_pairs": pairs,
        "n_anchors": len(pairs),
        "bandwidth": bandwidth,
        "run_rt_min": run_rt_min,
        "run_rt_max": run_rt_max,
        "ref_rt_min": ref_rt_min,
        "ref_rt_max": ref_rt_max,
        "identity_fallback": False,
        # Raw per-anchor drift data for plotting
        "anchor_ref_rt_abs": ref_rt_abs.tolist(),
        "anchor_ppm_errors": ppm_errors.tolist(),
        "anchor_im_deltas":  im_deltas.tolist(),
    }

    if len(pairs) < min_anchor_count:
        warnings.warn(
            f"Only {len(pairs)} anchor pairs found (need {min_anchor_count}); "
            "falling back to identity RT warp (no correction).",
            RuntimeWarning,
            stacklevel=3,
        )
        diag["identity_fallback"] = True
        return lambda rt: rt, _zero, _zero, diag

    warp_fn, window_centers, window_medians = _fit_rt_warp(
        pairs, ref_rt_min, ref_rt_max, run_rt_min, run_rt_max,
        window_width=bandwidth,
        sigma_clip=sigma_clip,
        clip_iters=clip_iters,
        spline_lam=spline_lam,
    )
    diag["warp_fn"] = warp_fn
    diag["window_centers"] = window_centers.tolist()
    diag["window_medians"] = window_medians.tolist()

    # Mass PPM drift: linear fit of (run_mass - ref_mass)/ref_mass*1e6 vs ref RT
    mass_corr_fn, mass_slope, mass_intercept = _fit_linear_drift(ref_rt_abs, ppm_errors)
    diag["mass_drift_slope_ppm_per_min"] = mass_slope
    diag["mass_drift_intercept_ppm"] = mass_intercept

    # IM drift: only fit when anchor pairs have non-zero IM on both sides
    im_valid = im_deltas != 0.0
    if im_valid.sum() >= 4:
        im_corr_fn, im_slope, im_intercept = _fit_linear_drift(
            ref_rt_abs[im_valid], im_deltas[im_valid]
        )
        diag["im_drift_slope_per_min"] = im_slope
        diag["im_drift_intercept"] = im_intercept
    else:
        im_corr_fn = _zero
        diag["im_drift_slope_per_min"] = 0.0
        diag["im_drift_intercept"] = 0.0

    return warp_fn, mass_corr_fn, im_corr_fn, diag


def _match_features(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    warp_fn: Callable[[np.ndarray], np.ndarray],
    mass_ppm: float,
    rt_window: float,
    im_tolerance: float,
    intensity_col: str,
    mass_corr_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    im_corr_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> tuple[pl.DataFrame, dict]:
    """Return (DataFrame with ref feature rows + warped intensity, match stats dict).

    Uses a sort + binary-search strategy: O(n log n) time, O(n+m) memory.
    No intermediate pair DataFrame is ever materialised.
    """
    n_ref = len(ref_df)
    n_run = len(run_df)

    # Pull numpy arrays for all columns we need — zero-copy via Arrow
    run_rt_abs    = run_df["rt_apex"].cast(pl.Float64).to_numpy()
    run_rt_warped = warp_fn(run_rt_abs)
    run_mass_raw  = run_df["mass"].cast(pl.Float64).to_numpy()
    run_charge    = run_df["charge"].to_numpy()
    run_im_raw    = run_df["im"].cast(pl.Float64).to_numpy()
    run_intensity = run_df[intensity_col].cast(pl.Float64).to_numpy()

    # Apply linear drift corrections derived from anchor pairs
    if mass_corr_fn is not None:
        # ppm_shift = (run - ref) / ref * 1e6  →  corrected = run / (1 + shift*1e-6)
        run_mass = run_mass_raw / (1.0 + mass_corr_fn(run_rt_abs) * 1e-6)
    else:
        run_mass = run_mass_raw
    if im_corr_fn is not None:
        run_im = np.where(run_im_raw != 0.0, run_im_raw - im_corr_fn(run_rt_abs), 0.0)
    else:
        run_im = run_im_raw

    ref_mass   = ref_df["mass"].cast(pl.Float64).to_numpy()
    ref_charge = ref_df["charge"].to_numpy()
    ref_rt     = ref_df["rt_apex"].cast(pl.Float64).to_numpy()
    ref_im     = ref_df["im"].cast(pl.Float64).to_numpy()

    # Build per-charge sorted index into the run array
    run_index: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for z in np.unique(run_charge):
        mask = run_charge == z
        idx = np.where(mask)[0]
        order = np.argsort(run_mass[idx])
        sorted_idx = idx[order]
        run_index[int(z)] = (run_mass[sorted_idx], sorted_idx)

    # Bhattacharyya score for tie-breaking (secondary sort after mass PPM)
    run_score = run_df["score"].cast(pl.Float64).to_numpy() if "score" in run_df.columns else np.zeros(n_run)
    ref_score = ref_df["score"].cast(pl.Float64).to_numpy() if "score" in ref_df.columns else np.zeros(n_ref)

    matched_intensity = np.zeros(n_ref, dtype=np.float64)
    runnerup_gaps = np.full(n_ref, np.inf)
    n_after_charge_mass_rt = 0
    n_after_im = 0

    # Greedy 1:1 matching: process ref features highest-score first so the most
    # confident features claim their best run counterpart before lower-scoring
    # ref features get a chance. Once a run feature is claimed it cannot be
    # matched again, preventing double-counting.
    claimed_run = np.zeros(n_run, dtype=bool)
    ref_order = np.argsort(-ref_score)

    for i in ref_order:
        z = int(ref_charge[i])
        if z not in run_index:
            continue

        sorted_masses, sorted_idx = run_index[z]

        # Binary search: find all run features within the mass PPM window
        m = ref_mass[i]
        delta = m * mass_ppm * 1e-6
        lo = int(np.searchsorted(sorted_masses, m - delta))
        hi = int(np.searchsorted(sorted_masses, m + delta, side="right"))
        if lo >= hi:
            continue

        cands = sorted_idx[lo:hi]

        # RT filter (warped run RT vs ref RT)
        cands = cands[np.abs(ref_rt[i] - run_rt_warped[cands]) <= rt_window]
        if len(cands) == 0:
            continue
        n_after_charge_mass_rt += 1

        # IM filter — skip when either side has no ion mobility (im == 0)
        rim = ref_im[i]
        cands = cands[
            (rim == 0.0)
            | (run_im[cands] == 0.0)
            | (np.abs(rim - run_im[cands]) <= im_tolerance)
        ]
        if len(cands) == 0:
            continue
        n_after_im += 1

        # Exclude run features already claimed by a higher-scoring ref feature
        cands = cands[~claimed_run[cands]]
        if len(cands) == 0:
            continue

        # Best candidate: primary sort by mass PPM, Bhattacharyya score breaks sub-ppm ties
        ppms = np.abs(m - run_mass[cands]) / m * 1e6
        sort_key = ppms - 1e-6 * run_score[cands]
        best_local = int(np.argmin(sort_key))
        best_run_idx = cands[best_local]
        claimed_run[best_run_idx] = True
        matched_intensity[i] = run_intensity[best_run_idx]

        # Runner-up gap: distance from best to second-best candidate (ppm)
        if len(cands) >= 2:
            sorted_ppms = np.sort(ppms)
            runnerup_gaps[i] = sorted_ppms[1] - sorted_ppms[0]

    n_matched = int((matched_intensity > 0).sum())
    result = ref_df.select(["mass", "charge", "rt_apex", "im"]).with_columns(
        pl.Series(intensity_col, matched_intensity, dtype=pl.Float64)
    )

    stats = {
        "n_ref": n_ref,
        "n_run": n_run,
        "after_charge_mass_rt": n_after_charge_mass_rt,
        "after_im": n_after_im,
        "after_dedup": n_matched,
        "n_matched": n_matched,
        "runnerup_gaps": runnerup_gaps,
    }
    return result, stats


def _make_decoy_ref(ref_df: pl.DataFrame, mass_shift_da: float = 100.0) -> pl.DataFrame:
    """Shift all reference masses by a fixed offset to create a null-model database.

    100 Da is ~10,000× larger than any real ppm tolerance window, ensuring
    decoy features cannot coincidentally match real run features.
    """
    return ref_df.with_columns([
        (pl.col("mass") + mass_shift_da).alias("mass"),
        (pl.col("mz") + mass_shift_da / pl.col("charge").cast(pl.Float64)).alias("mz"),
    ])


def estimate_match_fdr(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    warp_fn: Callable[[np.ndarray], np.ndarray],
    mass_ppm: float,
    rt_window: float,
    im_tolerance: float,
    intensity_col: str,
    mass_shift_da: float = 100.0,
) -> dict:
    """Estimate false discovery rate for feature matching using a decoy reference.

    Runs matching twice with identical tolerances: once against the real reference
    (target) and once against a mass-shifted reference (decoy). Because the decoy
    is 100 Da away from any real feature, every decoy match is a false positive —
    their count directly estimates the expected FP rate.

    FDR = n_decoy / n_target  (the standard non-competition formula; no factor of 2
    because target and decoy are matched in separate passes, not head-to-head).

    Returns
    -------
    dict with keys: n_target, n_decoy, fdr, mass_shift_da
    """
    _, target_stats = _match_features(ref_df, run_df, warp_fn, mass_ppm, rt_window, im_tolerance, intensity_col)
    decoy_ref = _make_decoy_ref(ref_df, mass_shift_da)
    _, decoy_stats = _match_features(decoy_ref, run_df, warp_fn, mass_ppm, rt_window, im_tolerance, intensity_col)

    n_t = max(target_stats["n_matched"], 1)
    n_d = decoy_stats["n_matched"]
    return {
        "n_target": target_stats["n_matched"],
        "n_decoy": n_d,
        "fdr": n_d / n_t,
        "mass_shift_da": mass_shift_da,
    }


def sweep_match_parameters(
    ref_df: pl.DataFrame,
    run_df: pl.DataFrame,
    warp_fn: Callable[[np.ndarray], np.ndarray],
    intensity_col: str,
    mass_ppm_values: list[float] | None = None,
    rt_window_values: list[float] | None = None,
    im_tolerance_values: list[float] | None = None,
    fdr_threshold: float = 0.01,
    mass_shift_da: float = 100.0,
) -> pl.DataFrame:
    """Grid search over (mass_ppm, rt_window, im_tolerance) combinations.

    For each combination runs both target and decoy matching to compute FDR.
    Returns a DataFrame sorted by passes_fdr DESC, n_target DESC, fdr ASC.
    The first row where passes_fdr=True is the recommended parameter set.
    """
    ppm_grid = mass_ppm_values or [5.0, 10.0, 15.0, 20.0, 30.0]
    rt_grid  = rt_window_values or [0.25, 0.5, 0.75, 1.0, 1.5]
    im_grid  = im_tolerance_values or [0.03, 0.05, 0.1, 0.5, 1.0]

    decoy_ref = _make_decoy_ref(ref_df, mass_shift_da)
    rows = []
    for ppm in ppm_grid:
        for rt in rt_grid:
            for im in im_grid:
                _, ts = _match_features(ref_df,   run_df, warp_fn, ppm, rt, im, intensity_col)
                _, ds = _match_features(decoy_ref, run_df, warp_fn, ppm, rt, im, intensity_col)
                n_t = max(ts["n_matched"], 1)
                n_d = ds["n_matched"]
                fdr = n_d / n_t
                rows.append({
                    "mass_ppm": ppm,
                    "rt_window": rt,
                    "im_tolerance": im,
                    "n_target": ts["n_matched"],
                    "n_decoy": n_d,
                    "fdr": fdr,
                    "passes_fdr": fdr <= fdr_threshold,
                })

    return (
        pl.DataFrame(rows)
        .sort(["passes_fdr", "n_target", "fdr"], descending=[True, True, False])
    )


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
