#!/usr/bin/env python3
"""
Interactive feature viewer for koth_ff output files.

Accepts either:
  • a features file (`features.tsv` / `features.parquet`)
  • a folder containing `features.{tsv,parquet}` AND `hills.{tsv,parquet}`

If a hills file is found, the elution panel overlays each constituent hill
(M, M+1, M+2, …) as its own coloured trace so co-elution can be assessed.
Otherwise it falls back to the summed `elution_profile` written into the
features file.

Usage:
    python scripts/feature_dashboard.py path/to/run_output/
    python scripts/feature_dashboard.py path/to/features.parquet --hills path/to/hills.parquet
    python scripts/feature_dashboard.py path/to/features.tsv --port 8052
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, dash_table, dcc, html, no_update
from plotly.subplots import make_subplots


# ──────────────────────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────────────────────

def load_features(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t")
    df = df.reset_index(drop=True)
    # `feature_id` is a synthetic row id used by the dashboard so callbacks can
    # map table rows back to the underlying DataFrame regardless of sorting.
    df["feature_id"] = df.index.astype("int64")
    # Pre-compute |ppm_error| for filtering convenience.
    if "ppm_error" in df.columns:
        df["abs_ppm_error"] = df["ppm_error"].abs()
    else:
        df["abs_ppm_error"] = 0.0
    return df


def load_hills(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t")
    df = df.reset_index(drop=True)
    if "hill_id" not in df.columns:
        df["hill_id"] = df.index.astype("int64")
    return df


def resolve_inputs(input_path: Path, hills_override: Path | None) -> tuple[Path, Path | None]:
    """Resolve a user-supplied path into (features_path, hills_path_or_None).

    - If `input_path` is a directory, look for features.{tsv,parquet} and
      hills.{tsv,parquet} inside.
    - If `input_path` is a file, use it as features; check for a sibling
      hills.{tsv,parquet} unless `hills_override` is given.
    """
    if input_path.is_dir():
        features = next(
            (input_path / name for name in ("features.parquet", "features.tsv")
             if (input_path / name).exists()),
            None,
        )
        if features is None:
            raise FileNotFoundError(
                f"No features.parquet or features.tsv found in {input_path}"
            )
        if hills_override is not None:
            hills = hills_override
        else:
            hills = next(
                (input_path / name for name in ("hills.parquet", "hills.tsv")
                 if (input_path / name).exists()),
                None,
            )
        return features, hills

    if hills_override is not None:
        return input_path, hills_override

    parent = input_path.parent
    hills = next(
        (parent / name for name in ("hills.parquet", "hills.tsv")
         if (parent / name).exists()),
        None,
    )
    return input_path, hills


# ──────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────────────────────────────────────

def _parse_json(raw) -> list:
    if raw is None:
        return []
    if isinstance(raw, float) and np.isnan(raw):
        return []
    if isinstance(raw, str):
        if not raw or raw == "[]":
            return []
        return json.loads(raw)
    if isinstance(raw, (list, tuple, np.ndarray)):
        return list(raw)
    return []


def score_histogram_figure(filtered: pd.DataFrame) -> go.Figure:
    """Three side-by-side histograms: isotope, chromato-cosine, combined."""
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("isotope_score", "cosine_score", "combined_score"),
        horizontal_spacing=0.07,
    )
    if len(filtered) == 0:
        fig.update_layout(
            height=260,
            margin=dict(l=40, r=20, t=40, b=40),
            annotations=[dict(
                text="No features match current filters",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=13, color="gray"),
            )],
        )
        return fig

    isotope = filtered["isotope_score"].to_numpy(dtype=np.float64)
    cosine = filtered["cosine_score"].to_numpy(dtype=np.float64)
    combined = filtered["combined_score"].to_numpy(dtype=np.float64)

    for col_idx, (data, color) in enumerate(
        [(isotope, "#1f77b4"), (cosine, "#2ca02c"), (combined, "#9467bd")],
        start=1,
    ):
        data = data[np.isfinite(data)]
        fig.add_trace(
            go.Histogram(
                x=data,
                marker_color=color,
                opacity=0.85,
                nbinsx=50,
                showlegend=False,
            ),
            row=1, col=col_idx,
        )

    fig.update_yaxes(title_text="count", row=1, col=1)
    fig.update_xaxes(title_text="isotope_score", row=1, col=1)
    fig.update_xaxes(title_text="cosine_score", row=1, col=2)
    fig.update_xaxes(title_text="combined_score", row=1, col=3)
    fig.update_layout(
        height=260,
        margin=dict(l=50, r=20, t=40, b=45),
        bargap=0.05,
    )
    return fig


def empty_figure(message: str, height: int = 560) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        annotations=[
            dict(
                text=message,
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=14, color="gray"),
            )
        ],
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        height=height,
    )
    return fig


ISOTOPE_PALETTE = [
    "#1f77b4",  # blue
    "#d62728",  # red
    "#2ca02c",  # green
    "#ff7f0e",  # orange
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # cyan
    "#bcbd22",  # olive
    "#7f7f7f",  # gray
]


def _isotope_color(i: int) -> str:
    return ISOTOPE_PALETTE[i % len(ISOTOPE_PALETTE)]


def _hill_anchors(hill_row: pd.Series) -> list[tuple[int, float]]:
    """(absolute_scan, rt) anchors contributed by one hill.

    Each hill provides up to three anchors: scan_start↔rt_start, scan_apex↔rt
    (apex), and scan_end↔rt_end. These are reliable points on the
    monotonically-increasing scan→RT mapping; we collect them across all
    constituent hills to build a single shared interpolation.
    """
    n = len(_parse_json(hill_row.get("intensity_profile")))
    if n == 0:
        return []
    anchors: list[tuple[int, float]] = []
    s_start = int(hill_row.get("scan_start", 0))
    s_end = int(hill_row.get("scan_end", s_start + n - 1))
    s_apex = int(hill_row.get("scan_apex", s_start))
    rt_start = float(hill_row.get("rt_start", 0.0))
    rt_end = float(hill_row.get("rt_end", 0.0))
    rt_apex = float(hill_row.get("rt", 0.0))
    if np.isfinite(rt_start):
        anchors.append((s_start, rt_start))
    if np.isfinite(rt_end) and s_end != s_start:
        anchors.append((s_end, rt_end))
    # rt may be 0.0 for older outputs that hit the apex-on-gap bug; skip those.
    if rt_apex > 0.0 and np.isfinite(rt_apex):
        anchors.append((s_apex, rt_apex))
    return anchors


def _build_scan_to_rt(hill_rows: list[pd.Series]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a single piecewise-linear scan→RT map shared by all constituent hills.

    Returns (all_scans, all_rts, _) where:
      • all_scans  = arange(min_scan_start, max_scan_end + 1) across all hills
      • all_rts    = scan_to_rt evaluated on all_scans, via np.interp on the
                     sorted/deduped (scan, rt) anchors from every hill

    Aligning by absolute scan number — instead of interpolating each hill in
    its own [rt_start, rt_end] window — guarantees that all isotope traces
    share an identical scan→RT mapping, which is the only honest way to
    visualise co-elution.
    """
    anchor_pairs: list[tuple[int, float]] = []
    scan_lo = None
    scan_hi = None
    for hr in hill_rows:
        if hr is None:
            continue
        n = len(_parse_json(hr.get("intensity_profile")))
        if n == 0:
            continue
        s_start = int(hr.get("scan_start", 0))
        s_end = int(hr.get("scan_end", s_start + n - 1))
        scan_lo = s_start if scan_lo is None else min(scan_lo, s_start)
        scan_hi = s_end if scan_hi is None else max(scan_hi, s_end)
        anchor_pairs.extend(_hill_anchors(hr))

    if scan_lo is None or scan_hi is None or not anchor_pairs:
        return np.array([], dtype=np.int64), np.array([]), np.array([])

    # Sort and dedup by scan (mean of RTs at duplicate scans).
    anchors_df = (
        pd.DataFrame(anchor_pairs, columns=["scan", "rt"])
        .groupby("scan", as_index=False)["rt"]
        .mean()
        .sort_values("scan")
    )
    anchor_scans = anchors_df["scan"].to_numpy(dtype=np.float64)
    anchor_rts = anchors_df["rt"].to_numpy(dtype=np.float64)

    all_scans = np.arange(int(scan_lo), int(scan_hi) + 1, dtype=np.int64)
    # np.interp clamps at the endpoints, so scans outside the anchor span get
    # the boundary RT (acceptable here: the global scan range is bounded by
    # anchor scans by construction, equal at the ends).
    all_rts = np.interp(all_scans.astype(np.float64), anchor_scans, anchor_rts)
    return all_scans, all_rts, anchor_scans


def feature_figure(row: pd.Series, hills_by_id: dict[int, pd.Series] | None = None) -> go.Figure:
    summed_elution = np.asarray(_parse_json(row.get("elution_profile")), dtype=np.float64)
    isotope = np.asarray(_parse_json(row.get("isotope_profile")), dtype=np.float64)
    theoretical = np.asarray(_parse_json(row.get("theoretical_pattern")), dtype=np.float64)
    hill_ids = [int(h) for h in _parse_json(row.get("hill_ids"))]

    have_per_hill = bool(hills_by_id) and bool(hill_ids)

    if len(summed_elution) == 0 and not have_per_hill and len(isotope) == 0 and len(theoretical) == 0:
        return empty_figure("Selected feature has no profile data.")

    rt_apex = float(row.get("rtApex", 0.0))

    elution_title = (
        "Per-isotope elution (M, M+1, M+2, …)"
        if have_per_hill
        else "Elution profile (linear, summed across isotopes)"
    )
    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=False,
        row_heights=[0.45, 0.25, 0.30],
        vertical_spacing=0.10,
        subplot_titles=(
            elution_title,
            "Elution profile (log₁₀)",
            "Isotope envelope (observed vs theoretical)",
        ),
    )

    # ── Elution panels ──
    if have_per_hill:
        # Step 1: build a single scan→RT map from the union of all hill
        # anchors. Aligning by absolute scan number ensures every isotope
        # trace lives on the same x-axis — required for honest co-elution
        # visualisation. We do this once per feature, not per-hill.
        hill_rows = [hills_by_id.get(hid) for hid in hill_ids]
        all_scans, all_rts, _ = _build_scan_to_rt(hill_rows)
        have_per_hill = all_scans.size > 0

    if have_per_hill:
        scan_lo = int(all_scans[0])
        # Step 2: each hill's intensity profile is placed on the shared scan
        # grid by absolute scan index; positions outside the hill's own range
        # are NaN so they don't draw spurious zero-segments.
        for iso_idx, (hid, hill_row) in enumerate(zip(hill_ids, hill_rows)):
            if hill_row is None:
                continue
            profile = np.asarray(_parse_json(hill_row.get("intensity_profile")), dtype=np.float64)
            n = len(profile)
            if n == 0:
                continue
            s_start = int(hill_row.get("scan_start", 0))
            offset = s_start - scan_lo
            if offset < 0 or offset + n > len(all_scans):
                # Defensive: shouldn't happen because all_scans spans every
                # hill's [scan_start, scan_end].
                continue

            hill_rts = all_rts[offset : offset + n]
            y_lin = np.full(all_scans.shape, np.nan, dtype=np.float64)
            y_lin[offset : offset + n] = profile

            color = _isotope_color(iso_idx)
            mz_val = float(hill_row.get("mz", 0.0))
            label = f"M+{iso_idx}  m/z {mz_val:.4f}  (hill {hid})"

            # Linear trace — connectgaps=False so NaN slots leave a real gap.
            fig.add_trace(
                go.Scatter(
                    x=all_rts, y=y_lin,
                    mode="lines+markers",
                    line=dict(color=color, width=2),
                    marker=dict(size=5, color=color),
                    name=label,
                    legendgroup=label,
                    connectgaps=False,
                    hovertemplate=(
                        f"M+{iso_idx}  m/z {mz_val:.4f}<br>"
                        "rt %{x:.4f} min<br>"
                        "intensity %{y:.3e}<extra></extra>"
                    ),
                ),
                row=1, col=1,
            )

            with np.errstate(divide="ignore"):
                y_log = np.where(profile > 0, np.log10(profile), np.nan)
            y_log_full = np.full(all_scans.shape, np.nan, dtype=np.float64)
            y_log_full[offset : offset + n] = y_log
            fig.add_trace(
                go.Scatter(
                    x=all_rts, y=y_log_full,
                    mode="lines+markers",
                    line=dict(color=color, width=2),
                    marker=dict(size=4, color=color),
                    name=label,
                    legendgroup=label,
                    showlegend=False,
                    connectgaps=False,
                    hovertemplate=(
                        f"M+{iso_idx}  m/z {mz_val:.4f}<br>"
                        "rt %{x:.4f} min<br>"
                        "log₁₀ %{y:.3f}<extra></extra>"
                    ),
                ),
                row=2, col=1,
            )

        fig.add_vline(
            x=rt_apex,
            line=dict(color="tomato", dash="dash", width=1.5),
            annotation_text=f"apex {rt_apex:.3f} min",
            annotation_position="top right",
            row=1, col=1,
        )
        fig.add_vline(
            x=rt_apex,
            line=dict(color="tomato", dash="dash", width=1.0),
            row=2, col=1,
        )
        # Pin x-range to the shared scan span (with small padding).
        rt_lo = float(all_rts[0])
        rt_hi = float(all_rts[-1])
        if rt_hi > rt_lo:
            pad = (rt_hi - rt_lo) * 0.03
            fig.update_xaxes(range=[rt_lo - pad, rt_hi + pad], row=1, col=1)
            fig.update_xaxes(range=[rt_lo - pad, rt_hi + pad], row=2, col=1)
    else:
        # Fallback: summed elution from features.tsv.
        n = len(summed_elution)
        rt_start = float(row.get("rtStart", 0.0))
        rt_end = float(row.get("rtEnd", 0.0))
        if rt_end > rt_start and n > 1:
            rts = np.linspace(rt_start, rt_end, n)
        else:
            rts = np.full(n, rt_apex)
        if n > 0:
            gap_mask = summed_elution == 0.0
            fig.add_trace(
                go.Scatter(
                    x=rts, y=summed_elution,
                    mode="lines+markers",
                    line=dict(color="steelblue", width=2),
                    marker=dict(
                        size=6,
                        color=["lightcoral" if g else "steelblue" for g in gap_mask],
                        line=dict(width=0),
                    ),
                    fill="tozeroy",
                    fillcolor="rgba(70,130,180,0.15)",
                    name="intensity (summed)",
                    hovertemplate=(
                        "rt %{x:.4f} min<br>"
                        "intensity %{y:.3e}<extra></extra>"
                    ),
                ),
                row=1, col=1,
            )
            with np.errstate(divide="ignore"):
                log_elution = np.where(summed_elution > 0, np.log10(summed_elution), np.nan)
            fig.add_trace(
                go.Scatter(
                    x=rts, y=log_elution,
                    mode="lines+markers",
                    line=dict(color="steelblue", width=2),
                    marker=dict(size=5),
                    fill="tozeroy",
                    fillcolor="rgba(70,130,180,0.15)",
                    name="log10(intensity)",
                    showlegend=False,
                    hovertemplate="rt %{x:.4f} min<br>log₁₀ %{y:.3f}<extra></extra>",
                ),
                row=2, col=1,
            )
            fig.add_vline(
                x=rt_apex,
                line=dict(color="tomato", dash="dash", width=1.5),
                annotation_text=f"apex {rt_apex:.3f} min",
                annotation_position="top right",
                row=1, col=1,
            )
            fig.add_vline(
                x=rt_apex,
                line=dict(color="tomato", dash="dash", width=1.0),
                row=2, col=1,
            )
            if rt_end > rt_start:
                pad = (rt_end - rt_start) * 0.03
                fig.update_xaxes(range=[rt_start - pad, rt_end + pad], row=1, col=1)
                fig.update_xaxes(range=[rt_start - pad, rt_end + pad], row=2, col=1)

    # ── Isotope envelope ──
    if len(isotope) > 0 or len(theoretical) > 0:
        obs_sum = float(isotope.sum()) if len(isotope) else 0.0
        obs_norm = isotope / obs_sum if obs_sum > 0 else isotope
        n_bars = max(len(obs_norm), len(theoretical))
        x_bars = np.arange(n_bars)

        if len(obs_norm) > 0:
            obs_colors = [_isotope_color(i) for i in range(len(obs_norm))]
            fig.add_trace(
                go.Bar(
                    x=x_bars[: len(obs_norm)],
                    y=obs_norm,
                    name="Observed",
                    marker_color=obs_colors,
                    opacity=0.85,
                    offsetgroup="obs",
                    showlegend=False,
                    hovertemplate="M+%{x}<br>rel %{y:.3f}<extra>obs</extra>",
                ),
                row=3, col=1,
            )
        if len(theoretical) > 0:
            fig.add_trace(
                go.Bar(
                    x=x_bars[: len(theoretical)],
                    y=theoretical,
                    name="Theoretical",
                    marker_color="#555555",
                    opacity=0.55,
                    offsetgroup="theo",
                    showlegend=False,
                    hovertemplate="M+%{x}<br>rel %{y:.3f}<extra>theo</extra>",
                ),
                row=3, col=1,
            )

    # ── Layout ──
    charge = int(row["charge"]) if not pd.isna(row.get("charge")) else 0
    charge_str = f"{charge}+" if charge > 0 else "?"
    mass_calib = row.get("massCalib")
    title = (
        f"Feature #{int(row['feature_id'])} — "
        f"m/z {row['mz']:.5f}  ·  z={charge_str}  ·  RT {rt_apex:.3f} min"
    )
    if pd.notna(mass_calib):
        title += f"  ·  mass {float(mass_calib):.4f} Da"
    if float(row.get("im", 0.0) or 0.0) != 0.0:
        title += f"  ·  IM {float(row['im']):.4f}"

    fig.update_layout(
        title=dict(text=title, font=dict(size=15)),
        height=720,
        margin=dict(l=60, r=30, t=80, b=50),
        showlegend=True,
        barmode="group",
        bargap=0.25,
        legend=dict(orientation="h", yanchor="bottom", y=-0.06, x=0.0),
    )
    fig.update_xaxes(title_text="Retention time (min)", row=2, col=1)
    fig.update_yaxes(title_text="Intensity", row=1, col=1, tickformat=".2e")
    fig.update_yaxes(title_text="log₁₀(intensity)", row=2, col=1)
    fig.update_xaxes(title_text="Isotope (M, M+1, …)", row=3, col=1,
                     tickmode="linear", tick0=0, dtick=1)
    fig.update_yaxes(title_text="Relative intensity", row=3, col=1)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Layout
# ──────────────────────────────────────────────────────────────────────────────

TABLE_COLUMNS = [
    {"name": "id", "id": "feature_id", "type": "numeric"},
    {"name": "m/z", "id": "mz", "type": "numeric", "format": {"specifier": ".5f"}},
    {"name": "mass", "id": "massCalib", "type": "numeric", "format": {"specifier": ".4f"}},
    {"name": "z", "id": "charge", "type": "numeric"},
    {"name": "RT (min)", "id": "rtApex", "type": "numeric", "format": {"specifier": ".4f"}},
    {"name": "n_iso", "id": "nIsotopes", "type": "numeric"},
    {"name": "n_scans", "id": "nScans", "type": "numeric"},
    {"name": "Σ int", "id": "intensitySum", "type": "numeric", "format": {"specifier": ".2e"}},
    {"name": "apex int", "id": "intensityApex", "type": "numeric", "format": {"specifier": ".2e"}},
    {"name": "iso", "id": "isotope_score", "type": "numeric", "format": {"specifier": ".3f"}},
    {"name": "cos", "id": "cosine_score", "type": "numeric", "format": {"specifier": ".3f"}},
    {"name": "comb", "id": "combined_score", "type": "numeric", "format": {"specifier": ".3f"}},
    {"name": "ppm", "id": "ppm_error", "type": "numeric", "format": {"specifier": ".2f"}},
    {"name": "Δn", "id": "neutron_offset", "type": "numeric"},
]


def build_app(
    df: pd.DataFrame,
    source_label: str,
    hills_by_id: dict[int, pd.Series] | None = None,
) -> Dash:
    mz_lo, mz_hi = float(df["mz"].min()), float(df["mz"].max())
    rt_lo, rt_hi = float(df["rtApex"].min()), float(df["rtApex"].max())
    score_lo, score_hi = float(df["combined_score"].min()), float(df["combined_score"].max())
    cos_lo, cos_hi = float(df["cosine_score"].min()), float(df["cosine_score"].max())
    ppm_max = float(df["abs_ppm_error"].max()) if "abs_ppm_error" in df else 50.0
    charges = sorted({int(c) for c in df["charge"].dropna().unique() if int(c) > 0})

    summary = (
        f"{len(df):,} features · m/z {mz_lo:.2f}–{mz_hi:.2f} · "
        f"rt {rt_lo:.2f}–{rt_hi:.2f} min · charges {','.join(str(c) for c in charges)}"
    )
    if hills_by_id is not None:
        summary += f" · {len(hills_by_id):,} hills loaded (per-isotope overlay)"
    else:
        summary += " · hills not loaded (summed elution only)"

    app = Dash(__name__, title="koth_ff feature viewer")

    filter_row = [
        html.Div([
            html.Label("m/z range", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.RangeSlider(
                id="mz-slider",
                min=mz_lo, max=mz_hi,
                value=[mz_lo, mz_hi],
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": False},
                marks=None,
            ),
        ], style={"flex": 2, "padding": "0 12px"}),
        html.Div([
            html.Label("RT range (min)", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.RangeSlider(
                id="rt-slider",
                min=rt_lo, max=rt_hi,
                value=[rt_lo, rt_hi],
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": False},
                marks=None,
            ),
        ], style={"flex": 2, "padding": "0 12px"}),
        html.Div([
            html.Label("charge", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Dropdown(
                id="charge-filter",
                options=[{"label": f"{c}+", "value": c} for c in charges],
                value=charges,
                multi=True,
                placeholder="All charges",
            ),
        ], style={"flex": 1, "padding": "0 12px", "minWidth": "140px"}),
        html.Div([
            html.Label("min n_iso", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="min-iso",
                type="number",
                min=1, step=1, value=1,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("min n_scans", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="min-scans",
                type="number",
                min=1, step=1, value=1,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("min combined", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="min-score",
                type="number",
                min=score_lo, max=score_hi, step=0.01,
                value=score_lo,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("min cosine", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="min-cosine",
                type="number",
                min=cos_lo, max=cos_hi, step=0.01,
                value=cos_lo,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("max |ppm|", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="max-ppm",
                type="number",
                min=0.0, step=0.1,
                value=ppm_max,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("Sort by", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Dropdown(
                id="sort-by",
                options=[
                    {"label": "intensitySum (desc)", "value": "intensitySum"},
                    {"label": "intensityApex (desc)", "value": "intensityApex"},
                    {"label": "combined_score (desc)", "value": "combined_score"},
                    {"label": "isotope_score (desc)", "value": "isotope_score"},
                    {"label": "cosine_score (desc)", "value": "cosine_score"},
                    {"label": "nIsotopes (desc)", "value": "nIsotopes"},
                    {"label": "nScans (desc)", "value": "nScans"},
                    {"label": "m/z (asc)", "value": "mz"},
                    {"label": "rtApex (asc)", "value": "rtApex"},
                    {"label": "charge (asc)", "value": "charge"},
                    {"label": "|ppm| (asc)", "value": "abs_ppm_error"},
                ],
                value="intensitySum",
                clearable=False,
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
    ]

    app.layout = html.Div([
        html.Div([
            html.H2("koth_ff feature viewer", style={"margin": "0", "fontFamily": "system-ui, sans-serif"}),
            html.Div(source_label, style={"fontFamily": "monospace", "fontSize": 12, "color": "#555"}),
            html.Div(summary, style={"fontSize": 13, "color": "#333", "marginTop": 6}),
        ], style={"padding": "10px 20px", "borderBottom": "1px solid #ddd", "background": "#fafafa"}),

        html.Div(
            filter_row,
            style={
                "display": "flex",
                "alignItems": "center",
                "flexWrap": "wrap",
                "padding": "10px 8px",
                "borderBottom": "1px solid #eee",
                "rowGap": "6px",
            },
        ),

        html.Div([
            dcc.Graph(
                id="score-histograms",
                figure=score_histogram_figure(df),
                config={"displayModeBar": False},
            ),
        ], style={"padding": "4px 10px 0", "borderBottom": "1px solid #eee"}),

        html.Div([
            html.Div([
                html.Div(id="row-count", style={"fontSize": 12, "color": "#666", "padding": "4px 4px 8px"}),
                dash_table.DataTable(
                    id="feature-table",
                    columns=TABLE_COLUMNS,
                    data=[],
                    page_size=25,
                    sort_action="native",
                    row_selectable="single",
                    selected_rows=[],
                    style_table={"height": "720px", "overflowY": "auto"},
                    style_cell={"fontFamily": "monospace", "fontSize": 12, "padding": "3px 6px"},
                    style_header={"fontWeight": "bold", "background": "#f3f3f3"},
                    style_data_conditional=[
                        {"if": {"state": "selected"}, "backgroundColor": "#dbeafe", "border": "1px solid #2563eb"},
                    ],
                ),
            ], style={"flex": 1, "padding": "10px", "minWidth": "560px"}),

            html.Div([
                dcc.Graph(id="feature-plot", figure=empty_figure("Select a feature from the table")),
                html.Div(id="feature-metadata", style={
                    "fontFamily": "monospace",
                    "fontSize": 12,
                    "padding": "10px",
                    "background": "#fafafa",
                    "border": "1px solid #eee",
                    "borderRadius": 6,
                    "marginTop": 8,
                    "whiteSpace": "pre",
                }),
            ], style={"flex": 1, "padding": "10px", "minWidth": "560px"}),
        ], style={"display": "flex", "alignItems": "stretch", "flexWrap": "wrap"}),
    ])

    # ── Callbacks ────────────────────────────────────────────────────────────
    @app.callback(
        Output("feature-table", "data"),
        Output("feature-table", "selected_rows"),
        Output("row-count", "children"),
        Output("score-histograms", "figure"),
        Input("mz-slider", "value"),
        Input("rt-slider", "value"),
        Input("charge-filter", "value"),
        Input("min-iso", "value"),
        Input("min-scans", "value"),
        Input("min-score", "value"),
        Input("min-cosine", "value"),
        Input("max-ppm", "value"),
        Input("sort-by", "value"),
    )
    def update_table(mz_range, rt_range, charge_sel, min_iso, min_scans,
                     min_score, min_cos, max_ppm, sort_by):
        mask = (
            df["mz"].between(mz_range[0], mz_range[1])
            & df["rtApex"].between(rt_range[0], rt_range[1])
            & (df["nIsotopes"] >= (min_iso or 1))
            & (df["nScans"] >= (min_scans or 1))
            & (df["combined_score"] >= (min_score if min_score is not None else -np.inf))
            & (df["cosine_score"] >= (min_cos if min_cos is not None else -np.inf))
            & (df["abs_ppm_error"] <= (max_ppm if max_ppm is not None else np.inf))
        )
        if charge_sel:
            mask &= df["charge"].isin(charge_sel)

        filtered = df.loc[mask]
        ascending = sort_by in ("mz", "rtApex", "charge", "abs_ppm_error")
        filtered = filtered.sort_values(sort_by, ascending=ascending)

        cap = 5000
        truncated = len(filtered) > cap
        view = filtered.head(cap)

        records = view[[c["id"] for c in TABLE_COLUMNS]].to_dict("records")

        msg = f"{len(filtered):,} features match filters"
        if truncated:
            msg += f" — showing first {cap:,}"
        hist_fig = score_histogram_figure(filtered)
        return records, [], msg, hist_fig

    @app.callback(
        Output("feature-plot", "figure"),
        Output("feature-metadata", "children"),
        Input("feature-table", "selected_rows"),
        State("feature-table", "data"),
    )
    def update_plot(selected_rows, table_data):
        if not selected_rows or not table_data:
            return empty_figure("Select a feature from the table"), ""
        idx = selected_rows[0]
        if idx >= len(table_data):
            return no_update, no_update
        feature_id = table_data[idx]["feature_id"]
        row = df.loc[df["feature_id"] == feature_id].iloc[0]
        fig = feature_figure(row, hills_by_id)
        meta = format_metadata(row)
        return fig, meta

    return app


def _fmt(val, spec: str, default: str = "—") -> str:
    if val is None:
        return default
    try:
        if pd.isna(val):
            return default
    except (TypeError, ValueError):
        pass
    try:
        return format(val, spec)
    except (TypeError, ValueError):
        return str(val)


def format_metadata(row: pd.Series) -> str:
    charge = int(row["charge"]) if not pd.isna(row.get("charge")) else 0
    hill_ids = _parse_json(row.get("hill_ids"))
    hill_ids_str = ", ".join(str(int(h)) for h in hill_ids) if hill_ids else "—"

    lines = [
        f"feature_id     {int(row['feature_id'])}",
        f"m/z            {_fmt(row.get('mz'), '.6f')}",
        f"massCalib      {_fmt(row.get('massCalib'), '.6f')} Da",
        f"charge         {charge}+",
        f"rt apex        {_fmt(row.get('rtApex'), '.4f')} min",
        f"rt range       {_fmt(row.get('rtStart'), '.4f')} – {_fmt(row.get('rtEnd'), '.4f')} min",
        f"n_isotopes     {_fmt(row.get('nIsotopes'), 'd')}",
        f"n_scans        {_fmt(row.get('nScans'), 'd')}",
        f"intensity sum  {_fmt(row.get('intensitySum'), '.4e')}",
        f"intensity apex {_fmt(row.get('intensityApex'), '.4e')}",
        f"isotope_score  {_fmt(row.get('isotope_score'), '.4f')}",
        f"cosine_score   {_fmt(row.get('cosine_score'), '.4f')}",
        f"combined_score {_fmt(row.get('combined_score'), '.4f')}",
        f"ppm error      {_fmt(row.get('ppm_error'), '.3f')}",
        f"neutron offset {_fmt(row.get('neutron_offset'), 'd')}",
    ]
    if float(row.get("im", 0.0) or 0.0) != 0.0:
        lines.append(f"ion mobility   {_fmt(row.get('im'), '.4f')}")
    lines.append(f"hill_ids       {hill_ids_str}")
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# Entrypoint
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "input",
        type=Path,
        help=(
            "Path to a features file (.tsv/.parquet) OR a folder containing "
            "features.{tsv,parquet} and (optionally) hills.{tsv,parquet}."
        ),
    )
    parser.add_argument(
        "--hills",
        type=Path,
        default=None,
        help="Override path to the hills file (otherwise auto-discovered).",
    )
    parser.add_argument("--port", type=int, default=8052, help="Port to serve the dashboard on (default 8052)")
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind to (default 127.0.0.1)")
    parser.add_argument("--debug", action="store_true", help="Enable Dash debug mode")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Path not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if args.hills is not None and not args.hills.exists():
        print(f"Hills file not found: {args.hills}", file=sys.stderr)
        sys.exit(1)

    features_path, hills_path = resolve_inputs(args.input, args.hills)

    print(f"Loading features from {features_path} …", file=sys.stderr)
    df = load_features(features_path)
    print(f"Loaded {len(df):,} features", file=sys.stderr)

    hills_by_id: dict[int, pd.Series] | None = None
    if hills_path is not None:
        print(f"Loading hills from {hills_path} …", file=sys.stderr)
        hills_df = load_hills(hills_path)
        # Indexing by hill_id gives O(1) row lookup; we keep the rows as Series.
        hills_by_id = {
            int(hid): row for hid, row in zip(hills_df["hill_id"], hills_df.to_dict("records"))
        }
        # `to_dict("records")` returns plain dicts; wrap them in pd.Series so
        # downstream code can use `.get()` consistently.
        hills_by_id = {hid: pd.Series(rec) for hid, rec in hills_by_id.items()}
        print(f"Loaded {len(hills_by_id):,} hills", file=sys.stderr)
    else:
        print("No hills file found — falling back to summed elution from features.", file=sys.stderr)

    source_label = str(features_path.resolve())
    if hills_path is not None:
        source_label += f"   +   {hills_path.resolve()}"

    app = build_app(df, source_label=source_label, hills_by_id=hills_by_id)
    print(f"Serving on http://{args.host}:{args.port}", file=sys.stderr)
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
