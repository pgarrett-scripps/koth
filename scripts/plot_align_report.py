"""HTML diagnostics report from koth_align's `align_report.json`.

Reads the JSON written by the Rust `koth_align` binary and produces a
self-contained interactive HTML file with:
  * Top-level summary cards (runs, consensus features, timing, reference)
  * Per-run alignment table (anchor counts + residual MAD before/after)
  * Per-run 4-panel diagnostic figure:
      - RT mapping: run_rt vs ref_rt + warp curve + identity line
      - RT correction: ref_rt - run_rt at anchors + fitted warp delta
      - Mass drift: ppm_error vs ref_rt_norm + linear fit
      - IM drift:   im_delta vs ref_rt_norm + linear fit
  * Residual distributions (before vs after) for ppm / IM / RT
  * Cross-run detection rate + q-value percentile chart

Usage:
    python scripts/plot_align_report.py <align_output_dir_or_json>
    python scripts/plot_align_report.py <path/to/align_report.json> -o report.html
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ── small helpers ────────────────────────────────────────────────────────────

ACTIVE_COLOR   = "#1f77b4"   # blue
INACTIVE_COLOR = "#d62728"   # red
FIT_COLOR      = "#2ca02c"   # green
REF_COLOR      = "#7f7f7f"   # grey


def _fmt_count(n: int) -> str:
    return f"{n:,}"


def _fmt_float(x: float, digits: int = 3) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "—"
    return f"{x:.{digits}f}"


def _fmt_pct(x: float) -> str:
    if x is None or math.isnan(x):
        return "—"
    return f"{x * 100:.1f}%"


def _warp_curve(knot_x: list[float], knot_y: list[float], n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Sample the piecewise-linear warp delta(run_norm) on [0, 1]."""
    xs = np.linspace(0.0, 1.0, n)
    ys = np.interp(xs, knot_x, knot_y)
    return xs, ys


def _has_im(anchors: list[dict]) -> bool:
    return any(a.get("ref_im", 0.0) != 0.0 and a.get("run_im", 0.0) != 0.0 for a in anchors)


# ── per-run 4-panel figure ───────────────────────────────────────────────────

def _run_alignment_figure(run: dict) -> go.Figure | None:
    align = run.get("alignment")
    if align is None:
        return None

    name = run["name"]
    anchors = align["anchors"]
    if not anchors:
        return None

    run_lo, run_hi = run["rt_range"]
    run_span = run_hi - run_lo or 1.0

    # Reference RT range is implicit in the absolute ref_rt values stored in anchors.
    ref_rt_arr = np.array([a["ref_rt"] for a in anchors])
    ref_lo = float(ref_rt_arr.min()) if ref_rt_arr.size else 0.0
    ref_hi = float(ref_rt_arr.max()) if ref_rt_arr.size else 1.0
    ref_span = ref_hi - ref_lo or 1.0

    has_im_data = _has_im(anchors)

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "RT mapping (run → reference)",
            "RT correction (ref_rt − run_rt)",
            "Mass drift vs normalised reference RT",
            "IM drift vs normalised reference RT" if has_im_data else "IM drift — no IM anchors",
        ),
        horizontal_spacing=0.10, vertical_spacing=0.13,
    )

    # Split anchors by each active mask. Build arrays once.
    arr_run_rt   = np.array([a["run_rt"]     for a in anchors])
    arr_ref_rt   = np.array([a["ref_rt"]     for a in anchors])
    arr_run_norm = np.array([a["run_rt_norm"] for a in anchors])
    arr_ref_norm = np.array([a["ref_rt_norm"] for a in anchors])
    arr_ppm      = np.array([a["ppm_error"]  for a in anchors])
    arr_im_d     = np.array([a["im_delta"]   for a in anchors])
    arr_ref_im   = np.array([a["ref_im"]     for a in anchors])
    arr_run_im   = np.array([a["run_im"]     for a in anchors])
    arr_ref_mz   = np.array([a["ref_mz"]     for a in anchors])
    arr_run_mz   = np.array([a["run_mz"]     for a in anchors])

    rt_active   = np.array([a["rt_active"]   for a in anchors], dtype=bool)
    mass_active = np.array([a["mass_active"] for a in anchors], dtype=bool)
    im_active   = np.array([a["im_active"]   for a in anchors], dtype=bool)

    # Anchor-level hover info that's the same across panels.
    hover_customdata = np.stack([arr_ref_mz, arr_ref_im, arr_run_rt, arr_ref_rt], axis=-1)

    # ── Panel (1,1): RT mapping ──
    hover_map = (
        "run_rt: %{x:.3f}<br>ref_rt: %{y:.3f}<br>"
        "ref_mz: %{customdata[0]:.4f}<br>ref_im: %{customdata[1]:.4f}<extra></extra>"
    )
    if rt_active.any():
        fig.add_trace(
            go.Scattergl(
                x=arr_run_rt[rt_active], y=arr_ref_rt[rt_active],
                mode="markers",
                marker=dict(size=4, color=ACTIVE_COLOR, opacity=0.55),
                name="active anchor", legendgroup="active",
                customdata=hover_customdata[rt_active],
                hovertemplate=hover_map,
            ),
            row=1, col=1,
        )
    if (~rt_active).any():
        fig.add_trace(
            go.Scattergl(
                x=arr_run_rt[~rt_active], y=arr_ref_rt[~rt_active],
                mode="markers",
                marker=dict(size=4, color=INACTIVE_COLOR, opacity=0.55, symbol="x"),
                name="clipped", legendgroup="clipped",
                customdata=hover_customdata[~rt_active],
                hovertemplate=hover_map,
            ),
            row=1, col=1,
        )

    # Warp curve in absolute minutes: x = run_lo + run_norm * run_span,
    # y = ref_lo + (run_norm + delta) * ref_span.
    knots = align["rt_warp"]
    xs_n, deltas = _warp_curve(knots["knot_x"], knots["knot_y"], n=400)
    warp_run_abs = run_lo + xs_n * run_span
    warp_ref_abs = ref_lo + (xs_n + deltas) * ref_span
    fig.add_trace(
        go.Scatter(
            x=warp_run_abs, y=warp_ref_abs, mode="lines",
            line=dict(color=FIT_COLOR, width=2),
            name="fitted warp", legendgroup="fit",
            hovertemplate="warp(run=%{x:.3f}) → ref=%{y:.3f}<extra></extra>",
        ),
        row=1, col=1,
    )

    # Identity line.
    id_lo = min(run_lo, ref_lo)
    id_hi = max(run_hi, ref_hi)
    fig.add_trace(
        go.Scatter(
            x=[id_lo, id_hi], y=[id_lo, id_hi], mode="lines",
            line=dict(color=REF_COLOR, width=1, dash="dash"),
            name="identity", legendgroup="identity",
            hoverinfo="skip",
        ),
        row=1, col=1,
    )

    # ── Panel (1,2): RT correction in minutes ──
    delta_abs = arr_ref_rt - arr_run_rt
    hover_d = (
        "run_rt: %{x:.3f}<br>Δ (ref−run): %{y:.3f} min<br>"
        "ref_mz: %{customdata[0]:.4f}<extra></extra>"
    )
    if rt_active.any():
        fig.add_trace(
            go.Scattergl(
                x=arr_run_rt[rt_active], y=delta_abs[rt_active],
                mode="markers",
                marker=dict(size=4, color=ACTIVE_COLOR, opacity=0.55),
                showlegend=False, legendgroup="active",
                customdata=hover_customdata[rt_active],
                hovertemplate=hover_d,
            ),
            row=1, col=2,
        )
    if (~rt_active).any():
        fig.add_trace(
            go.Scattergl(
                x=arr_run_rt[~rt_active], y=delta_abs[~rt_active],
                mode="markers",
                marker=dict(size=4, color=INACTIVE_COLOR, opacity=0.55, symbol="x"),
                showlegend=False, legendgroup="clipped",
                customdata=hover_customdata[~rt_active],
                hovertemplate=hover_d,
            ),
            row=1, col=2,
        )
    # Fitted correction curve in minutes:
    #   ref_abs - run_abs = (ref_lo + (xs_n + delta) * ref_span) - (run_lo + xs_n * run_span)
    correction_abs = (ref_lo + (xs_n + deltas) * ref_span) - (run_lo + xs_n * run_span)
    fig.add_trace(
        go.Scatter(
            x=warp_run_abs, y=correction_abs, mode="lines",
            line=dict(color=FIT_COLOR, width=2),
            showlegend=False, legendgroup="fit",
            hovertemplate="run=%{x:.3f}<br>correction=%{y:.3f} min<extra></extra>",
        ),
        row=1, col=2,
    )
    fig.add_hline(y=0.0, line=dict(color=REF_COLOR, dash="dash", width=1),
                  row=1, col=2)

    # ── Panel (2,1): Mass drift ──
    hover_m = (
        "ref_rt_norm: %{x:.3f}<br>ppm_error: %{y:.3f}<br>"
        "ref_mz: %{customdata[0]:.4f}<extra></extra>"
    )
    if mass_active.any():
        fig.add_trace(
            go.Scattergl(
                x=arr_ref_norm[mass_active], y=arr_ppm[mass_active],
                mode="markers",
                marker=dict(size=4, color=ACTIVE_COLOR, opacity=0.55),
                showlegend=False, legendgroup="active",
                customdata=hover_customdata[mass_active],
                hovertemplate=hover_m,
            ),
            row=2, col=1,
        )
    if (~mass_active).any():
        fig.add_trace(
            go.Scattergl(
                x=arr_ref_norm[~mass_active], y=arr_ppm[~mass_active],
                mode="markers",
                marker=dict(size=4, color=INACTIVE_COLOR, opacity=0.55, symbol="x"),
                showlegend=False, legendgroup="clipped",
                customdata=hover_customdata[~mass_active],
                hovertemplate=hover_m,
            ),
            row=2, col=1,
        )
    mass_fit = align["mass_drift"]
    xs_line = np.linspace(0.0, 1.0, 100)
    ys_line = mass_fit["intercept"] + mass_fit["slope"] * xs_line
    fig.add_trace(
        go.Scatter(
            x=xs_line, y=ys_line, mode="lines",
            line=dict(color=FIT_COLOR, width=2),
            showlegend=False, legendgroup="fit",
            hovertemplate=(
                f"fit: {mass_fit['intercept']:+.3f} + {mass_fit['slope']:+.3f}·x ppm"
                "<extra></extra>"
            ),
        ),
        row=2, col=1,
    )
    fig.add_hline(y=0.0, line=dict(color=REF_COLOR, dash="dash", width=1),
                  row=2, col=1)

    # ── Panel (2,2): IM drift ──
    if has_im_data:
        im_mask = (arr_ref_im != 0.0) & (arr_run_im != 0.0)
        active_in_mask = im_active & im_mask
        clipped_in_mask = (~im_active) & im_mask
        hover_im = (
            "ref_rt_norm: %{x:.3f}<br>im_delta: %{y:.4f}<br>"
            "ref_im: %{customdata[1]:.4f}<extra></extra>"
        )
        if active_in_mask.any():
            fig.add_trace(
                go.Scattergl(
                    x=arr_ref_norm[active_in_mask], y=arr_im_d[active_in_mask],
                    mode="markers",
                    marker=dict(size=4, color=ACTIVE_COLOR, opacity=0.55),
                    showlegend=False, legendgroup="active",
                    customdata=hover_customdata[active_in_mask],
                    hovertemplate=hover_im,
                ),
                row=2, col=2,
            )
        if clipped_in_mask.any():
            fig.add_trace(
                go.Scattergl(
                    x=arr_ref_norm[clipped_in_mask], y=arr_im_d[clipped_in_mask],
                    mode="markers",
                    marker=dict(size=4, color=INACTIVE_COLOR, opacity=0.55, symbol="x"),
                    showlegend=False, legendgroup="clipped",
                    customdata=hover_customdata[clipped_in_mask],
                    hovertemplate=hover_im,
                ),
                row=2, col=2,
            )
        im_fit = align["im_drift"]
        ys_im_line = im_fit["intercept"] + im_fit["slope"] * xs_line
        fig.add_trace(
            go.Scatter(
                x=xs_line, y=ys_im_line, mode="lines",
                line=dict(color=FIT_COLOR, width=2),
                showlegend=False, legendgroup="fit",
                hovertemplate=(
                    f"fit: {im_fit['intercept']:+.5f} + {im_fit['slope']:+.5f}·x"
                    "<extra></extra>"
                ),
            ),
            row=2, col=2,
        )
        fig.add_hline(y=0.0, line=dict(color=REF_COLOR, dash="dash", width=1),
                      row=2, col=2)

    # Axis labels.
    fig.update_xaxes(title_text="run RT (min)", row=1, col=1)
    fig.update_yaxes(title_text="reference RT (min)", row=1, col=1)
    fig.update_xaxes(title_text="run RT (min)", row=1, col=2)
    fig.update_yaxes(title_text="Δ RT (min)", row=1, col=2)
    fig.update_xaxes(title_text="ref RT (normalised)", row=2, col=1)
    fig.update_yaxes(title_text="ppm error", row=2, col=1)
    fig.update_xaxes(title_text="ref RT (normalised)", row=2, col=2)
    fig.update_yaxes(title_text="Δ 1/K₀", row=2, col=2)

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{name}</b>  ·  "
                f"{align['n_anchors']:,} anchors  ·  "
                f"RT active {align['n_anchors_rt_active']:,}  ·  "
                f"mass active {align['n_anchors_mass_active']:,}  ·  "
                f"IM active {align['n_anchors_im_active']:,}"
            ),
            font=dict(size=14),
        ),
        height=720,
        margin=dict(l=60, r=20, t=80, b=50),
        legend=dict(orientation="h", x=0.5, y=-0.08, xanchor="center"),
        plot_bgcolor="#fafafa",
    )

    return fig


# ── residual histograms ──────────────────────────────────────────────────────

def _residual_histogram_figure(run: dict) -> go.Figure | None:
    align = run.get("alignment")
    if align is None:
        return None

    anchors = align["anchors"]
    if not anchors:
        return None

    ppm_before = np.array([a["ppm_error"]                for a in anchors])
    ppm_after  = np.array([a["ppm_residual_after_fit"]   for a in anchors])
    rt_before  = np.array([a["ref_rt_norm"] - a["run_rt_norm"] for a in anchors])
    rt_after   = np.array([a["rt_residual_after_warp"]   for a in anchors])

    # IM residuals can be NaN where IM wasn't available.
    im_mask = np.array([(a["ref_im"] != 0.0 and a["run_im"] != 0.0) for a in anchors])
    im_before = np.array([a["im_delta"]                for a in anchors])[im_mask]
    im_after_all = np.array([a["im_residual_after_fit"] for a in anchors], dtype=float)
    im_after  = im_after_all[im_mask & np.isfinite(im_after_all)]

    has_im = im_mask.any()

    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("ppm residual", "Δ RT (normalised) residual",
                        "Δ IM residual" if has_im else "Δ IM — no IM anchors"),
        horizontal_spacing=0.08,
    )

    def _hist_pair(arr_before, arr_after, col, name_before, name_after):
        if arr_before.size:
            fig.add_trace(
                go.Histogram(
                    x=arr_before, name=name_before, marker_color=REF_COLOR,
                    opacity=0.55, nbinsx=60, legendgroup=name_before,
                ),
                row=1, col=col,
            )
        if arr_after.size:
            fig.add_trace(
                go.Histogram(
                    x=arr_after, name=name_after, marker_color=FIT_COLOR,
                    opacity=0.7, nbinsx=60, legendgroup=name_after,
                ),
                row=1, col=col,
            )

    _hist_pair(ppm_before, ppm_after, 1, "before", "after")
    _hist_pair(rt_before,  rt_after,  2, "before", "after")
    if has_im:
        _hist_pair(im_before, im_after, 3, "before", "after")

    fig.update_layout(
        barmode="overlay",
        title=dict(text=f"<b>{run['name']}</b> — residual distributions",
                   font=dict(size=14)),
        height=320,
        margin=dict(l=50, r=20, t=70, b=40),
        legend=dict(orientation="h", x=0.5, y=-0.18, xanchor="center"),
        plot_bgcolor="#fafafa",
    )
    fig.update_yaxes(title_text="count", row=1, col=1)
    return fig


# ── cross-run summary ────────────────────────────────────────────────────────

def _detection_rate_figure(runs: list[dict]) -> go.Figure:
    names = [r["name"] for r in runs]
    rates = [r["lfq"]["detection_rate"] * 100.0 for r in runs]
    n_det = [r["lfq"]["n_detected"] for r in runs]
    ref_mask = [r["is_reference"] for r in runs]

    colors = [FIT_COLOR if ref else ACTIVE_COLOR for ref in ref_mask]
    text = [f"{n_det[i]:,}" for i in range(len(names))]

    fig = go.Figure(
        data=[
            go.Bar(
                x=names, y=rates, text=text, textposition="outside",
                marker_color=colors,
                hovertemplate="%{x}<br>%{y:.1f}% detected<br>n=%{text}<extra></extra>",
            )
        ]
    )
    fig.update_layout(
        title=dict(text="Detection rate per run (green = reference)",
                   font=dict(size=14)),
        yaxis=dict(title="% detected", range=[0, 105]),
        xaxis=dict(title="run", tickangle=-30),
        height=380,
        margin=dict(l=50, r=20, t=60, b=120),
        plot_bgcolor="#fafafa",
    )
    return fig


def _qvalue_percentile_figure(runs: list[dict]) -> go.Figure:
    names = [r["name"] for r in runs]
    p25 = [r["lfq"]["qvalue_percentiles"]["p25"] for r in runs]
    p50 = [r["lfq"]["qvalue_percentiles"]["p50"] for r in runs]
    p75 = [r["lfq"]["qvalue_percentiles"]["p75"] for r in runs]
    p95 = [r["lfq"]["qvalue_percentiles"]["p95"] for r in runs]

    fig = go.Figure()
    for label, series, color in [
        ("p25", p25, "#74add1"),
        ("p50", p50, "#1f77b4"),
        ("p75", p75, "#f46d43"),
        ("p95", p95, "#d62728"),
    ]:
        fig.add_trace(go.Scatter(
            x=names, y=series, mode="lines+markers", name=label,
            line=dict(color=color, width=2), marker=dict(size=8),
        ))
    fig.update_layout(
        title=dict(text="Per-run q-value percentiles (detected features only)",
                   font=dict(size=14)),
        yaxis=dict(title="q-value", type="log", autorange=True),
        xaxis=dict(title="run", tickangle=-30),
        height=380,
        margin=dict(l=60, r=20, t=60, b=120),
        plot_bgcolor="#fafafa",
        legend=dict(orientation="h", x=0.5, y=-0.25, xanchor="center"),
    )
    return fig


# ── HTML assembly ────────────────────────────────────────────────────────────

CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       max-width: 1280px; margin: 32px auto; padding: 0 28px; color: #222;
       line-height: 1.5; background: #ffffff; }
h1 { color: #1a4a72; border-bottom: 3px solid #1a4a72; padding-bottom: 10px;
     margin-bottom: 8px; }
h2 { color: #1a4a72; margin-top: 44px; font-size: 1.2em;
     border-bottom: 1px solid #e0e0e0; padding-bottom: 4px; }
h3 { color: #244c70; font-size: 1.05em; margin-top: 28px; }
p.meta { color: #666; font-size: 13px; margin-top: -8px; }
.stat-grid { display: grid;
             grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
             gap: 16px; margin: 24px 0; }
.stat-box { background: #eef5fc; border-left: 5px solid #1a4a72;
            padding: 14px 18px; border-radius: 6px; }
.stat-box .val { font-size: 28px; font-weight: 700; color: #1a4a72; }
.stat-box .lbl { font-size: 12px; color: #555; margin-top: 2px; }
table { border-collapse: collapse; width: 100%; font-size: 13px;
        margin-top: 10px; }
th, td { padding: 7px 12px; text-align: left; border-bottom: 1px solid #ddd; }
th { background: #1a4a72; color: white; font-weight: 600; }
tr:nth-child(even) { background: #f7f9fc; }
.section { margin: 28px 0; }
details { margin: 20px 0; }
details summary { cursor: pointer; font-weight: 600; color: #1a4a72; }
pre.cfg { background: #f4f6fa; padding: 12px 16px; border-radius: 6px;
          font-size: 12px; overflow: auto; max-height: 320px; }
footer { margin-top: 60px; font-size: 11px; color: #aaa;
         border-top: 1px solid #eee; padding-top: 12px; }
"""


def _summary_cards(report: dict) -> str:
    runs = report["runs"]
    ref = report["reference_run"]
    timing = report["timing"]

    aligned_runs = [r for r in runs if not r["is_reference"]]
    total_anchors = sum(
        r["alignment"]["n_anchors"] for r in aligned_runs if r.get("alignment")
    )
    mean_detect = (
        sum(r["lfq"]["detection_rate"] for r in runs) / max(1, len(runs))
    )

    cards = [
        ("Runs aligned", _fmt_count(report["n_runs"])),
        ("Consensus features", _fmt_count(report["n_consensus_features"])),
        ("Reference run", f"<span style='font-size:14px'>{ref}</span>"),
        ("Total anchors", _fmt_count(total_anchors)),
        ("Mean detection rate", _fmt_pct(mean_detect)),
        ("Alignment time", f"{timing['alignment_sec']:.1f}s"),
        ("LFQ time", f"{timing['lfq_sec']:.1f}s"),
        ("Total time", f"{timing['total_sec']:.1f}s"),
    ]
    html = '<div class="stat-grid">'
    for lbl, val in cards:
        html += f'<div class="stat-box"><div class="val">{val}</div>'
        html += f'<div class="lbl">{lbl}</div></div>'
    html += "</div>"
    return html


def _per_run_table(runs: list[dict]) -> str:
    rows = []
    for r in runs:
        a = r.get("alignment")
        if a is None:
            n_anc = "(reference)"
            n_act = "—"
            ppm_b = ppm_a = im_b = im_a = rt_b = rt_a = "—"
        else:
            n_anc = _fmt_count(a["n_anchors"])
            n_act = (
                f"rt {a['n_anchors_rt_active']:,} · "
                f"mz {a['n_anchors_mass_active']:,} · "
                f"im {a['n_anchors_im_active']:,}"
            )
            rs = a["residual_stats"]
            ppm_b = _fmt_float(rs["ppm_before"]["mad"], 3)
            ppm_a = _fmt_float(rs["ppm_after"]["mad"],  3)
            im_b  = _fmt_float(rs["im_before"]["mad"],  5)
            im_a  = _fmt_float(rs["im_after"]["mad"],   5)
            rt_b  = _fmt_float(rs["rt_norm_before"]["mad"], 4)
            rt_a  = _fmt_float(rs["rt_norm_after"]["mad"],  4)

        rows.append(
            f"<tr>"
            f"<td>{r['name']}</td>"
            f"<td>{_fmt_count(r['n_features'])}</td>"
            f"<td>{_fmt_count(r['n_features_high_score'])}</td>"
            f"<td>{_fmt_count(r['n_hills'])}</td>"
            f"<td>{n_anc}</td>"
            f"<td style='font-size:11px'>{n_act}</td>"
            f"<td>{ppm_b} → {ppm_a}</td>"
            f"<td>{im_b} → {im_a}</td>"
            f"<td>{rt_b} → {rt_a}</td>"
            f"<td>{_fmt_count(r['lfq']['n_detected'])}</td>"
            f"<td>{_fmt_pct(r['lfq']['detection_rate'])}</td>"
            f"</tr>"
        )

    return (
        "<table><tr>"
        "<th>Run</th><th>Features</th><th>≥ min_anchor_combined_score</th><th>Hills</th>"
        "<th>Anchors</th><th>Active (rt · mz · im)</th>"
        "<th>ppm MAD (before → after)</th>"
        "<th>IM MAD (before → after)</th>"
        "<th>RTₙ MAD (before → after)</th>"
        "<th>Detected</th><th>%</th>"
        "</tr>" + "".join(rows) + "</table>"
    )


def build_html(report: dict) -> str:
    parts: list[str] = []

    # Header
    parts.append(
        f'<h1>koth_align — diagnostics</h1>'
        f'<p class="meta">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} · '
        f"schema v{report['schema_version']}</p>"
    )

    # Top summary
    parts.append(_summary_cards(report))

    # Per-run summary table
    parts.append('<div class="section"><h2>Per-run summary</h2>')
    parts.append(_per_run_table(report["runs"]))
    parts.append("</div>")

    # Cross-run charts
    parts.append('<div class="section"><h2>Cross-run quality</h2>')
    fig_det = _detection_rate_figure(report["runs"])
    fig_q   = _qvalue_percentile_figure(report["runs"])
    parts.append(fig_det.to_html(full_html=False, include_plotlyjs=False))
    parts.append(fig_q.to_html(full_html=False, include_plotlyjs=False))
    parts.append("</div>")

    # Per-run alignment figures
    parts.append('<div class="section"><h2>Per-run alignment diagnostics</h2>')
    parts.append(
        '<p style="font-size:12px;color:#666">Blue dots: anchors retained by '
        'the sigma-clip fit. Red ×: clipped anchors. Green line: fitted model. '
        'Grey dashed: identity / zero reference. Hover for per-anchor values.</p>'
    )
    for r in report["runs"]:
        if r["is_reference"]:
            continue
        fig = _run_alignment_figure(r)
        if fig is None:
            continue
        parts.append(f'<h3>{r["name"]}</h3>')
        parts.append(fig.to_html(full_html=False, include_plotlyjs=False))

        hist = _residual_histogram_figure(r)
        if hist is not None:
            parts.append(hist.to_html(full_html=False, include_plotlyjs=False))
    parts.append("</div>")

    # Config
    parts.append('<div class="section"><details><summary>Config used</summary>')
    parts.append(f'<pre class="cfg">{json.dumps(report["config"], indent=2)}</pre>')
    parts.append("</details></div>")

    parts.append(
        f"<footer>koth_align report · {_fmt_count(report['n_consensus_features'])} "
        f"consensus features · {report['n_runs']} runs · reference: "
        f"{report['reference_run']}</footer>"
    )

    body = "\n".join(parts)
    # Inline plotly.js once.
    plotlyjs = go.Figure().to_html(
        full_html=False, include_plotlyjs="inline"
    ).split("<div ")[0]

    return (
        '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
        '<title>koth_align report</title>'
        f"<style>{CSS}</style>"
        f"{plotlyjs}"
        f"</head><body>{body}</body></html>"
    )


# ── CLI ──────────────────────────────────────────────────────────────────────

def _resolve_json_path(arg: str) -> Path:
    p = Path(arg)
    if p.is_dir():
        cand = p / "align_report.json"
        if not cand.exists():
            raise SystemExit(
                f"No align_report.json in {p} — did you re-run koth_align "
                f"after the report change?"
            )
        return cand
    if not p.exists():
        raise SystemExit(f"Path not found: {p}")
    return p


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Generate an HTML alignment report from align_report.json"
    )
    ap.add_argument(
        "input",
        help="Path to align_report.json OR to the align_output directory containing it",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output HTML path (default: <input_dir>/align_report.html)",
    )
    args = ap.parse_args(argv)

    json_path = _resolve_json_path(args.input)
    with json_path.open() as fh:
        report: dict[str, Any] = json.load(fh)

    out_path = args.output or json_path.with_name("align_report.html")
    html = build_html(report)
    out_path.write_text(html, encoding="utf-8")
    print(f"Wrote {out_path}  ({len(html) / 1024:.1f} KiB)")


if __name__ == "__main__":
    main()
