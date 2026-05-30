#!/usr/bin/env python3
"""
Interactive hill viewer for koth_ff output files.

Works on both MS1 (`hills.tsv` / `hills.parquet`) and MS2
(`hills_ms2.tsv` / `hills_ms2.parquet`) files. MS2 inputs are auto-detected
from the presence of `iso_target_mz` columns and gain an isolation-window
filter.

Usage:
    python scripts/hill_dashboard.py path/to/hills.tsv
    python scripts/hill_dashboard.py path/to/hills_ms2.parquet --port 8051
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

def load_hills(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t")
    df = df.reset_index(drop=True)
    # `hill_id` is written by koth_ff >= the change that added hill provenance to
    # features.tsv. For older files, fall back to the row index so the dashboard
    # still works.
    if "hill_id" not in df.columns:
        df["hill_id"] = df.index.astype("int64")
    return df


def is_ms2(df: pd.DataFrame) -> bool:
    return "iso_target_mz" in df.columns and df["iso_target_mz"].notna().any()


def isolation_window_options(df: pd.DataFrame) -> list[dict]:
    if not is_ms2(df):
        return []
    iws = (
        df[["iso_target_mz", "iso_lower_mz", "iso_upper_mz"]]
        .dropna()
        .drop_duplicates()
        .sort_values("iso_target_mz")
        .to_dict("records")
    )
    options = [{"label": "All windows", "value": "all"}]
    for iw in iws:
        label = (
            f"{iw['iso_target_mz']:.3f}"
            f"  [{iw['iso_lower_mz']:.3f}–{iw['iso_upper_mz']:.3f}]"
        )
        # Encode as string key for the dropdown
        value = f"{iw['iso_target_mz']:.6f}"
        options.append({"label": label, "value": value})
    return options


# ──────────────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────────────

def _parse_profile(raw) -> np.ndarray:
    if isinstance(raw, str):
        return np.asarray(json.loads(raw), dtype=np.float64)
    if isinstance(raw, (list, tuple)):
        return np.asarray(raw, dtype=np.float64)
    if isinstance(raw, np.ndarray):
        return raw.astype(np.float64)
    return np.array([], dtype=np.float64)


def empty_figure(message: str) -> go.Figure:
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
        height=520,
    )
    return fig


def hill_figure(row: pd.Series) -> go.Figure:
    profile = _parse_profile(row["intensity_profile"])
    n = len(profile)
    if n == 0:
        return empty_figure("Selected hill has an empty intensity profile.")

    scan_start = int(row["scan_start"])
    scan_apex = int(row["scan_apex"])
    scans = np.arange(scan_start, scan_start + n)
    rt_start = float(row["rt_start"])
    rt_end = float(row["rt_end"])
    if rt_end > rt_start and n > 1:
        rts = np.linspace(rt_start, rt_end, n)
    else:
        rts = np.full(n, float(row["rt"]))

    gap_mask = profile == 0.0
    rel_apex = scan_apex - scan_start
    apex_rt = rts[rel_apex] if 0 <= rel_apex < n else float(row["rt"])

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.65, 0.35],
        vertical_spacing=0.07,
        subplot_titles=("Elution profile (linear)", "Elution profile (log₁₀)"),
    )

    # Linear trace
    fig.add_trace(
        go.Scatter(
            x=rts, y=profile,
            mode="lines+markers",
            line=dict(color="steelblue", width=2),
            marker=dict(
                size=6,
                color=["lightcoral" if g else "steelblue" for g in gap_mask],
                line=dict(width=0),
            ),
            fill="tozeroy",
            fillcolor="rgba(70,130,180,0.15)",
            name="intensity",
            customdata=np.stack([scans, profile], axis=-1),
            hovertemplate=(
                "scan %{customdata[0]}<br>"
                "rt %{x:.4f} min<br>"
                "intensity %{customdata[1]:.3e}<extra></extra>"
            ),
        ),
        row=1, col=1,
    )

    # Apex marker
    fig.add_vline(
        x=apex_rt,
        line=dict(color="tomato", dash="dash", width=1.5),
        annotation_text=f"apex (scan {scan_apex})",
        annotation_position="top right",
        row=1, col=1,
    )

    # Log trace
    with np.errstate(divide="ignore"):
        log_profile = np.where(profile > 0, np.log10(profile), np.nan)
    fig.add_trace(
        go.Scatter(
            x=rts, y=log_profile,
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
        x=apex_rt,
        line=dict(color="tomato", dash="dash", width=1.0),
        row=2, col=1,
    )

    title = f"Hill #{int(row['hill_id'])} — m/z {row['mz']:.5f} ± {row['mz_std']:.5f}"
    if "iso_target_mz" in row.index and pd.notna(row["iso_target_mz"]):
        title += (
            f"   |   isolation {row['iso_target_mz']:.3f}"
            f" [{row['iso_lower_mz']:.3f}–{row['iso_upper_mz']:.3f}]"
        )
    if row.get("im", 0.0) and float(row["im"]) != 0.0:
        title += f"   |   IM {float(row['im']):.4f}"

    fig.update_layout(
        title=dict(text=title, font=dict(size=15)),
        height=520,
        margin=dict(l=60, r=30, t=80, b=50),
        showlegend=False,
        hovermode="x unified",
    )
    fig.update_xaxes(title_text="Retention time (min)", row=2, col=1)
    fig.update_yaxes(title_text="Intensity", row=1, col=1, tickformat=".2e")
    fig.update_yaxes(title_text="log₁₀(intensity)", row=2, col=1)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Layout
# ──────────────────────────────────────────────────────────────────────────────

TABLE_COLUMNS_MS1 = [
    {"name": "id", "id": "hill_id", "type": "numeric"},
    {"name": "m/z", "id": "mz", "type": "numeric", "format": {"specifier": ".5f"}},
    {"name": "RT (min)", "id": "rt", "type": "numeric", "format": {"specifier": ".4f"}},
    {"name": "n_scans", "id": "n_scans", "type": "numeric"},
    {"name": "gaps", "id": "skipped_scans", "type": "numeric"},
    {"name": "Σ int", "id": "intensity_sum", "type": "numeric", "format": {"specifier": ".2e"}},
    {"name": "max int", "id": "intensity_max", "type": "numeric", "format": {"specifier": ".2e"}},
    {"name": "score", "id": "hill_score", "type": "numeric", "format": {"specifier": ".3f"}},
]
TABLE_COLUMNS_MS2 = TABLE_COLUMNS_MS1 + [
    {"name": "iso_target", "id": "iso_target_mz", "type": "numeric", "format": {"specifier": ".3f"}},
]


def build_app(df: pd.DataFrame, source_label: str) -> Dash:
    ms2 = is_ms2(df)
    iw_options = isolation_window_options(df)
    columns = TABLE_COLUMNS_MS2 if ms2 else TABLE_COLUMNS_MS1

    # Static stats
    mz_lo, mz_hi = float(df["mz"].min()), float(df["mz"].max())
    rt_lo, rt_hi = float(df["rt"].min()), float(df["rt"].max())

    summary = (
        f"{len(df):,} hills · m/z {mz_lo:.2f}–{mz_hi:.2f} · "
        f"rt {rt_lo:.2f}–{rt_hi:.2f} min"
    )
    if ms2:
        n_windows = df[["iso_target_mz", "iso_lower_mz", "iso_upper_mz"]].drop_duplicates().shape[0]
        summary += f" · {n_windows} isolation windows"

    app = Dash(__name__, title="koth_ff hill viewer")

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
            html.Label("min n_scans", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Input(
                id="min-scans",
                type="number",
                min=1, step=1, value=1,
                style={"width": "100%"},
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
        html.Div([
            html.Label("Sort by", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Dropdown(
                id="sort-by",
                options=[
                    {"label": "intensity_sum (desc)", "value": "intensity_sum"},
                    {"label": "intensity_max (desc)", "value": "intensity_max"},
                    {"label": "n_scans (desc)", "value": "n_scans"},
                    {"label": "hill_score (desc)", "value": "hill_score"},
                    {"label": "m/z (asc)", "value": "mz"},
                    {"label": "rt (asc)", "value": "rt"},
                ],
                value="intensity_sum",
                clearable=False,
            ),
        ], style={"flex": 1, "padding": "0 12px"}),
    ]
    if ms2:
        filter_row.append(html.Div([
            html.Label("Isolation window", style={"fontWeight": "bold", "fontSize": 12}),
            dcc.Dropdown(
                id="iw-filter",
                options=iw_options,
                value="all",
                clearable=False,
            ),
        ], style={"flex": 2, "padding": "0 12px"}))

    app.layout = html.Div([
        html.Div([
            html.H2("koth_ff hill viewer", style={"margin": "0", "fontFamily": "system-ui, sans-serif"}),
            html.Div(source_label, style={"fontFamily": "monospace", "fontSize": 12, "color": "#555"}),
            html.Div(summary, style={"fontSize": 13, "color": "#333", "marginTop": 6}),
        ], style={"padding": "10px 20px", "borderBottom": "1px solid #ddd", "background": "#fafafa"}),

        html.Div(
            filter_row,
            style={"display": "flex", "alignItems": "center", "padding": "10px 8px", "borderBottom": "1px solid #eee"},
        ),

        html.Div([
            html.Div([
                html.Div(id="row-count", style={"fontSize": 12, "color": "#666", "padding": "4px 4px 8px"}),
                dash_table.DataTable(
                    id="hill-table",
                    columns=columns,
                    data=[],
                    page_size=25,
                    sort_action="native",
                    row_selectable="single",
                    selected_rows=[],
                    style_table={"height": "520px", "overflowY": "auto"},
                    style_cell={"fontFamily": "monospace", "fontSize": 12, "padding": "3px 6px"},
                    style_header={"fontWeight": "bold", "background": "#f3f3f3"},
                    style_data_conditional=[
                        {"if": {"state": "selected"}, "backgroundColor": "#dbeafe", "border": "1px solid #2563eb"},
                    ],
                ),
            ], style={"flex": 1, "padding": "10px"}),

            html.Div([
                dcc.Graph(id="hill-plot", figure=empty_figure("Select a hill from the table")),
                html.Div(id="hill-metadata", style={
                    "fontFamily": "monospace",
                    "fontSize": 12,
                    "padding": "10px",
                    "background": "#fafafa",
                    "border": "1px solid #eee",
                    "borderRadius": 6,
                    "marginTop": 8,
                    "whiteSpace": "pre",
                }),
            ], style={"flex": 1, "padding": "10px"}),
        ], style={"display": "flex", "alignItems": "stretch"}),

        # Hidden state
        dcc.Store(id="is-ms2", data=ms2),
    ])

    # ── Callbacks ────────────────────────────────────────────────────────────
    @app.callback(
        Output("hill-table", "data"),
        Output("hill-table", "selected_rows"),
        Output("row-count", "children"),
        Input("mz-slider", "value"),
        Input("rt-slider", "value"),
        Input("min-scans", "value"),
        Input("sort-by", "value"),
        *([Input("iw-filter", "value")] if ms2 else []),
    )
    def update_table(mz_range, rt_range, min_scans, sort_by, *extras):
        mask = (
            df["mz"].between(mz_range[0], mz_range[1])
            & df["rt"].between(rt_range[0], rt_range[1])
            & (df["n_scans"] >= (min_scans or 1))
        )
        if ms2 and extras:
            iw_value = extras[0]
            if iw_value and iw_value != "all":
                target = float(iw_value)
                mask &= np.isclose(df["iso_target_mz"], target, atol=1e-4)

        filtered = df.loc[mask]
        ascending = sort_by in ("mz", "rt")
        filtered = filtered.sort_values(sort_by, ascending=ascending)

        cap = 5000
        truncated = len(filtered) > cap
        view = filtered.head(cap)

        records = view[[c["id"] for c in columns]].to_dict("records")

        msg = f"{len(filtered):,} hills match filters"
        if truncated:
            msg += f" — showing first {cap:,}"
        return records, [], msg

    @app.callback(
        Output("hill-plot", "figure"),
        Output("hill-metadata", "children"),
        Input("hill-table", "selected_rows"),
        State("hill-table", "data"),
    )
    def update_plot(selected_rows, table_data):
        if not selected_rows or not table_data:
            return empty_figure("Select a hill from the table"), ""
        idx = selected_rows[0]
        if idx >= len(table_data):
            return no_update, no_update
        hill_id = table_data[idx]["hill_id"]
        row = df.loc[df["hill_id"] == hill_id].iloc[0]
        fig = hill_figure(row)
        meta = format_metadata(row)
        return fig, meta

    return app


def format_metadata(row: pd.Series) -> str:
    lines = [
        f"hill_id        {int(row['hill_id'])}",
        f"m/z            {row['mz']:.6f}   (std {row['mz_std']:.6f})",
        f"rt apex        {row['rt']:.4f} min",
        f"rt range       {row['rt_start']:.4f} – {row['rt_end']:.4f} min  (width {row['rt_width']:.4f})",
        f"scan range     {int(row['scan_start'])} – {int(row['scan_end'])}  (apex {int(row['scan_apex'])})",
        f"n_scans        {int(row['n_scans'])}   (skipped {int(row['skipped_scans'])})",
        f"intensity sum  {row['intensity_sum']:.4e}",
        f"intensity max  {row['intensity_max']:.4e}",
        f"hill_score     {row['hill_score']:.4f}",
    ]
    if float(row.get("im", 0.0)) != 0.0:
        lines.append(f"ion mobility   {row['im']:.4f}   (std {row['im_std']:.4f})")
    if "iso_target_mz" in row.index and pd.notna(row["iso_target_mz"]):
        lines.append(
            f"isolation      target {row['iso_target_mz']:.4f} "
            f"[{row['iso_lower_mz']:.4f} – {row['iso_upper_mz']:.4f}]"
        )
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# Entrypoint
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hills", type=Path, help="Path to a hills file (.tsv or .parquet)")
    parser.add_argument("--port", type=int, default=8050, help="Port to serve the dashboard on (default 8050)")
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind to (default 127.0.0.1)")
    parser.add_argument("--debug", action="store_true", help="Enable Dash debug mode")
    args = parser.parse_args()

    if not args.hills.exists():
        print(f"File not found: {args.hills}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {args.hills} …", file=sys.stderr)
    df = load_hills(args.hills)
    print(f"Loaded {len(df):,} hills ({'MS2' if is_ms2(df) else 'MS1'})", file=sys.stderr)

    app = build_app(df, source_label=str(args.hills.resolve()))
    print(f"Serving on http://{args.host}:{args.port}", file=sys.stderr)
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
