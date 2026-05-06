"""
HTML report generator for koth_ff alignment results.

Called from run_align.py — not intended to be run directly.
"""

from __future__ import annotations

import base64
import io
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _img(b64: str, alt: str = "") -> str:
    if not b64:
        return ""
    return f'<img src="data:image/png;base64,{b64}" alt="{alt}" style="max-width:100%;margin:10px 0;">'


# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------

def _plot_detection_histogram(consensus_df: pl.DataFrame, n_runs: int) -> str:
    counts = (
        consensus_df.group_by("n_runs_detected")
        .agg(pl.len().alias("count"))
        .sort("n_runs_detected")
    )
    xs = counts["n_runs_detected"].to_list()
    ys = counts["count"].to_list()

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(xs, ys, color="#4c72b0", edgecolor="white", linewidth=0.6)
    ax.set_xlabel("Number of runs detected in", fontsize=11)
    ax.set_ylabel("Feature count", fontsize=11)
    ax.set_title("Feature detection across runs", fontsize=12)
    ax.set_xticks(range(1, n_runs + 1))
    for x, y in zip(xs, ys):
        ax.text(x, y + max(ys) * 0.01, f"{y:,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _kernel_warp_with_sd(
    anchor_pairs: list[tuple[float, float]],
    bandwidth: float,
    x_norm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (delta_norm, sd_norm) — kernel-weighted mean and weighted SD of anchor deltas.

    SD reflects the actual scatter of anchors around the smooth curve at each point,
    not the precision of the mean (SEM would be sd/sqrt(n) and would be invisible).
    """
    run_norms = np.array([p[1] for p in anchor_pairs])
    deltas    = np.array([p[0] - p[1] for p in anchor_pairs])

    # (n_query, n_anchors)
    w = np.exp(-0.5 * ((x_norm[:, np.newaxis] - run_norms[np.newaxis, :]) / bandwidth) ** 2)
    w_sum = w.sum(axis=1)
    safe  = w_sum > 1e-10

    delta_mean = np.where(safe, (w * deltas).sum(axis=1) / np.where(safe, w_sum, 1.0), 0.0)

    residuals = deltas[np.newaxis, :] - delta_mean[:, np.newaxis]
    var = (w * residuals ** 2).sum(axis=1) / np.where(safe, w_sum, 1.0)
    sd  = np.sqrt(var)

    return delta_mean, sd


def _plot_rt_warps(diagnostics: dict) -> str:
    run_diags = {
        k: v for k, v in diagnostics["runs"].items()
        if not v.get("identity_fallback", True)
    }
    if not run_diags:
        return ""

    n = len(run_diags)
    # Two panels per run: normalised warp  |  RT delta in minutes
    fig, axes = plt.subplots(n, 2, figsize=(12, 4 * n), squeeze=False)

    for row, (name, info) in enumerate(run_diags.items()):
        ax_warp  = axes[row][0]
        ax_delta = axes[row][1]

        all_pairs     = info.get("anchor_pairs", [])
        inlier_pairs  = info.get("inlier_pairs", all_pairs)
        bandwidth = info.get("bandwidth", 0.15)
        run_rt_min = info["run_rt_min"]
        run_rt_max = info["run_rt_max"]
        ref_rt_min = info["ref_rt_min"]
        ref_rt_max = info["ref_rt_max"]
        ref_span   = ref_rt_max - ref_rt_min or 1.0
        run_span   = run_rt_max - run_rt_min or 1.0

        x_norm = np.linspace(0, 1, 400)

        if inlier_pairs:
            run_norms = np.array([p[1] for p in inlier_pairs])
            ref_norms = np.array([p[0] for p in inlier_pairs])

            # Left: run RT → ref RT in absolute minutes (the raw mapping)
            anch_run_abs = run_norms * run_span + run_rt_min
            anch_ref_abs = ref_norms * ref_span + ref_rt_min
            ax_warp.scatter(anch_run_abs, anch_ref_abs, s=8, alpha=0.4,
                            color="#4c72b0", label=f"anchors (n={len(inlier_pairs)})", zorder=3)

            # Fitted warp curve in absolute minutes
            warp_fn_left = info.get("warp_fn")
            if warp_fn_left is not None:
                x_abs_left = x_norm * run_span + run_rt_min
                y_abs_left = warp_fn_left(x_abs_left)
                ax_warp.plot(x_abs_left, y_abs_left, color="#c0392b", linewidth=1.8,
                             label="kernel smooth", zorder=4)

            # Identity line in absolute space
            id_min = min(run_rt_min, ref_rt_min)
            id_max = max(run_rt_max, ref_rt_max)
            ax_warp.plot([id_min, id_max], [id_min, id_max], "k--",
                         linewidth=0.8, alpha=0.35, label="identity", zorder=2)

            delta_norm, sd_norm = _kernel_warp_with_sd(inlier_pairs, bandwidth, x_norm)

            # Delta panel: warp(run_rt) - run_rt gives the true absolute correction
            warp_fn = info.get("warp_fn")
            x_abs   = x_norm * run_span + run_rt_min
            if warp_fn is not None:
                correction_abs = warp_fn(x_abs) - x_abs
                sd_abs = sd_norm * ref_span

                # Anchor deltas behind the curve
                anch_run_abs   = run_norms * run_span + run_rt_min
                anch_ref_abs   = ref_norms * ref_span + ref_rt_min
                anch_delta_abs = anch_ref_abs - anch_run_abs
                ax_delta.scatter(anch_run_abs, anch_delta_abs, s=6, alpha=0.2,
                                 color="#4c72b0", zorder=2, label="anchor Δ")

                ax_delta.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.35)
                ax_delta.fill_between(x_abs, correction_abs - sd_abs, correction_abs + sd_abs,
                                      color="#c0392b", alpha=0.2, label="±1 SD")
                ax_delta.plot(x_abs, correction_abs - sd_abs, color="#c0392b",
                              linewidth=0.8, linestyle="--", alpha=0.7)
                ax_delta.plot(x_abs, correction_abs + sd_abs, color="#c0392b",
                              linewidth=0.8, linestyle="--", alpha=0.7)
                ax_delta.plot(x_abs, correction_abs, color="#c0392b", linewidth=2.0,
                              label="RT correction", zorder=4)

                # Sliding-window median control points (what the spline is fitted to).
                # wc is normalised run RT; wm is (ref_norm - run_norm).
                # Absolute correction = warp(run_rt_abs) - run_rt_abs
                #                     = (wc+wm)*ref_span + ref_rt_min  -  (wc*run_span + run_rt_min)
                # The naive wm*ref_span is only correct when gradients are identical.
                wc = np.array(info.get("window_centers", []))
                wm = np.array(info.get("window_medians", []))
                if len(wc) and len(wm):
                    wc_abs = wc * run_span + run_rt_min
                    wm_abs = (wc + wm) * ref_span + ref_rt_min - wc_abs
                    ax_delta.scatter(wc_abs, wm_abs, s=40, marker="D", color="#e67e22",
                                     edgecolors="white", linewidths=0.6, zorder=5,
                                     label=f"window medians (n={len(wc)})")

                # Y-axis tight around correction ± SD, not driven by outlier anchors
                pad = 1.0
                y_lo = float((correction_abs - sd_abs).min()) - pad
                y_hi = float((correction_abs + sd_abs).max()) + pad
                ax_delta.set_ylim(y_lo, y_hi)

        ax_warp.set_xlabel("Run RT (min)", fontsize=9)
        ax_warp.set_ylabel("Reference RT (min)", fontsize=9)
        ax_warp.set_title(f"{name[:40]}  — RT mapping (run → ref)", fontsize=9)
        ax_warp.legend(fontsize=7)

        ax_delta.set_xlabel("Run RT (min)", fontsize=9)
        ax_delta.set_ylabel("RT correction (min)", fontsize=9)
        ax_delta.set_title(f"{name[:40]}  — RT correction (min)", fontsize=9)
        ax_delta.legend(fontsize=7)

    fig.suptitle("RT warp curves (kernel smooth ± 1 SEM)", fontsize=11)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_correlation_heatmap(matrix_df: pl.DataFrame, sample_names: list[str]) -> str:
    intensity = matrix_df.select(sample_names).to_numpy().astype(float)
    detected_in_2 = (intensity > 0).sum(axis=1) >= 2
    intensity = intensity[detected_in_2]
    if intensity.shape[0] < 10:
        return ""

    log_int = np.log10(intensity + 1.0)
    n = len(sample_names)
    corr = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            mask = (intensity[:, i] > 0) & (intensity[:, j] > 0)
            if mask.sum() >= 5:
                corr[i, j] = pearsonr(log_int[mask, i], log_int[mask, j])[0]

    # Hierarchical clustering: reorder rows/columns by similarity so related
    # samples cluster together. Distance = 1 - |r|; NaN pairs treated as r=0.
    if n >= 3:
        corr_fill = np.nan_to_num(corr, nan=0.0)
        np.fill_diagonal(corr_fill, 1.0)
        dist = 1.0 - np.clip(corr_fill, 0.0, 1.0)
        np.fill_diagonal(dist, 0.0)
        Z = linkage(squareform(dist, checks=False), method="average")
        order = dendrogram(Z, no_plot=True)["leaves"]
    else:
        Z = None
        order = list(range(n))

    corr_ord = corr[np.ix_(order, order)]
    short    = [sample_names[i][:22] for i in order]

    size   = max(5, n * 0.75 + 1)
    dend_s = max(1.2, size * 0.18)   # dendrogram panel size

    if Z is not None:
        fig = plt.figure(figsize=(size + dend_s + 1.5, size + dend_s + 0.5))
        gs = fig.add_gridspec(
            2, 3,
            width_ratios=[dend_s, size, 0.35],
            height_ratios=[dend_s, size],
            hspace=0.01, wspace=0.01,
        )
        ax_corner    = fig.add_subplot(gs[0, 0])
        ax_top_dend  = fig.add_subplot(gs[0, 1])
        ax_left_dend = fig.add_subplot(gs[1, 0])
        ax_heatmap   = fig.add_subplot(gs[1, 1])
        ax_cbar      = fig.add_subplot(gs[1, 2])
        ax_corner.set_visible(False)

        _dend_kw = dict(no_labels=True, color_threshold=0, above_threshold_color="#666666")
        dendrogram(Z, ax=ax_top_dend,  orientation="top",  **_dend_kw)
        dendrogram(Z, ax=ax_left_dend, orientation="left", **_dend_kw)

        # Align dendrogram leaf coordinates with heatmap cell indices.
        # Scipy places leaf k at position 10*k+5, so n leaves span 0..10n.
        ax_top_dend.set_xlim(-5, 10 * n - 5)
        ax_top_dend.axis("off")
        ax_left_dend.set_ylim(10 * n - 5, -5)   # inverted to match heatmap top→bottom
        ax_left_dend.axis("off")

        im = ax_heatmap.imshow(corr_ord, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        ax_heatmap.set_xlim(-0.5, n - 0.5)
        ax_heatmap.set_ylim(n - 0.5, -0.5)
        ax_heatmap.set_xticks(range(n))
        ax_heatmap.set_yticks(range(n))
        ax_heatmap.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
        ax_heatmap.set_yticklabels(short, fontsize=7)
        for i in range(n):
            for j in range(n):
                if not np.isnan(corr_ord[i, j]):
                    ax_heatmap.text(j, i, f"{corr_ord[i, j]:.2f}",
                                    ha="center", va="center", fontsize=6, color="black")
        fig.colorbar(im, cax=ax_cbar).set_label("Pearson r", fontsize=8)
        ax_heatmap.set_title(
            "Pairwise Pearson r  (log₁₀ intensity, features ≥2 runs)\n"
            "ordered by hierarchical clustering  (average linkage, distance = 1−r)",
            fontsize=8,
        )
    else:
        fig, ax = plt.subplots(figsize=(size + 1, size))
        im = ax.imshow(corr_ord, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(short, fontsize=7)
        for i in range(n):
            for j in range(n):
                if not np.isnan(corr_ord[i, j]):
                    ax.text(j, i, f"{corr_ord[i, j]:.2f}", ha="center", va="center",
                            fontsize=6, color="black")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title("Pairwise Pearson r  (log₁₀ intensity, features in ≥2 runs)", fontsize=9)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_intensity_distributions(matrix_df: pl.DataFrame, sample_names: list[str]) -> str:
    data = [
        np.log10(matrix_df[name].filter(matrix_df[name] > 0).to_numpy() + 1)
        for name in sample_names
    ]

    fig, ax = plt.subplots(figsize=(max(8, len(sample_names) * 0.7 + 1), 5))
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for patch in bp["boxes"]:
        patch.set_facecolor("#4c72b0")
        patch.set_alpha(0.7)
    short = [s[:28] for s in sample_names]
    ax.set_xticks(range(1, len(sample_names) + 1))
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("log₁₀(intensity)", fontsize=10)
    ax.set_title("Intensity distribution per sample (detected features only)", fontsize=11)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_missing_rates(matrix_df: pl.DataFrame, sample_names: list[str]) -> str:
    n_consensus = len(matrix_df)
    rates = [
        100.0 * float((matrix_df[n] == 0).sum()) / n_consensus
        for n in sample_names
    ]
    short = [s[:28] for s in sample_names]

    fig, ax = plt.subplots(figsize=(max(7, len(sample_names) * 0.7 + 1), 4))
    bars = ax.bar(range(len(sample_names)), rates, color="#dd8047", edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Missing values (%)", fontsize=10)
    ax.set_ylim(0, 105)
    ax.set_title("Missing value rate per sample", fontsize=11)
    for i, r in enumerate(rates):
        ax.text(i, r + 1, f"{r:.1f}%", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Target / decoy FDR plots
# ---------------------------------------------------------------------------

def _plot_fdr_summary(diagnostics: dict) -> str:
    """Grouped bar chart: target vs decoy match counts per run, annotated with FDR %."""
    run_diags = {
        name: info for name, info in diagnostics["runs"].items()
        if "fdr" in info
    }
    if not run_diags:
        return ""

    names = list(run_diags.keys())
    n_targets = [run_diags[n]["fdr"]["n_target"] for n in names]
    n_decoys  = [run_diags[n]["fdr"]["n_decoy"]  for n in names]
    fdrs      = [run_diags[n]["fdr"]["fdr"] for n in names]

    x = np.arange(len(names))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 1.4 + 2), 5))
    bars_t = ax.bar(x - width / 2, n_targets, width, label="Target matches", color="#4c72b0", edgecolor="white")
    bars_d = ax.bar(x + width / 2, n_decoys,  width, label="Decoy matches",  color="#dd8047", edgecolor="white")

    max_y = max(n_targets + [1])
    for xi, (nt, nd, fdr) in enumerate(zip(n_targets, n_decoys, fdrs)):
        color = "#c0392b" if fdr > 0.01 else "#27ae60"
        ax.text(xi, nt + max_y * 0.01, f"FDR\n{fdr*100:.1f}%", ha="center", va="bottom",
                fontsize=7, color=color, fontweight="bold")

    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([n[:30] for n in names], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Feature matches", fontsize=10)
    ax.set_title("Target vs decoy matches per run  (competition FDR; red = > 1%)", fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_parameter_sweep(diagnostics: dict) -> str:
    """Heatmap of competition FDR over (mass_ppm, rt_window) per run and im_tolerance."""
    run_diags = {
        name: info for name, info in diagnostics["runs"].items()
        if "sweep" in info
    }
    if not run_diags:
        return ""

    figs = []
    for name, info in run_diags.items():
        sweep: pl.DataFrame = info["sweep"]
        im_vals = sorted(sweep["im_tolerance"].unique().to_list())
        ppm_vals = sorted(sweep["mass_ppm"].unique().to_list())
        rt_vals  = sorted(sweep["rt_window"].unique().to_list())

        n_im = len(im_vals)
        fig, axes = plt.subplots(1, n_im, figsize=(5 * n_im + 1, 4.5), squeeze=False)

        for col_idx, im_val in enumerate(im_vals):
            ax = axes[0][col_idx]
            sub = sweep.filter(pl.col("im_tolerance") == im_val)

            grid = np.full((len(ppm_vals), len(rt_vals)), np.nan)
            grid_n = np.full((len(ppm_vals), len(rt_vals)), 0)
            for row in sub.iter_rows(named=True):
                pi = ppm_vals.index(row["mass_ppm"])
                ri = rt_vals.index(row["rt_window"])
                grid[pi, ri] = row["fdr"]
                grid_n[pi, ri] = row["n_target"]

            im_img = ax.imshow(grid, vmin=0, vmax=0.05, cmap="RdYlGn_r", aspect="auto",
                               origin="lower")
            ax.set_xticks(range(len(rt_vals)))
            ax.set_yticks(range(len(ppm_vals)))
            ax.set_xticklabels([f"{v}" for v in rt_vals], fontsize=8)
            ax.set_yticklabels([f"{int(v)}" for v in ppm_vals], fontsize=8)
            ax.set_xlabel("RT window (min)", fontsize=8)
            ax.set_ylabel("Mass PPM", fontsize=8)
            ax.set_title(f"IM tol = {im_val}", fontsize=9)

            for pi in range(len(ppm_vals)):
                for ri in range(len(rt_vals)):
                    fdr_v = grid[pi, ri]
                    nt_v  = grid_n[pi, ri]
                    if not np.isnan(fdr_v):
                        ax.text(ri, pi, f"{fdr_v*100:.1f}%\n{nt_v:,}",
                                ha="center", va="center", fontsize=6)

            plt.colorbar(im_img, ax=ax, fraction=0.046, pad=0.04, label="FDR")

        fig.suptitle(f"Parameter sweep — {name[:50]}", fontsize=10)
        fig.tight_layout()
        figs.append(_fig_to_b64(fig))

    return "".join(_img(b, "Parameter sweep heatmap") for b in figs)


def _plot_runnerup_gap(diagnostics: dict) -> str:
    """Histogram of runner-up PPM gap across all matches (confidence proxy)."""
    all_gaps: list[np.ndarray] = []
    for info in diagnostics["runs"].values():
        gaps = info.get("match_stats", {}).get("runnerup_gaps")
        if gaps is not None:
            finite = gaps[np.isfinite(gaps)]
            if len(finite):
                all_gaps.append(finite)

    if not all_gaps:
        return ""

    combined = np.concatenate(all_gaps)
    capped = np.clip(combined, 0, np.percentile(combined, 95))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(capped, bins=60, color="#4c72b0", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Runner-up gap (PPM between best and 2nd-best candidate)", fontsize=10)
    ax.set_ylabel("Match count", fontsize=10)
    ax.set_title("Match ambiguity — runner-up PPM gap  (larger = more confident)", fontsize=11)
    ax.text(0.97, 0.96, f"n = {len(combined):,}\n(capped at 95th pct)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# m/z and IM drift plots
# ---------------------------------------------------------------------------

def _plot_mz_drift(diagnostics: dict) -> str:
    """Scatter of per-anchor PPM error vs RT with fitted linear trend, one panel per run."""
    run_diags = {
        k: v for k, v in diagnostics["runs"].items()
        if not v.get("identity_fallback", True) and v.get("anchor_ref_rt_abs")
    }
    if not run_diags:
        return ""

    n = len(run_diags)
    fig, axes = plt.subplots(n, 1, figsize=(9, 4 * n), squeeze=False)

    for row, (name, info) in enumerate(run_diags.items()):
        ax = axes[row][0]
        rt  = np.array(info["anchor_ref_rt_abs"])
        ppm = np.array(info["anchor_ppm_errors"])
        slope     = info.get("mass_drift_slope_ppm_per_min", 0.0)
        intercept = info.get("mass_drift_intercept_ppm", 0.0)

        ax.scatter(rt, ppm, s=8, alpha=0.45, color="#4c72b0",
                   label=f"anchors (n={len(rt):,})", zorder=3)
        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.35, zorder=2)

        if len(rt) >= 2:
            x_line = np.linspace(rt.min(), rt.max(), 300)
            ax.plot(x_line, slope * x_line + intercept, color="#c0392b",
                    linewidth=1.8, zorder=4,
                    label=f"fit: {slope:+.4f} ppm/min  (offset {intercept:+.3f} ppm)")

        # Sigma-clip y-axis so a few outliers don't crush the view
        if len(ppm) >= 4:
            p2, p98 = np.percentile(ppm, [2, 98])
            pad = max(0.5, (p98 - p2) * 0.25)
            ax.set_ylim(p2 - pad, p98 + pad)

        ax.set_xlabel("RT (min)", fontsize=9)
        ax.set_ylabel("Mass error (ppm)", fontsize=9)
        ax.set_title(f"{name[:50]}  — m/z drift vs RT", fontsize=9)
        ax.legend(fontsize=7)

    fig.suptitle("m/z drift correction (anchor pairs, linear fit)", fontsize=11)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_im_drift(diagnostics: dict) -> str:
    """Scatter of per-anchor IM delta vs RT with fitted linear trend, one panel per run."""
    run_diags = {
        k: v for k, v in diagnostics["runs"].items()
        if not v.get("identity_fallback", True) and v.get("anchor_ref_rt_abs")
    }
    # Only produce the plot if at least one run has IM anchor data
    has_im = any(
        np.any(np.array(v.get("anchor_im_deltas", [])) != 0.0)
        for v in run_diags.values()
    )
    if not run_diags or not has_im:
        return ""

    n = len(run_diags)
    fig, axes = plt.subplots(n, 1, figsize=(9, 4 * n), squeeze=False)

    for row, (name, info) in enumerate(run_diags.items()):
        ax = axes[row][0]
        rt_all  = np.array(info["anchor_ref_rt_abs"])
        im_all  = np.array(info["anchor_im_deltas"])
        mask    = im_all != 0.0
        slope     = info.get("im_drift_slope_per_min", 0.0)
        intercept = info.get("im_drift_intercept", 0.0)

        if mask.sum() == 0:
            ax.text(0.5, 0.5, "No IM data in anchors", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="#888")
            ax.set_title(f"{name[:50]}  — IM drift vs RT", fontsize=9)
            continue

        rt  = rt_all[mask]
        imd = im_all[mask]

        ax.scatter(rt, imd, s=8, alpha=0.45, color="#4c72b0",
                   label=f"anchors with IM (n={mask.sum():,})", zorder=3)
        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.35, zorder=2)

        if len(rt) >= 4:
            x_line = np.linspace(rt.min(), rt.max(), 300)
            ax.plot(x_line, slope * x_line + intercept, color="#c0392b",
                    linewidth=1.8, zorder=4,
                    label=f"fit: {slope:+.5f} 1/K₀/min  (offset {intercept:+.4f})")

        if len(imd) >= 4:
            p2, p98 = np.percentile(imd, [2, 98])
            pad = max(0.005, (p98 - p2) * 0.25)
            ax.set_ylim(p2 - pad, p98 + pad)

        ax.set_xlabel("RT (min)", fontsize=9)
        ax.set_ylabel("IM delta (1/K₀)", fontsize=9)
        ax.set_title(f"{name[:50]}  — IM drift vs RT", fontsize=9)
        ax.legend(fontsize=7)

    fig.suptitle("Ion-mobility drift correction (anchor pairs, linear fit)", fontsize=11)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def generate_report(
    consensus_df: pl.DataFrame,
    matrix_df: pl.DataFrame,
    runs: dict[str, pl.DataFrame],
    diagnostics: dict,
    out_path: Path,
) -> None:
    sample_names = [c for c in matrix_df.columns if c not in ("mass", "charge", "rt_apex")]
    n_runs = len(sample_names)
    ref_name = diagnostics["reference"]

    img_detection   = _plot_detection_histogram(consensus_df, n_runs)
    img_warps       = _plot_rt_warps(diagnostics)
    img_mz_drift    = _plot_mz_drift(diagnostics)
    img_im_drift    = _plot_im_drift(diagnostics)
    img_dist        = _plot_intensity_distributions(matrix_df, sample_names)
    img_corr        = _plot_correlation_heatmap(matrix_df, sample_names)
    img_missing     = _plot_missing_rates(matrix_df, sample_names)
    img_fdr         = _plot_fdr_summary(diagnostics)
    img_sweep       = _plot_parameter_sweep(diagnostics)
    img_runnerup    = _plot_runnerup_gap(diagnostics)

    # Per-sample summary table rows
    rows = []
    for name in sample_names:
        n_input    = len(runs.get(name, []))
        n_detected = int((matrix_df[name] > 0).sum())
        pct        = f"{100 * n_detected / len(consensus_df):.1f}%"
        med_int    = matrix_df[name].filter(matrix_df[name] > 0).median()
        n_anchors  = diagnostics["runs"].get(name, {}).get("n_anchors", "—")
        fallback   = diagnostics["runs"].get(name, {}).get("identity_fallback", False)
        warp_str   = "identity" if fallback else str(n_anchors)
        fdr_info   = diagnostics["runs"].get(name, {}).get("fdr")
        fdr_str    = f"{fdr_info['fdr']*100:.1f}%" if fdr_info else "—"
        rows.append((name, f"{n_input:,}", f"{n_detected:,}", pct,
                     f"{float(med_int or 0):.3e}", warp_str, fdr_str))

    def section(title: str, content: str) -> str:
        return f'<div class="section"><h2>{title}</h2>{content}</div>'

    table_rows = "".join(
        f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td>"
        f"<td>{r[3]}</td><td>{r[4]}</td><td>{r[5]}</td><td>{r[6]}</td></tr>"
        for r in rows
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>koth_ff Alignment Report</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    max-width: 1200px; margin: 40px auto; padding: 0 24px; color: #222;
    line-height: 1.5;
}}
h1   {{ color: #1a4a72; border-bottom: 3px solid #1a4a72; padding-bottom: 10px; }}
h2   {{ color: #1a4a72; margin-top: 40px; font-size: 1.2em; }}
p.meta {{ color: #666; font-size: 13px; margin-top: -8px; }}
.stat-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    gap: 16px; margin: 24px 0;
}}
.stat-box {{
    background: #eef5fc; border-left: 5px solid #1a4a72;
    padding: 14px 18px; border-radius: 6px;
}}
.stat-box .val {{ font-size: 30px; font-weight: 700; color: #1a4a72; }}
.stat-box .lbl {{ font-size: 12px; color: #555; margin-top: 2px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin-top: 10px; }}
th, td {{ padding: 7px 14px; text-align: left; border-bottom: 1px solid #ddd; }}
th {{ background: #1a4a72; color: white; font-weight: 600; }}
tr:nth-child(even) {{ background: #f7f9fc; }}
.section {{ margin: 36px 0; }}
.two-col {{ display: flex; flex-wrap: wrap; gap: 20px; }}
.two-col > div {{ flex: 1; min-width: 320px; }}
img {{ border-radius: 4px; box-shadow: 0 1px 4px rgba(0,0,0,0.12); }}
footer {{ margin-top: 60px; font-size: 11px; color: #aaa;
          border-top: 1px solid #eee; padding-top: 12px; }}
</style>
</head>
<body>

<h1>koth_ff — Alignment Report</h1>
<p class="meta">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")}
&nbsp;·&nbsp; Reference run: <strong>{ref_name}</strong>
&nbsp;·&nbsp; {n_runs} runs aligned</p>

<div class="stat-grid">
  <div class="stat-box">
    <div class="val">{n_runs}</div>
    <div class="lbl">Runs aligned</div>
  </div>
  <div class="stat-box">
    <div class="val">{len(consensus_df):,}</div>
    <div class="lbl">Consensus features</div>
  </div>
  <div class="stat-box">
    <div class="val">{int((consensus_df["n_runs_detected"] == n_runs).sum()):,}</div>
    <div class="lbl">Detected in all runs</div>
  </div>
  <div class="stat-box">
    <div class="val">{int((consensus_df["n_runs_detected"] >= 2).sum()):,}</div>
    <div class="lbl">Detected in ≥ 2 runs</div>
  </div>
  <div class="stat-box">
    <div class="val">{100 * int((consensus_df["n_runs_detected"] == n_runs).sum()) // max(1, len(consensus_df))}%</div>
    <div class="lbl">Complete detection rate</div>
  </div>
</div>

{section("Per-sample summary", f"""
<table>
<tr>
  <th>Sample</th><th>Input features</th><th>Matched</th>
  <th>% detected</th><th>Median intensity</th><th>RT anchors used</th><th>FDR (competition)</th>
</tr>
{table_rows}
</table>""")}

{section("Feature detection across runs", _img(img_detection, "Detection histogram"))}

{section("RT warp curves (kernel smooth ± 1 SEM)", _img(img_warps, "RT warp curves")) if img_warps else ""}

{section("m/z drift vs RT (linear correction)", _img(img_mz_drift, "m/z drift")) if img_mz_drift else ""}

{section("Ion-mobility drift vs RT (linear correction)", _img(img_im_drift, "IM drift")) if img_im_drift else ""}

{section("Target / decoy FDR", _img(img_fdr, "FDR summary")) if img_fdr else ""}

{section("Match ambiguity — runner-up PPM gap", _img(img_runnerup, "Runner-up gap")) if img_runnerup else ""}

{section("Parameter sweep — FDR vs tolerances", img_sweep) if img_sweep else ""}

{section("Intensity distributions &amp; missing values", f"""
<div class="two-col">
  <div>{_img(img_dist, "Intensity distributions")}</div>
  <div>{_img(img_missing, "Missing value rates")}</div>
</div>""")}

{section("Pairwise intensity correlation (hierarchically clustered)", _img(img_corr, "Correlation heatmap"))}

<footer>
  koth_ff alignment report &nbsp;·&nbsp;
  {len(consensus_df):,} consensus features &nbsp;·&nbsp;
  {n_runs} runs &nbsp;·&nbsp;
  reference: {ref_name}
</footer>

</body>
</html>"""

    out_path.write_text(html, encoding="utf-8")
    print(f"  report:   {out_path}")
