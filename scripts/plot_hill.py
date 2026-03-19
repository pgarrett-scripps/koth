#!/usr/bin/env python3
"""
Plot a single chromatographic hill from a koth_ff hills.tsv file.

Usage:
    python scripts/plot_hill.py hills.tsv 42
    python scripts/plot_hill.py hills.tsv 42 --save hill_42.png
"""

import argparse
import json
import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d


def smooth(profile: np.ndarray, sigma: float) -> np.ndarray:
    """Gaussian-smooth a 1-D profile, treating zeros as gaps (NaN during smooth)."""
    if sigma <= 0 or len(profile) < 3:
        return profile
    # Replace zeros with NaN so gaps don't bleed into neighbours
    work = np.where(profile > 0, profile, np.nan)
    # Fill NaNs with interpolated values just for the smoothing pass
    nans = np.isnan(work)
    if nans.all():
        return profile
    idx = np.arange(len(work))
    work[nans] = np.interp(idx[nans], idx[~nans], work[~nans])
    smoothed = gaussian_filter1d(work, sigma=sigma)
    # Restore zeros at original gap positions
    smoothed[profile == 0] = 0.0
    return smoothed


def load_hill(tsv_path: str, index: int) -> pd.Series:
    df = pd.read_csv(tsv_path, sep="\t", nrows=index + 1)
    if index >= len(df):
        print(f"Error: index {index} out of range (file has {len(df)} rows)", file=sys.stderr)
        sys.exit(1)
    return df.iloc[index]


def plot_hill(row: pd.Series, index: int, save_path: str | None = None, smooth_sigma: float = 2.0):
    profile = np.array(json.loads(row["intensity_profile"]), dtype=np.float64)
    n = len(profile)
    scans = np.arange(int(row["scan_start"]), int(row["scan_start"]) + n)
    apex_scan = int(row["scan_apex"])

    # Compute per-scan RT estimate (linear interpolation between rt_start and rt_end)
    rt_start = float(row["rt_start"])
    rt_end = float(row["rt_end"])
    if rt_end > rt_start and n > 1:
        rts = np.linspace(rt_start, rt_end, n)
    else:
        rts = np.full(n, float(row["rt"]))

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle(
        f"Hill #{index}  |  m/z {row['mz']:.5f} ± {row['mz_std']:.5f}"
        + (f"  |  IM {row['im']:.4f}" if row["im"] != 0.0 else ""),
        fontsize=13,
    )

    smoothed = smooth(profile, smooth_sigma)

    # --- Top: elution profile ---
    ax = axes[0]
    ax.fill_between(rts, profile, alpha=0.15, color="steelblue")
    ax.plot(rts, profile, color="steelblue", linewidth=0.8, alpha=0.5, label="raw")
    if smooth_sigma > 0:
        ax.plot(rts, smoothed, color="steelblue", linewidth=2.0, label=f"smoothed (σ={smooth_sigma})")

    # Mark gaps (zero intensity)
    gap_mask = profile == 0.0
    if gap_mask.any():
        ax.axvline(x=rts[gap_mask], color="lightcoral", linewidth=0.7, alpha=0.5, label="gap scans")

    # Mark apex
    rel_apex = apex_scan - int(row["scan_start"])
    if 0 <= rel_apex < n:
        ax.axvline(x=rts[rel_apex], color="tomato", linewidth=1.5, linestyle="--", label=f"apex (scan {apex_scan})")

    ax.set_xlabel("Retention time (min)")
    ax.set_ylabel("Intensity")
    ax.set_xlim(rts[0], rts[-1])
    ax.yaxis.get_major_formatter().set_scientific(True)
    ax.yaxis.get_major_formatter().set_powerlimits((-2, 4))
    ax.legend(fontsize=9)

    # Annotation box
    info = (
        f"RT apex: {row['rt']:.4f} min\n"
        f"RT range: {rt_start:.4f} – {rt_end:.4f} min\n"
        f"Scans: {int(row['scan_start'])} – {int(row['scan_end'])}  ({int(row['n_scans'])} total, {int(row['skipped_scans'])} gaps)\n"
        f"Intensity max: {row['intensity_max']:.3e}\n"
        f"Intensity sum: {row['intensity_sum']:.3e}"
    )
    ax.text(
        0.02, 0.97, info,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8),
    )

    # --- Bottom: log-scale profile ---
    ax2 = axes[1]
    with np.errstate(divide="ignore"):
        log_profile = np.where(profile > 0, np.log10(profile), np.nan)
        log_smoothed = np.where(smoothed > 0, np.log10(smoothed), np.nan)
    ax2.fill_between(rts, log_profile, alpha=0.15, color="steelblue")
    ax2.plot(rts, log_profile, color="steelblue", linewidth=0.8, alpha=0.5)
    if smooth_sigma > 0:
        ax2.plot(rts, log_smoothed, color="steelblue", linewidth=2.0)
    ax2.set_xlabel("Retention time (min)")
    ax2.set_ylabel("log₁₀(intensity)")
    ax2.set_xlim(rts[0], rts[-1])

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot a hill from hills.tsv")
    parser.add_argument("tsv", help="Path to hills.tsv")
    parser.add_argument("index", type=int, help="0-based row index of the hill to plot")
    parser.add_argument("--save", metavar="FILE", help="Save figure to file instead of displaying")
    parser.add_argument("--smooth", type=float, default=2.0, metavar="SIGMA",
                        help="Gaussian smoothing sigma (default 2.0; 0 to disable)")
    args = parser.parse_args()

    row = load_hill(args.tsv, args.index)
    plot_hill(row, args.index, save_path=args.save, smooth_sigma=args.smooth)


if __name__ == "__main__":
    main()
