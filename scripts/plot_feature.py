#!/usr/bin/env python3
"""
Plot a single isotope feature from a koth_ff features.tsv file.

Three panels:
  1. Elution profile — summed intensity across all isotopes vs RT scan
  2. Isotope envelope — observed vs theoretical averagine pattern at apex
  3. Per-isotope elution traces (one line per isotope hill)

Usage:
    python scripts/plot_feature.py features.tsv 0
    python scripts/plot_feature.py features.tsv 42 --save feature_42.png
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
    work = np.where(profile > 0, profile, np.nan)
    nans = np.isnan(work)
    if nans.all():
        return profile
    idx = np.arange(len(work))
    work[nans] = np.interp(idx[nans], idx[~nans], work[~nans])
    smoothed = gaussian_filter1d(work, sigma=sigma)
    smoothed[profile == 0] = 0.0
    return smoothed


def load_feature(tsv_path: str, index: int) -> pd.Series:
    df = pd.read_csv(tsv_path, sep="\t", nrows=index + 1)
    if index >= len(df):
        print(f"Error: index {index} out of range (file has {len(df)} rows)", file=sys.stderr)
        sys.exit(1)
    return df.iloc[index]


def parse_json_col(value) -> list:
    if pd.isna(value) or value == "" or value == "[]":
        return []
    return json.loads(value)


def plot_feature(row: pd.Series, index: int, save_path: str | None = None, smooth_sigma: float = 2.0):
    elution_profile = np.array(parse_json_col(row["elution_profile"]), dtype=np.float64)
    isotope_profile = np.array(parse_json_col(row["isotope_profile"]), dtype=np.float64)
    theoretical_pattern = np.array(parse_json_col(row["theoretical_pattern"]), dtype=np.float64)
    mono_scan_lists = parse_json_col(row["mono_hills_scan_lists"])
    mono_intensity_lists = parse_json_col(row["mono_hills_intensity_list"])

    charge = int(row["charge"]) if row["charge"] else 0
    n_isotopes = int(row["nIsotopes"]) if not pd.isna(row["nIsotopes"]) else len(isotope_profile)
    score = float(row["score"]) if not pd.isna(row["score"]) else 0.0
    cosine = float(row["cosine_similarity"]) if not pd.isna(row["cosine_similarity"]) else 0.0
    mz = float(row["mz"])
    mass = float(row["massCalib"]) if not pd.isna(row["massCalib"]) else 0.0
    rt_apex = float(row["rtApex"]) if not pd.isna(row["rtApex"]) else 0.0
    rt_start = float(row["rtStart"]) if not pd.isna(row["rtStart"]) else rt_apex
    rt_end = float(row["rtEnd"]) if not pd.isna(row["rtEnd"]) else rt_apex
    im = float(row["im"]) if not pd.isna(row["im"]) and row["im"] != "" else 0.0

    fig = plt.figure(figsize=(12, 9))
    gs = plt.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    ax_elution = fig.add_subplot(gs[0, :])   # top full-width: elution profile
    ax_isotope = fig.add_subplot(gs[1, 0])   # bottom-left: isotope envelope
    ax_traces = fig.add_subplot(gs[1, 1])    # bottom-right: per-isotope traces

    charge_str = f"{charge}+" if charge > 0 else "unknown"
    title = (
        f"Feature #{index}  |  m/z {mz:.5f}  |  z={charge_str}  |  "
        f"mass {mass:.4f} Da  |  RT {rt_apex:.3f} min"
    )
    if im != 0.0:
        title += f"  |  IM {im:.4f}"
    fig.suptitle(title, fontsize=11)

    # ---- Panel 1: elution profile ----
    if len(elution_profile) > 0:
        x = np.arange(len(elution_profile))
        # Estimate RT axis from rt_start/rt_end
        if rt_end > rt_start and len(elution_profile) > 1:
            rts = np.linspace(rt_start, rt_end, len(elution_profile))
        else:
            rts = x.astype(float)
        elution_smoothed = smooth(elution_profile, smooth_sigma)
        ax_elution.fill_between(rts, elution_profile, alpha=0.15, color="steelblue")
        ax_elution.plot(rts, elution_profile, color="steelblue", linewidth=0.8, alpha=0.5, label="raw")
        if smooth_sigma > 0:
            ax_elution.plot(rts, elution_smoothed, color="steelblue", linewidth=2.0,
                            label=f"smoothed (σ={smooth_sigma})")
        ax_elution.axvline(rt_apex, color="tomato", linewidth=1.5, linestyle="--", label=f"apex {rt_apex:.3f}")
        ax_elution.set_xlabel("Retention time (min)")
        ax_elution.set_ylabel("Summed intensity")
        ax_elution.set_title("Elution profile (all isotopes summed)")
        ax_elution.yaxis.get_major_formatter().set_scientific(True)
        ax_elution.yaxis.get_major_formatter().set_powerlimits((-2, 4))
        ax_elution.legend(fontsize=9)

        info = (
            f"score: {score:.4f}  |  cosine: {cosine:.4f}  |  "
            f"n_isotopes: {n_isotopes}  |  ppm error: {row['ppm_error']:.2f}  |  "
            f"neutron offset: {int(row['neutron_offset']) if not pd.isna(row['neutron_offset']) else 0}"
        )
        ax_elution.text(
            0.02, 0.97, info,
            transform=ax_elution.transAxes,
            fontsize=8, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8),
        )
    else:
        ax_elution.text(0.5, 0.5, "No elution profile data", ha="center", va="center")

    # ---- Panel 2: observed vs theoretical isotope envelope ----
    if len(isotope_profile) > 0:
        obs = np.array(isotope_profile, dtype=np.float64)
        obs_sum = obs.sum()
        obs_norm = obs / obs_sum if obs_sum > 0 else obs

        n_bars = max(len(obs_norm), len(theoretical_pattern))
        x_bars = np.arange(n_bars)
        width = 0.35

        ax_isotope.bar(x_bars[:len(obs_norm)] - width / 2, obs_norm, width,
                       label="Observed", color="steelblue", alpha=0.8)
        if len(theoretical_pattern) > 0:
            theo = np.array(theoretical_pattern, dtype=np.float64)
            ax_isotope.bar(x_bars[:len(theo)] + width / 2, theo, width,
                           label="Theoretical (averagine)", color="tomato", alpha=0.7)

        ax_isotope.set_xlabel("Isotope index (M, M+1, M+2, …)")
        ax_isotope.set_ylabel("Relative intensity")
        ax_isotope.set_title("Isotope envelope")
        ax_isotope.set_xticks(x_bars)
        ax_isotope.legend(fontsize=9)
    else:
        ax_isotope.text(0.5, 0.5, "No isotope profile data", ha="center", va="center")

    # ---- Panel 3: monoisotopic elution trace ----
    if mono_scan_lists and mono_intensity_lists:
        scans = np.array(mono_scan_lists)
        ints = np.array(mono_intensity_lists, dtype=np.float64)
        mono_smoothed = smooth(ints, smooth_sigma)
        ax_traces.fill_between(scans, ints, alpha=0.15, color="steelblue")
        ax_traces.plot(scans, ints, color="steelblue", linewidth=0.8, alpha=0.5, label="raw")
        if smooth_sigma > 0:
            ax_traces.plot(scans, mono_smoothed, color="steelblue", linewidth=2.0,
                           label=f"smoothed (σ={smooth_sigma})")
        ax_traces.set_xlabel("Scan index")
        ax_traces.set_ylabel("Intensity")
        ax_traces.set_title("Monoisotopic hill trace")
        ax_traces.yaxis.get_major_formatter().set_scientific(True)
        ax_traces.yaxis.get_major_formatter().set_powerlimits((-2, 4))
        ax_traces.legend(fontsize=9)
    else:
        ax_traces.text(0.5, 0.5, "No monoisotopic trace data", ha="center", va="center")

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot a feature from features.tsv")
    parser.add_argument("tsv", help="Path to features.tsv")
    parser.add_argument("index", type=int, help="0-based row index of the feature to plot")
    parser.add_argument("--save", metavar="FILE", help="Save figure to file instead of displaying")
    parser.add_argument("--smooth", type=float, default=2.0, metavar="SIGMA",
                        help="Gaussian smoothing sigma (default 2.0; 0 to disable)")
    args = parser.parse_args()

    row = load_feature(args.tsv, args.index)
    plot_feature(row, args.index, save_path=args.save, smooth_sigma=args.smooth)


if __name__ == "__main__":
    main()
