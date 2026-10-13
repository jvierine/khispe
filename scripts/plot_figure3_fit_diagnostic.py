#!/usr/bin/env python3
"""Render the figure-derived Figure 3 target against a nonlinear KHI run."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

BG = "#07131f"
PANEL = "#0d1d2b"
FG = "#eaf4f8"
MUTED = "#94adba"
CYAN = "#65d6e8"
GOLD = "#ffd166"
CORAL = "#ff6b6b"
GRID = "#274354"
DENSITY = LinearSegmentedColormap.from_list(
    "density", ["#09111b", "#173f5f", "#1d8a99", "#f2d06b", "#fff5dc"]
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", type=Path)
    parser.add_argument("simulation", type=Path)
    parser.add_argument("evaluation", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with h5py.File(args.target, "r") as h5:
        obs = h5["figure3/log10_ne_m3"][:]
        xo = h5["figure3/distance_km"][:]
        zo = h5["figure3/altitude_km"][:]
        t = h5["figure3/time_seconds_lt"][:].astype(float)
        time_text = [s.decode() for s in h5["figure3/time_text_lt"][:]]
    t -= t[0]

    with h5py.File(args.evaluation, "r") as h5:
        ridge_o = h5["observed_ridge_km"][:]
        ridge_p = h5["predicted_ridge_km"][:]
        sigma_o = h5["observed_sigma_km"][:]
        sigma_p = h5["predicted_sigma_km"][:]
        shift = int(h5.attrs["cyclic_shift_cells"])
        offset = float(h5.attrs["vertical_offset_km"])
        trend = float(h5.attrs["vertical_trend_km_per_600s"])
        train_rmse = float(h5.attrs["train_rmse_km"])
        held_rmse = float(h5.attrs["heldout_rmse_km"])
        growth_held = float(h5.attrs["sigma_growth_heldout_rmse_km"])

    chosen = [0, 5, 9, 14]
    model_images = []
    model_altitudes = []
    for index in chosen:
        path = args.simulation / f"frame_{index:05d}.h5"
        with h5py.File(path, "r") as h5:
            n = h5["ion_density_m-3"][:]
            xs = h5["x_m"][:] / 1000.0
            zs = h5["altitude_m"][:] / 1000.0
        n = np.roll(n, shift, axis=0)
        target_x = xo + 32.0
        registered = np.column_stack(
            [np.interp(target_x, xs, n[:, iz]) for iz in range(n.shape[1])]
        )
        model_images.append(np.log10(np.maximum(registered, 1.0)).T)
        model_altitudes.append(zs + offset + trend * t[index] / 600.0)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.facecolor": PANEL,
        "figure.facecolor": BG,
        "savefig.facecolor": BG,
        "text.color": FG,
        "axes.labelcolor": FG,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.edgecolor": GRID,
    })
    fig = plt.figure(figsize=(15.5, 9.2))
    grid = fig.add_gridspec(3, 4, height_ratios=[1, 1, 0.9])
    fig.subplots_adjust(left=0.055, right=0.985, top=0.875, bottom=0.075, wspace=0.065, hspace=0.19)

    for col, index in enumerate(chosen):
        ax = fig.add_subplot(grid[0, col])
        ax.pcolormesh(xo, zo, obs[index], cmap=DENSITY, vmin=10.95, vmax=11.08, shading="auto", rasterized=True)
        ax.plot(xo, ridge_o[index], color=GOLD, lw=1.5)
        ax.set_title(f"{time_text[index][:-3]} LT  •  panel {chr(97 + index)}", color=FG, fontsize=11, weight="bold")
        ax.set_ylim(113.0, 119.0)
        ax.set_xlim(-32, 54)
        if col == 0:
            ax.set_ylabel("Published figure\nAltitude [km]")
        else:
            ax.set_yticklabels([])
        ax.set_xticklabels([])

        ax = fig.add_subplot(grid[1, col])
        ax.pcolormesh(xo, model_altitudes[col], model_images[col], cmap=DENSITY, vmin=10.95, vmax=11.08, shading="auto", rasterized=True)
        ax.plot(xo, ridge_p[index], color=CYAN, lw=1.5)
        ax.set_ylim(113.0, 119.0)
        ax.set_xlim(-32, 54)
        ax.set_xlabel("Zonal distance [km]")
        if col == 0:
            ax.set_ylabel("Nonlinear run\nAltitude [km]")
        else:
            ax.set_yticklabels([])

    central = (xo >= -25.0) & (xo <= 48.0)
    axr = fig.add_subplot(grid[2, :2])
    ridge_error = np.nanmedian(np.abs(ridge_p[:, central] - ridge_o[:, central]), axis=1)
    axr.plot(t / 60.0, ridge_error, color=CYAN, marker="o", lw=2, ms=4)
    axr.axvspan(t[0] / 60.0, t[5] / 60.0, color=GOLD, alpha=0.10, label="fit: a–f")
    axr.axvspan(t[6] / 60.0, t[-1] / 60.0, color=CORAL, alpha=0.08, label="held out: g–o")
    axr.set_title("Ridge error stays bounded", loc="left", weight="bold")
    axr.set_xlabel("Elapsed time [min]")
    axr.set_ylabel("Median |ridge residual| [km]")
    axr.grid(color=GRID, alpha=0.45, lw=0.7)
    axr.legend(frameon=False, labelcolor=FG, ncols=2, loc="upper left")

    axw = fig.add_subplot(grid[2, 2:])
    obs_width = np.nanmedian(sigma_o[:, central], axis=1)
    pred_width = np.nanmedian(sigma_p[:, central], axis=1)
    axw.plot(t / 60.0, obs_width - obs_width[0], color=GOLD, marker="o", lw=2, ms=4, label="figure-derived core")
    axw.plot(t / 60.0, pred_width - pred_width[0], color=CYAN, marker="o", lw=2, ms=4, label="nonlinear run")
    axw.axhline(0, color=GRID, lw=0.8)
    axw.set_title("The modeled layer barely grows", loc="left", weight="bold")
    axw.set_xlabel("Elapsed time [min]")
    axw.set_ylabel("Median sigma change [km]")
    axw.grid(color=GRID, alpha=0.45, lw=0.7)
    axw.legend(frameon=False, labelcolor=FG, loc="upper left")

    fig.suptitle("SYISR Figure 3: nonlinear KHI run matches the ridge better, but not billow growth", fontsize=18, weight="bold", x=0.03, y=0.975, ha="left")
    fig.text(0.03, 0.925, "Figure-derived target • best full-resolution seed-amplitude case • panels g–o held out", color=MUTED, fontsize=11)
    fig.text(0.985, 0.925, f"ridge RMSE  {train_rmse:.3f} km fit  |  {held_rmse:.3f} km held out   •   width-growth RMSE  {growth_held:.3f} km held out", color=FG, fontsize=10, ha="right")
    fig.text(0.03, 0.018, "Source: Wang et al. (2025), DOI 10.1029/2025GL116932. Digitized publisher pixels; not a raw-data fit.", color=MUTED, fontsize=8.5)
    fig.text(0.985, 0.018, "Script: scripts/plot_figure3_fit_diagnostic.py", color=MUTED, fontsize=8.5, ha="right")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    print(args.output)


if __name__ == "__main__":
    main()
