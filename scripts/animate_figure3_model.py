#!/usr/bin/env python3
"""Animate the full-resolution 2-D KHI model state from Rust/MPI HDF5 frames."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

BG = "#07131f"
PANEL = "#0d1d2b"
FG = "#eaf4f8"
MUTED = "#91aab7"
CYAN = "#5fe0f0"
GOLD = "#ffd166"
GRID = "#294756"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("frames", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    paths = sorted(args.frames.glob("frame_*.h5"))
    if not paths:
        raise SystemExit("no frame_*.h5 files found")

    with h5py.File(paths[0], "r") as h5:
        config = json.loads(bytes(h5["config_json"][:]).decode())

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "figure.facecolor": BG,
        "axes.facecolor": PANEL,
        "text.color": FG,
        "axes.labelcolor": FG,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.edgecolor": GRID,
    })

    rendered: list[Image.Image] = []
    for frame_index, path in enumerate(paths):
        with h5py.File(path, "r") as h5:
            x = h5["x_m"][:] / 1000.0 - config["lx_m"] / 2000.0
            z = h5["altitude_m"][:] / 1000.0
            n = h5["ion_density_m-3"][:]
            u = h5["u_m_s"][:]
            w = h5["w_m_s"][:]
            omega = h5["omega_s-1"][:]
            time_s = float(h5.attrs["time_s"])

        # Remove the horizontally uniform background wind before drawing arrows.
        up = u - np.median(u, axis=0, keepdims=True)
        speed = np.hypot(up, w)
        qscale = max(float(np.nanpercentile(speed, 99)), 0.05)
        vmax_omega = max(float(np.nanpercentile(np.abs(omega), 99.5)), 1e-5)
        logn = np.log10(np.maximum(n, 1.0))

        fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.7), dpi=105)
        fig.subplots_adjust(left=0.062, right=0.965, top=0.82, bottom=0.12, wspace=0.17)
        fig.suptitle("Nonlinear 2-D KHI model", x=0.04, y=0.965, ha="left", fontsize=21, weight="bold")
        fig.text(0.04, 0.91, "Rust/MPI • full-resolution 512 × 192 state • density is a conservative tracer", color=MUTED, fontsize=11)
        fig.text(0.96, 0.95, f"t = {time_s/60:.1f} min   •   frame {frame_index+1:02d}/{len(paths):02d}", ha="right", color=FG, fontsize=13, weight="bold")
        pairing = config.get("seed_pairing_m_s", 0.0)
        fig.text(0.96, 0.91, f"Ri = {config['n2_s2'] * (config['shear_m']/config['wind_m_s'])**2:.3f}   •   seed λ = {config['seed_wavelength_m']/1000:.1f} km   •   primary/pairing w = {config['seed_m_s']:.1f}/{pairing:.1f} m s⁻¹", ha="right", color=GOLD, fontsize=10)

        ax = axes[0]
        im0 = ax.pcolormesh(x, z, logn.T, cmap="magma", vmin=10.0, vmax=11.0, shading="auto", rasterized=True)
        sx = slice(4, None, 26)
        sz = slice(4, None, 16)
        ax.quiver(
            x[sx], z[sz], up[sx, sz].T / qscale, w[sx, sz].T / qscale,
            color=CYAN, alpha=0.82, pivot="mid", scale=21, width=0.0032,
            headwidth=3.2, headlength=4.0,
        )
        ax.set_title("Tracer density + perturbation velocity", loc="left", fontsize=13, weight="bold")
        ax.set_xlabel("Zonal distance [km]")
        ax.set_ylabel("Altitude [km]")
        cb0 = fig.colorbar(im0, ax=ax, orientation="horizontal", pad=0.13, fraction=0.055)
        cb0.set_label("log₁₀ density [m⁻³]")

        ax = axes[1]
        im1 = ax.pcolormesh(x, z, omega.T, cmap="RdBu_r", vmin=-vmax_omega, vmax=vmax_omega, shading="auto", rasterized=True)
        levels = [10.6, 10.85, 10.97]
        ax.contour(x, z, logn.T, levels=levels, colors=[MUTED, GOLD, FG], linewidths=[0.7, 1.1, 1.5], alpha=0.95)
        ax.set_title("Neutral vorticity + density contours", loc="left", fontsize=13, weight="bold")
        ax.set_xlabel("Zonal distance [km]")
        ax.set_ylabel("Altitude [km]")
        cb1 = fig.colorbar(im1, ax=ax, orientation="horizontal", pad=0.13, fraction=0.055)
        cb1.set_label("Vorticity [s⁻¹] • frame-wise 99.5% range")

        for ax in axes:
            ax.set_xlim(x[0], x[-1])
            ax.set_ylim(z[0], z[-1])
            ax.grid(color=GRID, alpha=0.22, lw=0.6)

        footer = (
            "Symmetry-breaking pilot: primary cat's eyes can pair and merge."
            if pairing > 0
            else "Current result: the layer advects and ripples, but does not develop the observed billow thickness."
        )
        fig.text(0.04, 0.025, footer, color=GOLD, fontsize=10, weight="bold")
        fig.text(0.96, 0.025, "Script: scripts/animate_figure3_model.py", color=MUTED, fontsize=8.5, ha="right")
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=105, facecolor=BG)
        plt.close(fig)
        buffer.seek(0)
        rendered.append(Image.open(buffer).convert("RGB").copy())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    first = rendered[0].quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    gif_frames = [first] + [
        frame.quantize(palette=first, dither=Image.Dither.NONE) for frame in rendered[1:]
    ]
    gif_frames[0].save(
        args.output,
        save_all=True,
        append_images=gif_frames[1:],
        duration=420,
        loop=0,
        disposal=2,
        optimize=False,
    )
    rendered[-1].convert("RGB").save(args.output.with_suffix(".png"), quality=95)
    print(args.output)


if __name__ == "__main__":
    main()
