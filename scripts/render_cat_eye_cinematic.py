#!/usr/bin/env python3
"""Render the 2-D nonlinear KH run with a vortex-merger pseudocolor palette."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
from PIL import Image

BG = "#120909"
INK = "#fff4dc"
MUTED = "#dfb9a7"
CYAN = "#6ce9ed"
GOLD = "#ffd565"
GRID = "#7b3429"
CMAP = LinearSegmentedColormap.from_list(
    "vortex_merger",
    ["#07189f", "#0759ee", "#16bff4", "#6ee7bd", "#f4ef66", "#ffad3b", "#ef542f", "#8f2119", "#35100f"],
    N=512,
)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("frames", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--duration-ms", type=int, default=230)
    a = p.parse_args()
    paths = sorted(a.frames.glob("frame_*.h5"))
    if not paths:
        raise SystemExit("no HDF5 frames")
    with h5py.File(paths[0], "r") as h:
        cfg = json.loads(bytes(h["config_json"][:]).decode())
    # A single fixed physical range prevents color pumping between frames.
    fields = []
    for path in paths:
        with h5py.File(path, "r") as h:
            omega = h["omega_s-1"][:]
        z = np.linspace(-cfg["lz_m"]/2, cfg["lz_m"]/2, omega.shape[1])
        base = (cfg["wind_m_s"] / cfg["shear_m"]) / np.cosh(z / cfg["shear_m"])**2
        fields.append(omega - base[None, :])
    all_negative = np.concatenate([q[q < 0] for q in fields])
    vmin = float(np.nanpercentile(all_negative, 0.35))
    norm = Normalize(vmin=vmin, vmax=0.0, clip=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": BG, "axes.facecolor": BG, "text.color": INK})
    images: list[Image.Image] = []
    for i, (path, vort) in enumerate(zip(paths, fields)):
        with h5py.File(path, "r") as h:
            x = h["x_m"][:] / 1000 - cfg["lx_m"] / 2000
            z = h["altitude_m"][:] / 1000
            t = float(h.attrs["time_s"])

        fig = plt.figure(figsize=(12.8, 7.2), dpi=120)
        ax = fig.add_axes([0, 0, 1, 1])
        im = ax.pcolormesh(x, z, vort.T, cmap=CMAP, norm=norm, shading="gouraud", rasterized=True)
        ax.set_xlim(x[0], x[-1]); ax.set_ylim(z[0], z[-1])
        ax.set_axis_off()
        fig.text(.965, .93, f"{t/60:04.1f} min", ha="right", fontsize=18, weight="bold", color="#fff0c8",
                 bbox={"boxstyle":"round,pad=.35", "facecolor":"#160a09", "edgecolor":"none", "alpha":.55})
        bio = io.BytesIO(); fig.savefig(bio, format="png", facecolor=BG); plt.close(fig); bio.seek(0)
        images.append(Image.open(bio).convert("RGB").copy())

    a.output.parent.mkdir(parents=True, exist_ok=True)
    first = images[0].quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    gif_frames = [first] + [q.quantize(palette=first, dither=Image.Dither.FLOYDSTEINBERG) for q in images[1:]]
    gif_frames[0].save(a.output, save_all=True, append_images=gif_frames[1:], duration=a.duration_ms,
                       loop=0, disposal=2, optimize=False)
    images[-1].save(a.output.with_suffix(".png"), quality=96)
    print(f"{a.output}\nfixed vorticity range {vmin:.6g} to 0 s^-1")


if __name__ == "__main__":
    main()
