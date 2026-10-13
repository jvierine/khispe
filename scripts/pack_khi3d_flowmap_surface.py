#!/usr/bin/env python3
"""Extract advected Clebsch material sheets for the Three.js flow-map view.

The first Clebsch label is initialized as sin(2 pi x/Lx) and is then
materially advected by the 3-D solver.  Its level surfaces therefore form a
compact, reproducible flow-map visualization: they show how initially planar
material sheets are folded by the Kelvin--Helmholtz billows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import map_coordinates
from skimage.measure import marching_cubes


p = argparse.ArgumentParser()
p.add_argument("input", type=Path)
p.add_argument("stem", type=Path)
p.add_argument("--levels", type=float, nargs="+", default=[-0.72, -0.20, 0.32, 0.82])
p.add_argument("--vertical-exaggeration", type=float, default=4.0)
a = p.parse_args()

with h5py.File(a.input, "r") as f:
    alpha = f["clebsch_alpha1"][:]
    ox = f["omega_x_s-1"][:]
    oy = f["omega_y_s-1"][:]
    oz = f["omega_z_s-1"][:]
    x, y, z = (f[k][:] for k in ("x_m", "y_m", "altitude_m"))
    t = float(f.attrs["time_s"])

dx, dy, dz = (float(np.mean(np.diff(q))) for q in (x, y, z))
lx = float(np.ptp(x))
omag = np.sqrt(ox * ox + oy * oy + oz * oz)
oscale = max(float(np.percentile(omag, 99.5)), 1e-12)

all_vertices: list[np.ndarray] = []
all_faces: list[np.ndarray] = []
all_colors: list[np.ndarray] = []
offset = 0
palette = np.array(
    [[0.95, 0.57, 0.49], [0.87, 0.36, 0.29], [0.98, 0.70, 0.62], [0.76, 0.25, 0.20]],
    dtype=np.float32,
)

for number, level in enumerate(a.levels):
    if not alpha.min() < level < alpha.max():
        continue
    verts, faces, _, _ = marching_cubes(
        alpha, level=level, spacing=(dx, dy, dz), allow_degenerate=False
    )
    grid_index = np.column_stack((verts[:, 0] / dx, verts[:, 1] / dy, verts[:, 2] / dz)).T
    strength = np.clip(map_coordinates(omag, grid_index, order=1, mode="nearest") / oscale, 0, 1)
    base = palette[number % len(palette)]
    # Vorticity gives a restrained pearly modulation while keeping each
    # material sheet recognizably salmon/coral.
    shade = 0.72 + 0.28 * strength[:, None]
    colors = np.clip(base[None, :] * shade + 0.08 * (1 - shade), 0, 1)
    verts[:, 0] = (verts[:, 0] + x[0] - 0.5 * (x[0] + x[-1])) / lx
    verts[:, 1] = (verts[:, 1] + y[0] - 0.5 * (y[0] + y[-1])) / lx
    verts[:, 2] = a.vertical_exaggeration * (
        verts[:, 2] + z[0] - 0.5 * (z[0] + z[-1])
    ) / lx
    all_vertices.append(verts.astype("<f4"))
    all_faces.append((faces + offset).astype("<u4"))
    all_colors.append(colors.astype("<f4"))
    offset += len(verts)

vertices = np.concatenate(all_vertices)
faces = np.concatenate(all_faces)
colors = np.concatenate(all_colors)
a.stem.parent.mkdir(parents=True, exist_ok=True)
paths = {
    "positions": a.stem.with_suffix(".flow.pos.f32"),
    "colors": a.stem.with_suffix(".flow.color.f32"),
    "indices": a.stem.with_suffix(".flow.index.u32"),
}
vertices.tofile(paths["positions"])
colors.tofile(paths["colors"])
faces.ravel().tofile(paths["indices"])
meta = {
    "time_s": t,
    "surface_type": "advected Clebsch material sheets",
    "levels": a.levels,
    "vertices": len(vertices),
    "triangles": len(faces),
    "positions": paths["positions"].name,
    "colors": paths["colors"].name,
    "indices": paths["indices"].name,
    "source_hdf5": str(a.input),
    "vertical_exaggeration": a.vertical_exaggeration,
    "script": "scripts/pack_khi3d_flowmap_surface.py",
}
a.stem.with_suffix(".q.json").write_text(json.dumps(meta, indent=2) + "\n")
print(json.dumps({"vertices": len(vertices), "triangles": len(faces), "levels": a.levels}))
