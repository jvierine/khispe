#!/usr/bin/env python3
"""Pack evenly spaced meridional electron-density cuts for WebGL.

The HDF5 source is retained as the quantitative product.  This script only
quantizes a fixed log10(n_e) scale into compact 8-bit display textures.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


p = argparse.ArgumentParser()
p.add_argument("input", type=Path)
p.add_argument("stem", type=Path)
p.add_argument("--cuts", type=int, default=9)
p.add_argument("--log-min", type=float, default=8.0)
p.add_argument("--log-max", type=float, default=10.8)
p.add_argument("--vertical-exaggeration", type=float, default=5.0)
a = p.parse_args()

with h5py.File(a.input, "r") as f:
    ne = f["ion_density_m-3"][:]
    x = f["x_m"][:]
    y = f["y_m"][:]
    z = f["altitude_m"][:]
    t = float(f.attrs["time_s"])

indices = np.rint(np.linspace(0, len(y) - 1, a.cuts)).astype(int)
log_ne = np.log10(np.maximum(ne[:, indices, :], 1.0))
scaled = np.clip((log_ne - a.log_min) / (a.log_max - a.log_min), 0, 1)
# WebGL texture layout is cut, z(row), x(column).
textures = np.rint(255 * scaled.transpose(1, 2, 0)).astype(np.uint8)
a.stem.parent.mkdir(parents=True, exist_ok=True)
raw = a.stem.with_suffix(".necuts.u8")
textures.tofile(raw)
meta = {
    "time_s": t,
    "nx": len(x),
    "nz": len(z),
    "cuts": len(indices),
    "cut_indices": indices.tolist(),
    "cut_y_m": y[indices].tolist(),
    "domain_x_m": float(np.ptp(x)),
    "domain_y_m": float(np.ptp(y)),
    "domain_z_m": float(np.ptp(z)),
    "altitude_center_m": float(np.mean(z)),
    "log10_ne_min": a.log_min,
    "log10_ne_max": a.log_max,
    "vertical_exaggeration": a.vertical_exaggeration,
    "data": raw.name,
    "source_hdf5": str(a.input),
    "script": "scripts/pack_khi3d_ne_cuts.py",
}
a.stem.with_suffix(".necuts.json").write_text(json.dumps(meta, indent=2) + "\n")
print(json.dumps({"time_s": t, "cuts": len(indices), "bytes": textures.nbytes}))
