#!/usr/bin/env python3
"""Pack a KHI HDF5 checkpoint into one WebGL RGBA8 3-D texture.

This script performs data reduction only. Rendering is done entirely by the
WebGL2 ray marcher in webgl/khi3d-volume.html.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def robust_encode(a: np.ndarray, q: float = 99.5) -> tuple[np.ndarray, float]:
    scale = float(np.percentile(np.abs(a[np.isfinite(a)]), q))
    scale = max(scale, np.finfo(np.float64).tiny)
    out = np.clip(0.5 + 0.5 * a / scale, 0.0, 1.0)
    return np.rint(255.0 * out).astype(np.uint8), scale


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=Path)
    p.add_argument("output_stem", type=Path)
    args = p.parse_args()

    with h5py.File(args.input, "r") as f:
        wx = f["omega_x_s-1"][:]
        wy = f["omega_y_s-1"][:]
        wz = f["omega_z_s-1"][:]
        ne = f["ion_density_m-3"][:]
        w = f["w_m_s"][:]
        # Remove the horizontally uniform ion profile; retain the evolving sheet.
        ne0 = np.median(ne, axis=(0, 1), keepdims=True)
        dne = (ne - ne0) / np.maximum(ne0, 1.0)
        # Isolate instability vorticity from the horizontally uniform base
        # shear. Otherwise volume integration makes the whole shear slab opaque.
        pwx = wx - np.mean(wx, axis=(0, 1), keepdims=True)
        pwy = wy - np.mean(wy, axis=(0, 1), keepdims=True)
        pwz = wz - np.mean(wz, axis=(0, 1), keepdims=True)
        vort, vort_scale = robust_encode(pwy)
        ion, ion_scale = robust_encode(dne)
        vert, w_scale = robust_encode(w)
        a1, a1_scale = robust_encode(f["clebsch_alpha1"][:], 100.0)
        b1, b1_scale = robust_encode(f["clebsch_beta1_m2_s-1"][:], 100.0)
        a2, a2_scale = robust_encode(f["clebsch_alpha2"][:], 100.0)
        b2, b2_scale = robust_encode(f["clebsch_beta2_m2_s-1"][:], 100.0)
        # Vorticity magnitude is useful for opacity independent of sign.
        omag = np.sqrt(pwx**2 + pwy**2 + pwz**2)
        mag_scale = float(np.percentile(omag, 99.5))
        mag = np.rint(255 * np.clip(omag / max(mag_scale, 1e-30), 0, 1)).astype(np.uint8)
        x = f["x_m"][:]
        y = f["y_m"][:]
        z = f["altitude_m"][:]
        time_s = float(f.attrs["time_s"])

    # WebGL expects x as the fastest varying texture coordinate.
    # The primary texture is a vortex map: signed spanwise vorticity, two
    # materially advected Clebsch coordinates, and total-vorticity opacity.
    vortex_rgba = np.stack((vort, a1, b1, mag), axis=-1).transpose(2, 1, 0, 3).copy()
    # Auxiliary texture carries the observable layer and second Clebsch pair.
    aux_rgba = np.stack((ion, vert, a2, b2), axis=-1).transpose(2, 1, 0, 3).copy()
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    vortex_path = args.output_stem.with_suffix(".vortex.rgba8")
    aux_path = args.output_stem.with_suffix(".aux.rgba8")
    meta_path = args.output_stem.with_suffix(".json")
    vortex_rgba.tofile(vortex_path)
    aux_rgba.tofile(aux_path)
    meta = {
        "vortex_texture": vortex_path.name,
        "aux_texture": aux_path.name,
        "dims": [int(wy.shape[0]), int(wy.shape[1]), int(wy.shape[2])],
        "domain_km": [float(np.ptp(x) / 1e3), float(np.ptp(y) / 1e3), float(np.ptp(z) / 1e3)],
        "time_s": time_s,
        "scales": {
            "omega_y_s-1": vort_scale,
            "relative_ion_density": ion_scale,
            "w_m_s": w_scale,
            "omega_magnitude_s-1": mag_scale,
            "clebsch_alpha1": a1_scale,
            "clebsch_beta1_m2_s-1": b1_scale,
            "clebsch_alpha2": a2_scale,
            "clebsch_beta2_m2_s-1": b2_scale,
        },
        "source_hdf5": str(args.input),
        "layout": {
            "vortex": "RGBA = signed omega_y, Clebsch alpha1, Clebsch beta1, omega magnitude",
            "aux": "RGBA = signed ion anomaly, signed w, Clebsch alpha2, Clebsch beta2",
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(meta_path)


if __name__ == "__main__":
    main()
