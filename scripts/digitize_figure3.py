#!/usr/bin/env python3
"""Digitize the published, medium-resolution Figure 3 into an HDF5 target.

This is an image-derived target, not the released Wang et al. measurement
array.  It preserves the source pixels, physical-axis calibration, a colorbar
inversion, and conservative pixel/color uncertainties so image-space fitting
can proceed while the ScienceDB transfer is unavailable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from PIL import Image


TIMES_LT = (
    "06:41:38", "06:42:28", "06:43:19",
    "06:44:09", "06:45:00", "06:45:50",
    "06:46:40", "06:47:31", "06:48:21",
    "06:49:12", "06:50:02", "06:50:52",
    "06:51:43", "06:52:33", "06:53:24",
)

# Pixel bounds were read from the publisher's 500 x 434 PNG.  Bounds are
# inclusive at the left/top and exclusive at the right/bottom.
PANEL_X0 = (30, 180, 330)
PANEL_Y0 = (10, 92, 174, 257, 339)
PANEL_WIDTH = 140
PANEL_HEIGHT = 69


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    image = np.asarray(Image.open(args.input).convert("RGB"), dtype=np.uint8)
    if image.shape != (434, 500, 3):
        raise ValueError(f"expected 500x434 RGB source, got {image.shape}")

    # The visible colorbar is two pixels wide at x=476:478 and spans the five
    # panel rows.  Median the two columns and map its top/bottom to the printed
    # limits 11.4 and 10.7 log10(m^-3).
    cb_y0, cb_y1 = 93, 407
    palette_rgb = np.median(image[cb_y0:cb_y1, 476:478, :], axis=1)
    palette_value = np.linspace(11.4, 10.7, palette_rgb.shape[0])

    rgb = np.empty((15, PANEL_HEIGHT, PANEL_WIDTH, 3), dtype=np.uint8)
    density = np.full((15, PANEL_HEIGHT, PANEL_WIDTH), np.nan, dtype=np.float32)
    color_distance = np.full_like(density, np.nan)

    panel = 0
    for y0 in PANEL_Y0:
        for x0 in PANEL_X0:
            crop = image[y0:y0 + PANEL_HEIGHT, x0:x0 + PANEL_WIDTH, :]
            rgb[panel] = crop

            # Invert the published color scale.  Black contour/text pixels and
            # anti-aliased marks are rejected by their distance from the bar.
            flat = crop.reshape(-1, 3).astype(np.float32)
            delta = flat[:, None, :] - palette_rgb[None, :, :]
            distance = np.sqrt(np.sum(delta * delta, axis=2))
            nearest = np.argmin(distance, axis=1)
            dmin = distance[np.arange(distance.shape[0]), nearest]
            values = palette_value[nearest].astype(np.float32)
            values[dmin > 32.0] = np.nan
            density[panel] = values.reshape(PANEL_HEIGHT, PANEL_WIDTH)
            color_distance[panel] = dmin.reshape(PANEL_HEIGHT, PANEL_WIDTH)
            panel += 1

    distance_km = np.linspace(-32.0, 54.0, PANEL_WIDTH, dtype=np.float32)
    altitude_km = np.linspace(119.0, 113.0, PANEL_HEIGHT, dtype=np.float32)
    time_seconds = np.array(
        [
            int(t[0:2]) * 3600 + int(t[3:5]) * 60 + int(t[6:8])
            for t in TIMES_LT
        ],
        dtype=np.int32,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["product"] = "figure-derived observational target"
        h5.attrs["raw_measurement_data"] = False
        h5.attrs["article_doi"] = "10.1029/2025GL116932"
        h5.attrs["source_asset"] = "grl71225-fig-0003-m.png"
        h5.attrs["source_sha256"] = (
            "063c164fbfd7d81953a969763eb36091ede414f7df76783e23fe114fac4d25a0"
        )
        h5.attrs["digitizer"] = "scripts/digitize_figure3.py"
        h5.attrs["caveat"] = (
            "Colorbar-inverted publisher pixels; black contours, labels, and "
            "poor color matches are masked. Do not call this the released array."
        )

        g = h5.create_group("figure3")
        g.create_dataset("rgb", data=rgb, compression="gzip", shuffle=True)
        g.create_dataset("log10_ne_m3", data=density, compression="gzip", shuffle=True)
        g.create_dataset(
            "palette_rgb_distance", data=color_distance, compression="gzip", shuffle=True
        )
        g.create_dataset("distance_km", data=distance_km)
        g.create_dataset("altitude_km", data=altitude_km)
        g.create_dataset("time_seconds_lt", data=time_seconds)
        g.create_dataset("time_text_lt", data=np.asarray(TIMES_LT, dtype="S8"))
        g.create_dataset("panel_x0_pixel", data=np.asarray(PANEL_X0 * 5, dtype=np.int16))
        g.create_dataset(
            "panel_y0_pixel",
            data=np.repeat(np.asarray(PANEL_Y0, dtype=np.int16), 3),
        )

        u = h5.create_group("uncertainty")
        u.attrs["distance_pixel_km"] = float(86.0 / (PANEL_WIDTH - 1))
        u.attrs["altitude_pixel_km"] = float(6.0 / (PANEL_HEIGHT - 1))
        u.attrs["time_quantization_s"] = 1.0
        u.attrs["color_rejection_rgb_distance"] = 32.0
        u.attrs["colorbar_step_log10_ne"] = float(0.7 / (palette_rgb.shape[0] - 1))
        u.attrs["recommended_log10_ne_floor"] = 0.02
        u.attrs["axis_registration_pixels"] = 1.0

    finite = np.isfinite(density)
    print(
        f"wrote {args.output}: {finite.sum()} valid pixels / {finite.size} "
        f"({100.0 * finite.mean():.1f}%)"
    )


if __name__ == "__main__":
    main()
