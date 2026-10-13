#!/usr/bin/env python3
"""Render a smooth, strictly 2-D KHI vorticity or electron-density animation.

Intermediate display frames are formed on the fixed solver grid.  Before
blending adjacent checkpoints, each field is Fourier-shifted to the
intermediate time using the configured periodic-domain advection velocity.
This removes the dominant translational jump without inventing 3-D geometry.
The HDF5 checkpoints remain the quantitative states.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import fourier_shift


VORTICITY_COLORS = np.array(
    [
        [5, 20, 110],
        [8, 92, 238],
        [22, 190, 244],
        [104, 231, 189],
        [246, 239, 100],
        [255, 171, 58],
        [239, 78, 44],
        [136, 29, 22],
        [43, 10, 12],
    ],
    dtype=float,
)

NE_COLORS = np.array(
    [
        [3, 5, 31],
        [31, 14, 92],
        [104, 25, 135],
        [204, 48, 104],
        [247, 118, 59],
        [253, 190, 86],
        [255, 244, 174],
    ],
    dtype=float,
)


def translate_periodic(field: np.ndarray, cells: float) -> np.ndarray:
    spectrum = np.fft.fftn(field)
    return np.fft.ifftn(fourier_shift(spectrum, (cells, 0.0))).real


def palette(colors: np.ndarray) -> np.ndarray:
    x = np.linspace(0, 1, len(colors))
    q = np.linspace(0, 1, 240)
    rgb = np.column_stack([np.interp(q, x, colors[:, k]) for k in range(3)])
    extra = np.tile(rgb[-1], (16, 1))
    extra[-1] = [255, 244, 220]
    return np.vstack((rgb, extra)).clip(0, 255).astype(np.uint8)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("frames", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--inbetweens", type=int, default=6)
    p.add_argument("--duration-ms", type=int, default=42)
    p.add_argument("--width", type=int, default=1280)
    p.add_argument("--height", type=int, default=720)
    p.add_argument("--field", choices=("vorticity", "ne"), default="vorticity")
    a = p.parse_args()
    paths = sorted(a.frames.glob("frame_*.h5"))
    if len(paths) < 2:
        raise SystemExit("need at least two HDF5 checkpoints")

    fields, times = [], []
    with h5py.File(paths[0], "r") as h:
        cfg = json.loads(bytes(h["config_json"][:]).decode())
    with h5py.File(paths[0], "r") as h:
        nz = h["altitude_m"].shape[0]
    z0 = np.linspace(-cfg["lz_m"] / 2, cfg["lz_m"] / 2, nz)
    base = (cfg["wind_m_s"] / cfg["shear_m"]) / np.cosh(z0 / cfg["shear_m"]) ** 2
    for path in paths:
        with h5py.File(path, "r") as h:
            if a.field == "vorticity":
                fields.append(h["omega_s-1"][:] - base[None, :])
            else:
                fields.append(np.log10(np.maximum(h["ion_density_m-3"][:], 1.0)))
            times.append(float(h.attrs["time_s"]))
    if a.field == "vorticity":
        negative = np.concatenate([f[f < 0] for f in fields])
        vmin, vmax = float(np.percentile(negative, 0.35)), 0.0
        colors, colorbar_label = VORTICITY_COLORS, "Perturbation vorticity [s^-1]"
        tick_labels = (f"{vmin:.3f}", f"{0.5*vmin:.3f}", "0")
    else:
        vmin, vmax = 10.0, 11.0
        colors, colorbar_label = NE_COLORS, "log10 n_e [m^-3]"
        tick_labels = ("10", "10.5", "11")
    dx = cfg["lx_m"] / fields[0].shape[0]
    pal = palette(colors)
    pal_flat = pal.ravel().tolist()
    font_path = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
    font = ImageFont.truetype(font_path, 24) if Path(font_path).exists() else ImageFont.load_default()
    math_font_path = "/System/Library/Fonts/Supplemental/STIXTwoText.ttf"
    small_path = math_font_path if Path(math_font_path).exists() else font_path
    small = ImageFont.truetype(small_path, 17) if Path(small_path).exists() else ImageFont.load_default()
    script_font = ImageFont.truetype(small_path, 12) if Path(small_path).exists() else ImageFont.load_default()
    images: list[Image.Image] = []

    def render(field: np.ndarray, time_s: float) -> Image.Image:
        index = np.rint(239 * np.clip((field - vmin) / (vmax - vmin), 0, 1)).astype(np.uint8)
        image = Image.fromarray(np.flipud(index.T), mode="P")
        image.putpalette(pal_flat)
        image = image.resize((a.width, a.height), Image.Resampling.BICUBIC).convert("RGB")
        draw = ImageDraw.Draw(image)
        label = f"{time_s / 60:04.1f} min"
        box = draw.textbbox((0, 0), label, font=font)
        tw, th = box[2] - box[0], box[3] - box[1]
        draw.rounded_rectangle((a.width - tw - 43, 20, a.width - 19, 37 + th), radius=9, fill=(15, 8, 10))
        draw.text((a.width - tw - 31, 25), label, font=font, fill=(255, 244, 220))

        # Fixed physical colorbar embedded in every frame.
        bar_w, bar_h = 560, 14
        bar_x, bar_y = (a.width - bar_w) // 2, a.height - 57
        gradient = np.tile(np.arange(240, dtype=np.uint8), (bar_h, 1))
        gradient = Image.fromarray(gradient, mode="P"); gradient.putpalette(pal_flat)
        gradient = gradient.resize((bar_w, bar_h), Image.Resampling.BILINEAR).convert("RGB")
        draw.rounded_rectangle((bar_x - 14, bar_y - 27, bar_x + bar_w + 14, a.height - 8), radius=9, fill=(10, 7, 10))
        image.paste(gradient, (bar_x, bar_y))
        draw = ImageDraw.Draw(image)
        if a.field == "ne":
            # Build log_10 n_e [m^-3] from positioned runs so the scientific
            # notation remains correct even when a font lacks Unicode scripts.
            runs = (("log", small, 0), ("10", script_font, 7), (" n", small, 0),
                    ("e", script_font, 7), (" [m", small, 0), ("-3", script_font, -4),
                    ("]", small, 0))
            title_w = sum(draw.textlength(text, font=run_font) for text, run_font, _ in runs)
            cursor = (a.width - title_w) / 2
            for text, run_font, dy in runs:
                draw.text((cursor, bar_y - 23 + dy), text, font=run_font, fill=(255, 240, 220))
                cursor += draw.textlength(text, font=run_font)
        else:
            title_box = draw.textbbox((0, 0), colorbar_label, font=small)
            title_w = title_box[2] - title_box[0]
            draw.text(((a.width - title_w) // 2, bar_y - 23), colorbar_label, font=small, fill=(255, 240, 220))
        for xpos, text, anchor in ((bar_x, tick_labels[0], "la"), (bar_x + bar_w // 2, tick_labels[1], "ma"), (bar_x + bar_w, tick_labels[2], "ra")):
            draw.text((xpos, bar_y + 17), text, font=small, fill=(255, 240, 220), anchor=anchor)
        return image

    for k in range(len(fields) - 1):
        dt = times[k + 1] - times[k]
        shift = cfg["advection_m_s"] * dt / dx
        for j in range(a.inbetweens):
            f = j / a.inbetweens
            smooth = f * f * (3 - 2 * f)
            left = translate_periodic(fields[k], f * shift)
            right = translate_periodic(fields[k + 1], -(1 - f) * shift)
            field = (1 - smooth) * left + smooth * right
            time_s = (1 - smooth) * times[k] + smooth * times[k + 1]
            images.append(render(field, time_s))

    images.append(render(fields[-1], times[-1]))

    a.output.parent.mkdir(parents=True, exist_ok=True)
    master = Image.new("P", (1, 1)); master.putpalette(pal_flat)
    gif = [im.quantize(palette=master, dither=Image.Dither.FLOYDSTEINBERG) for im in images]
    gif[0].save(a.output, save_all=True, append_images=gif[1:], duration=a.duration_ms, loop=0, disposal=2, optimize=False)
    images[-1].save(a.output.with_suffix(".png"))
    provenance = {
        "source_directory": str(a.frames),
        "source_frames": len(paths),
        "display_frames": len(images),
        "interpolation": "periodic Fourier bulk-advection alignment plus smooth scalar blending",
        "inbetweens_per_interval": a.inbetweens,
        "field": a.field,
        "fixed_display_range": [vmin, vmax],
        "colorbar_label": colorbar_label,
        "script": "scripts/render_cat_eye_cinematic_smooth.py",
    }
    a.output.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
