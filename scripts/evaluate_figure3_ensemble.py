#!/usr/bin/env python3
"""Evaluate every coarse Figure 3 KHI ensemble member into one HDF5 summary."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

from evaluate_figure3_khi import observed_ridge, simulated_ridge


def evaluate(xo, to, zo, directory: Path):
    xs, _, zs = simulated_ridge(directory)
    central = (xo >= -25.0) & (xo <= 48.0)
    train_panels = np.arange(15) < 6
    valid = np.isfinite(zo) & central[None, :]
    target_x = xo + 32.0
    best = None
    for shift in range(xs.size):
        shifted = np.roll(zs, shift, axis=1)
        raw = np.vstack([np.interp(target_x, xs, row) for row in shifted])
        mask = valid & train_panels[:, None]
        design = np.column_stack(
            [np.ones(mask.sum()), np.broadcast_to(to[:, None] / 600.0, zo.shape)[mask]]
        )
        offset = np.linalg.lstsq(design, (zo - raw)[mask], rcond=None)[0]
        predicted = raw + offset[0] + offset[1] * to[:, None] / 600.0
        score = np.sqrt(np.mean((predicted[mask] - zo[mask]) ** 2))
        if best is None or score < best[0]:
            best = score, shift, offset, predicted, xs
    train, shift, offset, predicted, xs = best
    heldout = valid & (~train_panels[:, None])
    test = np.sqrt(np.mean((predicted[heldout] - zo[heldout]) ** 2))
    return train, test, shift, offset, predicted, valid, xs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", type=Path)
    parser.add_argument("ensemble", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    xo, to, zo = observed_ridge(args.target)
    cases = sorted(path for path in args.ensemble.iterdir() if path.is_dir())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["status"] = "figure-derived coarse parameter search"
        h5.attrs["article_doi"] = "10.1029/2025GL116932"
        h5.attrs["training_panels"] = "a-f"
        h5.attrs["heldout_panels"] = "g-o"
        h5.create_dataset("distance_km", data=xo)
        h5.create_dataset("time_s", data=to)
        h5.create_dataset("observed_ridge_km", data=zo, compression="gzip")
        ranking = []
        for case in cases:
            train, test, shift, offset, predicted, valid, xs = evaluate(xo, to, zo, case)
            group = h5.create_group(case.name)
            group.attrs["train_rmse_km"] = train
            group.attrs["heldout_rmse_km"] = test
            group.attrs["cyclic_shift_cells"] = shift
            group.attrs["cyclic_shift_km"] = shift * (xs[1] - xs[0])
            group.attrs["vertical_offset_km"] = offset[0]
            group.attrs["vertical_trend_km_per_600s"] = offset[1]
            group.create_dataset("predicted_ridge_km", data=predicted, compression="gzip")
            group.create_dataset("residual_km", data=predicted - zo, compression="gzip")
            group.create_dataset("valid_mask", data=valid.astype(np.uint8), compression="gzip")
            ranking.append((test, train, case.name))
        ranking.sort()
        text = "\n".join(f"{name} {train:.9f} {test:.9f}" for test, train, name in ranking)
        h5.create_dataset("ranking_text", data=np.frombuffer(text.encode(), dtype=np.uint8))
    for test, train, name in ranking:
        print(f"{name:12s} train={train:.6f} km heldout={test:.6f} km")


if __name__ == "__main__":
    main()
