#!/usr/bin/env python3
"""Apply the Figure 3 ridge observation operator to a Rust/MPI KHI run."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def weighted_vertical_moments(
    weight: np.ndarray, altitude_km: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the vertical centroid and standard deviation for each x column."""
    denominator = weight.sum(axis=1)
    ridge = np.divide(
        (weight * altitude_km[None, :, None]).sum(axis=1),
        denominator,
        out=np.full_like(denominator, np.nan, dtype=float),
        where=denominator > 0,
    )
    variance = np.divide(
        (
            weight
            * (altitude_km[None, :, None] - ridge[:, None, :]) ** 2
        ).sum(axis=1),
        denominator,
        out=np.full_like(denominator, np.nan, dtype=float),
        where=denominator > 0,
    )
    return ridge, np.sqrt(np.maximum(variance, 0.0))


def observed_moments(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        d = h5["figure3/log10_ne_m3"][:]
        z = h5["figure3/altitude_km"][:]
        x = h5["figure3/distance_km"][:]
        t = h5["figure3/time_seconds_lt"][:].astype(float)
    weight = np.where(np.isfinite(d), np.maximum(d - 10.95, 0.0), 0.0) ** 2
    ridge, sigma = weighted_vertical_moments(weight, z)
    return x, t - t[0], ridge, sigma


def simulated_moments(
    directory: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    paths = sorted(directory.glob("frame_*.h5"))
    if len(paths) != 15:
        raise ValueError(f"expected 15 simulation frames, found {len(paths)}")
    ridges = []
    sigmas = []
    times = []
    x = None
    for path in paths:
        with h5py.File(path, "r") as h5:
            n = h5["ion_density_m-3"][:]
            altitude = h5["altitude_m"][:] / 1000.0
            if x is None:
                x = h5["x_m"][:] / 1000.0
            times.append(float(h5.attrs["time_s"]))
        denominator = n.sum(axis=1)
        ridge = (n * altitude[None, :]).sum(axis=1) / denominator
        variance = (
            n * (altitude[None, :] - ridge[:, None]) ** 2
        ).sum(axis=1) / denominator
        ridges.append(ridge)
        sigmas.append(np.sqrt(np.maximum(variance, 0.0)))
    return np.asarray(x), np.asarray(times), np.asarray(ridges), np.asarray(sigmas)


def masked_rmse(predicted: np.ndarray, observed: np.ndarray, mask: np.ndarray) -> float:
    selected = mask & np.isfinite(predicted) & np.isfinite(observed)
    return float(np.sqrt(np.mean((predicted[selected] - observed[selected]) ** 2)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", type=Path)
    parser.add_argument("simulation", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    xo, to, zo, sigma_o = observed_moments(args.target)
    xs, ts, zs, sigma_s = simulated_moments(args.simulation)
    central = (xo >= -25.0) & (xo <= 48.0)
    train_panels = np.arange(15) < 6
    valid = np.isfinite(zo) & central[None, :]

    # The simulation is periodic over the same 86-km zonal interval. Optimize
    # one cyclic registration on panels a-f, then freeze it for panels g-o.
    target_x = xo + 32.0
    best = None
    for shift in range(xs.size):
        shifted = np.roll(zs, shift, axis=1)
        predicted = np.vstack([np.interp(target_x, xs, row) for row in shifted])
        m = valid & train_panels[:, None]
        design = np.column_stack(
            [np.ones(m.sum()), np.broadcast_to(to[:, None] / 600.0, zo.shape)[m]]
        )
        offset = np.linalg.lstsq(design, (zo - predicted)[m], rcond=None)[0]
        adjusted = predicted + offset[0] + offset[1] * to[:, None] / 600.0
        rmse = np.sqrt(np.mean((adjusted[m] - zo[m]) ** 2))
        if best is None or rmse < best[0]:
            best = (rmse, shift, offset, adjusted)

    train_rmse, shift, offset, predicted = best
    heldout = valid & (~train_panels[:, None])
    heldout_rmse = np.sqrt(np.mean((predicted[heldout] - zo[heldout]) ** 2))
    residual = predicted - zo

    # Apply the frozen spatial registration to layer width.  No scale or offset
    # is fitted: absolute sigma is a direct model prediction.  A second metric
    # compares growth relative to panel a, which is less sensitive to the
    # threshold bias of the image-derived density observation operator.
    shifted_sigma = np.roll(sigma_s, shift, axis=1)
    predicted_sigma = np.vstack(
        [np.interp(target_x, xs, row) for row in shifted_sigma]
    )
    sigma_valid = valid & np.isfinite(sigma_o) & np.isfinite(predicted_sigma)
    sigma_train = sigma_valid & train_panels[:, None]
    sigma_heldout = sigma_valid & (~train_panels[:, None])
    sigma_train_rmse = masked_rmse(predicted_sigma, sigma_o, sigma_train)
    sigma_heldout_rmse = masked_rmse(predicted_sigma, sigma_o, sigma_heldout)

    sigma_growth_o = sigma_o - sigma_o[0:1, :]
    sigma_growth_p = predicted_sigma - predicted_sigma[0:1, :]
    growth_valid = sigma_valid & np.isfinite(sigma_growth_o) & np.isfinite(sigma_growth_p)
    growth_train = growth_valid & train_panels[:, None]
    growth_heldout = growth_valid & (~train_panels[:, None])
    sigma_growth_train_rmse = masked_rmse(
        sigma_growth_p, sigma_growth_o, growth_train
    )
    sigma_growth_heldout_rmse = masked_rmse(
        sigma_growth_p, sigma_growth_o, growth_heldout
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["status"] = "figure-derived, not a raw-data fit"
        h5.attrs["article_doi"] = "10.1029/2025GL116932"
        h5.attrs["operator"] = (
            "density-weighted ridge and vertical sigma; frozen cyclic registration"
        )
        h5.attrs["observed_weight"] = "max(log10_ne_m3 - 10.95, 0)^2"
        h5.attrs["sigma_definition"] = "weighted vertical standard deviation"
        h5.attrs["sigma_growth_reference"] = "panel a at each horizontal sample"
        h5.attrs["training_panels"] = "a-f"
        h5.attrs["heldout_panels"] = "g-o"
        h5.attrs["train_rmse_km"] = train_rmse
        h5.attrs["heldout_rmse_km"] = heldout_rmse
        h5.attrs["sigma_train_rmse_km"] = sigma_train_rmse
        h5.attrs["sigma_heldout_rmse_km"] = sigma_heldout_rmse
        h5.attrs["sigma_growth_train_rmse_km"] = sigma_growth_train_rmse
        h5.attrs["sigma_growth_heldout_rmse_km"] = sigma_growth_heldout_rmse
        h5.attrs["cyclic_shift_cells"] = shift
        h5.attrs["cyclic_shift_km"] = shift * (xs[1] - xs[0])
        h5.attrs["vertical_offset_km"] = offset[0]
        h5.attrs["vertical_trend_km_per_600s"] = offset[1]
        h5.create_dataset("distance_km", data=xo)
        h5.create_dataset("time_s", data=to)
        h5.create_dataset("observed_ridge_km", data=zo, compression="gzip")
        h5.create_dataset("predicted_ridge_km", data=predicted, compression="gzip")
        h5.create_dataset("residual_km", data=residual, compression="gzip")
        h5.create_dataset("observed_sigma_km", data=sigma_o, compression="gzip")
        h5.create_dataset(
            "predicted_sigma_km", data=predicted_sigma, compression="gzip"
        )
        h5.create_dataset(
            "sigma_residual_km", data=predicted_sigma - sigma_o, compression="gzip"
        )
        h5.create_dataset(
            "observed_sigma_growth_km", data=sigma_growth_o, compression="gzip"
        )
        h5.create_dataset(
            "predicted_sigma_growth_km", data=sigma_growth_p, compression="gzip"
        )
        h5.create_dataset("valid_mask", data=valid.astype(np.uint8), compression="gzip")
    print(
        f"train RMSE={train_rmse:.6f} km; held-out RMSE={heldout_rmse:.6f} km; "
        f"sigma train/held-out={sigma_train_rmse:.6f}/{sigma_heldout_rmse:.6f} km; "
        f"sigma-growth train/held-out="
        f"{sigma_growth_train_rmse:.6f}/{sigma_growth_heldout_rmse:.6f} km; "
        f"shift={shift} cells"
    )


if __name__ == "__main__":
    main()
