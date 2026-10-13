#!/usr/bin/env python3
"""Create the coarse nonlinear Figure 3 parameter ensemble."""

from __future__ import annotations

import json
from pathlib import Path


BASE = {
    "nx": 256,
    "nz": 128,
    "lx_m": 86000.0,
    "lz_m": 6000.0,
    "altitude_m": 116000.0,
    "shear_m": 1060.0,
    "wind_m_s": 50.0,
    "advection_m_s": -100.0,
    "n2_s2": 0.0004,
    "viscosity_m2_s": 30.0,
    "diffusivity_m2_s": 5.0,
    "dt_s": 0.3,
    "end_s": 711.0,
    "output_s": 50.8,
    "seed_m_s": 0.5,
    "seed_wavelength_m": 21500.0,
    "ion_mode": "tracer",
    "b_t": 0.000035,
    "dip_deg": 20.0,
    "nu_in_s": 300.0,
    "ion_mass_amu": 55.845,
    "layer_sigma_m": 255.0,
    "peak_m3": 100000000000.0,
    "poisson_rtol": 1e-9,
    "validate_only": False,
}

CASES = {
    "ri012": {"wind_m_s": 61.2},
    "ri024": {"wind_m_s": 43.3},
    "seed020": {"seed_m_s": 0.2},
    "seed100": {"seed_m_s": 1.0},
    "shear0800": {"shear_m": 800.0, "wind_m_s": 37.7},
    "shear1400": {"shear_m": 1400.0, "wind_m_s": 66.0},
    "nu0010": {"viscosity_m2_s": 10.0},
    "nu0060": {"viscosity_m2_s": 60.0},
}


def main() -> None:
    directory = Path("configs/figure3-ensemble")
    directory.mkdir(parents=True, exist_ok=True)
    for name, changes in CASES.items():
        config = BASE | changes
        (directory / f"{name}.json").write_text(
            json.dumps(config, indent=2) + "\n", encoding="utf-8"
        )
    print(f"wrote {len(CASES)} configurations to {directory}")


if __name__ == "__main__":
    main()
