# Initial Revontuli checks, 12 October 2026

Canonical source: `/home/j/src/syisr-khi-rust` on `revontuli.uit.no`.
Rust 1.93.1, MPICH 4.3.2, HDF5 1.14.6. Rust tests and simulations ran there.

`scripts/validate-revontuli.sh`, `src/bin/check-validation.rs`:

- Three unit tests: flux upwinding, limiter, ion Lorentz/drag force balance.
- Distributed manufactured Poisson solution: max error about 1.1e-13
  (32×48) and 4.4e-13 (64×96) relative to the discrete reference.
- Continuous-reference Poisson error reduced by 3.9923 when doubling resolution,
  consistent with second-order spatial truncation.
- 1/2/4-rank five-second runs: max field differences below 3.7e-15 after
  normalization by each field's maximum amplitude.
- Maximum ion-mass relative error in those runs: 3.18e-15.
- Unseeded, inviscid, tracer control: no vertical wind; no tracer change.

The 128×128 four-rank pilot completed 600 simulated seconds in 77.9 wall seconds,
writing 61 HDF5 snapshots. Neutral vertical-wind RMS grew from 0.156 to 12.565 m/s
(factor 80.45). Pilot ion-mass relative error was -1.81e-15. The GIF renderer is
`src/bin/render-pilot.rs`; the visual is a diagnostic of this exact pilot.

These are engineering and initial instability checks. They do **not** establish
nonlinear grid convergence, correct billow growth rates, magnetic/electric
closure, or agreement with Wang's measurements. The model remains 2D with E=0
ion mobility; the next physical milestones are 3D turning shear, polarization
and current continuity, and the SYISR measurement operator. Sanya parameters
and an event-matched neutral shear remain to be constrained.
