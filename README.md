# SYISR KHI: Rust/MPI pilot


This first executable is a **2D x-z stratified Boussinesq shear-layer benchmark**,
not a replication of Wang's event or Hysell's 3D model. MPI divides
the east-west dimension into equal slabs; finite-difference and flux operators
exchange halos. A distributed conjugate-gradient pressure/streamfunction solve
uses collective residual reductions. Only snapshots are gathered to rank zero
for serial HDF5 output. Each timestep runs on all ranks; this is not a collection
of independent serial simulations.

## Equations and approximations

Coordinates: x east, y north (no spatial y dependence yet), z up. With
u=U(z)+dψ/dz, w=-dψ/dx, and ω=dw/dx-du/dz:

```
dω/dt + div(ω [u,w]) = db/dx + ν laplacian(ω)
db/dt + div(b [u,w]) = -N² w + κ laplacian(b)
-laplacian(ψ) = ω - ω_background
dn/dt + div(n [v_ix,v_iz]) = D laplacian(n)
```

The tanh base shear is not artificially restored: its viscous evolution enters
the perturbation vorticity equation. Buoyancy includes stable stratification.
Horizontal boundaries are periodic. The vertical streamfunction is Dirichlet
zero; scalar transport is zero-flux, diffusion zero-gradient. This is a bounded
layer benchmark; walls must be moved farther away in boundary sensitivity runs.
MUSCL/minmod conservative fluxes, SSP RK2, and advective/diffusive/buoyancy
timestep bounds are used. No density clipping hides mass errors.

`ion_mode=tracer` is an explicit control. `collisional` uses the steady
Lorentz/drag mobility response for Fe+ in a prescribed magnetic field and E=0:
ν_in(v_i-u)=Ω_i v_i×bhat. It is a **prescribed-field diagnostic approximation**:
the y velocity is included algebraically, but polarization, electron momentum,
current continuity, chemistry, ambipolar diffusion, and ion feedback are absent.
Therefore the ion snapshots cannot yet be interpreted as quantitative electron
density or SYISR observables. Next plasma milestone: solve current continuity
with magnetized electron/ion mobilities before claiming layer compression.

## Build and run on Revontuli

```sh
scripts/build-revontuli.sh
scripts/validate-revontuli.sh
mpirun -np 4 target/release/syisr-khi configs/pilot.json output/pilot
target/release/render-pilot output/pilot
```

Dependencies: Rust, MPI C development files, libclang, HDF5 development files.
`nx` must be divisible by the rank count. Snapshots store x/altitude in metres,
time in seconds, neutral velocity, buoyancy, perturbation vorticity, ion density,
and the complete configuration. `diagnostics.json` contains scalar diagnostics.
Spatial numerical data are HDF5, not CSV.

On Revontuli, matching MPICH headers were extracted from Ubuntu's
`libmpich-dev` package into `/home/j/.local/opt/mpich-dev`. The build script
uses `scripts/mpicc-local` with those headers and the existing system runtime;
no system MPI installation was changed. The GIF is rendered directly in Rust
by `src/bin/render-pilot.rs`, with a fixed colour scale across all frames.

Defaults are exploratory parameters: 12×6 km, 115 km center, 350 m shear scale,
±50 m/s shear, -30 m/s translation, N²=4e-4 s^-2. They are **not retrieved
SYISR winds**. The center Richardson number is approximately 0.0196; Ri<1/4
is a susceptibility indicator, not sufficient proof of instability.

## Literature and acceptance targets

- Wang et al. (2025), *Direct Kelvin–Helmholtz Instability Detection at the
  Ionospheric E Region Electron Density*, https://doi.org/10.1029/2025GL116932.
  Target: 23 September 2024 SYISR event, evolving Es near 115 km; compare layer
  morphology and advection only after applying the actual beam geometry,
  spatial response and temporal sampling. Wind shear at 115 km was not measured.
- Bernhardt (2002), *The modulation of sporadic-E layers by Kelvin–Helmholtz
  billows in the neutral atmosphere*, https://doi.org/10.1016/S1364-6826(02)00086-X.
- Hysell and Larsen (2021), https://doi.org/10.1029/2021JA029257.
- Hysell et al. (2024), https://doi.org/10.1029/2024JA033088.
  Reference distinction: incompressible 3D neutrals plus extended MHD, turning
  shear, and dynamically generated polarization fields. Its subauroral magnetic
  inclination and collisional parameters must not simply be copied to Sanya.

Before event fitting: manufactured Poisson convergence, mass conservation,
1/2/4-rank equivalence, seed-free control, resolution and wall-distance tests,
growth-rate/Richardson-number scans; then 3D turning shear and plasma closure.


## Current Ri=0.08 memo and rendering

The current memo is in the separate Overleaf project
The simulation and dedicated analysis source are kept in this private repository.
Configuration: `configs/ri-sweep-short-2h/ri008.json`.
Neutral and plasma model: `src/main.rs`. New paper and memo plots use
Rustplotlib (`plots/rustplotlib`); older plotting scripts are historical.
Script provenance belongs in LaTeX comments, not visible captions.
Large HDF5 outputs and animations are not included in Git.

The newly updated 3-D solver is a pilot under test; its configuration is
`configs/ri008-3d-meridional.json`. No validated 3-D result is claimed.
