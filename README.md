# KHI sporadic-E vortex merging

Rust/MPI simulations focused on reproducing **Kelvin–Helmholtz vortex pairing
and merging**, motivated by the [SYISR observations of Wang et al. (2025)](https://doi.org/10.1029/2025GL116932).

The current Ri=0.08 case combines a stratified neutral shear layer with
magnetic, collision-dominated ion transport. A primary disturbance and its
spatial subharmonic seed neighbouring rolls that can pair into larger structures.
This is an idealized mechanism study, not yet a fit to the measurements.

For a more extensive model with more realistic large-scale dynamics, turning
wind shears and plasma electrodynamics, see [Hysell et al. (2024)](https://doi.org/10.1029/2024JA033088).
Our primary goal here is to replicate the vortex merging.

Related work:

- [Bernhardt (2002)](https://doi.org/10.1016/S1364-6826(02)00086-X), *The modulation of sporadic-E layers by Kelvin–Helmholtz billows in the neutral atmosphere*.
- [Hysell and Larsen (2021)](https://doi.org/10.1029/2021JA029257), *VHF Imaging Radar Observations and Theory of Banded Midlatitude Sporadic E Ionization Layers*.
- [Hysell et al. (2024)](https://doi.org/10.1029/2024JA033088), *Observations and Model of Subauroral Sporadic E Layer Irregularities Driven by Turning Shears and Dynamic Instability*.

- Solver: [`src/main.rs`](src/main.rs)
- Selected run: [`configs/ri-sweep-short-2h/ri008.json`](configs/ri-sweep-short-2h/ri008.json)
- Figures and animations: [`plots/rustplotlib`](plots/rustplotlib), using Rustplotlib

Build with `cargo build --release` using Rust, MPI, libclang and HDF5 development
libraries. Run with `mpirun -np 128 target/release/syisr-khi configs/ri-sweep-short-2h/ri008.json output/ri008`.
Scientific outputs use HDF5; large outputs are kept outside Git.
