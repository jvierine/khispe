# Original Ri=0.08 simulation, rendered by Rustplotlib

Run `plots/rustplotlib/render.sh` from the project. The renderer uses a pinned Rustplotlib Git revision. HDF5_DIR selects the installed HDF5 runtime.

The HDF5 input includes the original density and neutral velocity snapshots.
The renderer computes the model's collisional magnetic ion vertical velocity
using precisely the solver's algebraic formula. It never substitutes neutral
vertical velocity for ion drift. Velocity is relative to the model's descending
phase frame; prescribed phase descent is not treated as an actual vertical wind.

Two MP4s use 721 actual 10-second checkpoints, without temporal interpolation.
Camera: mean zonal advection U0=-100 m/s. Coordinates x-U0*t and z-zc(t).
No spatial mask, background replacement, or endpoint clipping is used.

The scientific figures use real LaTeX labels at 10 pt and 178 mm physical width.
Animation page size is explicitly fixed across frames to avoid crop jitter.
Static source PDFs at t=0,60,120 minutes and HDF5/rendering provenance retained.

To reproduce the memo panels, use the included small simulated HDF5 fixture:

```sh
cargo run --release --manifest-path plots/rustplotlib/Cargo.toml --bin memo_figures -- OUTPUT_FIGURE_DIRECTORY
```

The panel sequence contains 60 seconds, 20 minutes, 40 minutes, and 80 minutes. White
arrows have a black outline. Provenance belongs in LaTeX comments, not rendered.
The fixture contains derived simulated fields, not raw SYISR observations.
