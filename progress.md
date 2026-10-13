# SYISR KHI analysis progress

Updated: 2026-10-12T11:37:00Z

## Current phase

The public Wang et al. SYISR dataset for DOI
`10.57760/sciencedb.24829` is still being retrieved. No data transfer is
currently active. In parallel, the open-access published Figures 2--4 from
DOI `10.1029/2025GL116932` are now the immediate observational target for a
clearly labeled figure-derived fit.

## Verified provenance

- ScienceDB dataset UUID: `138b232ffb464558aa0abee8f500a85f`
- Version: `V1`
- License: CC BY-NC-ND 4.0
- Total size: 1,338,826,524 bytes
- The public landing page lists five files with MD5 digests recorded in
  `download-status.json`.

## Completed engineering work

- The Rust/MPI 2-D neutral and collisional-ion pilot builds and passes its
  manufactured Poisson, rank-invariance, conservation, and seed-free-control
  checks on Revontuli.
- These checks establish implementation consistency only. They are not an
  observational fit or physical validation.
- The 15 Figure 3 panels have been digitized to HDF5 with calibrated axes,
  source pixels, colorbar inversion, rejection mask, and explicit pixel/color
  uncertainty. The retained color-valid fraction is 83.4%.
- An eight-rank Rust/MPI search fit panels a--f and evaluated panels g--o as a
  held-out sequence. The comparison is figure-derived and does not use the
  released measurement arrays.

## First figure-derived fit result

| Model | Training RMSE (km) | Held-out RMSE (km) | Training BIC |
|---|---:|---:|---:|
| Descending/sloped background | 0.236 | 0.382 | -2001 |
| Background + one polarized GW | 0.190 | 0.400 | -2262 |
| Background + two polarized GWs | 0.132 | 0.609 | -2737 |

The gravity-wave components reduce error on the six fitted frames, but neither
improves the later nine frames. The two-wave case fails most strongly out of
sample and requires approximately 2 km fitted component amplitudes, far larger
than the approximately 150 m early displacement reported in the paper. Both
searches also approach the 50 km upper bound on projected wavelength rather
than the approximately 22 km mature-billow spacing. Thus this first test does
not support two simultaneous oblique gravity waves as a predictive explanation
of the full sequence. It instead shows that a linear wave description
overfits the early ridge and that nonlinear KHI evolution is needed.

Outputs are `output/figure3-fit/fit.h5` and `fit.json`; the HDF5 file includes
the fitted ridge and physically linked displacement, horizontal velocity,
vertical velocity, and buoyancy fields for each model class.

## First nonlinear KHI comparison

An observation-scaled 512 x 192 neutral/tracer run now covers the full 86 km
Figure 3 zonal interval at approximately 168 m x 31 m resolution. It uses a
21.5 km periodic seed, Ri = 0.180 at the shear center, -100 m/s mean advection,
and 15 frames over 711 s. Eight MPI ranks completed the run in 469 s with a
relative tracer-mass error of 8.3e-15.

The frozen ridge observation operator fits one cyclic registration and an
affine vertical descent using panels a--f, then holds both fixed for panels
g--o. Its RMSE is 0.269 km on the training frames and 0.288 km on the held-out
frames. The held-out error is smaller than for the linear background (0.382
km), one-GW (0.400 km), and two-GW (0.609 km) descriptions. This is meaningful
predictive improvement, but the residual is still about 3.3 vertical pixels
and the training error is worse than the descriptive linear fits. The current
nonlinear run is therefore promising but not yet a satisfactory match.

The model frames remain in `output/figure3-khi/`; the compact, reproducible
comparison is `output/figure3-khi/evaluation.h5`.

Two full-resolution confirmations tested the best coarse-ensemble changes.
Reducing the seed amplitude from 0.5 to 0.2 m/s gave ridge RMSEs of 0.269 km
(training) and 0.281 km (held out). Raising the initial Richardson number from
0.18 to 0.24 gave 0.269 and 0.283 km. These are only small ridge improvements.

The evaluator now also measures the density-weighted vertical layer width
(`sigma`) and its change relative to panel a. This exposed a failure hidden by
the centroid-only score. The digitized high-density core has panel-median
sigma values of 0.137 km initially and roughly 0.19--0.23 km in several later
panels. All three nonlinear runs start at 0.255 km and thicken monotonically by
only about 0.014 km. Their held-out sigma-growth RMSE is 0.172 km, essentially
unchanged by the seed-amplitude and Richardson-number variations. The
threshold-derived core sigma is not the same quantity as the paper's 1--3 km
full billow envelope, so the absolute width offset must not be used to retune
the physical initial thickness blindly. The growth comparison is nevertheless
decisive: the present solutions do not reproduce the observed billow growth.

The 2-D fit is therefore not yet successful, and extension to 3-D remains on
hold. The next nonlinear search must use both ridge and width/envelope
observables rather than optimizing the ridge alone.

## Cat's-eye formation and vortex pairing diagnostic

The earlier four-wavelength initial condition had exact 21.5 km translational
symmetry. A deterministic periodic solver preserves that symmetry, so the four
billows could grow but could not pair or merge. The solver now accepts a 43 km
subharmonic pairing seed and an optional 86 km domain-scale seed, each with an
independent phase.

A morphology pilot with center Ri = 0.090, a 2.5 m/s primary seed, and a
1.5 m/s pairing seed forms four primary cat's eyes and then produces two
simultaneous two-to-one vortex mergers. By 25 minutes the tracer and neutral
vorticity contain two domain-filling paired structures. The 8-rank run retained
tracer mass to 1.1e-15 relative error. The animation is
`output/figure3-cat-eye-double-merger/figure3_cat_eye_double_merger.gif`.

This establishes that the nonlinear solver can form cat's eyes and vortex
pairing when the symmetry-breaking mode is present. It is a controlled
morphology result, not a fit to the Figure 3 timing: the mergers remain slower
than the approximately 12-minute published sequence. The next fit must infer
the primary and subharmonic amplitudes together with shear width and velocity
jump, while retaining the held-out ridge and width/envelope tests.

## Required next steps

1. Add a robust full-envelope observable (for example weighted 5--95 percent
   vertical bounds) and fit shear width, velocity jump, seed amplitude, and
   phase against both ridge and layer growth while retaining the frozen
   held-out comparison.
2. Test whether one or two physically polarized gravity waves are useful only
   as initial/background forcing of that nonlinear solution.
3. Digitize the Figure 4 latitude cuts for a north-south held-out comparison
   after the Figure 3 nonlinear fit is predictive.
4. Obtain stable public download URLs from ScienceDB and start resumable
   transfers on Revontuli.
5. Verify byte counts and MD5 digests before opening any measurement file.
6. Inspect the Figure 2/3 density, time, altitude, and geographic arrays and
   the meteor-wind file.
7. Replace or refine the figure-derived fit with the released arrays and use
   held-out measurements before extending the model to three dimensions.

## Figure-derived constraints

- Figure 3 spans 06:41--06:53 local time in 15 zonal-altitude panels.
- Early layer displacement is approximately +/-150 m.
- Mature billows have approximately 2 km vertical extent and 22 km horizontal
  separation.
- The structures advect westward at approximately 100 m/s and evolve from
  weak modulation to broken loops over approximately 10 minutes.
- Figure 4 supplies five latitude cuts for a held-out test of north-south
  consistency; observed widths are approximately 1--3 km.
- These values are constraints reported by the paper and read from its
  figures, not estimates from the unavailable released arrays.

The experimental 3-D solver draft remains unvalidated and is deliberately on
hold until the first measured two-dimensional cut is matched.

## 2026-10-13 quick descending 2-D pilot

Completed an hour on Avaruus, 512 x 96 grid, 361 physical frames at 10 s
cadence. 128 ranks benchmarked 0.205 s versus 0.415 s for 8 ranks over the
same 10 simulated seconds, with matching vertical-wind diagnostics. The main
run began with 8 ranks and resumed at 1650 s on 128 ranks; the continuation
completed in 67.85 s. The full initial-to-final density-mass drift checked in
the transferred product is -1.61e-16.

Configuration: configs/descending-2d-14-7min.json. Along-flow seeds are 42 and
84 km (7 and 14 min passage periods at the imposed -100 m/s mean flow), with
2.5 and 1.5 m/s velocity amplitudes. Initial layer sigma is 150 m, peak density
2e12 m^-3, center Richardson number 0.090. These are hypotheses, not fitted
measurements. The conservative tracer is long-lived; magnetic ion drift,
polarization electric fields, and convergent ion transport remain absent.

IMPORTANT: layer_descent_m_s=-5000/3600 prescribes the altitude of a vertically
translating shear coordinate. The solver dynamics evolve in that frame. This
is kinematic layer phase descent, not a calculated neutral vertical velocity
or a self-consistent migrating ion-convergence zero crossing. It must not be
presented as validation of the physical descent mechanism.

Remote full states: /home/j/src/syisr-khi-rust/output/descending-2d-14-7min.
Local compact density states: output/descending-2d-14-7min/density_frames.h5,
created by src/bin/pack-density.rs; all physical frames are retained.
GIF and MP4: output/descending-2d-14-7min/descending_2d.{gif,mp4}.
Renderer: scripts/render_descending_2d.py; fixed altitude axes, fixed
1e10--2e12 m^-3 logarithmic color range, spatial interpolation for display,
no inserted temporal frames. Final frame inspected. Broad paired structures
remain; this is not yet a measured density-shape fit.

## 2026-10-13 magnetic-ion correction in progress

The next run is `descending-2d-magnetic-ri018`, launched on Avaruus as
`khi2d-magnetic-ri018.service`, 128 MPI ranks, 512 x 192 cells, one hour,
10-second frames. Center Ri is 0.18, with signed shear amplitude -49.97 m/s.
The magnetic mobility uses Fe-like 55.845 amu, B=35 microtesla, dip=20 degrees,
and prescribed collision frequency 300 s^-1 (these are scenario assumptions,
not verified event wind/magnetic profiles). The background electric field in
the moving layer frame is E=-U_ref cross B, with U_ref=-100 m/s: east=0,
north=1.197 mV/m, up=3.289 mV/m. This makes the ion convergence zero crossing
coincide with the shear center. The negative dU/dz gives upward ion drift
below and downward drift above that center. The force-balance unit test verifies
this sign, stationary vertical drift at the center, and disappearance of the
magnetic convergence in the zero-B, zero-E limit. All four Rust unit tests
passed on Avaruus. No ad hoc ion_convergence_s term is enabled in this run.

This is a prescribed-background E/B collisional-ion model, not a
self-consistent polarization/current-continuity calculation. Layer descent
is still kinematically prescribed. At the selected constant nu, the analytic
background trapping time is 116.5 s; diffusion-only equilibrium sigma is
24.1 m, which is below the 31.1 m vertical grid spacing. Thin-sheet width
therefore needs a resolution check before being claimed quantitatively.
The next comparison config Ri=0.24 exists but has not been launched yet.

Public YouTube upload of the EARLIER passive-tracer animation completed:
https://youtu.be/HhtzZnahuhc . Title: KHI first light sporadic E simulation.
Description: KHI first light sporadic E simulation with electron density as a tracer.
This video is not the magnetic-ion run.

Magnetic Ri=0.18 run completed: 361 frames, 3600 s, 439.55 wall seconds,
128 ranks; native float64 tracer mass drift -1.13e-15. The packed float32
visualization product has about 1.1e-8 relative density-sum drift from rounding.
Median column 5--95% thickness changes from 0.493 to 0.175 km; final p90 is
0.655 km. Final peak is 1.535e13 m^-3, too high versus the supplied measured
scale. Convergence restores a thin sheet but does not establish a calibrated
match. Thin-sheet convergence remains below grid resolution. The density
color scale saturates at 2e12 m^-3; use native/packed data to inspect true peaks.
Vorticity/density comparison is core_density_diagnostic.{png,pdf,h5}; automated
local-vorticity maxima are candidates only, NOT tracked merger identities or
proof of two-to-one vortex merging. Renderer provenance is retained in the
HDF5 diagnostic and in magnetic_2d.json. RI=0.24 comparison remains unrun.

## 2026-10-13 requested Ri sweep, two-hour simulations

Launched sequential sweep on Avaruus: center Ri = 0.10, 0.12, 0.14, 0.16, 0.18.
Each run starts afresh, ends at 7200 s, and saves 721 frames at 10 s cadence;
512 x 192 grid, 128 MPI ranks per case. All ion mobility, E/B, diffusion,
initial density, two spatial modes, and phase/descent assumptions match the
magnetic Ri=0.18 run. Only signed tanh shear amplitude varies, with fixed
N2=4e-4 s^-2 and width=1060 m. As a consequence magnetic convergence changes
with the shear, rather than being artificially held at an independent rate.
Linear prescribed phase descent continues to 108 km at two hours; collision
frequency, B, E and background remain fixed scenario values. This longer run
is a morphology sensitivity experiment, not an altitude-calibrated atmosphere.

Local configs and manifest: configs/ri-sweep-2h/. Runner:
scripts/run_ri_sweep_2h.sh. Remote service:
khi2d-ri-sweep-2h-20261013.service. Remote outputs:
/home/j/src/syisr-khi-rust/output/ri-sweep-2h-20261013/.
The batch writes atomic status.json, per-case run.log and diagnostics.json,
then density_frames.h5 with density/vorticity and per-case complete.json.
The next case starts automatically. Existing runs are never overwritten.
The solver binary SHA256 and requested config are recorded. No 3D run is
started and other users remain uncapped.

Ri sweep COMPLETE, 2026-10-13 at 08:14:18 UTC (10:14 Oslo). Every case
reached 7200 s with 721 checkpoints and native float64 mass drift below
1.7e-15. Case wall times in ascending Ri: 1034, 977, 924, 885, 844 s.
Local viewing exports and deliverables: output/ri-sweep-2h-20261013/animations/.
ri010, ri012, ri014, ri016, ri018 each have MP4 (H.264, 1100x660, 25 fps,
721 actual frames, 28.84 s), final PNG, JSON rendering provenance, and
matching three-beam RTI PNG/PDF/HDF5. All use 103--122 km fixed video axes
and 1e10--2e12 m^-3 color limits. View HDF5 packs use log10 density quantized
at 1/4096 dex to limit transfer size; full scientific outputs remain on
Avaruus. View packs are visualization products, not replacements for the
native full-precision simulation output. RTI HDF5 records quantization.
Peaks in sampled RTIs exceed the displayed ceiling (roughly 2.2e13--4.3e13
m^-3); color saturation must not be interpreted as observed peak agreement.
Scripts: src/bin/pack-density.rs --visualization, scripts/render_descending_2d.py,
scripts/plot_2d_range_time_density.py. Two endpoint RTIs and all final video
frames visually inspected; all MP4s confirmed 721 frames by ffprobe.

## 2026-10-13 restored short wavelengths and extended Ri range

New batch: output/ri-sweep-short-2h-20261013 on Avaruus, service
khi2d-ri-sweep-short-2h-20261013.service. Restored primary 21.5 km and
pairing 43 km wavelengths. Domain is 172 km, not 168 km, so eight primary
and four pairing wavelengths fit exactly; otherwise the periodic seam would
introduce a spurious discontinuity. Grid 512 x 192, 128 ranks, 7200 s,
10-second output; same seed amplitudes, magnetic Fe-like mobility, E field,
and prescribed descending convergence zero crossing (118 to 108 km).

User extended sweep to Ri=0.08. All cases: 0.08,0.10,0.12,0.14,0.16,0.18.
Existing 0.10--0.18 batch was already active; additional Ri=0.08 is queued
under khi2d-ri008-short-2h-20261013.service and starts automatically after
successful completion of that batch. It does not compete for CPU with the
active 128-rank case. Its config and all six entries are in
configs/ri-sweep-short-2h/manifest.json. Separate existing outputs preserved.
Scripts: run_ri_sweep_short_2h.sh and run_ri008_short_2h.sh. Verify both
service/status files before reporting the extended sweep complete. The
0.18 case is retained as a stable comparison. No claim of merging yet.

Short-wavelength sweep COMPLETE (six cases), verified live on Avaruus.
Ri=0.08 completed at 09:44:04 UTC, 11:44:04 Oslo on 2026-10-13. All six
runs reached 7200 s with 721 frames and mass drift below 1.8e-15.
Local deliverables: output/ri-sweep-short-2h-20261013/animations/riXXX.mp4
and riXXX_rti.{png,pdf,h5}, for XXX=008,010,012,014,016,018.
MP4s: H.264, 1100x660, 25 fps, 721 actual frames, 28.84 s; all probed.
Six final-frame panels and six RTIs visually inspected as contact sheets.
Native full precision data/vorticity remain on Avaruus; local viewing packs
are log10 density quantized to 1/4096 dex and are explicitly visualization
products. Plot source provenance is in each rendering JSON and RTI HDF5.
Identical 1e10--2e12 m^-3 display limits saturate larger simulated peaks.
Reduced density-blob count at low Ri is not yet proof of neutral vortex merger;
retain that distinction until neutral core trajectories are examined.

## 2026-10-13 localized Ri=0.08 packet and finite synthetic radar

Started khi2d-localized-ri008.service on Avaruus with 128 ranks; not a new
3-D run. configs/localized-ri008.json tapers the streamfunction seed: central
129 km plateau (six primary wavelengths) with 10.75 km raised-cosine edges.
No density blobs planted. Unseeded magnetic-ion background control complete,
721 frames. Remote output/localized-ri008 is still progressing.

Local preview output/localized-ri008/localized_preview.mp4 contains 345 actual
10 s checkpoints, t=0--3440 s; this is not the full two-hour result. Camera
follows prescribed phase descent and unwrapped packet center. No measured
fit or confirmed three-vortex merger is claimed. Fixed density display limits
1e10--2e12 m^-3 saturate larger peaks. Full precision scientific outputs remain
on Avaruus; local preview is log-density quantized 1/4096 dex, float32 velocity.

scripts/synthetic_radar_2d.py records 1280 RK2 ion-material marker trajectories,
wrap counts/events and unwrapped packet bounds in radar_preview.h5. Markers
are not vortex-core identities. The observation footprint selects central 80%
(10--90 percentiles) plus 2 km margin; full edge trajectories retained. The
original all-marker and 90% footprints stretched to >=domain width and were
rejected; the 80% core is a declared observation cutoff, not a PDE boundary fix.
Outside footprint: separate unseeded background control; periodic copies are
excluded in observation, but periodic dynamics still affect the fluid.
Preview has 2633 material wrap events and 25/24/25 beam scans seeing core,
zero ambiguous-width scans. Synthetic +/-20 degree zonal zenith beams and
vertical beam, curved Earth, 0.9 degree angular FWHM, 37.5 m range FWHM,
0.5 s integration, staggered 50 s scan. Not actual reconstructed SYISR beams.
Artifacts: radar_preview.{png,pdf,h5}; cursor animation radar_preview.mp4.
Renderer scripts: render_localized_packet.py and animate_radar_rti.py.
Checks: Rust five tests passed remotely; Python periodic bilinear sampling,
magnetic reference-flow balance/convergence sign, compilation; preview MP4
ffprobe confirms 345 frames, 1100x660, duration13.8 s. Final preview frame and
finite-radar RTI visually inspected.

## Rustplotlib renders of the original Ri=0.08 run

Requested original ri008.mp4 data retained; new renderer
plots/rustplotlib/src/main.rs reads all721 packed HDF5 checkpoints including
u,w and density. Frame uses x-U0t and relative altitude z-zc(t), U0=-100 m/s;
no quantile mask or synthetic background replacement. Electron density
1e10--2e12 m^-3 fixed logarithmic scale. Second render colors upward ion drift
from the same collisional Lorentz + neutral drag + electric field formula as
Rust solver, not neutral w. Fixed +/-40 m/s scale from full max35.2804 m/s.
Velocity is in the modeled convergence frame: prescribed descent phase speed
is not inserted as a physical ion velocity. Right panel U(z)=U0+DeltaU*tanh
(z/a), same vertical coordinate; no flow shear fitting implied. White
background, true10pt LaTeX at178mm wide. Rustplotlib API builds actual figures,
not Python drawing. Input full floating scientific outputs remain on Avaruus;
local velocityfloat32 and logdensity1/4096dex packed copy is explicitlyviewing.
Build HDF5_DIR=/opt/anaconda3; maintained hdf5-metno reader supports1.14.6.
Reproduction script plots/rustplotlib/render.sh. Outputs under
output/ri-sweep-short-2h-20261013/animations/rustplotlib/.
Rustplotlib sequences completed:721frames each. Midpoint PDF figures visually
inspected: density and ion-upward velocity/background-wind panel, correct
LaTeX labels and stable page sizing. Full-run native maximum ion vertical drift
35.280421 m/s. MP4 encoding H26430fps,721frames,24.0333s. These are visualization
outputs, not an observational fit. Rustplotlib default tight-export/API crop
fixes and scientific HEIMDALL/Dimant figures are in separate rustplotlib repo.
