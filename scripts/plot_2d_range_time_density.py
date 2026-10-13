#!/usr/bin/env python3
"""Sample actual 2-D checkpoints at fixed vertical beams; preserve RTI in HDF5."""
import argparse
import json
from contextlib import nullcontext
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import PchipInterpolator


def main():
    p = argparse.ArgumentParser()
    p.add_argument('frames', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--x-km', nargs='+', type=float, default=[10, 43, 70])
    a = p.parse_args()
    packed = h5py.File(a.frames, 'r') if a.frames.is_file() else None
    paths = sorted(packed.keys()) if packed is not None else sorted(a.frames.glob('frame_*.h5'))
    def frame(path):
        return nullcontext(packed[path]) if packed is not None else h5py.File(path)
    if len(paths) < 2:
        raise ValueError('Need at least two physical checkpoints')
    with frame(paths[0]) as h:
        x = h['x_m'][:] / 1000
        z = h['altitude_m'][:] / 1000
        cfg = json.loads(bytes(h['config_json'][:]).decode())
    indices = [int(np.argmin(abs(x - v))) for v in a.x_km]
    times, samples, altitude_grids = [], [], []
    for path in paths:
        with frame(path) as h:
            if not np.array_equal(h['x_m'][:] / 1000, x):
                raise ValueError('Horizontal grid changed')
            altitude_grids.append(h['altitude_m'][:] / 1000)
            times.append(float(h.attrs['time_s']))
            if 'log10_ne_quantized' in h:
                samples.append(np.stack([10**(h['log10_ne_quantized'][i, :].astype(float)/4096) for i in indices]))
            else:
                samples.append(np.stack([h['ion_density_m-3'][i, :] for i in indices]))
    t = np.asarray(times)
    if np.any(np.diff(t) <= 0):
        raise ValueError('Checkpoint times must increase')
    original_ne = np.stack(samples, axis=1)  # beam,time,altitude
    altitude_grids = np.asarray(altitude_grids)
    z = np.linspace(altitude_grids.min(), altitude_grids.max(), 768)
    ne = np.stack([[np.interp(z, altitude_grids[j], original_ne[k,j], left=0, right=0)
                    for j in range(len(t))] for k in range(len(indices))])
    vmax = cfg['peak_m3']
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(a.output.with_suffix('.h5'), 'w') as h:
        h.attrs['description'] = 'Synthetic fixed vertical beams from 2-D passive ion tracer'
        h.attrs['observational_fit'] = False
        h.attrs['script'] = 'scripts/plot_2d_range_time_density.py'
        h.attrs['source_directory'] = str(a.frames.resolve())
        h.attrs['temporal_interpolation'] = 'none; physical checkpoints only'
        h.attrs['figure_display_interpolation'] = 'shape-preserving cubic interpolation of log10 density in time and altitude; display only'
        h.attrs['time_origin'] = 'simulation start; no local-time registration'
        h.attrs['visualization_log10_quantum'] = 1/4096 if packed is not None and 'log10_ne_quantized' in packed[paths[0]] else 0.
        h.create_dataset('time_s', data=t)
        h.create_dataset('altitude_km', data=z)
        h.create_dataset('beam_x_km', data=x[indices])
        h.create_dataset('ne_m3', data=ne, compression='gzip', shuffle=True)
        h.create_dataset('native_altitude_km', data=altitude_grids)
        h.create_dataset('native_ne_m3', data=original_ne, compression='gzip', shuffle=True)
        h.attrs['altitude_remapping'] = 'linear density interpolation to fixed altitude; zero outside moving domain'
        h.create_dataset('source_files', data=np.array([str(v) for v in paths], dtype=h5py.string_dtype()))
        h.create_dataset('initial_config_json', data=json.dumps(cfg).encode())
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(len(indices), 1, figsize=(11, 6.2), sharex=True,
                             sharey=True, layout='constrained', squeeze=False)
    ri = cfg['n2_s2'] / (cfg['wind_m_s']/cfg['shear_m'])**2
    fig.suptitle(f'Ri = {ri:.2f} | primary wavelength {cfg["seed_wavelength_m"]/1000:g} km', fontsize=12)
    display_t = np.linspace(t[0], t[-1], 1601)
    display_z = np.linspace(z[0], z[-1], 512)
    for k, ax in enumerate(axes[:, 0]):
        log_n = np.log10(np.maximum(ne[k], 1))
        smooth_t = PchipInterpolator(t, log_n, axis=0)(display_t)
        smooth = PchipInterpolator(z, smooth_t, axis=1)(display_z)
        im = ax.imshow(10**smooth.T, origin='lower', aspect='auto',
                       extent=(t[0]/60, t[-1]/60, z[0], z[-1]),
                       cmap='magma', norm=LogNorm(1e10, vmax), interpolation='bilinear')
        ax.set_ylabel('Altitude [km]')
        ax.set_title(f'x = {x[indices[k]]:.1f} km', loc='left', fontsize=11)
        ax.set_ylim(z.min(), z.max())
        ax.set_xlim(t[0]/60, t[-1]/60)
    axes[-1, 0].set_xlabel('Elapsed simulation time [min]')
    ticks = [1e10, 1e11, vmax] if vmax > 1e11 else [1e10, 10**10.5, 1e11]
    cb = fig.colorbar(im, ax=axes[:, 0].tolist(), pad=0.025, fraction=0.035, ticks=ticks)
    cb.set_ticklabels([r'$10^{10}$', r'$10^{11}$', f'{vmax:.0e}'] if vmax > 1e11
                     else [r'$10^{10}$', r'$10^{10.5}$', r'$10^{11}$'])
    cb.minorticks_off()
    cb.set_label(r'Synthetic $n_e$ [m$^{-3}$]')
    fig.savefig(a.output.with_suffix('.png'), dpi=180, facecolor='white')
    fig.savefig(a.output.with_suffix('.pdf'), facecolor='white')
    plt.close(fig)
    print(json.dumps({'checkpoints': len(paths), 'time_range_s': [t[0], t[-1]],
                      'beam_x_km': x[indices].tolist(), 'peak_m3': float(ne.max()),
                      'output': str(a.output)}, indent=2))


if __name__ == '__main__':
    main()
