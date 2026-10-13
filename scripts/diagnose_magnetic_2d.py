#!/usr/bin/env python3
"""Inspect density thickness and candidate neutral-vorticity cores together."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from scipy.ndimage import gaussian_filter,maximum_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
with h5py.File(a.input) as h:
 keys=sorted(h.keys()); times=np.array([h[k].attrs['time_s'] for k in keys])
 cfg=json.loads(bytes(h[keys[0]]['config_json'][:]).decode())
 x=h[keys[0]]['x_m'][:]/1000
 summary=[]; snapshots=[]
 for k in keys:
  g=h[k];t=float(g.attrs['time_s']);n=g['ion_density_m-3'][:].astype(float);z=g['altitude_m'][:]/1000
  cum=np.cumsum(n,axis=1);cum/=cum[:,-1,None]
  low=np.array([np.interp(.05,row,z) for row in cum]);high=np.array([np.interp(.95,row,z) for row in cum])
  center=cfg['altitude_m']/1000+cfg['layer_descent_m_s']*t/1000
  rel=(z-center)*1000
  total=g['omega_s-1'][:].astype(float)-cfg['wind_m_s']/cfg['shear_m']/np.cosh(rel/cfg['shear_m'])**2
  smooth=gaussian_filter(total,(1.5,1.5),mode=('wrap','nearest'))
  maxima=(smooth==maximum_filter(smooth,size=(33,17),mode=('wrap','nearest')))&(smooth>.5*smooth.max())
  candidates=np.argwhere(maxima)
  summary.append((t,np.median(high-low),np.percentile(high-low,90),n.max(),n.sum(),len(candidates)))
  if any(abs(t-v)<1e-6 for v in [600,1500,2400,3600]):
   snapshots.append((t,n,z,total,candidates))
summary=np.asarray(summary)
with h5py.File(a.output.with_suffix('.h5'),'w') as h:
 h.attrs['script']='scripts/diagnose_magnetic_2d.py'
 h.attrs['core_definition']='candidate maxima of smoothed total neutral vorticity above 50% frame maximum; not a merger proof'
 h.create_dataset('time_s',data=summary[:,0]);h.create_dataset('median_5_95_width_km',data=summary[:,1]);h.create_dataset('p90_5_95_width_km',data=summary[:,2]);h.create_dataset('peak_density_m3',data=summary[:,3]);h.create_dataset('density_sum',data=summary[:,4]);h.create_dataset('candidate_core_count',data=summary[:,5])
fig,axs=plt.subplots(2,len(snapshots),figsize=(14,5.8),layout='constrained',sharex=True)
for j,(t,n,z,total,candidates) in enumerate(snapshots):
 center=cfg['altitude_m']/1000+cfg['layer_descent_m_s']*t/1000
 zr=z-center
 v=axs[0,j].imshow(total.T,origin='lower',extent=(x[0],x[-1],zr[0],zr[-1]),aspect='auto',cmap='viridis',vmin=0,vmax=.06)
 q=axs[1,j].imshow(np.maximum(n.T,1),origin='lower',extent=(x[0],x[-1],zr[0],zr[-1]),aspect='auto',cmap='magma',norm=LogNorm(1e10,2e12))
 axs[0,j].set_title(f'{t/60:g} min');axs[1,j].set_xlabel('Zonal distance [km]')
 for ax in axs[:,j]:ax.set_ylim(-3,3)
axs[0,0].set_ylabel('Height from layer center [km]');axs[1,0].set_ylabel('Height from layer center [km]')
fig.colorbar(v,ax=axs[0,:].tolist(),label=r'Neutral vorticity [s$^{-1}$]',shrink=.85)
cb=fig.colorbar(q,ax=axs[1,:].tolist(),label=r'$n_e$ [m$^{-3}$]',shrink=.85);cb.minorticks_off()
fig.savefig(a.output.with_suffix('.png'),dpi=150);fig.savefig(a.output.with_suffix('.pdf'));plt.close(fig)
print(json.dumps({'initial_median_width_km':summary[0,1],'final_median_width_km':summary[-1,1],'final_p90_width_km':summary[-1,2],'final_peak_m3':summary[-1,3],'relative_mass_drift':summary[-1,4]/summary[0,4]-1,'candidate_core_counts':summary[::60,5].tolist()},indent=2))
