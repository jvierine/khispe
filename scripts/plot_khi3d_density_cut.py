#!/usr/bin/env python3
"""Plot the meridional-midplane electron-density cut from a 3-D KHI frame."""
from __future__ import annotations
import argparse
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
import numpy as np

p=argparse.ArgumentParser(); p.add_argument('input',type=Path); p.add_argument('output',type=Path); p.add_argument('--lat',type=float,default=78.15); p.add_argument('--lon0',type=float,default=16.02); a=p.parse_args()
with h5py.File(a.input,'r') as f:
    ne=f['ion_density_m-3'][:]; x=f['x_m'][:]; z=f['altitude_m'][:]; t=float(f.attrs['time_s'])
j=ne.shape[1]//2; cut=ne[:,j,:].T
x0=.5*(x[0]+x[-1]); lon=a.lon0+(x-x0)/(111_320*np.cos(np.deg2rad(a.lat)))
alt=z/1e3
fig,ax=plt.subplots(figsize=(10.8,4.8),layout='constrained')
levels=np.linspace(9.0,11.0,101); logne=np.log10(np.maximum(cut,1e9))
im=ax.contourf(lon,alt,logne,levels=levels,cmap='turbo',extend='both')
cs=ax.contour(lon,alt,cut/1e11,levels=[.2,.4,.6,.8],colors='white',linewidths=.65,alpha=.8)
ax.clabel(cs,fmt=lambda v:f'{v:.1f}',fontsize=7,inline=True)
ax.set(xlabel='Longitude [°E]',ylabel='Altitude [km]',title=f'Electron density at the meridional midplane, t = {t/60:.1f} min')
ax.tick_params(direction='out'); ax.spines[['top','right']].set_visible(False)
cb=fig.colorbar(im,ax=ax,pad=.02); cb.set_label(r'$\log_{10} n_e$ [m$^{-3}$]')
a.output.parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.output,dpi=220,facecolor='white'); fig.savefig(a.output.with_suffix('.pdf'),facecolor='white'); print(a.output)
