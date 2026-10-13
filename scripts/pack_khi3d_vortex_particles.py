#!/usr/bin/env python3
"""Extract the strongest perturbation-vorticity parcels for GPU rendering."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import h5py
import numpy as np

p=argparse.ArgumentParser(); p.add_argument('input',type=Path); p.add_argument('stem',type=Path); p.add_argument('--percentile',type=float,default=98.0); a=p.parse_args()
with h5py.File(a.input,'r') as f:
    om=[]
    for name in ('omega_x_s-1','omega_y_s-1','omega_z_s-1'):
        q=f[name][:]; om.append(q-q.mean(axis=(0,1),keepdims=True))
    mag=np.sqrt(sum(q*q for q in om)); cut=float(np.percentile(mag,a.percentile)); scale=float(np.percentile(mag,99.8))
    ii,jj,kk=np.nonzero(mag>=cut)
    nx,ny,nz=mag.shape
    x=ii/(nx-1)-.5; y=(jj/(ny-1)-.5)*float(np.ptp(f['y_m'][:])/np.ptp(f['x_m'][:])); z=(kk/(nz-1)-.5)*float(4*np.ptp(f['altitude_m'][:])/np.ptp(f['x_m'][:]))
    wy=om[1][ii,jj,kk]/max(scale,1e-30)
    c1=f['clebsch_alpha1'][:][ii,jj,kk]; c1/=max(np.max(np.abs(c1)),1e-30)
    c2=f['clebsch_beta1_m2_s-1'][:][ii,jj,kk]; c2/=max(np.max(np.abs(c2)),1e-30)
    strength=np.clip(mag[ii,jj,kk]/max(scale,1e-30),0,1.5)
    pts=np.column_stack((x,y,z,wy,c1,c2,strength)).astype('<f4')
    t=float(f.attrs['time_s'])
a.stem.parent.mkdir(parents=True,exist_ok=True); raw=a.stem.with_suffix('.particles.f32'); pts.tofile(raw)
meta={'particles':raw.name,'count':len(pts),'stride_floats':7,'time_s':t,'threshold_percentile':a.percentile,'threshold_s-1':cut,'scale_s-1':scale,'source_hdf5':str(a.input)}
a.stem.with_suffix('.particles.json').write_text(json.dumps(meta,indent=2)+'\n'); print(len(pts))
