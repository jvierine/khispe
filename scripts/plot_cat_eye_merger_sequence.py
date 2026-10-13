#!/usr/bin/env python3
"""Plot selected states from the controlled KHI double-merger pilot."""
from __future__ import annotations
import argparse
from pathlib import Path
import h5py, matplotlib.pyplot as plt, numpy as np
BG='#ffffff'; PANEL='#ffffff'; FG='#16242c'; MUTED='#51656f'; GOLD='#b56a00'; GRID='#cbd6dc'
def main():
 p=argparse.ArgumentParser(); p.add_argument('frames',type=Path); p.add_argument('output',type=Path); a=p.parse_args()
 paths=sorted(a.frames.glob('frame_*.h5')); targets=np.array([0,300,600,900,1200,1500.])
 alltimes=[]
 for q in paths:
  with h5py.File(q) as h: alltimes.append(float(h.attrs['time_s']))
 selected=[paths[int(np.argmin(abs(np.asarray(alltimes)-t)))] for t in targets]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.facecolor':BG,'axes.facecolor':PANEL,'text.color':FG,'axes.labelcolor':FG,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':GRID})
 fig,axes=plt.subplots(2,3,figsize=(12.5,6.4)); fig.subplots_adjust(left=.07,right=.98,top=.95,bottom=.10,wspace=.09,hspace=.23)
 for i,(ax,path) in enumerate(zip(axes.flat,selected)):
  with h5py.File(path) as h:
   x=h['x_m'][:]/1000-43; z=h['altitude_m'][:]/1000; n=h['ion_density_m-3'][:]; om=h['omega_s-1'][:]; t=float(h.attrs['time_s'])
  logn=np.log10(np.maximum(n,1)); ax.pcolormesh(x,z,logn.T,cmap='magma',vmin=10,vmax=11,shading='auto',rasterized=True)
  vmax=max(np.percentile(abs(om),99.5),1e-5); ax.contour(x,z,om.T,levels=[-.5*vmax,.5*vmax],colors=['#007f8b','#e6ad00'],linewidths=.8,alpha=.9)
  ax.set(xlim=(x[0],x[-1]),ylim=(z[0],z[-1]),title=f'{t/60:.0f} min')
  if i//3==1: ax.set_xlabel('Zonal distance [km]')
  else: ax.set_xticklabels([])
  if i%3==0: ax.set_ylabel('Altitude [km]')
  else: ax.set_yticklabels([])
 a.output.parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.output,dpi=190,facecolor=BG); fig.savefig(a.output.with_suffix('.pdf'),facecolor=BG); print(a.output)
if __name__=='__main__': main()
