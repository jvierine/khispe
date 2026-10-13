#!/usr/bin/env python3
"""Render a minimal 3-D KHI vorticity/tracer isosurface view from HDF5."""
from __future__ import annotations
import argparse
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from skimage.measure import marching_cubes


def surface(ax, field, level, spacing, origin, color, alpha):
    verts, faces, _, _ = marching_cubes(field, level=level, spacing=spacing, allow_degenerate=False)
    verts = verts / 1000.0 + np.asarray(origin)[None, :]
    mesh = Poly3DCollection(verts[faces], linewidths=0, alpha=alpha)
    mesh.set_facecolor(color)
    ax.add_collection3d(mesh)


def main():
    p=argparse.ArgumentParser(); p.add_argument('frame',type=Path); p.add_argument('output',type=Path); a=p.parse_args()
    with h5py.File(a.frame,'r') as h:
        x=h['x_m'][:]; y=h['y_m'][:]; z=h['altitude_m'][:]
        ox=h['omega_x_s-1'][:]; oy=h['omega_y_s-1'][:]; oz=h['omega_z_s-1'][:]
        n=h['ion_density_m-3'][:]; t=float(h.attrs['time_s'])
    # Halve x only; preserve the meridional and vertical structure.
    x=x[::2]; ox=ox[::2]; oy=oy[::2]; oz=oz[::2]; n=n[::2]
    om=np.sqrt(ox*ox+oy*oy+oz*oz)
    dx=np.median(np.diff(x)); dy=np.median(np.diff(y)); dz=np.median(np.diff(z)); spacing=(dx,dy,dz)
    origin=(x[0]/1000-43., y[0]/1000, z[0]/1000-116.)
    fig=plt.figure(figsize=(12.8,7.2),dpi=150,facecolor='#070b0d'); ax=fig.add_axes([0,0,1,1],projection='3d',facecolor='#070b0d')
    # Coral material sheet; two pale/coral vorticity levels reveal the rolled surface.
    surface(ax,n,.42*np.nanmax(n),spacing,origin,'#ff6b4a',.18)
    surface(ax,om,.36*np.nanmax(om),spacing,origin,'#ed725e',.34)
    surface(ax,om,.67*np.nanmax(om),spacing,origin,'#fff0dc',.74)
    ax.view_init(elev=25,azim=-58); ax.set_xlim(-43,43); ax.set_ylim(y[0]/1000,y[-1]/1000); ax.set_zlim((z[0]-116e3)/1000,(z[-1]-116e3)/1000)
    ax.set_box_aspect((86,32,18)); ax.set_axis_off()
    fig.text(.035,.93,'3-D KHI',color='#fff5e7',fontsize=22,weight='bold')
    fig.text(.035,.885,'vorticity surfaces + ion-layer sheet',color='#d39a88',fontsize=11)
    fig.text(.965,.93,f'{t/60:.1f} min',ha='right',color='#ffd7a3',fontsize=16,weight='bold')
    fig.text(.965,.055,'3× vertical',ha='right',color='#8e7169',fontsize=9)
    a.output.parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.output,facecolor=fig.get_facecolor(),dpi=150); plt.close(fig); print(a.output)
if __name__=='__main__': main()
