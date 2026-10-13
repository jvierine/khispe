#!/usr/bin/env python3
"""Extract a smooth Q-criterion vortex surface for the Three.js renderer."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import h5py,numpy as np
from scipy.ndimage import map_coordinates
from skimage.measure import marching_cubes

p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('stem',type=Path);p.add_argument('--positive-percentile',type=float,default=90.0);a=p.parse_args()
with h5py.File(a.input,'r') as f:
 u,v,w=[f[k][:] for k in ('u_m_s','v_m_s','w_m_s')]; x=f['x_m'][:];y=f['y_m'][:];z=f['altitude_m'][:];wy=f['omega_y_s-1'][:];cl=f['clebsch_alpha1'][:];ne=f['ion_density_m-3'][:];t=float(f.attrs['time_s'])
dx,dy,dz=np.mean(np.diff(x)),np.mean(np.diff(y)),np.mean(np.diff(z))
vel=(u,v,w); grad=[[np.gradient(vel[i],(dx,dy,dz)[j],axis=j,edge_order=2) for j in range(3)] for i in range(3)]
ss=np.zeros_like(u);oo=np.zeros_like(u)
for i in range(3):
 for j in range(3):
  sij=.5*(grad[i][j]+grad[j][i]);oij=.5*(grad[i][j]-grad[j][i]);ss+=sij*sij;oo+=oij*oij
Q=.5*(oo-ss);positive=Q[Q>0];level=float(np.percentile(positive,a.positive_percentile))
verts,faces,normals,_=marching_cubes(Q,level=level,spacing=(dx,dy,dz),allow_degenerate=False,step_size=1)
# Center, normalize by the zonal domain, and exaggerate altitude fourfold.
lx=float(np.ptp(x)); verts[:,0]=(verts[:,0]+x[0]-.5*(x[0]+x[-1]))/lx;verts[:,1]=(verts[:,1]+y[0]-.5*(y[0]+y[-1]))/lx;verts[:,2]=4*(verts[:,2]+z[0]-.5*(z[0]+z[-1]))/lx
idx=np.column_stack(((verts[:,0]*lx+.5*(x[0]+x[-1])-x[0])/dx,(verts[:,1]*lx+.5*(y[0]+y[-1])-y[0])/dy,(verts[:,2]*lx/4+.5*(z[0]+z[-1])-z[0])/dz)).T
swy=map_coordinates(wy,idx,order=1,mode='nearest');scl=map_coordinates(cl,idx,order=1,mode='nearest');ws=max(float(np.percentile(np.abs(swy),99)),1e-30)
# A restrained coral/ivory range: vorticity sets tone and Clebsch labels add material variation.
r=np.clip(.72+.23*np.tanh(swy/ws)+.05*scl,.30,1);g=np.clip(.24+.25*(1-np.abs(swy/ws))+.06*scl,.12,.70);b=np.clip(.17+.22*(1-np.abs(swy/ws))+.04*scl,.10,.62);colors=np.column_stack((r,g,b)).astype('<f4')
stem=a.stem;stem.parent.mkdir(parents=True,exist_ok=True)
paths={'positions':stem.with_suffix('.q.pos.f32'),'colors':stem.with_suffix('.q.color.f32'),'indices':stem.with_suffix('.q.index.u32')}
verts.astype('<f4').tofile(paths['positions']);colors.tofile(paths['colors']);faces.astype('<u4').ravel().tofile(paths['indices'])
iv,ifa,_,_=marching_cubes(ne,level=.18*float(ne.max()),spacing=(dx,dy,dz),allow_degenerate=False,step_size=1)
iv[:,0]=(iv[:,0]+x[0]-.5*(x[0]+x[-1]))/lx;iv[:,1]=(iv[:,1]+y[0]-.5*(y[0]+y[-1]))/lx;iv[:,2]=4*(iv[:,2]+z[0]-.5*(z[0]+z[-1]))/lx
ionpaths={'ion_positions':stem.with_suffix('.ion.pos.f32'),'ion_indices':stem.with_suffix('.ion.index.u32')};iv.astype('<f4').tofile(ionpaths['ion_positions']);ifa.astype('<u4').ravel().tofile(ionpaths['ion_indices'])
meta={'time_s':t,'vertices':len(verts),'triangles':len(faces),'q_level_s-2':level,'positive_percentile':a.positive_percentile,'positions':paths['positions'].name,'colors':paths['colors'].name,'indices':paths['indices'].name,'ion_positions':ionpaths['ion_positions'].name,'ion_indices':ionpaths['ion_indices'].name,'ion_vertices':len(iv),'ion_triangles':len(ifa),'source_hdf5':str(a.input),'vertical_exaggeration':4.0}
stem.with_suffix('.q.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps({'vertices':len(verts),'triangles':len(faces),'q_level':level}))
