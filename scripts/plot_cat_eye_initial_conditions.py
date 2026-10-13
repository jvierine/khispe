#!/usr/bin/env python3
"""Write and plot the exact spatial seed used by the double-merger KHI pilot."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

BG="#ffffff"; PANEL="#ffffff"; FG="#16242c"; MUTED="#51656f"; CYAN="#087f8c"; GOLD="#b56a00"; CORAL="#c53f39"; GRID="#cbd6dc"
MERGER = LinearSegmentedColormap.from_list(
 "merger", ["#0823a8", "#0878f4", "#32d7e8", "#f5ed63", "#ff9a3c", "#cc3a24", "#541313"]
)

def main():
 p=argparse.ArgumentParser(); p.add_argument('config',type=Path); p.add_argument('h5',type=Path); p.add_argument('figure',type=Path); a=p.parse_args()
 c=json.loads(a.config.read_text()); nx,nz=c['nx'],c['nz']; lx,lz=c['lx_m'],c['lz_m']; dz=lz/(nz+1)
 x=(np.arange(nx)+0.5)*lx/nx; z=(np.arange(nz)+1)*dz-lz/2
 k1=2*np.pi/c['seed_wavelength_m']; k2=k1/2
 wall=np.sin(np.pi*(np.arange(nz)+1)/(nz+1)); env=np.exp(-0.5*(z/(2*c['shear_m']))**2)*wall
 psi1=(c['seed_m_s']/k1)*np.sin(k1*x)[:,None]*env[None,:]
 psi2=(c['seed_pairing_m_s']/k2)*np.sin(k2*x+c['seed_pairing_phase_rad'])[:,None]*env[None,:]
 kd=2*np.pi/lx
 psid=(c.get('seed_domain_m_s',0.)/kd)*np.sin(kd*x+c.get('seed_domain_phase_rad',0.))[:,None]*env[None,:]
 w1=-c['seed_m_s']*np.cos(k1*x)[:,None]*env[None,:]
 w2=-c['seed_pairing_m_s']*np.cos(k2*x+c['seed_pairing_phase_rad'])[:,None]*env[None,:]
 wd=-c.get('seed_domain_m_s',0.)*np.cos(kd*x+c.get('seed_domain_phase_rad',0.))[:,None]*env[None,:]
 u1=np.gradient(psi1,z,axis=1); u2=np.gradient(psi2,z,axis=1); ud=np.gradient(psid,z,axis=1)
 wt=w1+w2+wd; u=c['advection_m_s']+c['wind_m_s']*np.tanh(z/c['shear_m'])
 n=c['peak_m3']*np.exp(-0.5*(z/c['layer_sigma_m'])**2)
 a.h5.parent.mkdir(parents=True,exist_ok=True)
 with h5py.File(a.h5,'w') as h:
  h.attrs['status']='initial conditions for controlled morphology pilot, not an observational fit'; h.attrs['config_path']=str(a.config); h.attrs['generator']='scripts/plot_cat_eye_initial_conditions.py'
  for key,val in c.items():
   if isinstance(val,(int,float,str,bool)): h.attrs[key]=val
  h['x_m']=x; h['z_m']=z; h['base_wind_m_s']=u; h['ion_density_m-3']=n; h['primary_u_m_s']=u1; h['primary_w_m_s']=w1; h['pairing_u_m_s']=u2; h['pairing_w_m_s']=w2; h['domain_u_m_s']=ud; h['domain_w_m_s']=wd; h['total_u_m_s']=u1+u2+ud; h['total_w_m_s']=wt
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.facecolor':BG,'axes.facecolor':PANEL,'text.color':FG,'axes.labelcolor':FG,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':GRID})
 fig=plt.figure(figsize=(12.5,6.8)); gs=fig.add_gridspec(2,3,height_ratios=[0.72,1],left=.07,right=.97,top=.95,bottom=.14,wspace=.18,hspace=.34)
 ax=fig.add_subplot(gs[0,:2]); iz=np.argmin(abs(z)); xx=x/1000-lx/2000
 ax.plot(xx,w1[:,iz],color=CYAN,lw=2,label=r'primary: $\lambda_1=21.5$ km')
 ax.plot(xx,w2[:,iz],color=GOLD,lw=2,label=r'pairing: $\lambda_2=43$ km')
 ax.plot(xx,wt[:,iz],color=FG,lw=2.5,label='sum')
 ax.axhline(0,color=GRID,lw=.8); ax.grid(color=GRID,alpha=.35); ax.set(xlabel='Zonal distance [km]',ylabel="Initial $w'$ [m s$^{-1}$]",title='Horizontal seed at the shear center'); ax.legend(frameon=False,labelcolor=FG,ncols=3,loc='upper right')
 ax=fig.add_subplot(gs[0,2]); ax.plot(u,z/1000,color=CORAL,lw=2.5); ax.axhline(0,color=GRID,lw=.8); ax.grid(color=GRID,alpha=.35); ax.set(xlabel='Base wind [m s$^{-1}$]',ylabel='Height from center [km]',title='Tanh shear + flow direction')
 iq=np.linspace(10,nz-11,5,dtype=int); ax.quiver(np.zeros_like(iq,dtype=float),z[iq]/1000,u[iq],np.zeros_like(iq),angles='xy',scale_units='xy',scale=1,color=CYAN,width=.009,headwidth=4.2,headlength=5.2,alpha=.9)
 titles=['Primary mode: four equal billows','Pairing subharmonic: two pairs','Combined initial vertical velocity']
 fields=[w1,w2,wt]; ufields=[u1,u2,u1+u2+ud]
 vmax=max(np.max(abs(q)) for q in fields)
 for j,(title,q,uq) in enumerate(zip(titles,fields,ufields)):
  ax=fig.add_subplot(gs[1,j]); im=ax.pcolormesh(xx,z/1000,q.T,cmap=MERGER,vmin=-vmax,vmax=vmax,shading='auto',rasterized=True); ax.contour(xx,z/1000,q.T,levels=[0],colors=[FG],linewidths=.7,alpha=.65)
  sx=slice(24,None,64); sz=slice(20,None,42); qs=max(np.nanpercentile(np.hypot(uq,w1 if j==0 else w2 if j==1 else wt),98),.1)
  ax.quiver(xx[sx],z[sz]/1000,(uq[sx,sz]/qs).T,(q[sx,sz]/qs).T,color=FG,alpha=.72,pivot='mid',scale=18,width=.0042,headwidth=3.4,headlength=4.4)
  ax.set(xlabel='Zonal distance [km]',title=title); 
  if j==0: ax.set_ylabel('Height from shear center [km]')
  else: ax.set_yticklabels([])
 cb=fig.colorbar(im,ax=fig.axes[2:5],orientation='horizontal',fraction=.035,pad=.16,aspect=55); cb.set_label("Initial vertical velocity [m s$^{-1}$]; arrows show perturbation velocity")
 a.figure.parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.figure,dpi=190,facecolor=BG); fig.savefig(a.figure.with_suffix('.pdf'),facecolor=BG); print(a.figure)
if __name__=='__main__': main()
