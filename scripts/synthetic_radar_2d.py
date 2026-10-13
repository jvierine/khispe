#!/usr/bin/env python3
"""Finite synthetic radar of one material packet image, with unwrapped tracking.

The periodic PDE is not made nonperiodic by this observation operator. We
select a single unwrapped packet footprint, use the unseeded background control
outside it, and reject ambiguous footprints spanning a full periodic domain.
Beam angles, scan offsets, and Gaussian responses are scenario assumptions.
"""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def bilinear(field,x,z,lx,z0,dz):
    nx,nz=field.shape
    xi=np.asarray(x)/lx*nx-.5
    zi=(np.asarray(z)-z0)/dz
    i=np.floor(xi).astype(int);j=np.floor(zi).astype(int)
    fx=xi-i;fz=zi-j
    j=np.clip(j,0,nz-2);fz=np.clip(fz,0,1)
    return ((1-fx)*(1-fz)*field[i%nx,j]+fx*(1-fz)*field[(i+1)%nx,j]
            +(1-fx)*fz*field[i%nx,j+1]+fx*fz*field[(i+1)%nx,j+1])

def ion_velocity(u,w,c):
    mu=1.602176634e-19/(c['ion_mass_amu']*1.66053906660e-27*c['nu_in_s'])
    beta=mu*c['b_t'];by=np.cos(np.deg2rad(c['dip_deg']));bz=-np.sin(np.deg2rad(c['dip_deg']))
    fx=u+mu*c['electric_east_v_m'];fy=mu*c['electric_north_v_m'];fz=w+mu*c['electric_up_v_m']
    return ((fx+beta*(fy*bz-fz*by))/(1+beta**2),
            (fz+beta*fx*by+beta**2*bz*(fy*by+fz*bz))/(1+beta**2))

def density(g):
    if 'log10_ne_quantized' in g:return 10**(g['log10_ne_quantized'][:].astype(float)/4096)
    return g['ion_density_m-3'][:].astype(float)

def main():
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('background',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--station-x-km',type=float,default=-150.)
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(a.input) as h,h5py.File(a.background) as bg:
        keys=sorted(h.keys());times=np.array([h[k].attrs['time_s'] for k in keys])
        c=json.loads(bytes(h[keys[0]]['config_json'][:]).decode());lx=c['lx_m']
        zgrid=h[keys[0]]['altitude_m'][:]-c['altitude_m'];z0=zgrid[0];dz=zgrid[1]-zgrid[0]
        # These label ion material within the seeded packet, not vortex cores.
        xx=np.linspace(c['seed_packet_center_m']-c['seed_packet_half_width_m']-c['seed_packet_taper_m'],
                       c['seed_packet_center_m']+c['seed_packet_half_width_m']+c['seed_packet_taper_m'],256)
        px,pz=np.meshgrid(xx,[-300.,-150.,0.,150.,300.]);px=px.ravel();pz=pz.ravel()
        weights=np.exp(-.5*(pz/c['layer_sigma_m'])**2)
        track_x=[px.copy()];track_z=[pz.copy()];events=[]
        u0=h[keys[0]]['u_m_s'][:];w0=h[keys[0]]['w_m_s'][:]
        for k in range(1,len(keys)):
            u1=h[keys[k]]['u_m_s'][:];w1=h[keys[k]]['w_m_s'][:]
            dt=times[k]-times[k-1];steps=int(np.ceil(dt));dd=dt/steps
            for j in range(steps):
                f=j/steps;fm=(j+.5)/steps
                uu=(1-f)*bilinear(u0,px,pz,lx,z0,dz)+f*bilinear(u1,px,pz,lx,z0,dz)
                ww=(1-f)*bilinear(w0,px,pz,lx,z0,dz)+f*bilinear(w1,px,pz,lx,z0,dz)
                vx,vz=ion_velocity(uu,ww,c);mx=px+.5*dd*vx;mz=pz+.5*dd*vz
                uu=(1-fm)*bilinear(u0,mx,mz,lx,z0,dz)+fm*bilinear(u1,mx,mz,lx,z0,dz)
                ww=(1-fm)*bilinear(w0,mx,mz,lx,z0,dz)+fm*bilinear(w1,mx,mz,lx,z0,dz)
                vx,vz=ion_velocity(uu,ww,c)
                before=np.floor(px/lx).astype(int);px=px+dd*vx;pz=np.clip(pz+dd*vz,zgrid[0],zgrid[-1])
                after=np.floor(px/lx).astype(int)
                for q in np.where(before!=after)[0]:events.append((times[k-1]+(j+1)*dd,q,before[q],after[q]))
            track_x.append(px.copy());track_z.append(pz.copy());u0=u1;w0=w1
        track_x=np.asarray(track_x);track_z=np.asarray(track_z)
        bounds=np.quantile(track_x,[.1,.9],axis=1).T+np.array([-2000.,2000.])
        center=np.average(track_x,axis=1,weights=weights)
        # Finite slanted rays on a curved Earth in the zonal plane.
        ranges=np.arange(100000.,140000.1,37.5);scan=np.arange(0.,times[-1]-34,50.)
        offsets=np.array([0.,50/3,100/3]);angles=np.deg2rad([20.,0.,-20.])
        ray_offsets=np.deg2rad(np.linspace(-1.2,1.2,9));aw=np.exp(-.5*(np.rad2deg(ray_offsets)/(.9/2.355))**2);aw/=aw.sum()
        obs=np.zeros((3,len(scan),len(ranges)));mask=np.zeros_like(obs,dtype=bool)
        ambiguous=np.zeros((3,len(scan)),dtype=bool)
        bkeys=sorted(bg.keys());bt=np.array([bg[k].attrs['time_s'] for k in bkeys]);bn=np.stack([density(bg[k]).mean(axis=0) for k in bkeys])
        R=6371000.;cached={}
        def state_at(t):
            q=int(np.clip(np.searchsorted(times,t)-1,0,len(times)-2));f=(t-times[q])/(times[q+1]-times[q])
            for v in [q,q+1]:
                if v not in cached:cached[v]=density(h[keys[v]])
            for v in list(cached):
                if v not in [q,q+1]:del cached[v]
            return q,float(np.clip(f,0,1))
        for b,theta in enumerate(angles):
            for j,t0 in enumerate(scan+offsets[b]):
                for time_offset in [-.25,0.,.25]:
                    t=float(np.clip(t0+time_offset,times[0],times[-1]));q,f=state_at(t)
                    low=(1-f)*bounds[q,0]+f*bounds[q+1,0];high=(1-f)*bounds[q,1]+f*bounds[q+1,1]
                    ambiguous[b,j]=high-low>=lx
                    bq=int(np.clip(np.searchsorted(bt,t)-1,0,len(bt)-2));bf=(t-bt[bq])/(bt[bq+1]-bt[bq]);bg_n=(1-bf)*bn[bq]+bf*bn[bq+1]
                    for dth,wgt in zip(ray_offsets,aw):
                        th=theta+dth
                        alt=np.sqrt((R+ranges*np.cos(th))**2+(ranges*np.sin(th))**2)-R
                        world_x=a.station_x_km*1000+R*np.arctan2(ranges*np.sin(th),R+ranges*np.cos(th))
                        rel_z=alt-(c['altitude_m']+c['layer_descent_m_s']*t)
                        original=(world_x>=low)&(world_x<=high)&(~ambiguous[b,j])
                        inside=(rel_z>=zgrid[0])&(rel_z<=zgrid[-1]);original &= inside
                        # Grid periodicity is accessed only for the selected
                        # unwrapped material image; elsewhere use the control.
                        vals=(1-f)*bilinear(cached[q],world_x,rel_z,lx,z0,dz)+f*bilinear(cached[q+1],world_x,rel_z,lx,z0,dz)
                        background=np.interp(rel_z,zgrid,bg_n,left=0.,right=0.)
                        obs[b,j]+=wgt/3*np.where(original,vals,background)
                        mask[b,j] |= original
                obs[b,j]=gaussian_filter1d(obs[b,j],1/2.355,mode='nearest')
        with h5py.File(a.output.with_suffix('.h5'),'w') as out:
            out.attrs['script']='scripts/synthetic_radar_2d.py';out.attrs['source']=str(a.input.resolve());out.attrs['background_source']=str(a.background.resolve())
            out.attrs['geometry']='Synthetic east/vertical/west zonal-plane beams, not reconstructed SYISR pointings'
            out.attrs['outside_packet']='Unseeded magnetic-ion background control; periodic copies excluded'
            out.attrs['packet_footprint']='Central 80 percent of ion material markers (10th to 90th percentiles), with 2 km margin; all edge trajectories retained'
            out.attrs['range_resolution_fwhm_m']=37.5;out.attrs['beamwidth_fwhm_deg']=.9;out.attrs['integration_s']=.5;out.attrs['scan_cycle_s']=50.
            out.attrs['tracking']='RK2 ion material trajectories; linear interpolation between 10 s velocity snapshots; not vortex-core identities'
            out.attrs['wrap_event_time_resolution_s']=10.;out.attrs['station_x_m']=a.station_x_km*1000
            out.create_dataset('config_json',data=json.dumps(c).encode());out.create_dataset('time_s',data=times)
            out.create_dataset('particle_x_unwrapped_m',data=track_x,compression='gzip');out.create_dataset('particle_z_relative_m',data=track_z,compression='gzip')
            out.create_dataset('particle_wrap_count',data=np.floor(track_x/lx).astype(np.int16),compression='gzip')
            out.create_dataset('wrap_events_time_particle_old_new',data=np.asarray(events).reshape(-1,4),compression='gzip')
            out.create_dataset('packet_bounds_unwrapped_m',data=bounds);out.create_dataset('packet_center_unwrapped_m',data=center)
            out.create_dataset('range_m',data=ranges);out.create_dataset('scan_time_s',data=scan);out.create_dataset('beam_time_offset_s',data=offsets);out.create_dataset('beam_zenith_angle_rad',data=angles)
            out.create_dataset('synthetic_ne_m3',data=obs,compression='gzip');out.create_dataset('original_packet_mask',data=mask,compression='gzip');out.create_dataset('ambiguous_packet_image',data=ambiguous)
        fig,axs=plt.subplots(3,1,figsize=(11,6.5),sharex=True,layout='constrained')
        for b,ax in enumerate(axs):
            alt=(np.sqrt((R+ranges*np.cos(angles[b]))**2+(ranges*np.sin(angles[b]))**2)-R)/1000
            im=ax.pcolormesh((scan+offsets[b])/60,alt,np.maximum(obs[b].T,1),norm=LogNorm(1e10,2e12),cmap='magma',shading='gouraud',rasterized=True)
            ax.set_ylim(103,122);ax.set_ylabel('Altitude [km]');ax.set_title(['East: elevation 70°','Vertical','West: elevation 70°'][b],loc='left',fontsize=11)
        axs[-1].set_xlabel('Elapsed simulation time [min]');axs[-1].set_xlim(0,times[-1]/60)
        cb=fig.colorbar(im,ax=axs,pad=.02,fraction=.035,ticks=[1e10,1e11,2e12]);cb.set_ticklabels([r'$10^{10}$',r'$10^{11}$',r'$2\times10^{12}$']);cb.minorticks_off();cb.set_label(r'Synthetic $n_e$ [m$^{-3}$]')
        fig.savefig(a.output.with_suffix('.png'),dpi=160);fig.savefig(a.output.with_suffix('.pdf'))
        fig.canvas.draw();width,height=fig.get_size_inches()*160
        boxes=np.array([[ax.get_position().x0*width,(1-ax.get_position().y1)*height,ax.get_position().x1*width,(1-ax.get_position().y0)*height] for ax in axs])
        with h5py.File(a.output.with_suffix('.h5'),'a') as out:out.attrs['rti_axes_bounds_pixel']=boxes
        plt.close(fig)
        print(json.dumps({'frames':len(keys),'wrap_events':len(events),'center_initial_final_m':[float(center[0]),float(center[-1])],'packet_beam_scans':mask.any(axis=2).sum(axis=1).tolist(),'ambiguous_scans':ambiguous.sum(axis=1).tolist()},indent=2))

if __name__=='__main__':main()
