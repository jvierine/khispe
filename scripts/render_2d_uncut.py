#!/usr/bin/env python3
"""Unmasked periodic solver field, in a frame moving with mean wind/descent.
No particle quantile clipping, background replacement, or temporal morphing.
"""
import argparse,json,subprocess
from pathlib import Path
import h5py,numpy as np
from PIL import Image,ImageDraw,ImageFont
from synthetic_radar_2d import bilinear
from render_cat_eye_cinematic_smooth import NE_COLORS,VORTICITY_COLORS
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path);p.add_argument('--field',choices=['ne','vorticity'],default='ne');a=p.parse_args()
W,H=1400,800;left,top,pw,ph=80,55,1280,640
font=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',22);small=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',18)
colors=NE_COLORS if a.field=='ne' else VORTICITY_COLORS
pal=np.column_stack([np.interp(np.linspace(0,1,256),np.linspace(0,1,len(colors)),colors[:,k]) for k in range(3)]).astype(np.uint8)
a.output.parent.mkdir(parents=True,exist_ok=True)
enc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pixel_format','rgb24','-video_size',f'{W}x{H}','-framerate','30','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)],stdin=subprocess.PIPE)
with h5py.File(a.source) as h:
 keys=sorted(h);times=[]
 c=json.loads(bytes(h[keys[0]]['config_json'][:]).decode());lx=c['lx_m']
 if a.field=='vorticity':
  limits=[]
  for k in keys[::24]:
   g=h[k];zz=g['altitude_m'][:]-c['altitude_m']-c['layer_descent_m_s']*g.attrs['time_s'];v=g['omega_s-1'][:]-c['wind_m_s']/c['shear_m']/np.cosh(zz/c['shear_m'])**2
   limits.append(np.quantile(v,[.001,.999]))
  vmin,vmax=np.min(limits),np.max(limits)
 else:vmin,vmax=10,np.log10(2e12)
 for j,k in enumerate(keys):
  g=h[k];t=float(g.attrs['time_s']);times.append(t);z=g['altitude_m'][:];zc=c['altitude_m']+c['layer_descent_m_s']*t
  center=c.get('seed_packet_center_m',lx/2)+c['advection_m_s']*t;x=np.linspace(center-lx/2,center+lx/2,pw,endpoint=False)
  zz=np.linspace(z[-1],z[0],ph);xx,zmesh=np.meshgrid(x,zz)
  if a.field=='ne':field=g['log10_ne_quantized'][:].astype(float)/4096 if 'log10_ne_quantized' in g else np.log10(np.maximum(g['ion_density_m-3'][:],1))
  else:field=g['omega_s-1'][:]-c['wind_m_s']/c['shear_m']/np.cosh((z-zc)/c['shear_m'])**2
  vals=bilinear(field,xx,zmesh,lx,z[0],z[1]-z[0]);idx=np.rint(np.clip((vals-vmin)/(vmax-vmin),0,1)*255).astype(np.uint8)
  canvas=Image.new('RGB',(W,H),'#07121b');canvas.paste(Image.fromarray(pal[idx]),(left,top));d=ImageDraw.Draw(canvas)
  d.text((left,15),'Electron density' if a.field=='ne' else 'Neutral spanwise vorticity',font=font,fill='white');d.text((W-200,15),f'{t/60:.1f} min',font=font,fill='#ffd565')
  for ztick in np.arange(np.ceil(z[0]/1000),z[-1]/1000+.01):
   yy=top+(z[-1]-ztick*1000)/(z[-1]-z[0])*ph;d.text((left-55,yy-10),f'{ztick:g}',font=small,fill='white')
  d.text((5,top-25),'km',font=small,fill='white')
  for xtick in [-80,-40,0,40,80]:
   xp=left+(xtick*1000+lx/2)/lx*pw;d.text((xp-18,top+ph+8),str(xtick),font=small,fill='white')
  d.text((left+pw//2-165,top+ph+35),'Distance in mean-flow frame [km]',font=small,fill='white')
  bx,by,bw=460,765,480;canvas.paste(Image.fromarray(np.tile(pal[None,:,:],(12,1,1))).resize((bw,12)),(bx,by))
  lo,hi=('10¹⁰','2 × 10¹² m⁻³') if a.field=='ne' else (f'{vmin:.3f}',f'{vmax:.3f} s⁻¹')
  d.text((bx-65,by-4),lo,font=small,fill='white');d.text((bx+bw+10,by-4),hi,font=small,fill='white');enc.stdin.write(canvas.tobytes())
  if j in [0,180,360,540,720]:canvas.save(a.output.with_name(a.output.stem+f'_frame{j:04d}.png'))
 enc.stdin.close();assert enc.wait()==0;canvas.save(a.output.with_suffix('.png'))
a.output.with_suffix('.json').write_text(json.dumps({'script':'scripts/render_2d_uncut.py','source':str(a.source.resolve()),'field':a.field,'physical_frames':len(keys),'time_range_s':[times[0],times[-1]],'mask':False,'background_replacement':False,'coordinate_frame':'mean wind comoving; prescribed descending phase altitude','temporal_interpolation':False,'periodic_fluid_domain':True,'observational_fit':False,'display_limits':[float(vmin),float(vmax)]},indent=2))
print(a.output)
