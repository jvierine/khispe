#!/usr/bin/env python3
"""Render a single unwrapped packet and finite radar, from physical checkpoints."""
import argparse,json,subprocess
from pathlib import Path
import h5py,numpy as np
from PIL import Image,ImageDraw,ImageFont
from synthetic_radar_2d import bilinear
from render_cat_eye_cinematic_smooth import NE_COLORS
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('tracking',type=Path);p.add_argument('background',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.parent.mkdir(parents=True,exist_ok=True)
W,H=1100,660;left,top,pw,ph=80,65,980,500
font=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',20);small=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',16)
pal=np.column_stack([np.interp(np.linspace(0,1,256),np.linspace(0,1,len(NE_COLORS)),NE_COLORS[:,k]) for k in range(3)]).astype(np.uint8)
enc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pixel_format','rgb24','-video_size',f'{W}x{H}','-framerate','25','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','22','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)],stdin=subprocess.PIPE)
with h5py.File(a.source) as h,h5py.File(a.tracking) as tr,h5py.File(a.background) as bg:
 keys=sorted(h);bkeys=sorted(bg);bt=np.array([bg[k].attrs['time_s'] for k in bkeys]);center=tr['packet_center_unwrapped_m'][:];bounds=tr['packet_bounds_unwrapped_m'][:];times=tr['time_s'][:]
 for j,k in enumerate(keys):
  g=h[k];c=json.loads(bytes(g['config_json'][:]).decode());lx=c['lx_m'];t=times[j];z=g['altitude_m'][:];zc=c['altitude_m']+c['layer_descent_m_s']*t
  x=np.linspace(center[j]-lx/2,center[j]+lx/2,pw);zz=np.linspace(zc+3000,zc-3000,ph);xx,zmesh=np.meshgrid(x,zz)
  ln=bilinear(g['log10_ne_quantized'][:].astype(float)/4096,xx,zmesh,lx,z[0],z[1]-z[0]);bq=int(np.argmin(abs(bt-t)));bprof=bg[bkeys[bq]]['log10_ne_quantized'][:].astype(float).mean(axis=0)/4096
  background=np.interp(zz,z,bprof,left=0,right=0)[:,None]
  selected=(xx>=bounds[j,0])&(xx<=bounds[j,1]);ambiguous=bounds[j,1]-bounds[j,0]>=lx
  ln=np.where(selected & ~np.bool_(ambiguous),ln,background);idx=np.rint(np.clip((ln-10)/(np.log10(2e12)-10),0,1)*255).astype(np.uint8)
  canvas=Image.new('RGB',(W,H),'#07121b');canvas.paste(Image.fromarray(pal[idx]),(left,top));d=ImageDraw.Draw(canvas)
  d.text((left,20),'Localized Ri = 0.08 | single unwrapped packet',font=font,fill='white');d.text((W-190,20),f'{t/60:.1f} min',font=font,fill='#ffd565')
  for ztick in np.arange(np.ceil((zc-3000)/1000),(zc+3000)/1000+.01):
   yy=top+(zc+3000-ztick*1000)/6000*ph;d.text((left-55,yy-8),f'{ztick:g}',font=small,fill='white')
  d.text((3,top-25),'km',font=small,fill='white')
  for xtick in np.arange(np.ceil(x[0]/40000)*40,x[-1]/1000+.01,40):
   xp=left+(xtick*1000-x[0])/lx*pw;d.text((xp-20,top+ph+8),f'{xtick:g}',font=small,fill='white')
  d.text((left+pw//2-125,top+ph+30),'Unwrapped zonal distance [km]',font=small,fill='white')
  station=float(tr.attrs['station_x_m']);xp=left+(station-x[0])/lx*pw
  if left<=xp<=left+pw:d.line((xp,top,xp,top+ph),fill='#eeeeee',width=1);d.text((xp+5,top+8),'Radar',font=small,fill='white')
  bx,by,bw=300,625,480;canvas.paste(Image.fromarray(np.tile(pal[None,:,:],(12,1,1))).resize((bw,12)),(bx,by));d.text((bx-55,by-2),'1e10',font=small,fill='white');d.text((bx+bw+8,by-2),'2e12 m^-3',font=small,fill='white');enc.stdin.write(canvas.tobytes())
 enc.stdin.close();assert enc.wait()==0;canvas.save(a.output.with_suffix('.png'))
a.output.with_suffix('.json').write_text(json.dumps({'script':'scripts/render_localized_packet.py','source':str(a.source.resolve()),'tracking':str(a.tracking.resolve()),'physical_frames':len(keys),'time_range_s':[float(times[0]),float(times[-1])],'vertical_view':'follows prescribed descending phase frame','observational_fit':False},indent=2))
print(a.output)
