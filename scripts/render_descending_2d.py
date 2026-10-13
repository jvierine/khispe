#!/usr/bin/env python3
"""Render actual checkpoints on fixed altitude axes; prescribed phase descent."""
import argparse
import json
import subprocess
from contextlib import nullcontext
from pathlib import Path
import h5py
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import map_coordinates
from render_cat_eye_cinematic_smooth import NE_COLORS

p=argparse.ArgumentParser()
p.add_argument('frames',type=Path)
p.add_argument('output',type=Path)
p.add_argument('--zmin',type=float,default=109.)
p.add_argument('--zmax',type=float,default=122.)
a=p.parse_args()
packed=h5py.File(a.frames,'r') if a.frames.is_file() else None
paths=sorted(packed.keys()) if packed is not None else sorted(a.frames.glob('frame_*.h5'))
font=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',20)
small=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',16)
W,H=1100,660
left,top,pw,ph=80,65,980,500
zmin,zmax=a.zmin,a.zmax
pal=np.column_stack([np.interp(np.linspace(0,1,256),np.linspace(0,1,len(NE_COLORS)),NE_COLORS[:,k]) for k in range(3)]).astype(np.uint8)
images=[]
times=[]
a.output.parent.mkdir(parents=True,exist_ok=True)
encoder=None
if a.output.suffix=='.mp4':
 encoder=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pixel_format','rgb24','-video_size',f'{W}x{H}','-framerate','25','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','23','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)],stdin=subprocess.PIPE)
for path in paths:
 with (nullcontext(packed[path]) if packed is not None else h5py.File(path)) as h:
  n=None if 'log10_ne_quantized' in h else h['ion_density_m-3'][:]
  native_ln=h['log10_ne_quantized'][:].astype(float)/4096 if n is None else np.log10(np.maximum(n,1))
  z=h['altitude_m'][:]/1000
  t=float(h.attrs['time_s'])
  cfg=json.loads(bytes(h['config_json'][:]).decode())
 xi=np.linspace(0,native_ln.shape[0]-1,pw)
 zi=(np.linspace(zmax,zmin,ph)-z[0])/(z[1]-z[0])
 xx,zz=np.meshgrid(xi,zi)
 ln=map_coordinates(native_ln,[xx,zz],order=1,mode='constant',cval=0)
 idx=np.rint(np.clip((ln-10)/(np.log10(2e12)-10),0,1)*255).astype(np.uint8)
 canvas=Image.new('RGB',(W,H),'#07121b')
 canvas.paste(Image.fromarray(pal[idx]),(left,top))
 d=ImageDraw.Draw(canvas)
 label = '2-D magnetic ion drift | prescribed layer descent' if cfg['ion_mode']=='collisional' else '2-D ion-layer density | prescribed descending shear frame'
 d.text((left,20),label,font=font,fill='white')
 ri=cfg['n2_s2']/(cfg['wind_m_s']/cfg['shear_m'])**2
 d.text((left+660,20),f'Ri = {ri:.2f}',font=font,fill='white')
 d.text((W-185,20),f'{t/60:.1f} min',font=font,fill='#ffd565')
 for zk in range(int(np.ceil(zmin/2)*2),int(zmax)+1,2):
  yy=top+(zmax-zk)/(zmax-zmin)*ph
  d.text((left-45,yy-8),str(zk),font=small,fill='white')
 d.text((3,top-25),'km',font=small,fill='white')
 for xk in [0,42,84,126,168]:
  xp=left+xk/(cfg['lx_m']/1000)*pw
  d.text((xp-12,top+ph+8),str(xk),font=small,fill='white')
 d.text((left+pw//2-95,top+ph+30),'Zonal distance [km]',font=small,fill='white')
 bx,by,bw=300,625,480
 grad=Image.fromarray(np.tile(pal[None,:,:],(12,1,1))).resize((bw,12))
 canvas.paste(grad,(bx,by))
 d.text((bx-55,by-2),'1e10',font=small,fill='white')
 d.text((bx+bw+8,by-2),'2e12 m^-3',font=small,fill='white')
 if encoder is not None:
  encoder.stdin.write(canvas.tobytes())
 else:
  images.append(canvas)
 times.append(t)
a.output.parent.mkdir(parents=True,exist_ok=True)
# One frame per actual 10 s checkpoint; no time interpolation.
if encoder is not None:
 encoder.stdin.close()
 if encoder.wait()!=0: raise RuntimeError('ffmpeg failed')
 canvas.save(a.output.with_suffix('.png'))
else:
 images[0].save(a.output,save_all=True,append_images=images[1:],duration=40,loop=0)
 images[-1].save(a.output.with_suffix('.png'))
a.output.with_suffix('.json').write_text(json.dumps({'physical_frames':len(paths),'time_range_s':[times[0],times[-1]],'color_limits_m3':[1e10,2e12],'prescribed_layer_descent_m_s':cfg['layer_descent_m_s'],'ion_mode':cfg['ion_mode'],'observational_fit':False,'script':'scripts/render_descending_2d.py'},indent=2)+'\n')
print(a.output)
