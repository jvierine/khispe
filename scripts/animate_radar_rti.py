#!/usr/bin/env python3
"""Animate a time cursor over the saved finite synthetic radar RTI."""
import argparse,subprocess,json
from pathlib import Path
import h5py,numpy as np
from PIL import Image,ImageDraw,ImageFont
p=argparse.ArgumentParser();p.add_argument('rti',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
base=Image.open(a.rti.with_suffix('.png')).convert('RGB');W,H=base.size;W-=W%2;H-=H%2;base=base.crop((0,0,W,H))
with h5py.File(a.rti.with_suffix('.h5')) as h:times=h['time_s'][:];offset=h['beam_time_offset_s'][:]/60;axes=h.attrs['rti_axes_bounds_pixel']
enc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pixel_format','rgb24','-video_size',f'{W}x{H}','-framerate','25','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','23','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)],stdin=subprocess.PIPE)
for t in times:
 im=base.copy();d=ImageDraw.Draw(im)
 for b,(x0,y0,x1,y1) in enumerate(axes):
  x=x0+(t/times[-1])*(x1-x0);d.line((x,y0,x,y1),fill='#00a8aa',width=3)
 enc.stdin.write(im.tobytes())
enc.stdin.close();assert enc.wait()==0
print(a.output)
