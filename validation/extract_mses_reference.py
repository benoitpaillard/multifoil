#!/usr/bin/env python3
"""Reproduce ICAS2006 Figure6 MSES curve samples; default checks, --write saves.

The bundled figure is an unchanged embedded image extracted from PDF page6.
Axis calibration and curve selection were visually audited. Samples represent
the published graphical curve, not recovered native MSES solution files.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

DIRECTORY = Path(__file__).resolve().parent/'data/nlr7301_mses'


def extract(directory=DIRECTORY):
    cfg=json.loads((directory/'digitization.json').read_text())
    image=directory/'source_figure6.png'
    if hashlib.sha256(image.read_bytes()).hexdigest()!=cfg['image_sha256']:
        raise ValueError('Source figure SHA-256 mismatch')
    pixels=np.asarray(Image.open(image).convert('RGB'),dtype=int)
    if pixels.shape[:2]!=(cfg['image_height'],cfg['image_width']):
        raise ValueError('Source figure dimensions changed')
    green=(pixels[:,:,1]>cfg['green_min']) & (
        pixels[:,:,1]>pixels[:,:,0]+cfg['green_channel_margin']) & (
        pixels[:,:,1]>pixels[:,:,2]+cfg['green_channel_margin'])
    (x0,a0),(x1,a1)=cfg['alpha_calibration']
    (y0,c0),(y1,c1)=cfg['cl_calibration']
    rows=[];audit=[]
    for alpha in cfg['sampling_alpha_deg']:
        x=x0+(x1-x0)*(alpha-a0)/(a1-a0)
        mask=green & (abs(np.arange(pixels.shape[1])-x)<=cfg['column_half_width_px'])[None,:]
        ys,_=np.where(mask)
        ys=ys[(ys>=cfg['curve_y_limits_px'][0])&(ys<cfg['curve_y_limits_px'][1])]
        if len(ys)<3 or np.ptp(ys)>8:
            raise ValueError(f'Ambiguous graphical sample at alpha={alpha}; inspect source')
        y=float(np.median(ys));cl=c0+(c1-c0)*(y-y0)/(y1-y0)
        cl=round(cl,3);allowance=cfg['cl_reading_allowance']
        rows.append(dict(alpha_deg=alpha,CL=cl,CL_low=round(cl-allowance,3),
            CL_high=round(cl+allowance,3),alpha_reading_allowance_deg=cfg['alpha_reading_allowance_deg'],
            reference_kind='published_solver',reference_solver='MSES'))
        audit.append(dict(alpha_deg=alpha,pixel_x=round(x,6),pixel_y=y,
            green_pixels=len(ys),pixel_y_min=int(min(ys)),pixel_y_max=int(max(ys))))
    return rows,audit


def csv_text(rows):
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return stream.getvalue()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    rows,audit=extract()
    for name,values in [('polar.csv',rows),('digitization_pixels.csv',audit)]:
        target=DIRECTORY/name;contents=csv_text(values)
        if args.write:
            target.write_text(contents)
        elif target.read_text()!=contents:
            raise SystemExit(f'Reference differs from reproducible extraction: {target}')
    if args.write:
        image=Image.open(DIRECTORY/'source_figure6.png').convert('RGB').resize((988,912))
        draw=ImageDraw.Draw(image)
        for r in audit:
            x,y=2*r['pixel_x'],2*r['pixel_y']
            draw.ellipse((x-8,y-8,x+8,y+8),outline='#c000c0',width=2)
        image.save(DIRECTORY.parent.parent/'audit/nlr7301_mses_digitization.png')
    print(f'{len(rows)} published MSES CL samples '+('written' if args.write else 'reproduced exactly'))


if __name__=='__main__':main()
