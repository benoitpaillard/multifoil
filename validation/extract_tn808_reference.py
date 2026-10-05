"""Auditable manual digitization of TN808 Figure8, fixed gap0.02c, CL1.5.

The clicks locate published marker centres. An affine grid calibration and
bracketed linear interpolation are explicit. No aerodynamic model is fitted.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
DATA=ROOT/'validation/data/naca23012_tn808'
# x-pixel,y-pixel,CL,CD in tn808_fig8_drag.png. The image is a 5x rendering
# of PDF18 cropped to fitz.Rect(123,228,530,506); pixel origin is top left.
CALIBRATION=np.array([[835.,1233.5,1.2,.02],[1033.,1232.5,1.6,.02],
                      [834.,1040.,1.2,.06],[1032.,1039.,1.6,.06]])
CLICKS={10:[[896.,1181.],[1010.,1152.],[1121.,1111.]],
        20:[[833.,1203.],[959.,1177.],[1069.,1150.]],
        30:[[840.,1159.],[953.,1131.],[1057.,1107.]]}


def extract():
    transform=np.linalg.lstsq(np.c_[CALIBRATION[:,:2],np.ones(4)],CALIBRATION[:,2:],rcond=None)[0]
    rows,targets=[],[]
    for angle,pixels in CLICKS.items():
        pixels=np.array(pixels)
        values=np.c_[pixels,np.ones(len(pixels))]@transform
        for point,value in zip(pixels,values):
            rows.append(dict(deflection_deg=angle,pixel_x=point[0],pixel_y=point[1],
                             CL=float(value[0]),CD=float(value[1])))
        target=1.5
        j=int(np.searchsorted(values[:,0],target))
        if not 0<j<len(values):raise ValueError('Target lift is not bracketed')
        cd=float(np.interp(target,values[:,0],values[:,1]))
        targets.append(dict(deflection_deg=angle,CL=target,CD=cd,
            lower_CL=float(values[j-1,0]),upper_CL=float(values[j,0]),
            assumed_measurement_allowance=.002,digitization_allowance=.0012,
            total_cd_allowance=.0032))
    return rows,targets,transform


def main():
    import matplotlib.pyplot as plt
    rows,targets,transform=extract()
    for filename,data in [('digitized_markers.csv',rows),('fixed_lift_reference.csv',targets)]:
        with (DATA/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,list(data[0]));writer.writeheader();writer.writerows(data)
    audit=dict(source='NACA TN808 Figure8, PDF18',reference_kind='experiment',
        image='tn808_fig8_drag.png',image_sha256=hashlib.sha256((DATA/'tn808_fig8_drag.png').read_bytes()).hexdigest(),
        calibration=CALIBRATION.tolist(),affine_pixel_to_cl_cd=transform.tolist(),
        marker_clicks=CLICKS,target_cl=1.5,interpolation='linear between recorded bounding markers; no extrapolation',
        source_precision='CD limits: minimum ±0.0003, CL1 ±0.0006, CL2.5 ±0.002. No value stated at CL1.5.',
        allowance_qualification='Use ±0.002 at CL1.5 as a conservative working assumption, '
            'plus ±0.0012 for approximately3-pixel clicks, calibration, lift mapping and interpolation. '
            'These are comparison allowances, not a measured total-uncertainty bound or confidence interval. '
            'Unknown configuration-dependent fitting bias is not bounded.',
        extraction_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DATA/'digitization.json').write_text(json.dumps(audit,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(11,6),layout='constrained')
    ax.imshow(plt.imread(DATA/'tn808_fig8_drag.png'))
    for (angle,points),color in zip(CLICKS.items(),('tab:blue','tab:orange','tab:green')):
        points=np.array(points)
        ax.scatter(points[:,0],points[:,1],s=100,facecolors='none',edgecolors=color,
                   linewidths=2,label=f'{angle} deg markers')
    ax.scatter(CALIBRATION[:,0],CALIBRATION[:,1],s=55,c='magenta',marker='+',label='Grid calibration')
    ax.set_xlim(770,1150);ax.set_ylim(1250,1010)
    ax.set_xlabel('Pixel x in archived crop');ax.set_ylabel('Pixel y')
    ax.set_title('NACA TN808 Fig.8: digitization audit (source image with marker clicks)')
    ax.legend(loc='upper left',fontsize=9)
    fig.savefig(DATA/'digitization_audit.png',dpi=180)
    print(json.dumps(targets,indent=2))


if __name__=='__main__':main()
