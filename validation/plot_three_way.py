#!/usr/bin/env python3
"""Plot archived MultiFoil, published MSES and experimental results.

Run from any directory. No new flow solves or fitted corrections.
Optional --pdf-source checks the new experimental digitization against the
original COBEM2005 PDF and requires PyMuPDF. See THREE_WAY.md for provenance.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from validate import verify_data, read_csv, number

VALIDATION=ROOT/'validation'
DATA=VALIDATION/'data'
RESULTS=VALIDATION/'results'
BLUE='#0875b5'
ORANGE='#d47713'
INK='#252b32'
RED='#bd2d39'


def check_experiment(pdf=None):
    """Check all nine experimental graph vertices and the optional source PDF."""
    directory=DATA/'naca0012_mses_drag'
    cfg=json.loads((directory/'experiment_digitization.json').read_text())
    segments=cfg['experiment_segments_pt']
    if pdf is not None:
        import fitz
        if hashlib.sha256(pdf.read_bytes()).hexdigest()!=cfg['source_pdf_sha256']:
            raise ValueError('Source PDF hash mismatch')
        recovered=[]
        with fitz.open(pdf) as doc:
            for d in doc[cfg['pdf_page_index']].get_drawings():
                c=d['color']
                if c and c[1]>.8 and c[0]<.01 and c[2]<.01 and \
                        fitz.Rect(*cfg['curve_selection_rect_pt']).contains(d['rect']):
                    recovered += [[list(i[1]),list(i[2])] for i in d['items'] if i[0]=='l']
        if recovered!=segments:
            raise ValueError('Experimental vectors differ from original PDF')
    points=[segments[0][0]]
    for start,end in segments:
        if start!=points[-1] or end[0]<=start[0]:
            raise ValueError('Experimental curve is discontinuous or reverses incidence')
        points.append(end)
    rows=read_csv(directory/'experiment.csv')
    if len(rows)!=len(points):
        raise ValueError('Experimental sample count mismatch')
    (x0,a0),(x1,a1)=cfg['alpha_calibration_pt']
    (y0,c0),(y1,c1)=cfg['cd_calibration_pt']
    for r,(x,y) in zip(rows,points):
        alpha=round(a0+(x-x0)*(a1-a0)/(x1-x0),6)
        cd=round(c0+(y-y0)*(c1-c0)/(y1-y0),6)
        allowance=cfg['cd_reading_allowance']
        expected={'alpha_deg':alpha,'CD':cd,'CD_low':round(cd-allowance,6),
                  'CD_high':round(cd+allowance,6)}
        if any(float(r[k])!=v for k,v in expected.items()) or r['reference_kind']!='experiment':
            raise ValueError('Experimental CSV differs from calibrated vectors')
    return rows


def model(name):
    folder=RESULTS/name
    rows=read_csv(folder/'polar.csv')
    meta=json.loads((folder/'metadata.json').read_text())
    for row,diagnostic in zip(rows,meta['diagnostics']):
        row['supersonic_flag']=any(diagnostic.get('local_supersonic_detected',{}).values())
    return rows


def series(ax,rows,quantity,kind,limit=None):
    selected=[r for r in rows if number(r.get(quantity)) is not None and
              (limit is None or float(r['alpha_deg'])<=limit)]
    x=np.array([float(r['alpha_deg']) for r in selected])
    y=np.array([float(r[quantity]) if r.get('converged','1')=='1' else np.nan for r in selected])
    if kind=='MultiFoil':
        ax.plot(x,y,'o-',color=BLUE,lw=2,ms=4,label=kind,zorder=3)
        flagged=np.array([r.get('supersonic_flag',False) for r in selected],dtype=bool)
        if flagged.any():
            ax.plot(x[flagged],y[flagged],'x',color=RED,ms=8,mew=1.6,
                    label='Local supersonic flag',zorder=6)
    else:
        color=ORANGE if kind.startswith('MSES') else INK
        marker='D' if kind=='MSES official' else ('s' if kind=='MSES' else 'o')
        style='--' if kind=='MSES' else 'none'
        label='MSES: MIT example (qualified)' if kind=='MSES official' else ('Published MSES' if kind=='MSES' else 'Experiment')
        ax.plot(x,y,marker=marker,ls=style,lw=1.5,ms=4.8,mfc='white',mec=color,
                color=color,label=label,zorder=4)
        if selected and all(number(r.get(quantity+'_low')) is not None for r in selected):
            low=np.array([float(r[quantity+'_low']) for r in selected])
            high=np.array([float(r[quantity+'_high']) for r in selected])
            ax.errorbar(x,y,yerr=[y-low,high-y],fmt='none',ecolor=color,
                        alpha=.55,lw=.9,capsize=2,zorder=2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=RESULTS/'three_way')
    parser.add_argument('--pdf-source',type=Path)
    args=parser.parse_args()
    verify_data()
    naca_exp=check_experiment(args.pdf_source)
    nlr_exp=read_csv(DATA/'nlr7301/polar.csv')
    nlr_mses=read_csv(DATA/'nlr7301_mses/polar.csv')
    nlr_official=read_csv(DATA/'nlr7301_mses_official/polar.csv')
    naca_mses=read_csv(DATA/'naca0012_mses_drag/polar.csv')
    nlr_model=model('nlr7301_mses')
    naca_model=model('naca0012_mses_drag')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,
        'axes.labelsize':11,'axes.titlesize':12,'axes.spines.top':False,
        'axes.spines.right':False,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.7))
    panels=[
        ('A','NLR7301 — lift','Re = 2.51 million  |  M = 0.185  |  nominal gap = 2.6%',
         'CL',nlr_model,nlr_mses,nlr_exp,None),
        ('B','NLR7301 — drag','Re = 2.51 million  |  M = 0.185  |  nominal gap = 2.6%',
         'CD',nlr_model,[],nlr_exp,None),
        ('C','NACA0012 — drag','Re = 6 million  |  M = 0.30',
         'CD',naca_model,naca_mses,naca_exp,None),
        ('D','NACA0012 — drag, low-angle detail','Same case as C  |  0–10°',
         'CD',naca_model,naca_mses,naca_exp,10)]
    export=[]
    for ax,(letter,title,conditions,q,computed,mses,experiment,limit) in zip(axes.flat,panels):
        sources=[('MultiFoil',computed),('MSES',mses),('Experiment',experiment)]
        if letter in 'AB':sources.append(('MSES official',nlr_official))
        for kind,rows in sources:
            if not rows:continue
            series(ax,rows,q,kind,limit)
            for row in rows:
                if number(row.get(q)) is None or (limit is not None and float(row['alpha_deg'])>limit):continue
                export.append(dict(panel=letter,case='nlr7301' if letter in 'AB' else 'naca0012_mses_drag',
                    series=kind,quantity=q,alpha_deg=float(row['alpha_deg']),value=float(row[q]),
                    reading_low=number(row.get(q+'_low')),reading_high=number(row.get(q+'_high')),
                    converged=row.get('converged',''),local_supersonic_flag=row.get('supersonic_flag',''),
                    qualification=row.get('qualification','')))
        ax.set_title(f'{letter}  {title}',loc='left',fontweight='semibold',pad=27)
        ax.text(0,1.025,conditions,transform=ax.transAxes,color='#5b636c',fontsize=9.2)
        ax.set_xlabel('Angle of attack, α [deg]')
        ax.set_ylabel('$C_L$' if q=='CL' else '$C_D$')
        ax.grid(color='#dde2e7',lw=.65,alpha=.8)
        ax.set_axisbelow(True)
        ax.legend(loc='upper left',frameon=True,framealpha=.94,edgecolor='none',fontsize=9)
        if letter=='A':ax.set(xlim=(-.5,16.6),ylim=(1.45,3.4))
        if letter=='B':
            ax.set(xlim=(-.5,15.6),ylim=(0,.123))
            ax.text(.035,.67,'One printed MSES point at 6°\nExact input geometry unconfirmed',transform=ax.transAxes,
                    fontsize=9,color='#6a5250',bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
        if letter=='C':ax.set(xlim=(-.5,18.6),ylim=(0,.194))
        if letter=='D':ax.set(xlim=(-.4,10.4),ylim=(.003,.0142),xticks=[0,2,4,6,8,10])
    fig.suptitle('MultiFoil vs MSES vs experiment',x=.075,y=.987,ha='left',
                 fontsize=20,fontweight='semibold',color='#1e2935')
    fig.text(.075,.951,'Archived MultiFoil solutions • published reference curves • nominal conditions matched',
             fontsize=11,color='#5b636c')
    fig.subplots_adjust(left=.075,right=.985,bottom=.125,top=.866,hspace=.48,wspace=.22)
    fig.text(.075,.067,'Bars show graphical reading allowances: MSES lift ±0.015; NACA drag ±0.0012 (12 counts).',fontsize=9.2,color='#4c5560')
    fig.text(.075,.044,'MultiFoil assumes Ncrit = 9. Red crosses flag local supersonic flow; the model does not resolve shocks.',fontsize=9.2,color='#4c5560')
    fig.text(.075,.021,'Sources: Neto et al. (ICAS 2006); Drela / MIT MSES example; Ranneberg / NLR; Lima e Silva et al. (COBEM 2005), Fig. 3.',fontsize=8.5,color='#64707b')
    for suffix in ['png','pdf']:
        fig.savefig(out/f'multifoil_mses_experiment.{suffix}',dpi=180)
    plt.close(fig)
    with (out/'plot_data.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(export[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(export)
    print(f'Saved PNG, PDF and {len(export)} plotted samples in {out}')
    print('Nine experimental NACA drag samples reproduced from source vectors; no new solver runs.')


if __name__=='__main__':main()
