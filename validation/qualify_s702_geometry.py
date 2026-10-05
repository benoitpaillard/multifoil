#!/usr/bin/env python3
"""Compare Figure 1 with four independently positioned pressure-figure outlines.

PyMuPDF is needed only to reproduce extraction. No aerodynamic solve or fitting
to pressure/drag is performed. The published x axes set every chordwise datum;
one common vertical translation registers each drawing to the fore element.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from extract_s702 import DATA, PDF_SHA256, URL, vertices, write_csv
from extract_s702_flap_pressure import CASES
from s702_geometry import ROOT, geometry

OUT = ROOT/'validation/results/s702_geometry_reconciliation_20261002'
DRAWINGS = {43:(1100,1094),44:(1133,1127),48:(1103,1097),49:(1135,1129)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register_outline(raw_fore, raw_aft, x0, scale, baseline_fore):
    """Keep the independently calibrated x and isotropic scale; translate z only."""
    if scale <= 0 or raw_fore.shape != baseline_fore.shape:
        raise ValueError('Invalid axis calibration or changed source vertices')
    offset = float(np.mean(baseline_fore[:,1]+raw_fore[:,1]/scale))
    return [np.column_stack(((raw[:,0]-x0)/scale,offset-raw[:,1]/scale))
            for raw in (raw_fore,raw_aft)], offset


def extract(out):
    import fitz
    pdf=DATA/'NASA-CR-20220014964.pdf'
    if sha(pdf)!=PDF_SHA256:raise ValueError('Source PDF changed')
    out.mkdir(parents=True,exist_ok=True)
    doc=fitz.open(pdf);ds=doc[26].get_drawings()
    chord=ds[2]['items'][0];origin=max(chord[1].y,chord[2].y)
    scale=abs(chord[1].y-chord[2].y);z0=chord[1].x
    inputs={str(pdf.relative_to(ROOT)):sha(pdf),str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)}
    # Recheck the original transform directly, without changing any old data.
    for name,i in [('fore',24),('aft_d0',0),('aft_dm2',3),('aft_d5',10),('aft_d10',17)]:
        raw=np.array(vertices(ds[i]));xy=np.column_stack(((origin-raw[:,1])/scale,(z0-raw[:,0])/scale))
        stored=np.loadtxt(DATA/(name+'.dat'),skiprows=1)
        if xy.shape!=stored.shape or np.max(abs(xy-stored))>2e-15:raise ValueError('Original reconstruction changed')
        inputs[str((DATA/(name+'.dat')).relative_to(ROOT))]=sha(DATA/(name+'.dat'))
    for name in ('validation/extract_s702.py','validation/extract_s702_flap_pressure.py','validation/s702_geometry.py'):
        inputs[name]=sha(ROOT/name)
    rows=[];sources=[];contours={}
    for c in CASES:
        page=doc[c['page']-1];ds=page.get_drawings();ids=DRAWINGS[c['page']]
        x0=ds[0]['items'][0][1].x;x1=ds[c['x1_path']]['items'][0][2].x
        baseline=dict(geometry(c['flap']))
        raw=[np.array(vertices(ds[i])) for i in ids]
        transformed,offset=register_outline(*raw,x0,x1-x0,baseline['fore'])
        contours[c['name']]=transformed
        for (name,xy),i,points in zip(zip(('fore','aft'),transformed),ids,raw):
            old=baseline[name]
            if xy.shape!=old.shape:raise ValueError('Contour ordering/vertex counts changed')
            difference=xy-old
            if np.max(abs(difference))>2e-5:raise ValueError('Drawn contours differ: inspect source before using geometry')
            rows.append(dict(case=c['name'],flap=c['flap'],element=name,pdf_page=c['page'],drawing_index=i,
                vertex_count=len(xy),minimum_x=float(xy[:,0].min()),max_abs_dx=float(abs(difference[:,0]).max()),
                max_abs_dz=float(abs(difference[:,1]).max()),rms_distance=float(np.sqrt(np.mean(np.sum(difference**2,axis=1)))),
                foremost_aft_tap_x=.69941 if name=='aft' else None,
                aft_nose_minus_tap_x=float(xy[:,0].min()-.69941) if name=='aft' else None))
            np.savetxt(out/(c['name']+'_'+name+'.dat'),xy,header='S702 pressure-figure outline; common x axis, fore-registered z',comments='',fmt='%.16g')
        sources.append(dict(case=c['name'],pdf_page=c['page'],drawing_indices=ids,x_axis_zero=x0,
            points_per_chord=x1-x0,z_offset=offset,pdf_vertices=[p.tolist() for p in raw],
            transform='x/c=(PDF_x-axis_zero)/points_per_chord; z/c=z_offset-PDF_y/points_per_chord',
            registration='One common z translation fitted to fore vertices only; no x translation, scale fit, rotation or aft adjustment'))
        if c['page'] in (43,48):page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(out/f'original_page_{c["page"]}.png')
    write_csv(out/'outline_comparison.csv',rows)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    for ax,flap in zip(axes,(5,10)):
        aft=dict(geometry(flap))['aft'];ax.plot(*aft.T,color='black',lw=2,label='Figure 1 input')
        for c in CASES:
            if c['flap']==flap:ax.plot(*contours[c['name']][1].T,'--',label=f'Pressure figure, page {c["page"]}')
        ax.axvline(.69941,color='#c44e52',label='Table II first tap x')
        ax.set(xlim=(.6988,.704),ylim=(-.033,-.019),xlabel='x / c',ylabel='z / c',title=f'Aft nose, flap {flap}°')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Independent figure axes reproduce the same contour/tap inconsistency\nNo geometry correction is established by these drawings')
    for ext in ('png','svg','pdf'):fig.savefig(out/('source_geometry_comparison.'+ext),dpi=180)
    plt.close(fig)
    summary=dict(original_figure1_reextraction_matches=True,pressure_outline_comparisons=len(rows),
        maximum_chordwise_difference=max(r['max_abs_dx'] for r in rows),
        maximum_ordinate_difference_after_fore_registration=max(r['max_abs_dz'] for r in rows),
        aft_nose_minus_tap_range=[min(r['aft_nose_minus_tap_x'] for r in rows if r['element']=='aft'),max(r['aft_nose_minus_tap_x'] for r in rows if r['element']=='aft')],
        correction_justified=False,exact_tested_coordinates_recovered=False,new_flow_solves=0,
        interpretation='Reconstruction follows repeated published outlines; the common source contour/tap inconsistency remains. It does not establish the as-built shape.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest=dict(source_url=URL,product='NASA/CR-20220014964',pdf_sha256=PDF_SHA256,access_mode='public anonymous, existing verified PDF',
        download_status='existing source verified',scope='Figure 1 and pressure-figure geometry on PDF pages 43,44,48,49; Table II(b) foremost tap',
        pymupdf=fitz.VersionBind,input_sha256=inputs,drawings=sources,
        files_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name not in ('source_manifest.json','FINDINGS.md','public_access_manifest.json')})
    (out/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=OUT)
    extract(parser.parse_args().out)
