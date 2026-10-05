#!/usr/bin/env python3
"""Plot audited S702 geometry and matched-lift MultiFoil/MSES/experiment drag."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT/'validation/results/s702_panel_convergence_analysis'
STUDY = ROOT/'validation/results/s702_panel_convergence_20260922'
REFERENCE = ROOT/'validation/data/s702/re700k_holdout'
OUT = ROOT/'validation/results/s702_cl12_geometry_mses_20260924'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    files = {}
    def track(path, expected=None):
        digest = sha(path)
        if expected is not None and digest != expected:
            raise ValueError('Changed source: '+str(path))
        files[str(path.relative_to(ROOT))] = digest
        return path
    def read_json(path, expected=None):
        return json.loads(track(path,expected).read_text())
    def rows(path, expected=None):
        with track(path,expected).open(newline='') as stream:
            return list(csv.DictReader(stream))

    receipt = read_json(AUDIT/'verification.json')
    def audited_json(name):
        return read_json(AUDIT/name,receipt['outputs_sha256'][name])
    def audited_rows(name):
        return rows(AUDIT/name,receipt['outputs_sha256'][name])
    pinned = audited_json('input_manifest.json')['files']
    summary = audited_json('summary.json')
    if not summary['full_state_audit_pass']:
        raise ValueError('Selected-state audit did not pass')
    paired = next(r for r in audited_rows('flap_deltas.csv') if float(r['target_cl']) == 1.2)
    forces = {int(r['flap']):float(r['CD']) for r in audited_rows('forces.csv')
              if float(r['target_cl']) == 1.2 and r['grid'] == 'finer'}
    model = (forces[10]-forces[5])*1e4
    assert abs(model-float(paired['finer_flap_delta_CD_counts'])) < 1e-10

    source = read_json(REFERENCE/'source.json')
    assert (source['reynolds'],source['mach'],source['mses_ncrit'],source['transition']) == (700000.,.07,9.,'free')
    values, interpolation = {}, []
    for series in ('mses','experiment'):
        cd = {}
        for flap in (5,10):
            name = f'd+{flap}/{series}_drag.csv'
            table = rows(REFERENCE/name,source['sha256'][name])
            cl = np.array([float(r['CL']) for r in table])
            drag = np.array([float(r['CD']) for r in table])
            if not (np.isfinite(cl).all() and np.isfinite(drag).all()
                    and np.all(np.diff(cl)>0) and cl[0] < 1.2 < cl[-1]):
                raise ValueError('Expected finite, ordered reference spanning CL=1.2')
            cd[flap] = float(np.interp(1.2,cl,drag))
            i = int(np.searchsorted(cl,1.2))
            interpolation.append(dict(series=series,flap=flap,CD=cd[flap],
                bracketing_CL=cl[i-1:i+1].tolist(),bracketing_CD=drag[i-1:i+1].tolist()))
        values[series] = (cd[10]-cd[5])*1e4
    comparison_path = ROOT/'validation/results/s702_drag_mechanism/s702_flap_5_to_10.csv'
    historic = next(r for r in rows(comparison_path,pinned[str(comparison_path.relative_to(ROOT))])
                    if (float(r['reynolds']),float(r['mach']),float(r['target_cl'])) == (700000.,.07,1.2))
    assert abs(values['mses']-float(historic['mses_n9_delta_10_minus_5_counts'])) < 1e-10
    assert abs(values['experiment']-float(paired['experiment_flap_delta_counts'])) < 1e-10

    geometry = {}
    for flap in (5,10):
        for element in ('fore','aft'):
            path = STUDY/f'd{flap}_cl1p2_finer'/f'{element}_state.npz'
            with np.load(track(path,pinned[str(path.relative_to(ROOT))]),allow_pickle=False) as state:
                geometry[flap,element] = state['xy'].T.copy()
    assert np.array_equal(geometry[5,'fore'],geometry[10,'fore'])
    assert all(np.isfinite(xy).all() for xy in geometry.values())

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
    fig = plt.figure(figsize=(12,9.2),facecolor='white')
    fig.text(.045,.955,'S702: geometry and drag response',fontsize=22,weight='bold',color='#142c3d')
    fig.text(.045,.909,'CL = 1.2  •  flap 5° → 10°  •  Re = 700,000  •  Mach = 0.07',fontsize=13,color='#4a5967')
    fig.text(.045,.858,'Geometry used by MultiFoil',fontsize=13,weight='bold',color='#142c3d')
    fig.text(.045,.830,'Reconstructed contours; common reference axes, equal x/y scale.',fontsize=10,color='#4a5967')
    full = fig.add_axes([.06,.64,.54,.18])
    zoom = fig.add_axes([.70,.63,.265,.225])
    fore = geometry[5,'fore']
    for ax in (full,zoom):
        ax.fill(fore[:,0],fore[:,1],fc='#e8edf0',ec='#606e78',lw=1.4,label='Fore element')
        for flap,color,style in ((5,'#c58516','-'),(10,'#b6424b','--')):
            xy = geometry[flap,'aft']; xy = np.vstack((xy,xy[0]))
            ax.plot(xy[:,0],xy[:,1],color=color,ls=style,lw=1.8,label=f'Flap {flap}°')
        ax.set_aspect('equal',adjustable='box')
        ax.set_xlabel('x / reference chord',fontsize=9)
        ax.set_ylabel('y / reference chord',fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=.17)
        for edge in ('top','right'): ax.spines[edge].set_visible(False)
    full.set(xlim=(-.025,1.04),ylim=(-.065,.10))
    full.add_patch(Rectangle((.86,-.045),.16,.08,fill=False,ec='#7c8c98',ls=':',lw=1))
    full.legend(loc='upper center',bbox_to_anchor=(.5,-.5),ncol=3,frameon=False,fontsize=9)
    zoom.set(xlim=(.86,1.02),ylim=(-.045,.035))
    zoom.set_title('Flap close-up (dotted box)',fontsize=11,pad=10)

    ax = fig.add_axes([.245,.245,.71,.30])
    labels = ['MultiFoil\nNcrit = 6.25 · 440/320 panels','Published MSES\nNcrit = 9','Experiment\ndigitized data']
    numbers = [model,values['mses'],values['experiment']]
    colors = ['#316ea5','#8554a3','#168778']
    ax.barh([2,1,0],numbers,height=.48,color=colors,zorder=3)
    ax.axvline(0,color='#3f4e5a',lw=1.2,zorder=2)
    ax.set_yticks([2,1,0],labels)
    ax.tick_params(axis='y',length=0,pad=12,labelsize=11)
    ax.set(xlim=(-12.5,14.5),ylim=(-.55,2.6),xticks=[-10,-5,0,5,10])
    ax.set_xlabel('Drag change, CD(10°) − CD(5°)  [counts]',labelpad=12)
    ax.grid(axis='x',color='#e2e8ed',lw=.8,zorder=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    for y,value,color in zip([2,1,0],numbers,colors):
        ax.text(value+(.35 if value>0 else -.35),y,
                f'{value:+.2f}\n'+('more drag' if value>0 else 'less drag'),
                va='center',ha='left' if value>0 else 'right',color=color,fontsize=13,weight='bold')
    fig.text(.045,.135,'Positive = drag increase. Negative = drag reduction.  1 count = 0.0001 CD.',fontsize=11,color='#4a5967')
    fig.text(.045,.090,'MSES and experiment: published NASA/CR–20220014964, Figs. 25(a) and 26(b); interpolated in CL.',fontsize=10,color='#4a5967')
    fig.text(.045,.057,'MSES was not run here. Transition settings differ. MultiFoil numerical convergence is not established.',fontsize=10,color='#4a5967')
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ('png','svg','pdf'):
        fig.savefig(OUT/('cl12_geometry_mses.'+ext),dpi=200,facecolor='white')
    plt.close(fig)
    for (flap,element),xy in geometry.items():
        np.savetxt(OUT/f'geometry_d{flap}_{element}.csv',xy,delimiter=',',header='x_over_cref,y_over_cref',comments='')
    track(Path(__file__).resolve())
    data = dict(target_cl=1.2,reynolds=700000,mach=.07,flap_degrees=[5,10],
        model_ncrit=6.25,published_mses_ncrit=9.,model_panels=[440,320],
        drag_changes_counts=dict(multifoil=model,published_mses=values['mses'],experiment=values['experiment']),
        model_minus_experiment_counts=model-values['experiment'],
        mses_minus_experiment_counts=values['mses']-values['experiment'],
        reference_interpolation=interpolation,source_url=source['url'],reference_figures=['25(a)','26(b)'],
        geometry='Audited MultiFoil 440/320 solver contours from reconstructed coordinates; not an original MSES geometry deck.',
        geometry_axes='Common reference frame, not rotated to each solution incidence; equal x/y scales.',
        qualification='Published MSES only. Ncrit differs between methods. No uncertainty intervals inferred. MultiFoil grid-convergence screen failed.',
        reproduce='python3 validation/plot_s702_cl12_comparison.py',input_sha256=files)
    (OUT/'data.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(drag_changes_counts=data['drag_changes_counts'],output=str(OUT)),indent=2))


if __name__ == '__main__':
    main()
