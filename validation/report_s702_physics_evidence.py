#!/usr/bin/env python3
"""Quantitative evidence for a physical diagnosis; uses saved states, no CFD."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from report_s702_ncrit_flap_error import Inputs, ROOT, MECH, json_write, write_csv, select, require, sha


def interp(points, target, x='CL', y='CD'):
    data = sorted((float(r[x]), float(r[y])) for r in points)
    require(data[0][0] <= target <= data[-1][0], 'Reference extrapolation prohibited')
    return float(np.interp(target, *np.asarray(data).T))


def report(out):
    inputs = Inputs()
    comparison = inputs.read(MECH + 'comparison_points.csv')
    states = inputs.read(MECH + 'selected_states.csv')
    surfaces = inputs.read(MECH + 'surface_summary.csv')
    separation = inputs.read(MECH + 'separation_intervals.csv')
    deltas = inputs.read(MECH + 'element_drag_deltas.csv')
    curves, angles, neighborhoods, offsets, minima = {}, [], [], [], []
    for dataset, base in (('calibration', 'validation/data/s702/'),
                          ('reserved', 'validation/data/s702/re700k_holdout/')):
        for flap in (5, 10):
            for method in ('experiment', 'mses'):
                curve = inputs.read(base + f'd+{flap}/{method}_drag.csv')
                curves[dataset, flap, method] = curve
                value = interp(curve, .8)
                original = select(comparison, dataset=dataset, target_cl='0.8', flap_deg=str(float(flap)),
                                  series='experiment' if method == 'experiment' else 'mses_n9')
                require(abs(value-float(original['CD'])) < 1e-12, 'Independent polar interpolation differs')
            lift = inputs.read(base + f'd+{flap}/experiment_lift.csv')
            branch = [r for r in lift if float(r['alpha_deg']) < 8]
            for cl in (.8, 1.2):
                state = select(states, dataset=dataset, target_cl=str(cl), ncrit='6.25', flap=str(flap))
                raw = inputs.read(state['result_path'], 'json')
                require(raw['feasible'] and abs(raw['CD']-float(state['CD'])) < 1e-12, 'Saved model state differs')
                ef = raw['element_forces']
                angles.append(dict(dataset=dataset, flap=flap, target_cl=cl,
                    experiment_alpha_deg=interp(branch, cl, y='alpha_deg'),
                    model_alpha_deg=float(state['alpha_deg']), model_fore_CL=ef['fore']['CL'],
                    model_aft_CL=ef['aft']['CL'], source_state=state['state_id']))
        for cl in (.7, .75, .8, .85, .9, 1., 1.2):
            cd5 = interp(curves[dataset, 5, 'experiment'], cl)
            cd10 = interp(curves[dataset, 10, 'experiment'], cl)
            neighborhoods.append(dict(dataset=dataset, target_cl=cl, experiment_cd5=cd5,
                                      experiment_cd10=cd10, experiment_delta_counts=(cd10-cd5)*1e4))
        e5 = interp(curves[dataset, 5, 'experiment'], .8)
        xy = sorted((float(r['CL']), float(r['CD'])) for r in curves[dataset, 5, 'mses'])
        crossings = [a[0]+(e5-a[1])/(b[1]-a[1])*(b[0]-a[0]) for a,b in zip(xy, xy[1:])
                     if .6<a[0]<1.2 and (a[1]-e5)*(b[1]-e5)<=0 and a[1]!=b[1]]
        require(len(crossings)==1, 'Expected one local MSES flap-5 upper-bucket crossing')
        offsets.append(dict(dataset=dataset, experiment_target_cl=.8, matched_cd=e5,
                            mses_flap5_cl_for_same_cd=crossings[0], cl_offset=crossings[0]-.8))
        # Include interval endpoints to calculate the minimum of the piecewise-linear curve on this support.
        restricted = [(float(r['CL']), float(r['CD'])) for r in curves[dataset,10,'mses'] if .5<=float(r['CL'])<=1.2]
        restricted += [(cl, interp(curves[dataset,10,'mses'], cl)) for cl in (.5, 1.2)]
        location, value = min(restricted, key=lambda r:r[1])
        observed = interp(curves[dataset,10,'experiment'], .8)
        minima.append(dict(dataset=dataset, support_cl_min=.5, support_cl_max=1.2,
            mses_flap10_min_cl=location, mses_flap10_min_cd=value,
            experiment_flap10_cl08_cd=observed, minimum_gap_counts=(value-observed)*1e4))
    focus_ids = {r['source_state'] for r in angles if r['dataset']=='reserved'}
    selected_surfaces = [r for r in surfaces if r['state_id'] in focus_ids]
    selected_separation = [r for r in separation if r['state_id'] in focus_ids]
    require(not any(r['reaches_trailing_edge']=='True' for r in selected_separation if r['element']=='aft'),
            'Diagnosis requires revisiting: a selected aft separated interval reaches the TE')
    out.mkdir(parents=True, exist_ok=True)
    for filename, rows in [('matched_lift_angles.csv',angles), ('experimental_neighborhood.csv',neighborhoods),
                           ('bucket_position.csv',offsets), ('mses_flap10_minimum.csv',minima),
                           ('selected_surfaces.csv',selected_surfaces), ('selected_separation.csv',selected_separation)]:
        write_csv(out/filename, rows)
    fig, axes = plt.subplots(1,2,figsize=(12,4.9),layout='constrained')
    for flap,color in ((5,'#3278a5'),(10,'#d18323')):
        for method in ('experiment','mses'):
            rows = sorted((float(r['CL']),float(r['CD'])*1e4) for r in curves['reserved',flap,method]
                          if .45<=float(r['CL'])<=1.35)
            axes[0].plot(*np.asarray(rows).T, color=color,
                         linestyle='-' if method=='experiment' else '--', marker='o' if method=='experiment' else None,
                         ms=3,label=f'{flap}° '+('experiment' if method=='experiment' else 'published MSES N=9'))
    axes[0].axvline(.8,color='#777',lw=.8)
    axes[0].set(xlabel='CL',ylabel='CD, drag counts',title='Two different polar discrepancies')
    axes[0].legend(fontsize=8)
    for j,cl in enumerate((.8,1.2)):
        row = select(deltas,dataset='reserved',ncrit='6.25',target_cl=str(cl))
        for offset,key,color,label in [(-.22,'fore_delta_wake_counts','#3278a5','Model fore contribution'),
                                      (0.,'aft_delta_wake_counts','#d18323','Model aft contribution')]:
            axes[1].bar(j+offset,float(row[key]),width=.19,color=color,label=label if j==0 else None)
        axes[1].plot(j+.22,float(row['delta_CD_counts']),'o',color='#333',label='Model total' if j==0 else None)
        axes[1].plot(j+.22,float(row['experiment_delta_counts']),'D',color='#a13940',label='Measured total' if j==0 else None)
    axes[1].axhline(0,color='#888',lw=.8)
    axes[1].set_xticks([0,1],['CL = 0.8','CL = 1.2'])
    axes[1].set(ylabel='CD(10°) − CD(5°), drag counts',title='Modeled fore benefit is offset by aft losses')
    axes[1].legend(fontsize=8,loc='lower right')
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle('S702 · Re 700,000 · M 0.07\nRight: MultiFoil Ncrit 6.25, original 260/180 panels; element contributions are model outputs.',fontsize=11)
    for ext in ('png','svg','pdf'):fig.savefig(out/('physics_evidence.'+ext),dpi=180)
    plt.close(fig)
    json_write(out/'summary.json',dict(bucket_offsets=offsets,mses_flap10_minima=minima,
        new_flow_solves=0,independent_experimental_element_drag_available=False,
        selected_aft_separation_reaches_te=False,consumed_inputs=len(inputs.consumed)))
    json_write(out/'input_manifest.json',dict(files=inputs.consumed,script_sha256=sha(Path(__file__)),
        helper_sha256=sha(ROOT/'validation/report_s702_ncrit_flap_error.py'),
        scope='Existing digitized curves and archived selected states; no new experimental measurements'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/s702_physics_diagnosis_20261001')
    args=parser.parse_args();report(args.out)
