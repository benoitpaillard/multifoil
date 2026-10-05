#!/usr/bin/env python3
"""Plot the audited S702 CL=1.2 flap-drag comparison without running solves."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'validation/results/s702_panel_convergence_analysis'
OUT = ROOT/'validation/results/s702_cl12_drag_change_20260924'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    receipt = json.loads((SOURCE/'verification.json').read_text())
    for name in ('flap_deltas.csv', 'forces.csv', 'summary.json'):
        if digest(SOURCE/name) != receipt['outputs_sha256'][name]:
            raise ValueError('Audited input changed: '+name)
    with (SOURCE/'flap_deltas.csv').open(newline='') as stream:
        row = next(r for r in csv.DictReader(stream) if float(r['target_cl']) == 1.2)
    with (SOURCE/'forces.csv').open(newline='') as stream:
        forces = {int(r['flap']):float(r['CD']) for r in csv.DictReader(stream)
                  if float(r['target_cl']) == 1.2 and r['grid'] == 'finer'}
    model = float(row['finer_flap_delta_CD_counts'])
    experiment = float(row['experiment_flap_delta_counts'])
    gap = model-experiment
    assert abs(model-(forces[10]-forces[5])*1e4) < 1e-10
    assert abs(gap-float(row['finer_flap_delta_error_counts'])) < 1e-10
    summary = json.loads((SOURCE/'summary.json').read_text())
    assert summary['full_state_audit_pass'] and not summary['decision']['numerical_screen_pass']

    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':12,
                         'svg.fonttype':'none', 'axes.unicode_minus':True})
    fig, ax = plt.subplots(figsize=(10.5,5.8), facecolor='white')
    fig.subplots_adjust(left=.23, right=.97, top=.70, bottom=.28)
    fig.text(.04,.92,'Increasing flap deflection: prediction vs measurement',
             fontsize=19, weight='bold', color='#142c3d')
    fig.text(.04,.85,'S702  •  CL = 1.2  •  flap 5° → 10°  •  Re = 700,000  •  Mach = 0.07',
             fontsize=12, color='#4a5967')
    colors = ('#316ea5','#168778')
    ax.barh([1,0],[model,experiment],height=.43,color=colors,zorder=3)
    ax.axvline(0,color='#3f4e5a',lw=1.2,zorder=2)
    ax.set_yticks([1,0],['MultiFoil\n440/320 panels','Experiment\ndigitized data'])
    ax.tick_params(axis='y',length=0,pad=12)
    ax.set_xlim(-12.5,8.5); ax.set_ylim(-.5,1.65)
    ax.set_xticks([-10,-5,0,5])
    ax.set_xlabel('Drag change, CD(10°) − CD(5°)  [counts]',labelpad=12)
    ax.grid(axis='x',color='#e2e8ed',lw=.8,zorder=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.text(model+.32,1,f'+{model:.2f}\nmore drag',va='center',ha='left',
            color=colors[0],fontsize=13,weight='bold')
    ax.text(experiment-.32,0,f'{experiment:.2f}\nless drag',va='center',ha='right',
            color=colors[1],fontsize=13,weight='bold')
    ax.annotate('',xy=(experiment,1.46),xytext=(model,1.46),
                arrowprops=dict(arrowstyle='<->',color='#667582',lw=1.1))
    ax.text((experiment+model)/2,1.52,f'Difference: {gap:.2f} counts',
            ha='center',va='bottom',fontsize=11,color='#465664')
    fig.text(.04,.13,'Positive = drag increase. Negative = drag reduction.  1 count = 0.0001 CD.',
             fontsize=10.5,color='#4a5967')
    fig.text(.04,.075,'Model: Ncrit = 6.25, finest tested grid. Numerical convergence is not established.',
             fontsize=10,color='#4a5967')
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ('png','svg','pdf'):
        fig.savefig(OUT/('cl12_drag_change.'+ext),dpi=200,facecolor='white')
    plt.close(fig)
    (OUT/'data.json').write_text(json.dumps(dict(
        target_cl=1.2,flap_change_degrees=[5,10],reynolds=700000,mach=.07,ncrit_model=6.25,
        model_panels=[440,320],model_drag_change_counts=model,
        experiment_drag_change_counts=experiment,model_minus_experiment_counts=gap,
        drag_count_in_CD=1e-4,experimental_data='Existing audited digitization',
        qualification='No uncertainty interval is inferred. Model grid-convergence screen failed.',
        reproduce='python3 validation/plot_s702_cl12_drag_change.py',
        input_sha256={str(p.relative_to(ROOT)):digest(p) for p in
                      [SOURCE/name for name in ('flap_deltas.csv','forces.csv','summary.json','verification.json')]+[Path(__file__).resolve()]}
    ),indent=2)+'\n')
    print(json.dumps(dict(model=model,experiment=experiment,gap=gap,output=str(OUT)),indent=2))


if __name__ == '__main__':
    main()
