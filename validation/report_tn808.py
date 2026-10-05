"""Plot and report the qualified TN808 challenge, retaining every study attempt."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tn808_geometry import geometry,geometry_audit,MAIN_TABLE,FLAP_TABLE,SLOT_TABLE

ROOT=Path(__file__).resolve().parent.parent
DATA=ROOT/'validation/data/naca23012_tn808'
LEGACY=ROOT/'validation/audit/before_paneling_guard'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def load_study(path):
    path=path.resolve()
    spec=json.loads((path/'study.json').read_text())
    if 'retry_of' in spec:
        original=ROOT/spec['retry_of']['path']
        if sha(original)!=spec['retry_of']['sha256']:
            raise ValueError('Parent failure record changed after retry')
    provenance={}
    for name,digest in spec['sha256'].items():
        current=ROOT/name
        archived=LEGACY/name
        if current.is_file() and sha(current)==digest:provenance[name]=str(current.relative_to(ROOT))
        elif archived.is_file() and sha(archived)==digest:provenance[name]=str(archived.relative_to(ROOT))
        else:raise ValueError(f'No matching current or archived source for {path.name}: {name}')
    records=[]
    for job in spec['jobs']:
        run=path/f"d{job['deflection']:02d}__{job['setting']}"
        record=json.loads((run/'result.json').read_text())
        record['study']=path.name
        record['path']=str(run.relative_to(ROOT))
        if record['feasible']:
            if not np.isfinite([record['CL'],record['CD']]).all() or abs(record['CL']-1.5)>2e-5:
                raise ValueError('Invalid fixed-lift record')
            if record['local_supersonic']:raise ValueError('Inadmissible compressibility flag')
        records.append(record)
    return records,provenance


def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',action='append',type=Path)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/tn808_comparison')
    args=parser.parse_args()
    paths=args.study or [ROOT/'validation/results'/name for name in
                         ('tn808_base','tn808_guarded','tn808_sensitivity',
                          'tn808_retry_base30','tn808_retry_coarse10','tn808_retry_base30_warm',
                          'tn808_retry_ncrit6_20','tn808_retry_datum30')
                         if (ROOT/'validation/results'/name).exists()]
    attempts=[];provenance={}
    for path in paths:
        records,sources=load_study(path)
        attempts.extend(records);provenance[path.name]=sources
    with (DATA/'fixed_lift_reference.csv').open() as f:
        reference={int(row['deflection_deg']):{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)}
    # Each successful setting/angle is retained separately. If a retry duplicates
    # a successful pair it must be reviewed explicitly, not overwrite it here.
    by_setting={}
    for row in attempts:
        if not row['feasible']:continue
        values=by_setting.setdefault(row['setting'],{})
        if row['deflection'] in values:raise ValueError('Duplicate successful setting/angle; review retry explicitly')
        values[row['deflection']]=row
    primary_setting='base'
    if len(by_setting.get('base',{}))<3 and len(by_setting.get('coarse',{}))==3:
        primary_setting='coarse'
    nominal=by_setting.get(primary_setting,{})
    comparisons=[];bl_diagnostics=[]
    for row in attempts:
        if not row['feasible']:continue
        metadata=json.loads((ROOT/row['path']/'flow/metadata.json').read_text())
        for element,sides in metadata['diagnostics'][-1]['boundary_layers'].items():
            for side,values in sides.items():
                bl_diagnostics.append(dict(study=row['study'],deflection_deg=row['deflection'],
                    setting=row['setting'],element=element,side=side,
                    separated_arclength_fraction=values['separated_arclength_fraction'],
                    max_Hk=values['max_Hk'],min_cf=values['min_cf']))
    for angle,r in nominal.items():
        ref=reference[angle]
        comparisons.append(dict(deflection_deg=angle,setting=primary_setting,CL=r['CL'],alpha_deg=r['alpha_deg'],
            model_CD=r['CD'],experiment_CD=ref['CD'],
            delta_drag_counts=(r['CD']-ref['CD'])*1e4,
            relative_difference_percent=(r['CD']/ref['CD']-1)*100,
            experimental_working_allowance=ref['total_cd_allowance']))
    pairs=[]
    for a,b in itertools.combinations((10,20,30),2):
        delta=reference[b]['CD']-reference[a]['CD']
        allowance=reference[a]['total_cd_allowance']+reference[b]['total_cd_allowance']
        sign=0 if abs(delta)<=allowance else int(np.sign(delta))
        for setting,values in by_setting.items():
            model_delta=(values[b]['CD']-values[a]['CD']) if a in values and b in values else None
            pairs.append(dict(a_deg=a,b_deg=b,setting=setting,experiment_delta_counts=delta*1e4,
                experimental_difference_allowance_counts=allowance*1e4,
                model_delta_counts=None if model_delta is None else model_delta*1e4,
                interpretation=('incomplete' if model_delta is None else 'reference_unresolved' if sign==0
                    else 'same_direction_under_working_allowance' if np.sign(model_delta)==sign else
                         'opposite_direction_under_working_allowance')))
    args.out.mkdir(parents=True,exist_ok=True)
    write_csv(args.out/'absolute_drag.csv',comparisons)
    write_csv(args.out/'pairwise_drag_changes.csv',pairs)
    write_csv(args.out/'boundary_layer_diagnostics.csv',bl_diagnostics)
    (args.out/'attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
    (args.out/'source_verification.json').write_text(json.dumps(provenance,indent=2)+'\n')
    source_ordinates=dict(units='original wing chord',source='NACA TN808 Table I',
        main_columns=['x','upper_y','lower_y'],main=MAIN_TABLE.tolist(),
        flap_columns=['x','upper_y','lower_y'],flap=FLAP_TABLE.tolist(),
        slot_columns=['x','y'],slot=SLOT_TABLE.tolist(),audit=geometry_audit())
    (DATA/'ordinates.json').write_text(json.dumps(source_ordinates,indent=2)+'\n')
    for angle in (10,20,30):
        raw,audit=geometry(angle)
        bundle=dict(deflection=angle,audit=audit,elements={name:xy.tolist() for name,xy in raw})
        (DATA/f'geometry_d{angle}.json').write_text(json.dumps(bundle,indent=2)+'\n')
    summary=dict(attempted=len(attempts),converged=sum(r['feasible'] for r in attempts),
        requested_flow_evaluations=sum(r['flow_evaluations'] for r in attempts),
        independent_optimization_accuracy_established=False,primary_setting=primary_setting,
        comparison_kind='qualified reconstruction vs digitized experiment',
        note='Working allowances and observed sensitivity ranges are not total uncertainty bounds. '
             'No original MSES calculation or contour-shape gradient validation.',
        absolute_drag=comparisons,pairwise_drag_changes=pairs)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig=plt.figure(figsize=(12,8),layout='constrained')
    grid=fig.add_gridspec(2,2,height_ratios=[.8,1.5])
    axg=fig.add_subplot(grid[0,:]);axc=fig.add_subplot(grid[1,0]);axd=fig.add_subplot(grid[1,1])
    colors={10:'#377eb8',20:'#ff7f00',30:'#4daf4a'}
    for angle in (10,20,30):
        raw,_=geometry(angle)
        if angle==10:axg.plot(*raw[0][1].T,color='#333333',lw=1.3,label='Main with tabulated slot')
        axg.plot(*raw[1][1].T,color=colors[angle],lw=1.4,label=f'Flap {angle} deg')
    axg.set_aspect('equal');axg.set_xlabel('x / original wing chord');axg.set_ylabel('y / c')
    axg.legend(ncol=4,loc='lower left',fontsize=9)
    axg.set_title('Reconstructed slotted flap 1-b; vertical gap 0.02 c')
    angles=np.array([10,20,30]);ref=np.array([reference[d]['CD'] for d in angles])
    axc.errorbar(angles,ref,yerr=.0032,fmt='ko',capsize=4,label='Experiment + working allowance',zorder=5)
    axd.errorbar(angles,(ref-ref[1])*1e4,yerr=[64,0,64],fmt='ko',capsize=4,zorder=5)
    styles={'base':('tab:blue','o-','MultiFoil base'),
            'coarse':('tab:orange','s--','Coarse panels'),
            'ncrit6':('tab:green','^:','Ncrit 6'),
            'table_axis':('tab:purple','d-.','Alternate flap datum')}
    for setting,values in by_setting.items():
        color,fmt,label=styles.get(setting,('gray','x-',setting))
        cd=np.array([values[d]['CD'] if d in values else np.nan for d in angles])
        axc.plot(angles,cd,fmt,color=color,ms=4,label=label)
        axd.plot(angles,(cd-cd[1])*1e4,fmt,color=color,ms=4)
    axc.set_ylabel('Section profile drag coefficient CD');axc.set_ylim(bottom=0)
    axd.axhline(0,color='gray',lw=.8);axd.set_ylabel('CD change from 20 deg (drag counts)')
    axc.set_title('Absolute drag at CL = 1.5');axd.set_title('Change between configurations')
    for ax in (axc,axd):
        ax.set_xlabel('Flap deflection (deg)');ax.set_xticks(angles);ax.grid(alpha=.2)
    axc.legend(loc='lower left',fontsize=8)
    fig.suptitle('MultiFoil vs experiment: NACA TN808 fixed-lift challenge',fontsize=15)
    fig.get_layout_engine().set(rect=(0,.075,1,.925))
    fig.text(.5,.032,'Re = 2.19 million; Mach 0.105 assumed. Figure8 digitized; geometry and transition qualified.\n'
             'Reference bars: comparison allowance, not total uncertainty. Lines connect tested flap settings.',
             ha='center',fontsize=9)
    fig.savefig(args.out/'tn808_comparison.png',dpi=180)
    fig.savefig(args.out/'tn808_comparison.pdf')
    plt.close(fig)

    lines=['# TN808 executed comparison','',
        'MultiFoil now runs a qualified experimental flap-angle challenge. The results do not establish '
        'drag-optimization accuracy. Geometry reconstruction, transition assumptions and unquantified '
        'fitting drag remain qualifications. No source convention or physical coefficient was fitted to drag.','',
        '[Primary source: Lowry (1941), NACA TN808](https://ntrs.nasa.gov/api/citations/19930085051/downloads/19930085051.pdf). '
        'Figure8, flap1-b, gap 0.02c, physical Re 2.19 million; CL 1.5.','',
        f'Primary comparison setting: **{primary_setting}**. One setting is used for this whole table; '
        'failed points are not replaced individually with another resolution. All attempts and successful '
        'settings remain in the files and plot.','',
        '| Flap | MultiFoil CD | Experiment CD | Difference |',
        '|---|---:|---:|---:|']
    for row in sorted(comparisons,key=lambda r:r['deflection_deg']):
        lines.append(f"| {row['deflection_deg']} deg | {row['model_CD']:.5f} | {row['experiment_CD']:.5f} | "
                     f"{row['relative_difference_percent']:+.1f}% |")
    lines+=['','These are discrepancies for the reconstructed case, not a proven decomposition of solver '
            'error versus geometry, tunnel and fitting effects. The reference working allowance is ±0.0032 '
            'per point. It is not a measured total-uncertainty bound.','',
            '## Changes between geometries','',
            '| Pair (B minus A) | Setting | Experiment (counts) | MultiFoil (counts) | Interpretation |',
            '|---|---|---:|---:|---|']
    for row in pairs:
        value='missing' if row['model_delta_counts'] is None else f"{row['model_delta_counts']:+.1f}"
        lines.append(f"| {row['b_deg']} - {row['a_deg']} deg | {row['setting']} | "
                     f"{row['experiment_delta_counts']:+.1f} | {value} | {row['interpretation'].replace('_',' ')} |")
    lines+=['','The reference difference allowance is ±64 counts. A shared fitting-bias cancellation is '
        'not assumed. Agreement in a resolved direction would support only coarse screening within the tested '
        'conditions; it does not establish gradient magnitude or a continuous optimum.','',
        '## Flow-model limitations exposed by the challenge','',
        'The converged base 10/20 cases predict separation over about 20/21% of the main lower-side '
        'arc length around the slot entrance, with maximum Hk about 79/33. These are model diagnostics, '
        'not experimental boundary-layer measurements. They suggest that the sharp cove strongly '
        'strains the thin-layer approximation; they do not uniquely identify the cause of the drag discrepancy. '
        'Full per-side values are in boundary_layer_diagnostics.csv.','',
        'The initial coarse10 search also produced two converged values at incidences differing by '
        'only 1.16e-7 deg, with a CD change of about 1.01 drag counts after different continuation histories. '
        'The retained continuation_sensitivity.json records the two rows. This is numerical path '
        'sensitivity, not proof of physical hysteresis, and matters when optimizing sub-count changes.','',
        '## Execution and reproducibility','',
        f"{summary['attempted']} attempts; {summary['converged']} fixed-lift solutions. "
        f"{summary['requested_flow_evaluations']} requested flow evaluations archived, with internal continuation diagnostics. "
        'The original three attempts failed the geometry check before any flow solve and are retained.','',
        'The original cubic inverse-spacing fold was corrected by a conditional monotone inverse. '
        'All 46 tests passed, including the new regression and existing aerodynamic checks. '
        'Source files matching the pre-fix studies are retained in the audit overlay; '
        'source_verification.json records the actual matching path for each study hash.','',
        'See [TN808_BENCHMARK.md](../../TN808_BENCHMARK.md) for the reconstruction, '
        'experimental allowances and commands. The broad and narrower incidence procedures are distinct '
        'and must not be interpreted as a proof of branch uniqueness.','',
        '![TN808 comparison](tn808_comparison.png)','']
    (args.out/'TN808_RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
