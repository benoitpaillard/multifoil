#!/usr/bin/env python3
"""Compare frozen calibration with Reynolds, mesh, Mach and pressure checks."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from report_s702_ncrit_sweep import (ROOT, RESULTS, DATA, read, write, sha,
    load_study, csv_write, flow_directory, score, pairwise)
from report_s702_discrepancy import numerical_resolution, final_flow


def pressure_check(folder, out):
    _, points, details = load_study(folder)
    if len(details) != 1 or points[0]['flap'] != 0 or points[0]['target_cl'] != 1.04:
        raise ValueError('Pressure check must be flap 0, CL=1.04')
    current = details[0]
    if not current['feasible']:
        return dict(complete=False, error=current['error'])
    ncrit = points[0]['ncrit']
    sources = []
    for n in (9., 6.):
        target = RESULTS/f's702_pressure_check/ncrit{n:g}'
        result = read(target/'result.json')
        sources.append((f'MultiFoil N={n:g}', target/'flow', result, ['#5684b4', '#a0a7ad'][int(n==6.)]))
    sources.append((f'MultiFoil N={ncrit:g} selected', flow_directory(current), current, '#c76b28'))
    errors = []
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.8), layout='constrained')
    for ax, (element, side) in zip(axes.flat, [('fore','upper'),('aft','upper'),('fore','lower'),('aft','lower')]):
        label = element+'_'+side
        exp = np.genfromtxt(DATA/f'pressure_cl1.04/experiment_{label}.csv', delimiter=',', names=True)
        if side == 'lower':
            exp = exp[1:-1]  # Upper endpoint joins are not independent lower taps.
        ax.plot(exp['x_over_cref'], exp['Cp'], 'o', ms=3.5, mfc='white', mec='#222', label='Experiment', zorder=5)
        curves = []
        mses = np.genfromtxt(DATA/f'pressure_cl1.04/mses_{label}.csv', delimiter=',', names=True)
        curves.append(('Published MSES N=9', mses['x_over_cref'], mses['Cp'], '#555', '--'))
        for name, flow, result, color in sources:
            polar = np.genfromtxt(flow/'polar.csv', delimiter=',', names=True)
            index = np.flatnonzero((abs(polar['CL']-result['CL']) < 1e-10)&(abs(polar['CD']-result['CD']) < 1e-12))[-1]
            surface = np.genfromtxt(flow/f'surface_{index:03d}_{element}.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
            le = np.argmin(surface['x_over_cref'])
            part = surface[:le+1] if side == 'upper' else surface[le:]
            curves.append((name, part['x_over_cref'], part['Cp'], color, '-'))
        for name, x, cp, color, style in curves:
            order = np.argsort(x)
            x, cp = x[order], cp[order]
            x, unique = np.unique(x, return_index=True)
            cp = cp[unique]
            mask = (exp['x_over_cref'] >= x.min())&(exp['x_over_cref'] <= x.max())
            points = exp[mask]
            residual = np.interp(points['x_over_cref'], x, cp)-points['Cp']
            rms = np.sqrt(trapezoid(residual**2, points['x_over_cref'])/np.ptp(points['x_over_cref']))
            errors.append(dict(source=name, element=element, side=side, samples=len(points),
                               excluded_outside_overlap=int((~mask).sum()), rms_chord_weighted=float(rms)))
            ax.plot(x, cp, color=color, ls=style, lw=1.5, label=name)
        ax.invert_yaxis()
        ax.set(title=element.title()+' '+side, xlabel='x / assembly chord', ylabel='$C_p$')
        ax.grid(alpha=.18)
    axes[0,0].legend(fontsize=7)
    fig.suptitle('Pressure check excluded from calibration · flap 0° · $C_L=1.04$', fontsize=13)
    for ext in ('png','svg'):
        fig.savefig(out/f'ncrit_pressure.{ext}', dpi=180)
    plt.close(fig)
    csv_write(out/'pressure_errors.csv', errors)
    return dict(complete=True, ncrit=ncrit, errors=errors)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis', type=Path, default=RESULTS/'s702_ncrit_analysis')
    ap.add_argument('--holdout', type=Path, nargs='+', required=True)
    ap.add_argument('--fine', type=Path, nargs='+', required=True)
    ap.add_argument('--mach0', type=Path, required=True)
    ap.add_argument('--pressure', type=Path, required=True)
    ap.add_argument('--control', type=Path)
    ap.add_argument('--retry', type=Path, nargs='*', default=[])
    args = ap.parse_args()
    out = args.analysis.resolve()
    frozen = read(out/'frozen_selection.json')
    ncrit = frozen['selected_ncrit']
    for filename, field in [('points.csv', 'calibration_points_sha256'),
                            ('summary.json', 'calibration_summary_sha256')]:
        if sha(out/filename) != frozen[field]:
            raise ValueError(f'Calibration changed after selection was frozen: {filename}')
    if sha(DATA/'re700k_holdout/source.json') != frozen['holdout_reference_sha256']:
        raise ValueError('Reserved reference changed after selection was frozen')
    if frozen['holdout_used_in_selection'] or read(out/'summary.json')['best']['ncrit'] != ncrit:
        raise ValueError('Frozen selection does not match the calibration-only minimum')
    calibration = list(__import__('csv').DictReader((out/'points.csv').open()))
    lookup = {(float(r['ncrit']), int(r['flap']), float(r['target_cl'])):r for r in calibration}
    summary = dict(frozen_selection=frozen, source_manifests={})
    controls = {}
    if args.control:
        folder = args.control.resolve()
        _, control_rows, _ = load_study(folder)
        for row in control_rows:
            key = (row['ncrit'], row['flap'], row['target_cl'])
            old = lookup[key]
            row['protocol_minus_legacy_counts'] = (row['CD']-float(old['CD']))*1e4 if row['feasible'] else None
            controls[key] = row
        summary['control'] = control_rows
        csv_write(out/'continuation_control.csv', control_rows)
        summary['source_manifests'][str(folder.relative_to(ROOT))] = sha(folder/'study.json')
    holdout = []
    interrupted_holdout = []
    for folder in args.holdout:
        folder = folder.resolve()
        spec, rows, _ = load_study(folder, DATA/'re700k_holdout')
        if spec['ncrit_values'] not in ([9.], [ncrit]):
            raise ValueError('Reserved validation contains an unfrozen candidate')
        holdout += rows
        interrupted_holdout += [entry for path in folder.glob('resume_*.json')
                                for entry in read(path)['interrupted_attempts']]
        summary['source_manifests'][str(folder.relative_to(ROOT))] = sha(folder/'study.json')
    required = [(f, cl) for f in [0, 5, 10] for cl in [.5,.8,1.2]]
    scores = [dict(ncrit=n, **score([r for r in holdout if r['ncrit']==n], required)) for n in (9., ncrit)]
    mses = [dict(r, CD=r['CD_mses_n9'], feasible=True) for r in holdout if r['ncrit']==9.]
    summary['holdout'] = dict(scores=scores, published_mses_n9=score(mses, required),
        interrupted_attempts=interrupted_holdout,
        pairs=[dict(ncrit=n, **p) for n in (9.,ncrit)
               for p in pairwise([r for r in holdout if r['ncrit']==n], [0,5,10], [.5,.8,1.2])])
    common = [key for key in required if all(any(r['feasible'] and r['ncrit']==n
              and (r['flap'],r['target_cl'])==key for r in holdout) for n in (9.,ncrit))]
    if common and len(common) != len(required):
        summary['holdout']['common_subset'] = dict(points=common,
            excluded_points=[key for key in required if key not in common],
            qualification='Diagnostic on the identical converged subset only; not a complete reserved-validation score',
            scores=[dict(ncrit=n, **score([r for r in holdout if r['ncrit']==n], common)) for n in (9.,ncrit)],
            published_mses_n9=score(mses,common))
    csv_write(out/'holdout_points.csv', holdout)
    after = [dict(r, local_retry_of='') for r in holdout]
    retry_rows, retried_keys = [], set()
    for folder in args.retry:
        folder = folder.resolve()
        spec, rows, _ = load_study(folder, DATA/'re700k_holdout')
        original = spec.get('local_retry_of')
        if len(rows)!=1 or not original or rows[0]['ncrit']!=ncrit:
            raise ValueError('Unexpected local retry selection')
        source_spec = read((ROOT/original).parents[2]/'study.json')
        for dep in ('multifoil.py','viscous.py','paneling.py','optimization.py','vendor/mfoil.py'):
            if spec['sha256'][dep] != source_spec['sha256'][dep]:
                raise ValueError('Local retry changed the numerical core')
        index = next((i for i,r in enumerate(holdout) if r['result_path']==original), None)
        key = (rows[0]['ncrit'],rows[0]['flap'],rows[0]['target_cl'])
        if index is None or holdout[index]['feasible']:
            raise ValueError('Local retry must address one original failed reserved target')
        if (after[index]['ncrit'],after[index]['flap'],after[index]['target_cl']) != key:
            raise ValueError('Local retry target changed')
        replacement = dict(rows[0],local_retry_of=original)
        used = not after[index]['feasible']
        if used:
            after[index] = replacement
            for previous in retry_rows:
                if previous['local_retry_of']==original:
                    previous['selected_for_score'] = False
        retry_rows.append(dict(replacement,
            initialization='warm replay' if spec.get('warm_seed_alphas') else 'cold local seed',
            selected_for_score=used and replacement['feasible']))
        retried_keys.add(key)
        summary['source_manifests'][str(folder.relative_to(ROOT))] = sha(folder/'study.json')
    if retry_rows:
        summary['holdout']['after_local_retries'] = dict(
            qualification='Separate numerical follow-up with frozen Ncrit and physics; all original failures retained',
            selection_rule='First successful initialization in the supplied order; reproduction uses cold then warm; no drag-based selection',
            retries=retry_rows,
            scores=[dict(ncrit=n,**score([r for r in after if r['ncrit']==n],required)) for n in (9.,ncrit)])
        csv_write(out/'local_retry_points.csv', retry_rows)
        csv_write(out/'holdout_after_retries.csv', after)
    fine = []
    for folder in args.fine:
        folder = folder.resolve()
        spec, rows, details = load_study(folder)
        if spec['panels'] != [340,240]:
            raise ValueError('Unexpected refinement panel count')
        summary['source_manifests'][str(folder.relative_to(ROOT))] = sha(folder/'study.json')
        for row, detail in zip(rows, details):
            key = (row['ncrit'], row['flap'], row['target_cl'])
            base = controls.get(key, lookup[key])
            row['base_CD'] = float(base['CD']) if base['CD'] is not None else None
            row['base_reference'] = 'matched-protocol control' if key in controls else 'calibration sweep'
            row['fine_minus_base_counts'] = (row['CD']-row['base_CD'])*1e4 if row['feasible'] and row['base_CD'] is not None else None
            row['resolution'] = numerical_resolution((ROOT/row['result_path']).parent, detail) if row['feasible'] else []
            fine.append(row)
    summary['fine'] = fine
    csv_write(out/'fine_points.csv', [{k:v for k,v in r.items() if k!='resolution'} for r in fine])
    folder = args.mach0.resolve()
    _, mach, _ = load_study(folder)
    summary['source_manifests'][str(folder.relative_to(ROOT))] = sha(folder/'study.json')
    for row in mach:
        base = lookup[(row['ncrit'], row['flap'], row['target_cl'])]
        row['mach0_minus_mach01_counts'] = (row['CD']-float(base['CD']))*1e4 if row['feasible'] else None
    summary['mach_zero'] = mach
    csv_write(out/'mach_zero_points.csv', mach)
    summary['pressure'] = pressure_check(args.pressure.resolve(), out)
    summary['source_manifests'][str(args.pressure.resolve().relative_to(ROOT))] = sha(args.pressure.resolve()/'study.json')
    summary['holdout_reference_sha256'] = sha(DATA/'re700k_holdout/source.json')
    write(out/'transfer_summary.json', summary)
    fig, axes = plt.subplots(1,3,figsize=(11.5,4.3),layout='constrained')
    for ax, cl in zip(axes,[.5,.8,1.2]):
        x = np.arange(3)
        sources = [('MultiFoil N=9',9.,'#5684b4'),(f'MultiFoil N={ncrit:g}',ncrit,'#c76b28')]
        for i,(name,n,color) in enumerate(sources):
            pts = [next(r for r in holdout if r['ncrit']==n and r['flap']==f and r['target_cl']==cl) for f in (0,5,10)]
            ax.bar(x+(i-.5)*.26,[r['CD']*1e4 if r['feasible'] else np.nan for r in pts],.25,color=color,label=name)
            for position, point in zip(x+(i-.5)*.26, pts):
                if not point['feasible']:
                    retried = next(r for r in after if r['ncrit']==n and r['flap']==point['flap'] and r['target_cl']==cl)
                    if retried['feasible']:
                        ax.bar(position,retried['CD']*1e4,.25,color=color,alpha=.5,hatch='///',edgecolor='#333')
                    else:
                        ax.annotate('no convergence', xy=(position,.04),
                                    xycoords=('data','axes fraction'), color=color,
                                    rotation=90, ha='center', va='bottom', fontsize=7)
        ax.plot(x,[r['CD_experiment']*1e4 for r in pts],'ko',ms=5,label='Experiment')
        ax.plot(x,[r['CD_mses_n9']*1e4 for r in pts],'D',color='#777',ms=4,label='Published MSES N=9')
        ax.set(xticks=x,xticklabels=['0°','5°','10°'],xlabel='Flap',ylabel='$C_D$ (counts)',title=f'$C_L={cl:g}$')
        ax.grid(axis='y',alpha=.18)
    handles, labels = axes[0].get_legend_handles_labels()
    if any(r['feasible'] for r in retry_rows):
        from matplotlib.patches import Patch
        handles.append(Patch(facecolor='#c76b28',alpha=.5,hatch='///',edgecolor='#333'))
        labels.append('Local initialization retry')
    axes[0].legend(handles,labels,fontsize=7)
    fig.suptitle('Reserved Reynolds-number check · Re = 700,000 · M = 0.07',fontsize=13)
    for ext in ('png','svg'):
        fig.savefig(out/f'ncrit_holdout.{ext}',dpi=180)
    plt.close(fig)
    print(json.dumps(dict(ncrit=ncrit, holdout=scores, fine_points=len(fine), mach_zero_points=len(mach)),indent=2))


if __name__ == '__main__':
    main()
