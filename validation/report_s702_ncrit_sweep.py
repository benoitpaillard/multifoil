#!/usr/bin/env python3
"""Report a completed transition sweep, preserving failed and reused cases."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from s702_ncrit_sweep import ROOT, read, write, sha, tag
from s702_sweep_metrics import score, select, leave_one_flap_out, pairwise

RESULTS = ROOT/'validation/results'
DATA = ROOT/'validation/data/s702'
FLAPS, TARGETS = [-2, 0, 5, 10], [.5, .8, 1.2]


def flow_directory(result):
    # The first coarse-run script labeled the old N=9 export as "flow".
    # Its authoritative source_path and CL identify the original per-CL export.
    # New runs write the correct path; retain original coarse bytes/provenance.
    if result.get('reused') and result['ncrit'] == 9.:
        return (ROOT/result['source_path']).parent/f"cl{result['target_cl']:.4f}"
    return ROOT/result['flow_directory']


def csv_write(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reference(folder, flap, cl, kind):
    values = np.genfromtxt(folder/f'd{flap:+d}/{kind}_drag.csv', delimiter=',', names=True)
    order = np.argsort(values['CL'])
    x, y = values['CL'][order], values['CD'][order]
    if not x.min() <= cl <= x.max() or not np.isfinite(y).all() or y.min() <= 0:
        raise ValueError('Reference would require extrapolation or has invalid drag')
    return float(np.interp(cl, x, y))


def load_study(folder, reference_folder=DATA):
    folder = Path(folder).resolve()
    spec = read(folder/'study.json')
    for path, expected in spec['sha256'].items():
        if sha(ROOT/path) != expected:
            snapshot = folder/'source_snapshot'/path
            if not snapshot.exists() or sha(snapshot) != expected:
                raise ValueError('Source differs from executed study: '+path)
    expected_conditions = (700000., .07) if reference_folder.name == 're700k_holdout' else (1e6, .1)
    matched_reference = (spec['reynolds'], spec['mach']) == expected_conditions
    if not matched_reference and spec['mach'] != 0.:
        raise ValueError('Reference conditions are not matched')
    rows, details = [], []
    for job in spec['jobs']:
        for cl in spec['targets']:
            path = folder/job['name']/f'cl{tag(cl)}'/'result.json'
            result = read(path)
            if result['feasible']:
                if abs(result['CL']-cl) > spec['cl_tolerance'] or result['CD'] <= 0:
                    raise ValueError('Invalid converged target')
                if result.get('reused'):
                    if sha(ROOT/result['source_path']) != result['source_sha256']:
                        raise ValueError('Reused source changed')
                else:
                    selected = result['records'][result['selected_flow_index']]
                    if (not selected['admissible'] or selected['residual'] > spec['residual_tolerance'] or
                        abs(selected['CD']-result['CD']) > 1e-12 or abs(selected['CL']-result['CL']) > 1e-12):
                        raise ValueError('Selected state does not match converged target')
            exp = reference(reference_folder, job['flap'], cl, 'experiment')
            mses = reference(reference_folder, job['flap'], cl, 'mses')
            rows.append(dict(ncrit=job['ncrit'], flap=job['flap'], target_cl=cl,
                feasible=result['feasible'], reused=result['reused'],
                CL=result['CL'], CD=result['CD'], alpha_deg=result['alpha_deg'],
                CD_experiment=exp, CD_mses_n9=mses,
                reference_conditions_matched=matched_reference,
                error_counts=(result['CD']-exp)*1e4 if result['feasible'] and matched_reference else None,
                error_percent=(result['CD']/exp-1)*100 if result['feasible'] and matched_reference else None,
                result_path=str(path.relative_to(ROOT))))
            details.append(result)
    return spec, rows, details


def analyze(rows, candidates, flaps=FLAPS, targets=TARGETS):
    required = [(f, cl) for f in flaps for cl in targets]
    scores, pairs = [], []
    for n in candidates:
        subset = [r for r in rows if r['ncrit'] == n]
        result = score(subset, required)
        changes = pairwise(subset, flaps, targets)
        result.update(pairwise_delta_rmse_counts=float(np.sqrt(np.mean([r['delta_error_counts']**2 for r in changes]))) if result['complete'] else None,
                      matching_pair_signs=sum(r['same_sign'] for r in changes) if result['complete'] else None,
                      pair_comparisons=len(changes))
        scores.append(dict(ncrit=n, **result))
        pairs.extend(dict(ncrit=n, **r) for r in changes)
    return scores, pairs, select(rows, required, candidates)


def figures(rows, scores, pairs, best, out):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.7), layout='constrained')
    complete = [s for s in scores if s['complete']]
    ns = [s['ncrit'] for s in complete]
    axes[0].plot(ns, [s['mape_percent'] for s in complete], 'o-', color='#142d47', lw=2, label='All 12 points')
    for flap, color in zip(FLAPS, ('#8e6aa6', '#409783', '#df8a30', '#497ab4')):
        values = [score([r for r in rows if r['ncrit'] == n], [(flap, cl) for cl in TARGETS]) for n in ns]
        axes[0].plot(ns, [v['mape_percent'] for v in values], '.--', lw=1, alpha=.8,
                     color=color, label=f'Flap {flap}°')
    axes[0].axvline(best['ncrit'], color='#142d47', ls=':', lw=1,
                    label=f"Selected N={best['ncrit']:g}")
    incomplete = [s for s in scores if not s['complete']]
    if incomplete:
        axes[0].text(.02, .98, 'Unscored: '+', '.join(f"N={s['ncrit']:g} ({s['evaluated']}/{s['required']})" for s in incomplete),
                     transform=axes[0].transAxes, va='top', fontsize=8,
                     bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
        axes[0].set_xlim(min(s['ncrit'] for s in scores)-.15, max(s['ncrit'] for s in scores)+.15)
    axes[0].set(xlabel='$N_{crit}$, common to both elements', ylabel='Mean absolute relative drag error (%)',
                title='Calibration at matched lift')
    axes[0].legend(fontsize=8)
    colors = ['#409783', '#df8a30', '#497ab4']
    sampled_ns = sorted({r['ncrit'] for r in rows})
    for cl, color in zip(TARGETS, colors):
        points = sorted([p for p in pairs if p['first_flap'] == 5 and p['second_flap'] == 10
                         and p['target_cl'] == cl], key=lambda p:p['ncrit'])
        deltas = {p['ncrit']: p['predicted_delta_counts'] for p in points}
        axes[1].plot(sampled_ns, [deltas.get(n, np.nan) for n in sampled_ns],
                     'o-', color=color, ms=4, label=f'MultiFoil, $C_L={cl:g}$')
        if points:
            axes[1].axhline(points[0]['measured_delta_counts'], color=color, ls='--', lw=1)
    axes[1].axhline(0, color='#aaa', lw=.7)
    axes[1].plot([], [], color='#555', ls='--', label='Experiment, same color')
    axes[1].set(xlabel='$N_{crit}$', ylabel='$C_D$(10°) − $C_D$(5°), drag counts', title='Does the flap ranking transfer?')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.18)
    fig.suptitle('S702 transition sweep · Re = 1 million · M = 0.1', fontsize=13)
    for ext in ('png', 'svg'):
        fig.savefig(out/f'ncrit_sweep.{ext}', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(3, 4, figsize=(12, 8.2), layout='constrained')
    for i, cl in enumerate(TARGETS):
        for j, flap in enumerate(FLAPS):
            ax = axes[i,j]
            pts = sorted([r for r in rows if r['flap'] == flap and r['target_cl'] == cl], key=lambda r:r['ncrit'])
            ax.plot([r['ncrit'] for r in pts], [r['CD']*1e4 if r['feasible'] else np.nan for r in pts],
                    'o-', ms=3, label='MultiFoil')
            ax.axhline(pts[0]['CD_experiment']*1e4, color='#222', ls='--', label='Experiment')
            ax.plot(9., pts[0]['CD_mses_n9']*1e4, 's', color='#df8a30', ms=5, label='Published MSES N=9')
            ax.axvline(best['ncrit'], color='#777', ls=':', lw=.8,
                       label=f"Selected N={best['ncrit']:g}")
            ax.set(title=f'Flap {flap}°, $C_L={cl:g}$', xlabel='$N_{crit}$', ylabel='$C_D$ (counts)')
            ax.grid(alpha=.18)
    axes[0,0].legend(fontsize=7)
    fig.suptitle('Every calibration point; gaps retain failed solves', fontsize=13)
    for ext in ('png', 'svg'):
        fig.savefig(out/f'ncrit_pointwise.{ext}', dpi=180)
    plt.close(fig)


def transition_figure(rows, out):
    transitions = []
    for row in rows:
        if row['target_cl'] != .8 or row['flap'] not in (5,10) or not row['feasible']:
            continue
        result = read(ROOT/row['result_path'])
        flow = flow_directory(result)
        polar = list(csv.DictReader((flow/'polar.csv').open()))
        indices = [i for i,p in enumerate(polar) if abs(float(p['CL'])-row['CL']) < 1e-10
                   and abs(float(p['CD'])-row['CD']) < 1e-12]
        if not indices:
            raise ValueError('Transition output cannot be tied to the selected force result')
        layers = read(flow/'metadata.json')['diagnostics'][indices[-1]]['boundary_layers']
        for element in ('fore','aft'):
            for side in ('upper','lower'):
                tr = layers[element][side]['transition']
                transitions.append(dict(ncrit=row['ncrit'], flap=row['flap'], target_cl=.8,
                    element=element, side=side, natural_transition_before_te=tr is not None,
                    x_over_local_chord=tr['x_over_local_chord'] if tr else None,
                    x_over_assembly_chord=tr['xy_over_cref'][0] if tr else None))
    csv_write(out/'transition_locations.csv', transitions)
    fig, axes = plt.subplots(2,2,figsize=(10,6.2),layout='constrained')
    for i,element in enumerate(('fore','aft')):
        for j,flap in enumerate((5,10)):
            ax = axes[i,j]
            for side,color in [('upper','#236da5'),('lower','#c48432')]:
                pts = sorted([r for r in transitions if r['element']==element and r['flap']==flap
                              and r['side']==side],key=lambda r:r['ncrit'])
                n = np.array([r['ncrit'] for r in pts])
                x = np.array([r['x_over_local_chord'] if r['natural_transition_before_te'] else 1. for r in pts])
                natural = np.array([r['natural_transition_before_te'] for r in pts])
                ax.plot(n,x,'-',color=color,label=side.title())
                ax.plot(n[natural],x[natural],'o',ms=4,color=color)
                ax.plot(n[~natural],x[~natural],'s',ms=4,mfc='white',mec=color)
            ax.set(title=f'{element.title()} element, flap {flap}°',xlabel='$N_{crit}$',
                   ylabel='Transition x / local chord',ylim=(0,1.055))
            ax.grid(alpha=.18)
    axes[0,0].plot([],[],'s',mfc='white',mec='#555',label='Laminar to trailing edge')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Modeled transition at matched $C_L=0.8$ · no measured transition positions fitted',fontsize=12)
    for ext in ('png','svg'):
        fig.savefig(out/f'ncrit_transition.{ext}',dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--coarse', type=Path, default=RESULTS/'s702_ncrit_coarse')
    ap.add_argument('--refined', type=Path, nargs='+')
    ap.add_argument('--out', type=Path, default=RESULTS/'s702_ncrit_analysis')
    args = ap.parse_args()
    args.coarse, args.out = args.coarse.resolve(), args.out.resolve()
    if args.refined:
        args.refined = [p.resolve() for p in args.refined]
    spec, coarse, detail = load_study(args.coarse)
    rows = list(coarse)
    for folder in args.refined or []:
        _, refined, extra = load_study(folder)
        if set(r['ncrit'] for r in refined) & set(r['ncrit'] for r in rows):
            raise ValueError('Refined and coarse grid values must be distinct')
        rows += refined
        detail += extra
    candidates = sorted({r['ncrit'] for r in rows})
    scores, pairs, best = analyze(rows, candidates)
    metric_optima = {key:min((s for s in scores if s['complete']), key=lambda s:(s[key],abs(s['ncrit']-9.)))
                     for key in ('mape_percent','rmse_counts','pairwise_delta_rmse_counts')}
    folds, cv = leave_one_flap_out(coarse, FLAPS, TARGETS, spec['ncrit_values'])
    args.out.mkdir(parents=True, exist_ok=True)
    csv_write(args.out/'points.csv', rows)
    csv_write(args.out/'scores.csv', scores)
    csv_write(args.out/'pairwise.csv', pairs)
    csv_write(args.out/'cross_validation.csv', folds)
    near = [s['ncrit'] for s in scores if s['complete'] and s['mape_percent'] <= best['mape_percent']+1.]
    interruptions = [entry for folder in [args.coarse]+(args.refined or [])
                     for path in folder.glob('resume_*.json')
                     for entry in read(path)['interrupted_attempts']]
    summary = dict(best=best, near_optimal_sampled_values=near,
        metric_optima=metric_optima,
        coarse_ncrit_values=spec['ncrit_values'],
        refined_ncrit_values=sorted(set(candidates)-set(spec['ncrit_values'])),
        near_optimal_definition='MAPE within 1 percentage point of best; engineering tolerance, not a confidence interval',
        coarse_grid_cross_validation=folds, cross_validated_score=cv,
        scored_points=len(rows), feasible_points=sum(r['feasible'] for r in rows),
        reused_points=sum(r['reused'] for r in rows),
        recorded_requested_evaluations=sum(len(r['records']) for r in detail if not r['reused']),
        failed_requested_evaluations=sum(not e['admissible'] for r in detail if not r['reused'] for e in r['records']),
        interrupted_attempts=interruptions,
        calibration='Equal-point MAPE at Re=1e6, M=.1; same Ncrit on both elements',
        qualification='Flap data were previously inspected; grouped checks are not a blind experiment. Refined grid excluded from CV.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                      [args.coarse/'study.json']+[p/'study.json' for p in args.refined or []]+
                      sorted(DATA.glob('d*/experiment_drag.csv'))+sorted(DATA.glob('d*/mses_drag.csv'))})
    write(args.out/'summary.json', summary)
    figures(rows, scores, pairs, best, args.out)
    transition_figure(rows, args.out)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
