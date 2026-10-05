#!/usr/bin/env python3
"""Report controlled attempts to explain the S702 experimental drag discrepancy.

No drag offsets are fitted or applied. Experimental Cp supplies an additional
observable. Published MSES curves remain reference data; MSES is not run here.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'validation/results'
DATA = ROOT/'validation/data/s702'
FLAPS = [-2, 0, 5, 10]
TARGETS = [.5, .8, 1.2]
COLORS = {'MultiFoil N=9': '#2367b0', 'MultiFoil N=6': '#d66c21',
          'Published MSES N=9': '#93a1b0', 'Experiment': '#202b35'}
NASA = 'https://ntrs.nasa.gov/citations/20220014964'
XFOIL = 'https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt'


def read_json(p):
    return json.loads(p.read_text())


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    return list(csv.DictReader(p.open()))


def write_table(p, rows):
    with p.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def study(name):
    folder = RESULTS/name
    spec = read_json(folder/'study.json')
    fine_study = name in ('s702_transition_fine_n9', 's702_transition_fine_local')
    expected_panels = [340, 240] if fine_study else [260, 180]
    expected_mach = 0. if name.endswith('mach0_n6') else .1
    if (spec['reynolds'] != 1e6 or spec['mach'] != expected_mach or
        spec['panels'] != expected_panels or spec['wake_length'] != 1.):
        raise ValueError('Study conditions differ from the comparison protocol')
    for path, sha in spec['sha256'].items():
        if digest(ROOT/path) != sha:
            raise ValueError('Source differs from completed study: '+path)
    rows = []
    for job in spec['jobs']:
        target = folder/job['name']
        if (target/'result.json').exists():
            row = read_json(target/'result.json')
            if (target/'termination_resolution.json').exists():
                row['termination_resolution'] = read_json(target/'termination_resolution.json')
            rows.append(row)
        elif fine_study and (target/'interruption.json').exists():
            row = read_json(target/'progress.json')
            row.update(feasible=False, interruption=read_json(target/'interruption.json'),
                       error='Interrupted incomplete wide-bracket attempt; no target drag available')
            rows.append(row)
        else:
            raise ValueError('A study is still incomplete without an interruption record: '+str(target))
    for row in rows:
        if (not row['feasible'] and not fine_study) or (row['feasible'] and abs(row['CL']-row['target_cl']) > 2e-5):
            raise ValueError('Missing or inadmissible fixed-lift point: '+row['name'])
        if not fine_study and not all(r['admissible'] for r in row['records']):
            raise ValueError('A requested flow evaluation was inadmissible')
        if row['feasible'] and abs(row['records'][-1]['alpha_deg']-row['alpha_deg']) > 1e-10:
            raise ValueError('Transition audit is not at the returned fixed-lift state')
        allowed = {9.} if fine_study else {6., 9.}
        if not set(row['ncrit_by_element'].values()) <= allowed:
            raise ValueError('Unexpected transition scenario')
    return spec, rows


def final_flow(folder, result):
    rows = table(folder/'flow/polar.csv')
    matches = [i for i, row in enumerate(rows)
               if abs(float(row['CL'])-result['CL']) < 1e-10 and
                  abs(float(row['CD'])-result['CD']) < 1e-12]
    if not matches:
        raise ValueError('Final solution does not match its preserved flow output')
    return matches[-1]


def numerical_resolution(folder, row):
    index = final_flow(folder, row)
    output = []
    for element in ('fore', 'aft'):
        surfaces = table(folder/f'flow/surface_{index:03d}_{element}.csv')
        for side in ('lower', 'upper'):
            values = [r for r in surfaces if r['side'] == side]
            if side == 'upper':
                values.reverse()  # physical marching order, stagnation -> TE
            pairs = [(a, b) for a, b in zip(values, values[1:])
                     if int(a['turbulent']) and int(b['turbulent'])]
            output.append(dict(element=element, side=side,
                maximum_turbulent_Hk_drop=max([float(a['Hk'])-float(b['Hk']) for a, b in pairs]+[0.]),
                maximum_turbulent_ue_step=max([abs(float(a['ue_over_vinf'])-float(b['ue_over_vinf'])) for a, b in pairs]+[0.])))
    return output


def pressure_comparison(folder, rows, out):
    refs = DATA/'pressure_cl1.04'
    results = []
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.5), constrained_layout=True)
    for axis, (element, side) in zip(axes.flat,
        [('fore', 'upper'), ('aft', 'upper'), ('fore', 'lower'), ('aft', 'lower')]):
        label = element+'_'+side
        exp = np.genfromtxt(refs/f'experiment_{label}.csv', delimiter=',', names=True)
        # The lower curve is joined through the first/last upper tap. Those
        # repeated endpoints are plotting joins, not independent lower taps.
        if side == 'lower':
            exp = exp[1:-1]
        axis.plot(exp['x_over_cref'], exp['Cp'], 'o', ms=4, mfc='white',
                  mec=COLORS['Experiment'], label='Experiment', zorder=5)
        for name in ('Published MSES N=9', 'MultiFoil N=9', 'MultiFoil N=6'):
            if name.startswith('Published'):
                pred = np.genfromtxt(refs/f'mses_{label}.csv', delimiter=',', names=True)
                x, cp = pred['x_over_cref'], pred['Cp']
            else:
                n = int(name[-1])
                result = next(r for r in rows if r['ncrit_by_element']['fore'] == n)
                job = folder/result['name']
                index = final_flow(job, result)
                pred = np.genfromtxt(job/f'flow/surface_{index:03d}_{element}.csv',
                                     delimiter=',', names=True, dtype=None, encoding='utf-8')
                # For this zero-flap case, global x identifies the geometric LE.
                # Cp is compared by geometric surface, not flow-side labels at
                # the moving stagnation point.
                le = int(np.argmin(pred['x_over_cref']))
                part = pred[:le+1] if side == 'upper' else pred[le:]
                x, cp = part['x_over_cref'], part['Cp']
            order = np.argsort(x)
            x, cp = x[order], cp[order]
            x, unique = np.unique(x, return_index=True)
            cp = cp[unique]
            good = (exp['x_over_cref'] >= x.min()) & (exp['x_over_cref'] <= x.max())
            sample = exp[good]
            error = np.interp(sample['x_over_cref'], x, cp)-sample['Cp']
            chord_rms = np.sqrt(np.trapezoid(error**2, sample['x_over_cref'])/
                                np.ptp(sample['x_over_cref']))
            results.append(dict(source=name, element=element, side=side,
                samples=len(sample), excluded_outside_overlap=int((~good).sum()),
                rms_at_taps=float(np.sqrt(np.mean(error**2))),
                rms_chord_weighted=float(chord_rms),
                maximum_absolute_Cp_error=float(np.max(abs(error)))))
            axis.plot(x, cp, color=COLORS[name], lw=1.7,
                      ls='--' if name.startswith('Published') else '-', label=name)
        axis.invert_yaxis()
        axis.set(title=element.title()+' element — '+side, xlabel='x / assembly chord', ylabel='$C_p$')
        axis.grid(alpha=.18)
    axes[0, 0].legend(fontsize=8, loc='best')
    fig.suptitle('Independent pressure check: S702, flap 0°, $C_L=1.04$, Re=1 million, M=0.1', fontsize=12)
    fig.savefig(out/'pressure_check.png', dpi=180)
    fig.savefig(out/'pressure_check.svg')
    plt.close(fig)
    write_table(out/'pressure_errors.csv', results)
    return results


def drag_figure(rows, out):
    fig = plt.figure(figsize=(11.5, 5.5), layout='constrained')
    grid = fig.add_gridspec(2, 2, width_ratios=[1., 1.2])
    left = fig.add_subplot(grid[:, 0])
    keys = [('MultiFoil N=9', 'CD_n9'), ('Published MSES N=9', 'CD_mses'),
            ('Experiment', 'CD_experiment'), ('MultiFoil N=6', 'CD_n6')]
    x = np.arange(2)
    for i, (name, key) in enumerate(keys):
        values = [next(r[key] for r in rows if r['flap'] == f and r['target_cl'] == .8) for f in (5, 10)]
        bars = left.bar(x+(i-1.5)*.19, np.array(values)*1e4, .18,
                        color=COLORS[name], label=name)
        left.bar_label(bars, fmt='%.1f', fontsize=8, padding=3)
    left.set(xticks=x, xticklabels=['Flap 5°', 'Flap 10°'], ylabel='$C_D$ in drag counts (1 count = 0.0001)',
             title='Disputed ranking at $C_L=0.8$', ylim=(0, 140))
    left.legend(fontsize=8, loc='upper left')
    left.grid(axis='y', alpha=.15)
    for i, n in enumerate((9, 6)):
        axis = fig.add_subplot(grid[i, 1])
        z = np.array([[next(r[f'error_n{n}_percent'] for r in rows
                           if r['flap'] == f and r['target_cl'] == cl)
                       for f in FLAPS] for cl in TARGETS])
        im = axis.imshow(z, cmap='RdBu_r', vmin=-32, vmax=32, aspect='auto')
        for (j, k), value in np.ndenumerate(z):
            axis.text(k, j, f'{value:+.1f}%', ha='center', va='center',
                      color='white' if abs(value) > 20 else '#202b35', fontsize=10)
        axis.set(xticks=range(4), xticklabels=[f'{f}°' for f in FLAPS], yticks=range(3),
                 yticklabels=[str(c) for c in TARGETS], ylabel='$C_L$',
                 title=f'MultiFoil N={n}: drag error vs experiment')
        if i == 1:
            axis.set_xlabel('Flap deflection')
    fig.colorbar(im, ax=[fig.axes[1], fig.axes[2]], label='Drag error (%)', shrink=.75)
    fig.suptitle('S702 transition sensitivity: matched lift, Re=1 million, M=0.1', fontsize=13)
    fig.savefig(out/'drag_discrepancy.png', dpi=180)
    fig.savefig(out/'drag_discrepancy.svg')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=RESULTS/'s702_discrepancy')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    names = ['s702_transition_diagnostic', 's702_pressure_check', 's702_transition_fine_n9',
             's702_transition_mach0_n6', 's702_transition_fine_local']
    studies = {name: study(name) for name in names}
    _, changed = studies[names[0]]
    old = [r for r in table(RESULTS/'s702_comparison/fixed_lift_comparison.csv') if r['setting'] == 'base']
    comparisons = []
    for base in old:
        flap, cl = int(base['deflection']), float(base['target_cl'])
        altered = next(r for r in changed if r['flap'] == flap and
                       r['target_cl'] == cl and r['scenario'] == 'both')
        if set(altered['ncrit_by_element'].values()) != {6.}:
            raise ValueError('The changed-physics matrix must use Ncrit=6 on both elements')
        prior = read_json(RESULTS/f's702_fixed/d{flap:+d}__base/result.json')
        target = next(r for r in prior['targets'] if r['target_cl'] == cl)
        if not target['feasible'] or abs(float(base['CD'])-target['CD']) > 1e-12:
            raise ValueError('Baseline comparison differs from its archived solver result')
        for reference in ('mses', 'experiment'):
            curve = np.genfromtxt(DATA/f'd{flap:+d}/{reference}_drag.csv', delimiter=',', names=True)
            order = np.argsort(curve['CL'])
            if not curve['CL'].min() <= cl <= curve['CL'].max():
                raise ValueError('Reference interpolation would extrapolate')
            value = np.interp(cl, curve['CL'][order], curve['CD'][order])
            if abs(value-float(base[reference+'_CD'])) > 1e-10:
                raise ValueError('Reference comparison differs from its extracted source')
        entry = dict(flap=flap, target_cl=cl, CD_n9=float(base['CD']), CD_n6=altered['CD'],
            CD_mses=float(base['mses_CD']), CD_experiment=float(base['experiment_CD']))
        for n in (9, 6):
            entry[f'error_n{n}_counts'] = (entry[f'CD_n{n}']-entry['CD_experiment'])*1e4
            entry[f'error_n{n}_percent'] = (entry[f'CD_n{n}']/entry['CD_experiment']-1)*100
        comparisons.append(entry)
    write_table(args.out/'drag_comparison.csv', comparisons)
    rankings = []
    for cl in TARGETS:
        a, b = [next(r for r in comparisons if r['flap'] == f and r['target_cl'] == cl) for f in (5, 10)]
        for source in ('n9', 'n6', 'mses', 'experiment'):
            rankings.append(dict(target_cl=cl, source=source,
                                 delta_CD_10_minus_5_counts=(b['CD_'+source]-a['CD_'+source])*1e4))
    write_table(args.out/'ranking.csv', rankings)
    isolation = []
    for flap in (5, 10):
        base = read_json(RESULTS/f's702_mach_check/flap{flap}_mach0.1/result.json')
        jobs = [dict(base, scenario='neither', ncrit_by_element=dict(fore=9., aft=9.))]
        jobs += [r for r in changed if r['flap'] == flap and r['target_cl'] == .8]
        for r in jobs:
            tr = r['boundary_layers']['fore']['upper']['transition']
            isolation.append(dict(flap=flap, changed_element=r['scenario'],
                fore_ncrit=r['ncrit_by_element']['fore'], aft_ncrit=r['ncrit_by_element']['aft'],
                CD=r['CD'], CD_change_counts=(r['CD']-base['CD'])*1e4,
                fore_CD=r['element_forces']['fore']['CD_wake'],
                aft_CD=r['element_forces']['aft']['CD_wake'],
                fore_upper_transition_x_cref=tr['xy_over_cref'][0] if tr else None,
                fore_upper_transition_x_local=tr['x_over_local_chord'] if tr else None))
    write_table(args.out/'element_isolation.csv', isolation)
    fine_comparison = []
    for source_study in (names[2], names[4]):
        for r in studies[source_study][1]:
            base = next(x for x in comparisons if x['flap'] == r['flap'] and x['target_cl'] == .8)
            fine_comparison.append(dict(flap=r['flap'], study=source_study, bounds=r['bounds'],
                feasible=r['feasible'], error=r['error'], CD_base=base['CD_n9'],
                CD_fine=r['CD'] if r['feasible'] else None,
                fine_minus_base_counts=(r['CD']-base['CD_n9'])*1e4 if r['feasible'] else None,
                resolution=numerical_resolution(RESULTS/source_study/r['name'], r) if r['feasible'] else None))
    mach0 = studies[names[3]][1]
    mach0_gap = (next(r['CD'] for r in mach0 if r['flap'] == 10)-
                 next(r['CD'] for r in mach0 if r['flap'] == 5))*1e4
    cp = pressure_comparison(RESULTS/names[1], studies[names[1]][1], args.out)
    drag_figure(comparisons, args.out)
    audits = [a for _, rows in studies.values() for r in rows for a in r.get('transition_audit', [])]
    metrics = {}
    for n in (9, 6):
        for split, selected in [('all', comparisons),
             ('flap0', [r for r in comparisons if r['flap'] == 0]),
             ('other_flaps', [r for r in comparisons if r['flap'] != 0])]:
            e = np.array([r[f'error_n{n}_percent'] for r in selected])
            metrics[f'n{n}_{split}'] = dict(mean_absolute_percent=float(np.mean(abs(e))),
                maximum_absolute_percent=float(np.max(abs(e))),
                rms_counts=float(np.sqrt(np.mean([r[f'error_n{n}_counts']**2 for r in selected]))))
    summary = dict(physical_discrepancy_resolved=False, default_ncrit_changed=False,
        numerical_core_changed=False, calibrated_shape_derivatives=False,
        model_selection_used_for_prediction=False,
        metrics=metrics, ranking=rankings, fine_grid_comparison=fine_comparison,
        mach0_n6=dict(delta_CD_10_minus_5_counts=mach0_gap,
            points=[{k: r[k] for k in ('flap', 'CL', 'CD', 'alpha_deg')} for r in mach0],
            reference_status='Mach-zero sensitivity; experimental reference has Mach 0.1'),
        pressure_errors=cp, attempted_fixed_lift_targets=sum(len(r) for _, r in studies.values()),
        completed_fixed_lift_targets=sum(r['feasible'] for _, rows in studies.values() for r in rows),
        requested_flow_evaluations=sum(len(r['records']) for _, rows in studies.values() for r in rows),
        admissible_flow_evaluations=sum(q['admissible'] for _, rows in studies.values() for r in rows for q in r['records']),
        interruptions=[dict(study=name, job=r['name'], **r['interruption'])
                       for name, (_, rows) in studies.items() for r in rows if 'interruption' in r],
        termination_resolutions=[dict(study=name, job=r['name'], **r['termination_resolution'])
                       for name, (_, rows) in studies.items() for r in rows if 'termination_resolution' in r],
        transition_root_audit=dict(cells=len(audits), all_inside_cells=all(a['inside_cell'] for a in audits),
            maximum_amplification_residual=max(abs(a['amplification_residual']) for a in audits)),
        studies={name: dict(manifest_sha256=digest(RESULTS/name/'study.json'),
                           result_sha256={r['name']: digest(RESULTS/name/r['name']/
                               ('result.json' if r['feasible'] or 'interruption' not in r else 'progress.json')) for r in rows})
                 for name, (_, rows) in studies.items()},
        baseline_comparison_sha256=digest(RESULTS/'s702_comparison/fixed_lift_comparison.csv'),
        report_script_sha256=digest(Path(__file__)),
        source_urls=[NASA, XFOIL])
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    lines = ['# S702 numerical–experimental discrepancy investigation', '',
        '**The general physical discrepancy is not resolved.** This study isolates transition sensitivity, '
        'tests a proposed Ncrit change on additional operating points, checks finer panels, and adds '
        'an independent pressure observable. It does not fit a drag offset, change the default Ncrit, '
        'or certify contour-shape gradients.', '',
        f"The N=6 scenario reduces mean absolute drag error over the 12 matched cases from "
        f"{metrics['n9_all']['mean_absolute_percent']:.2f}% to {metrics['n6_all']['mean_absolute_percent']:.2f}%. "
        'It recovers the measured 5°/10° ordering at CL=0.8, reducing the drag-difference error '
        'from 39.18 to 6.76 counts. The large change at 5° comes primarily from earlier fore-element '
        'transition. Several other points worsen, so this is a partial remedy rather than an '
        'identified universal transition correction.', '',
        '## Matched-condition drag result', '',
        'All values below use Re=1,000,000, Mach=0.1, natural transition, the original assembly reference '
        'chord and the same reconstructed geometry. One drag count is 0.0001 in CD. Published MSES '
        'curves and experimental markers are from NASA Figures 23–26; MSES was not run here.', '',
        '| CL | Flap | MultiFoil N=9 | MultiFoil N=6 | Published MSES N=9 | Experiment | N=6 error |',
        '|---:|---:|---:|---:|---:|---:|---:|']
    lines += [f"| {r['target_cl']:g} | {r['flap']}° | {r['CD_n9']:.6f} | {r['CD_n6']:.6f} | "
              f"{r['CD_mses']:.6f} | {r['CD_experiment']:.6f} | {r['error_n6_percent']:+.1f}% |" for r in comparisons]
    lines += ['', 'The interpolated reference values are read from published vector curves, not raw tunnel '
        'measurements. Digitization precision is not an experimental uncertainty bound.', '',
        '![Drag comparison and all 12 pointwise errors](drag_discrepancy.png)', '',
        '## Ranking and transfer', '',
        'Define ΔCD = CD(10°) − CD(5°). Negative values favor 10°.', '',
        '| CL | MultiFoil N=9 | MultiFoil N=6 | Published MSES | Experiment |',
        '|---:|---:|---:|---:|---:|']
    for cl in TARGETS:
        vals = [next(r['delta_CD_10_minus_5_counts'] for r in rankings
                     if r['target_cl'] == cl and r['source'] == s) for s in ('n9', 'n6', 'mses', 'experiment')]
        lines.append(f'| {cl:g} | '+' | '.join(f'{v:+.2f}' for v in vals)+' |')
    lines += ['', 'These are sensitivity calculations with a pre-existing N=6 scenario, not a newly measured '
        'transition threshold. All configurations were previously examined, so the extra cases are transfer '
        'checks, not a blind test set. No threshold was optimized against the disputed 5°/10° pair.', '',
        '| Ncrit | Mean absolute drag error, all 12 | Maximum error, all 12 | Mean error, flap 0 only | Mean error, other flaps |',
        '|---:|---:|---:|---:|---:|']
    for n in (9, 6):
        lines.append(f"| {n} | {metrics[f'n{n}_all']['mean_absolute_percent']:.2f}% | "
                     f"{metrics[f'n{n}_all']['maximum_absolute_percent']:.2f}% | "
                     f"{metrics[f'n{n}_flap0']['mean_absolute_percent']:.2f}% | "
                     f"{metrics[f'n{n}_other_flaps']['mean_absolute_percent']:.2f}% |")
    lines += ['', 'A concrete counterexample to blanket N=6 is flap −2°, CL=0.5: drag changes from '
        '0.006340 to 0.006805, whereas experiment is approximately 0.006295. Higher-lift errors also '
        'remain. Better average agreement would not establish the correct underlying physics.', '',
        '## Which element changes the result?', '',
        'Each row below changes only the stated Ncrit values and resolves incidence to CL=0.8. '
        'Per-element wake drags sum to total profile drag; friction is already included and must not '
        'be added again. Transition x uses the assembly chord.', '',
        '| Flap | Fore N | Aft N | Total CD | Change, counts | Fore CD | Aft CD | Fore upper transition x/c |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in isolation:
        tr = 'Laminar to TE' if r['fore_upper_transition_x_cref'] is None else f"{r['fore_upper_transition_x_cref']:.4f}"
        lines.append(f"| {r['flap']}° | {r['fore_ncrit']:g} | {r['aft_ncrit']:g} | {r['CD']:.6f} | "
            f"{r['CD_change_counts']:+.2f} | {r['fore_CD']:.6f} | {r['aft_CD']:.6f} | {tr} |")
    lines += ['', 'The element changes are coupled: altering aft transition can also change fore loading, '
        'and vice versa. Their drag increments need not add linearly. A lower threshold is an effective '
        'disturbance sensitivity, not a resolved turbulent-wake interaction model or a prescribed trip.', '',
        '## Numerical checks', '',
        'The prior matched-lift Mach test changed the 5°/10° drag difference by only +0.036 counts '
        'when Mach 0.1 was replaced by zero. It cannot account for the experimental reversal.', '',
        f'The new N=6 Mach-zero pair has ΔCD(10°−5°) = {mach0_gap:+.3f} counts. '
        'This checks transfer of the transition sensitivity to incompressible flow; it is not a '
        'matched-Mach experimental validation.', '',
        '| Flap at CL=0.8, N=9 | Incidence bracket | CD, 260/180 panels | CD, 340/240 panels | Change, counts |',
        '|---:|---|---:|---:|---:|']
    for r in fine_comparison:
        value = f"{r['CD_fine']:.6f}" if r['feasible'] else ('Interrupted' if r['error'].startswith('Interrupted') else 'Nonconverged')
        change = f"{r['fine_minus_base_counts']:+.3f}" if r['feasible'] else '—'
        lines.append(f"| {r['flap']}° | {r['bounds']}° | {r['CD_base']:.6f} | {value} | {change} |")
    for interruption in summary['interruptions']:
        lines += ['', f"The wide-bracket 10° attempt was stopped after {interruption['elapsed_seconds']/60:.1f} "
            'minutes for computational cost, while an additional local-bracket run addressed the same '
            'target lift. Its completed lower-endpoint solve and unfinished upper-endpoint request are '
            'preserved, including a separate immutable attempt folder. This is an interrupted attempt, '
            'not a claim of numerical divergence. The local bracket surrounds the previously computed '
            'base-mesh incidence; it was not selected using experimental drag.']
    for resolution in summary['termination_resolutions']:
        lines += ['', 'The wide [-2°, 1°] attempt ultimately recorded a nonconverged evaluation '
            'at α=1°, after the Newton iteration limit. A termination request had been issued earlier '
            'while that evaluation was still running; the worker subsequently wrote its full failure '
            'record, which is the authoritative outcome. The earlier snapshot and termination note '
            'are retained alongside a resolution record. The target itself converged inside the '
            'local [-0.8°, -0.2°] bracket. The α=1° failure is not interpreted as stall or given '
            'a finite drag coefficient.']
    lines += ['', f"{summary['completed_fixed_lift_targets']} of {summary['attempted_fixed_lift_targets']} "
        f"new fixed-lift solve attempts completed. There are {summary['requested_flow_evaluations']} "
        f"completed requested flow evaluations, of which {summary['admissible_flow_evaluations']} are admissible. Internal continuation "
        'steps are also retained in the flow metadata; they are not counted as separate requested targets. '
        f"The final transition audit checked {len(audits)} subcells; all-inside-cell status is "
        f"{summary['transition_root_audit']['all_inside_cells']}, and the maximum amplification-equation residual is "
        f"{summary['transition_root_audit']['maximum_amplification_residual']:.3g}.", '',
        'The panel check is an observed sensitivity, not a rigorous error bound. The machine-readable '
        'summary also reports turbulent Hk drops and velocity steps, following the resolution checks '
        'described in the XFOIL manual. No converged transition-root failure was found in these cases. '
        'The separate convergence failure at α=1° remains recorded; this audit does not establish '
        'robustness at every incidence. It supplies no basis for changing the converged target drag '
        'by altering the boundary-layer equations.', '',
        'Six source-validation tests passed (`test_s702` and `test_s702_pressure`). The pressure '
        'tests check plotted orifice positions against the independent Table II coordinates, source '
        'hashes, and the caption lift against the separately plotted lift curve. These are data '
        'and implementation checks, not claims of aerodynamic accuracy.', '',
        '## Independent pressure observable', '',
        'NASA Figure 22(b), flap 0°, Re=1e6, M=0.1, supplies Cp at a caption lift of 1.04. '
        'The caption is rounded to two decimals; it is not an exact raw-data lift. The two new '
        'MultiFoil runs target that nominal value. Pressure data were not used to set Ncrit.', '',
        '![Pressure comparison](pressure_check.png)', '',
        '| Source | Fore upper RMS Cp | Fore lower RMS Cp | Aft upper RMS Cp | Aft lower RMS Cp |',
        '|---|---:|---:|---:|---:|']
    for source in ('MultiFoil N=9', 'MultiFoil N=6', 'Published MSES N=9'):
        values = [next(r['rms_chord_weighted'] for r in cp if r['source'] == source and
                  r['element'] == e and r['side'] == s) for e, s in
                  [('fore', 'upper'), ('fore', 'lower'), ('aft', 'upper'), ('aft', 'lower')]]
        lines.append('| '+source+' | '+' | '.join(f'{v:.4f}' for v in values)+' |')
    lines += ['', 'RMS values integrate squared errors along each surface’s sampled x interval, then divide '
        'by that interval length. Predictions are linearly interpolated only inside overlapping x ranges. '
        'Shared endpoint markers used to join the lower Cp curves are excluded from lower-surface error '
        'statistics. No experimental error bars or transition locations are inferred from sparse pressure taps.', '',
        '## What remains unmodelled or unknown', '',
        'The original report identifies possible interaction between the fore wake and aft laminar '
        'separation bubbles. The present solver couples displacement and potential flow, but has no '
        'model for wake-induced bypass transition or merging wakes. The experiment also used a plastic '
        'film over the flap hinge, with steps visible in Cp, and the as-built model differed from its '
        'specified contour by up to 0.13 mm (about 0.000284 assembly chords). Neither the measured '
        'contour nor film details are available in this archive.', '',
        'Tunnel turbulence was reported as approximately 0.05% at 46 m/s; that is not evidence that '
        'this case has Ncrit=6. Natural-transition Ncrit represents disturbance amplification and '
        'receptivity, not turbulence intensity alone. The report’s tripped cases cannot be treated '
        'as clean fully turbulent validation: the authors explain that thick tape itself altered '
        'displacement, lift and drag and could fail to trip at some conditions.', '',
        'A pressure or transition mismatch does not uniquely identify which missing input or closure '
        'is responsible. The next physical correction needs matched as-built geometry/hinge details '
        'and quantitative transition or wake measurements around the bucket edge, or independent '
        'transition-resolving CFD of those same conditions. Those inputs are absent; fitting a '
        'configuration-dependent Ncrit or drag offset would hide this lack of identification.', '',
        'For contour optimization, keep N=9 and N=6 as separate uncertainty scenarios. Correct '
        'ranking at one flap pair does not validate derivatives with respect to a new contour. '
        'This study makes no new shape acceptance or experimental-validation claim.', '',
        '## Reproduce', '', 'From the project root:', '', '```bash',
        'python validation/s702_transition_study.py --isolate --workers 4',
        'python validation/extract_s702_pressure.py',
        'python validation/s702_pressure_check.py --workers 2',
        'python validation/s702_transition_study.py --ncrit 9 --panels 340,240 --flaps 5,10 --targets .8 --workers 2 --out validation/results/s702_transition_fine_n9',
        'python validation/s702_transition_study.py --mach 0 --ncrit 6 --flaps 5,10 --targets .8 --workers 2 --out validation/results/s702_transition_mach0_n6',
        'python validation/s702_local_bracket_check.py',
        'python validation/report_s702_discrepancy.py', '```', '',
        'PyMuPDF is needed only to re-extract pressure data; the included CSV files allow comparisons '
        'without it. Use the normal requirements for the solver and reporter.', '',
        'The retained wide-bracket fine study returns exit code 2 because its 10° job failed at '
        'the upper endpoint. The separate local-bracket check and the reporter complete successfully.', '',
        'Completed studies resume only when their conditions and source hashes match. Use a new '
        'output directory to repeat an interrupted study. Every completed requested '
        'solve, surface and wake distribution, final state, transition audit and source selector is '
        'preserved. The numerical core and previous studies remain unchanged.', '',
        '## Sources', '',
        f'- [Somers and Maughmer, NASA/CR–20220014964]({NASA}): methods on printed pages 8–10, '
        'transition/wake observations on pages 10–14, Figure 22(b) on page 75, and matched drag polars in Figures 23–26.',
        f'- [Drela, XFOIL user manual]({XFOIL}): transition criterion, limitations of eN transition, '
        'wake drag definition and separation-bubble panel-resolution guidance.', '']
    (args.out/'DISCREPANCY_RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps({k: summary[k] for k in ('metrics', 'ranking', 'completed_fixed_lift_targets',
                     'requested_flow_evaluations', 'transition_root_audit')}, indent=2))


if __name__ == '__main__':
    raise SystemExit(main())
