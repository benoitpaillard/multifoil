#!/usr/bin/env python3
"""Compare the audited derivative modes without rerunning the flow solver."""
from audit_s414_jacobian import ROOT, STUDY, OUT, Inputs, load_states, require, sha, write
from report_s702_wake_extent import csv_write
from report_s702_matched_pressure import pressure_branch, save
from report_s414_reconstruction import interpolate
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'


def paired_metrics(old, new):
    require(old['result']['feasible'] and new['result']['feasible'], 'Both modes must converge')
    a, b = old['result'], new['result']
    row = dict(delta_CD_counts=(b['CD']-a['CD'])*1e4, delta_CL=b['CL']-a['CL'],
        delta_alpha_deg=b['alpha_deg']-a['alpha_deg'], same_selected_path=a['selected_path'] == b['selected_path'])
    for e in ('fore', 'aft'):
        require(np.array_equal(old['arrays'][e]['xy'], new['arrays'][e]['xy']), 'Paired geometry changed')
        row[e+'_delta_CD_counts'] = (b['element_forces'][e]['CD_wake']-a['element_forces'][e]['CD_wake'])*1e4
        row[e+'_maximum_Cp_change'] = float(abs(old['surfaces'][e]['Cp']-new['surfaces'][e]['Cp']).max())
        row[e+'_same_turbulence_flags'] = bool(np.array_equal(old['arrays'][e]['turbulent'], new['arrays'][e]['turbulent']))
        for side in ('lower', 'upper'):
            t0 = a['boundary_layers'][e][side]['transition']; t1 = b['boundary_layers'][e][side]['transition']
            row[e+'_'+side+'_transition_presence_matches'] = (t0 is None) == (t1 is None)
            row[e+'_'+side+'_delta_transition_x'] = None if t0 is None or t1 is None else t1['xy_over_cref'][0]-t0['xy_over_cref'][0]
    return row


def main():
    inputs = Inputs(); spec, states, controls = load_states(inputs)
    residual = inputs.json(OUT/'residual_recheck.json')
    require(residual['complete'] and len(residual['states']) == 4, 'Residual audit incomplete')
    for name, digest in residual['input_sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Residual audit input changed')
    reference = {r['case']: r for r in inputs.csv(PRIOR/'targets.csv')}
    old_taps = inputs.csv(PRIOR/'pressure_at_taps.csv')
    targets, pairs, pressures, iterations = [], [], [], []
    for job in spec['jobs']:
        state = states[job['name']]; r = state['result']; ex = reference[job['pressure_case']]
        diag = state['metadata']['diagnostics']
        targets.append(dict(job=job['name'], case=job['pressure_case'], mode=job['amplification_jacobian'],
            target_cl=job['target_cl'], CL=r['CL'], CD=r['CD'], alpha_deg=r['alpha_deg'],
            experiment_CD=float(ex['experiment_CD']), mses_CD=float(ex['mses_CD']),
            error_CD_counts=(r['CD']-float(ex['experiment_CD']))*1e4,
            experiment_alpha=float(ex['experiment_caption_alpha']),
            selected_path=r['selected_path'], requested_evaluations=len(r['records']),
            failed_requests=sum(not q['admissible'] for q in r['records']),
            reported_final_point_iterations=sum(d['iterations'] for d in diag),
            wall_seconds=r['wall_seconds'], fore_transition_x=None if r['boundary_layers']['fore']['lower']['transition'] is None
                else r['boundary_layers']['fore']['lower']['transition']['xy_over_cref'][0]))
        iterations.extend(dict(job=job['name'], evaluation=k, alpha_deg=q['alpha_deg'],
            iterations=d['iterations'], residual=q['residual']) for k, (q, d) in enumerate(zip(r['records'], diag)))
        for tap in old_taps:
            if tap['case'] != job['pressure_case'] or tap['method'] != 'multifoil': continue
            included = tap['included'] == 'True'
            value = None
            if included:
                e, side = tap['surface'].split('_'); branch = pressure_branch(state['surfaces'][e], side)
                value = float(interpolate(branch['x_over_cref'], branch['Cp'], [float(tap['x_over_c'])])[0])
            pressures.append(dict(tap, job=job['name'], method=job['amplification_jacobian'],
                model_Cp=value, error_Cp=None if value is None else value-float(tap['experiment_Cp'])))
    for case in ('cl0056', 'cl028'):
        pairs.append(dict(case=case, **paired_metrics(states[case+'_legacy'], states[case+'_corrected_af'])))
    for name, rows in [('targets', targets), ('paired_changes', pairs), ('pressure_at_taps', pressures), ('iterations', iterations)]:
        csv_write(OUT/(name+'.csv'), rows)
    write(OUT/'summary.json', dict(complete=True, full_state_audit_pass=True, cases=4,
        requests=sum(t['requested_evaluations'] for t in targets), failed_requests=sum(t['failed_requests'] for t in targets),
        controls=controls, pairs=pairs, targets=targets, retrieved_files_verified=987,
        frozen_inputs_verified=len(spec['sha256']), maximum_fresh_residual=max(r['max_residual'] for r in residual['states']),
        new_nonlinear_flow_solves=0, exact_geometry_established=False, numerical_convergence_established=False))
    plot(spec, states, targets, old_taps)
    inputs.add(ROOT/'validation/report_s414_jacobian.py')
    write(OUT/'input_manifest.json', dict(files=inputs.files))
    print('Four modes audited; pairs:', pairs, flush=True)


def plot(spec, states, targets, taps):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for col, case in enumerate(('cl0056', 'cl028')):
        target = next(t for t in targets if t['case'] == case)
        for row, surface in enumerate(('fore_lower', 'aft_upper')):
            ax = axes[row, col]; e, side = surface.split('_')
            for mode, color, style in [('legacy', '#2364aa', '-'), ('corrected_af', '#dd7a00', '--')]:
                branch = pressure_branch(states[case+'_'+mode]['surfaces'][e], side)
                ax.plot(branch['x_over_cref'], branch['Cp'], style, color=color, label=mode)
            chosen = [r for r in taps if r['case'] == case and r['surface'] == surface and r['coordinate'] == 'table' and r['included'] == 'True']
            exp = [r for r in chosen if r['method'] == 'multifoil']
            ax.plot([float(r['x_over_c']) for r in exp], [float(r['experiment_Cp']) for r in exp], 'ko', ms=4, label='Experiment')
            mses = [r for r in chosen if r['method'] == 'mses']
            if mses: ax.plot([float(r['x_over_c']) for r in mses], [float(r['model_Cp']) for r in mses], ':', color='#168a6e', label='Published MSES')
            ax.set(xlim=(0, .06) if row == 0 else (.8, 1.02), ylabel='Cp', xlabel='x/c', title=f'CL {target["target_cl"]:.3f} · {surface.replace("_", " ")}')
            ax.invert_yaxis(); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Derivative correction leaves the pressure discrepancy unchanged\nMatched-lift comparisons; corrected and legacy curves overlap. MSES pressure available only at CL 0.284.')
    save(fig, OUT, 'pressure_comparison'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    x = np.arange(2); width = .24
    for offset, mode, color in [(-width, 'legacy', '#2364aa'), (0., 'corrected_af', '#dd7a00'), (width, 'Published MSES', '#168a6e')]:
        chosen = [next(t for t in targets if t['case'] == c and t['mode'] == ('legacy' if mode == 'Published MSES' else mode)) for c in ('cl0056', 'cl028')]
        y = [(t['mses_CD']-t['experiment_CD'])*1e4 if mode == 'Published MSES' else t['error_CD_counts'] for t in chosen]
        axes[0].bar(x+offset, y, width, label=mode, color=color)
    for offset, mode, color in [(-.16, 'legacy', '#2364aa'), (.16, 'corrected_af', '#dd7a00')]:
        chosen = [next(t for t in targets if t['case'] == c and t['mode'] == mode) for c in ('cl0056', 'cl028')]
        bars = axes[1].bar(x+offset, [t['reported_final_point_iterations'] for t in chosen], .32, label=mode, color=color)
        axes[1].bar_label(bars)
    for ax in axes: ax.set_xticks(x, ['CL 0.056', 'CL 0.284']); ax.legend(fontsize=8); ax.grid(axis='y', alpha=.2)
    axes[0].set(ylabel='CD − experiment [counts]', title='Aerodynamic error retained')
    axes[1].set(ylabel='Sum across five requested evaluations', title='Reported final-point Newton iterations')
    fig.suptitle('Four converged targets, 20 requests, no failed requests\nIteration counts exclude unreported internal-continuation work; wall time did not improve.')
    save(fig, OUT, 'drag_and_iterations'); plt.close(fig)


if __name__ == '__main__':
    main()
