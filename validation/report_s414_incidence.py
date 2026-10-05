#!/usr/bin/env python3
"""Compare audited measured-incidence and matched-lift states without conflation."""
from audit_s414_incidence import ROOT, STUDY, OUT, Inputs, load_states, require, sha, write
from report_s702_wake_extent import csv_write
from report_s702_matched_pressure import pressure_branch, save
from report_s414_reconstruction import interpolate, pressure_metric, recovery_intervals
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
SURFACES = ('fore_upper', 'fore_lower', 'aft_upper', 'aft_lower')
COLORS = {'fixed_lift': '#2364aa', 'fixed_alpha': '#d86a00'}
LABELS = {'fixed_lift': 'MultiFoil at matched lift', 'fixed_alpha': 'MultiFoil at measured angle'}


def condition_metrics(result, job, reference):
    """Keep a same-angle result's lift error visible and its matched-lift score absent."""
    require(result['feasible'], 'Failed state cannot be scored')
    basis = job['comparison_basis']
    require(basis in ('fixed_lift', 'fixed_alpha'), 'Unknown comparison basis')
    ref_cl = float(reference['target_cl']); ref_cd = float(reference['experiment_CD'])
    if basis == 'fixed_lift':
        require(job['target_cl'] == ref_cl, 'Lift association changed')
    else:
        require('target_cl' not in job and job['reference_experimental_cl'] == ref_cl, 'Fixed-angle case has a lift target')
        require(job['target_alpha'] == float(reference['experiment_caption_alpha']), 'Measured angle changed')
    return dict(comparison_basis=basis, reference_caption_CL=ref_cl,
        reference_caption_alpha=float(reference['experiment_caption_alpha']),
        reference_appendix_CL=float(reference['appendix_CL']), reference_appendix_alpha=float(reference['appendix_alpha']),
        reference_appendix_CD=float(reference['appendix_CD']),
        CL=result['CL'], CD_context=result['CD'], alpha_deg=result['alpha_deg'],
        CL_minus_caption=result['CL']-ref_cl, CL_minus_appendix=result['CL']-float(reference['appendix_CL']),
        alpha_minus_caption=result['alpha_deg']-float(reference['experiment_caption_alpha']),
        matched_lift_experiment_CD=ref_cd,
        matched_lift_CD_error_counts=(result['CD']-ref_cd)*1e4 if basis == 'fixed_lift' else None)


def main():
    inputs = Inputs(); spec, states, controls = load_states(inputs)
    residual = inputs.json(OUT/'residual_recheck.json')
    require(residual['complete'] and len(residual['states']) == 4, 'Full residual audit incomplete')
    for name, digest in residual['input_sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Residual receipt changed: '+name)
    checks = {r['job']: r for r in residual['states']}
    refs = {r['case']: r for r in inputs.csv(PRIOR/'targets.csv')}
    source = inputs.csv(PRIOR/'pressure_at_taps.csv')
    targets, pressures, metrics, peaks, layers, separations, forces, amplification = [], [], [], [], [], [], [], []
    for job in spec['jobs']:
        c = job['pressure_case']; name = job['name']; state = states[name]; result = state['result']
        target = dict(job=name, case=c, **condition_metrics(result, job, refs[c]),
            selected_path=result['selected_path'], requested_evaluations=len(result['records']),
            failed_requests=sum(not r['admissible'] for r in result['records']),
            residual=checks[name]['max_residual'], wall_seconds=result['wall_seconds'])
        target['fore_lower_transition_x'] = None
        for e in ('fore', 'aft'):
            forces.append(dict(job=name, case=c, comparison_basis=job['comparison_basis'], element=e, **result['element_forces'][e]))
            for side in ('lower', 'upper'):
                b = result['boundary_layers'][e][side]; tr = b['transition']
                layers.append(dict(job=name, case=c, comparison_basis=job['comparison_basis'], element=e, side=side,
                    transition_x=None if tr is None else tr['xy_over_cref'][0], min_cf=b['min_cf'], max_Hk=b['max_Hk'],
                    separated_fraction=b['separated_arclength_fraction'], separation_regions=len(b['separation_regions'])))
                for index, region in enumerate(b['separation_regions']):
                    separations.append(dict(job=name, element=e, side=side, interval=index,
                        start_x=region['start']['xy_over_cref'][0], end_x=region['end']['xy_over_cref'][0],
                        reaches_trailing_edge=region['reaches_trailing_edge']))
                if e == 'fore' and side == 'lower':
                    target['fore_lower_transition_x'] = None if tr is None else tr['xy_over_cref'][0]
                    a = state['arrays'][e]
                    for k, j in enumerate(a['lower']):
                        lam = a['turbulent'][j] == 0
                        amplification.append(dict(job=name, case=c, comparison_basis=job['comparison_basis'], flow_index=k,
                            x=float(a['xy'][0, j]), xi=float(a['xi'][j]), laminar=bool(lam),
                            N=float(a['U'][2,j]) if lam else None, shear=None if lam else float(a['U'][2,j])))
        targets.append(target)
        for tap in source:
            if tap['case'] != c or tap['method'] != 'multifoil': continue
            included = tap['included'] == 'True'; cp = None
            if included:
                e, side = tap['surface'].split('_'); b = pressure_branch(state['surfaces'][e], side)
                cp = float(interpolate(b['x_over_cref'], b['Cp'], [float(tap['x_over_c'])])[0])
            pressures.append(dict(job=name, case=c, comparison_basis=job['comparison_basis'],
                surface=tap['surface'], coordinate=tap['coordinate'], table_index=tap['table_index'],
                plot_index=tap['plot_index'], x_over_c=None if tap['x_over_c']=='' else float(tap['x_over_c']),
                experiment_Cp=float(tap['experiment_Cp']), model_Cp=cp,
                error_Cp=None if cp is None else cp-float(tap['experiment_Cp']), included=included, exclusion=tap['exclusion']))
        for coord in ('table', 'plotted'):
            group = [r for r in pressures if r['job'] == name and r['coordinate'] == coord]
            for surface in (*SURFACES, 'all_surfaces'):
                rows = [r for r in group if surface == 'all_surfaces' or r['surface'] == surface]
                for region, chosen in [('all', rows), ('nose_x_le_006', [r for r in rows if r['surface'].startswith('fore_') and r['x_over_c'] is not None and r['x_over_c'] <= .06])]:
                    if not chosen: continue
                    metrics.append(dict(job=name, case=c, comparison_basis=job['comparison_basis'], coordinate=coord,
                        surface=surface, region=region, **pressure_metric(chosen)))
            aft = [r for r in group if r['surface'] == 'aft_upper' and r['included'] and .8 <= r['x_over_c'] <= .98]
            x = [r['x_over_c'] for r in aft]
            predicted = recovery_intervals(x, [r['model_Cp'] for r in aft])
            measured = recovery_intervals(x, [r['experiment_Cp'] for r in aft])
            for statistic, field in [('largest_rise', 'delta_Cp'), ('steepest_slope', 'slope')]:
                a, b = max(predicted, key=lambda r:r[field]), max(measured, key=lambda r:r[field])
                peaks.append(dict(job=name, case=c, coordinate=coord, statistic=statistic,
                    model_left_x=a['left_x'], model_right_x=a['right_x'], model_value=a[field],
                    experiment_left_x=b['left_x'], experiment_right_x=b['right_x'], experiment_value=b[field]))
    for c in ('cl0056', 'cl028'):
        a, b = states[c+'_control'], states[c+'_measured_alpha']
        require(all(np.array_equal(a['arrays'][e]['xy'], b['arrays'][e]['xy']) for e in ('fore', 'aft')), 'Geometry changed with incidence')
    for name, rows in [('targets', targets), ('pressure_at_taps', pressures), ('pressure_metrics', metrics),
                       ('recovery_peaks', peaks), ('boundary_layers', layers), ('separation_intervals', separations),
                       ('element_forces', forces), ('fore_amplification', amplification)]:
        csv_write(OUT/(name+'.csv'), rows)
    physical = {(r['case'], r['surface'], r['table_index']) for r in pressures if r['included']}
    write(OUT/'summary.json', dict(complete=True, full_state_audit_pass=True, cases=4,
        requests=sum(t['requested_evaluations'] for t in targets), failed_requests=sum(t['failed_requests'] for t in targets),
        controls=controls, targets=targets, maximum_fresh_residual=max(r['max_residual'] for r in residual['states']),
        retrieved_files_verified=1977, frozen_inputs_verified=len(spec['sha256']),
        unique_supported_pressure_observations=len(physical), pressure_comparison_rows=len(pressures),
        scored_pressure_rows=sum(r['included'] for r in pressures),
        mses_pressure_cases=['cl028'], mses_pressure_basis='Existing matched-lift reference only',
        new_nonlinear_flow_solves=0, exact_geometry_established=False, numerical_convergence_established=False,
        solution_uniqueness_established=False))
    plot(states, targets, source, amplification)
    inputs.add(ROOT/'validation/report_s414_incidence.py')
    write(OUT/'input_manifest.json', dict(files=inputs.files))
    print('Four selected states compared; supported experimental taps:', len(physical), flush=True)


def plot(states, targets, source, amplification):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    for full in (False, True):
        surfaces = SURFACES if full else ('fore_lower', 'aft_upper')
        fig, axes = plt.subplots(len(surfaces), 2, figsize=(12, 13 if full else 8), layout='constrained')
        for col, c in enumerate(('cl0056', 'cl028')):
            t = next(t for t in targets if t['case']==c and t['comparison_basis']=='fixed_alpha')
            for row, surface in enumerate(surfaces):
                ax = axes[row, col]; e, side = surface.split('_')
                for suffix, basis in [('control', 'fixed_lift'), ('measured_alpha', 'fixed_alpha')]:
                    b = pressure_branch(states[c+'_'+suffix]['surfaces'][e], side)
                    ax.plot(b['x_over_cref'], b['Cp'], color=COLORS[basis], label=LABELS[basis])
                raw = [r for r in source if r['case']==c and r['surface']==surface and r['coordinate']=='table']
                supported = [r for r in raw if r['method']=='multifoil' and r['included']=='True']
                excluded = [r for r in raw if r['method']=='multifoil' and r['included']!='True' and r['table_index']!='']
                ax.plot([float(r['x_over_c']) for r in supported], [float(r['experiment_Cp']) for r in supported], 'ko', ms=3.5, label='Experiment: supported taps')
                if excluded:
                    ax.plot([float(r['x_over_c']) for r in excluded], [float(r['experiment_Cp']) for r in excluded], 'x', color='#888888', label='Excluded source taps')
                mses = [r for r in raw if r['method']=='mses' and r['included']=='True']
                if mses: ax.plot([float(r['x_over_c']) for r in mses], [float(r['model_Cp']) for r in mses], ':', color='#168a6e', label='Published MSES at matched lift')
                if not full: ax.set_xlim((0,.06) if row==0 else (.8,1.02))
                ax.set(xlabel='x/c', ylabel='Cp', title=f'Experiment: α {t["reference_caption_alpha"]:.2f}°, CL {t["reference_caption_CL"]:.3f} · {surface.replace("_", " ")}')
                ax.invert_yaxis(); ax.grid(alpha=.2); ax.legend(fontsize=7)
        fig.suptitle('S414 pressure: matched lift versus measured incidence\nOrange model lifts are 0.1072 / 0.3592; experimental lifts are 0.056 / 0.284. No MSES Cp at the first case.')
        save(fig, OUT, 'full_pressure_comparison' if full else 'pressure_diagnostic'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout='constrained')
    for basis in ('fixed_lift','fixed_alpha'):
        rows = [t for t in targets if t['comparison_basis']==basis]
        axes[0].plot([t['reference_caption_CL'] for t in rows], [t['CL'] for t in rows], 'o', color=COLORS[basis], label=LABELS[basis])
        for t in rows:
            if basis=='fixed_alpha': axes[0].annotate(f'ΔCL +{t["CL_minus_caption"]:.4f}', (t['reference_caption_CL'], t['CL']), xytext=(6,4), textcoords='offset points', fontsize=8)
    axes[0].plot([0,.41],[0,.41],'k:',label='Equal lift'); axes[0].set(xlim=(0,.41), ylim=(0,.42), xlabel='Experimental CL', ylabel='Computed CL', title='Lift remains too high at the measured angle')
    for c, style in [('cl0056','-'),('cl028','--')]:
        for suffix, basis in [('control','fixed_lift'),('measured_alpha','fixed_alpha')]:
            name=c+'_'+suffix; rows=[r for r in amplification if r['job']==name]
            axes[1].plot([r['x'] for r in rows], [r['N'] if r['N'] is not None else np.nan for r in rows], style, color=COLORS[basis], label=f'{c}: '+('matched lift' if basis=='fixed_lift' else 'measured angle'))
            t=next(t for t in targets if t['job']==name)
            if t['fore_lower_transition_x'] is not None: axes[1].plot(t['fore_lower_transition_x'],9,'o',color=COLORS[basis])
    axes[1].axhline(9,color='k',ls=':',label='Ncrit 9'); axes[1].set(xlim=(0,.06),ylim=(0,10),xlabel='x/c',ylabel='Laminar amplification N',title='Early fore-lower transition persists')
    for ax in axes: ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Nominal geometry, Re 1m, Mach 0.10, Ncrit 9 and corrected derivative unchanged\nDifferent initialization histories remain a limitation; turbulent shear is not plotted as N.')
    save(fig,OUT,'lift_and_amplification'); plt.close(fig)


if __name__ == '__main__': main()
