#!/usr/bin/env python3
"""Compare the fully audited S414 lower-bucket sequence; no flow solves or fits."""
import platform
import importlib.metadata
from decimal import Decimal
from audit_s414_lower_bucket import ROOT, STUDY, OUT, STUDY_HASH, Inputs, load_states, require, sha, write
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_s702_wake_extent import csv_write, clearance, side_data
from report_s702_matched_pressure import pressure_branch, save
from report_s414_reconstruction import interpolate, common_support, pressure_metric, recovery_intervals

SOURCE = ROOT/'validation/results/s414_trip_qualification_20261003'
GEOMETRY = ROOT/'validation/results/s414_geometry_qualification_20261002'
SURFACES = ('fore_upper', 'fore_lower', 'aft_upper', 'aft_lower')
REGIONS = {'fore_lower': [('fore_front', -np.inf, .60), ('fore_slot_entry', .60, np.inf)],
    'aft_upper': [('aft_nose', -np.inf, .80), ('aft_front_recovery', .80, .86),
                  ('aft_bubble_recovery', .86, .94), ('aft_tail', .94, np.inf)]}
COLORS = {'multifoil': '#2364aa', 'mses': '#168a6e', 'experiment': '#222222'}


def reference_at(rows, lift):
    """Explicit missing coverage, never clamp or extrapolate a reference."""
    x = np.array([float(r['CL']) for r in rows])
    y = np.array([float(r['CD']) for r in rows])
    require(len(x) >= 2 and np.all(np.diff(x) > 0) and np.isfinite(y).all(), 'Invalid drag branch')
    if not x[0] <= lift <= x[-1]:
        return None
    return float(interpolate(x, y, [lift])[0])


def pressure_taps(rows, jobs):
    """Preserve physical table indices independently of duplicate plotting endpoints."""
    taps = []
    for job in jobs:
        for surface in SURFACES:
            selected = [r for r in rows if r['figure'] == job['pressure_figure'] and r['surface'] == surface]
            require(selected, 'Missing pressure surface')
            k = 0
            for r in selected:
                duplicate = r['endpoint_duplicate'] == 'True'
                require(float(r['CL']) == job['target_cl'], 'Pressure/lift association changed')
                taps.append(dict(r, case=job['pressure_case'], endpoint_duplicate=int(duplicate),
                    table_index=None if duplicate else k,
                    table_x=None if duplicate else float(r['table_x']), x_over_c=float(r['plotted_x'])))
                if not duplicate:
                    k += 1
    return taps


def report():
    OUT.mkdir(exist_ok=True)
    inputs = Inputs()
    spec, states, controls = load_states(inputs)
    receipt = inputs.json(OUT/'residual_recheck.json')
    require(receipt['complete'] and len(receipt['states']) == 4, 'Full residual audit required')
    for name, digest in receipt['input_sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Residual receipt invalidated: '+name)
    checks = {r['job']: r for r in receipt['states']}
    jobs = sorted(spec['jobs'], key=lambda j: j['target_cl'])
    tables = inputs.csv(SOURCE/'experimental_tables.csv')
    measured = [r for r in tables if r['configuration'] == 'free' and r['reynolds'] == '1000000'
                and r['prestall_branch'] == 'True' and r['drag_available'] == 'True']
    theory = [r for r in inputs.csv(SOURCE/'published_mses_drag.csv')
              if r['configuration'] == 'free' and r['reynolds'] == '1000000']
    taps = pressure_taps(inputs.csv(SOURCE/'free_boundary_pressure.csv'), jobs)
    mses_pressure = inputs.csv(GEOMETRY/'pressure/mses_2010_envelope.csv')
    unsupported = {(r['element'], r['surface'], int(r['index']))
        for r in inputs.csv(GEOMETRY/'pressure_tap_geometry_support.csv') if int(r['unsupported_center'])}
    require(len(unsupported) == 3, 'Unsupported source tap set changed')
    targets, layers, separations, forces, wake_rows, histories = [], [], [], [], [], []
    for job in jobs:
        name, c = job['name'], job['pressure_case']
        state = states[name]
        result = state['result']
        cl = job['target_cl']
        require(checks[name]['max_residual'] <= spec['residual_tolerance'], 'Residual failure')
        context = {}
        for method, rows in [('experiment', measured), ('mses', theory)]:
            cd = reference_at(rows, cl)
            require(cd is not None or method == 'mses' and not job['mses_drag_available'], 'Required drag reference missing')
            rounding = [reference_at(rows, cl+d) for d in (-.0005, 0, .0005)]
            context.update({method+'_CD': cd, method+'_CD_at_actual_model_CL': reference_at(rows, result['CL']),
                method+'_rounding_CD_low': rounding[0], method+'_rounding_CD_high': rounding[2],
                method+'_rounding_fully_supported': all(v is not None for v in rounding)})
        near = min(measured, key=lambda r: abs(float(r['CL'])-cl))
        require(abs(Decimal(near['CL'])-Decimal(str(cl))) <= Decimal('.0005'),
                'Caption/table row identity changed')
        target = dict(job=name, case=c, figure=job['pressure_figure'], target_cl=cl,
            CL=result['CL'], CD=result['CD'], alpha_deg=result['alpha_deg'],
            experiment_caption_alpha=job['experiment_alpha'],
            appendix_row_index=int(near['index']), appendix_CL=float(near['CL']), appendix_CD=float(near['CD']),
            appendix_alpha=float(near['alpha_deg']), caption_minus_appendix_CL=cl-float(near['CL']),
            error_CD_counts=(result['CD']-context['experiment_CD'])*1e4,
            mses_error_CD_counts=None if context['mses_CD'] is None else (context['mses_CD']-context['experiment_CD'])*1e4,
            lift_error=result['CL']-cl, residual=state['residual'], recomputed_residual=checks[name]['max_residual'],
            selected_path=result['selected_path'], requested_evaluations=len(result['records']),
            failed_requests=sum(not r['admissible'] for r in result['records']),
            failed_internal_points=sum(not p['converged'] for d in state['metadata']['diagnostics'] for p in d['continuation']),
            wall_seconds=result['wall_seconds'], mses_pressure_available=job['mses_pressure_available'], **context)
        targets.append(target)
        histories.append(dict(job=name, records=result['records'], attempts=result['attempts']))
        for element in ('fore', 'aft'):
            gap = clearance(state['arrays'][element]['wake_xy'].T,
                state['arrays']['aft' if element == 'fore' else 'fore']['xy'].T)
            require(gap > 0, 'Wake intersects other element')
            forces.append(dict(case=c, target_cl=cl, element=element, other_body_clearance_c=gap,
                               **result['element_forces'][element]))
            for side in ('upper', 'lower'):
                bl = result['boundary_layers'][element][side]
                flow = side_data(state['surfaces'][element], side)
                tr = bl['transition']
                layers.append(dict(case=c, target_cl=cl, element=element, side=side,
                    transition_x=None if tr is None else tr['xy_over_cref'][0],
                    min_cf=bl['min_cf'], max_Hk=bl['max_Hk'], separation_regions=len(bl['separation_regions']),
                    separated_fraction=bl['separated_arclength_fraction'],
                    TE_theta=float(flow['theta'][-1]), TE_delta_star=float(flow['delta_star'][-1]),
                    TE_cf=float(flow['cf'][-1]), TE_turbulent=bool(flow['turbulent'][-1])))
                for k, region in enumerate(bl['separation_regions']):
                    separations.append(dict(case=c, target_cl=cl, element=element, side=side, interval=k,
                        start_x=region['start']['xy_over_cref'][0], end_x=region['end']['xy_over_cref'][0],
                        reaches_TE=region['reaches_trailing_edge']))
            for k, row in enumerate(state['wakes'][element]):
                wake_rows.append(dict(case=c, target_cl=cl, element=element, index=k,
                    **{field: float(row[field]) for field in row.dtype.names}))

    pressures, metrics, intervals, peaks, excluded = [], [], [], [], []
    for job in jobs:
        c, state = job['pressure_case'], states[job['name']]
        for surface in SURFACES:
            element, side = surface.split('_')
            branch = pressure_branch(state['surfaces'][element], side)
            branches = {'multifoil': (branch['x_over_cref'], branch['Cp'])}
            if job['mses_pressure_available']:
                rows = [r for r in mses_pressure if r['case'] == c and r['surface'] == surface]
                require(len(rows) > 1, 'Missing declared MSES pressure curve')
                branches['mses'] = (np.array([float(r['x_over_c']) for r in rows]), np.array([float(r['Cp']) for r in rows]))
            for tap in [r for r in taps if r['case'] == c and r['surface'] == surface]:
                reason = common_support(tap, branches, unsupported)
                if reason:
                    excluded.append(dict(case=c, surface=surface, plot_index=int(tap['index']),
                        table_index=tap['table_index'], table_x=tap['table_x'], plotted_x=tap['x_over_c'], reason=reason))
                for coordinate in ('table', 'plotted'):
                    x = tap['table_x'] if coordinate == 'table' else tap['x_over_c']
                    for method, (bx, by) in branches.items():
                        cp = None if reason else float(interpolate(bx, by, [x])[0])
                        ex = float(tap['Cp'])
                        pressures.append(dict(case=c, target_cl=job['target_cl'], surface=surface,
                            plot_index=int(tap['index']), table_index=tap['table_index'], coordinate=coordinate,
                            x_over_c=x, method=method, experiment_Cp=ex, model_Cp=cp,
                            error_Cp=None if cp is None else cp-ex, included=not reason, exclusion=reason))
            for coordinate in ('table', 'plotted'):
                for method in branches:
                    rows = [r for r in pressures if r['case'] == c and r['surface'] == surface
                            and r['coordinate'] == coordinate and r['method'] == method]
                    for region, lo, hi in [('all', -np.inf, np.inf)]+REGIONS.get(surface, []):
                        chosen = [r for r in rows if r['x_over_c'] is not None and lo <= r['x_over_c'] < hi]
                        metrics.append(dict(case=c, target_cl=job['target_cl'], surface=surface,
                            coordinate=coordinate, method=method, region=region, **pressure_metric(chosen)))
                    if surface != 'aft_upper':
                        continue
                    selected = [r for r in rows if r['included'] and .80 <= r['x_over_c'] <= .98]
                    x = [r['x_over_c'] for r in selected]
                    actual = recovery_intervals(x, [r['model_Cp'] for r in selected])
                    observed = recovery_intervals(x, [r['experiment_Cp'] for r in selected])
                    for model, exp in zip(actual, observed):
                        intervals.append(dict(case=c, target_cl=job['target_cl'], coordinate=coordinate,
                            method=method, **model, experiment_delta_Cp=exp['delta_Cp'], experiment_slope=exp['slope'],
                            error_delta_Cp=model['delta_Cp']-exp['delta_Cp']))
                    for statistic, field in [('largest_rise', 'delta_Cp'), ('steepest_slope', 'slope')]:
                        mp, ep = max(actual, key=lambda r: r[field]), max(observed, key=lambda r: r[field])
                        peaks.append(dict(case=c, target_cl=job['target_cl'], coordinate=coordinate, method=method,
                            statistic=statistic, model_left_x=mp['left_x'], model_right_x=mp['right_x'], model_value=mp[field],
                            experiment_left_x=ep['left_x'], experiment_right_x=ep['right_x'], experiment_value=ep[field]))
    exports = [('targets', targets), ('boundary_layers', layers), ('separation_intervals', separations),
        ('element_forces', forces), ('wake_profiles', wake_rows), ('pressure_at_taps', pressures),
        ('pressure_metrics', metrics), ('recovery_intervals', intervals), ('recovery_peaks', peaks), ('excluded_taps', excluded)]
    for filename, rows in exports:
        csv_write(OUT/(filename+'.csv'), rows)
    summary = dict(required_targets=4, converged_targets=4, full_state_audit_pass=True,
        full_residual_recomputation_pass=True, new_flow_solves=0, exact_repeat_pass=controls['passed'],
        exact_saved_array_comparisons=len(controls['array_comparisons']),
        retrieved_files_verified=238, frozen_inputs_verified=len(spec['sha256']),
        maximum_lift_error=max(abs(r['lift_error']) for r in targets),
        maximum_recomputed_residual=max(r['recomputed_residual'] for r in targets),
        requested_evaluations=sum(r['requested_evaluations'] for r in targets),
        failed_requests=sum(r['failed_requests'] for r in targets),
        failed_internal_points=sum(r['failed_internal_points'] for r in targets),
        selected_paths=sorted({r['selected_path'] for r in targets}),
        minimum_wake_body_clearance_c=min(r['other_body_clearance_c'] for r in forces),
        raw_pressure_markers=len(taps), supported_physical_pressure_observations=len(taps)-len(excluded),
        duplicated_endpoint_observations=sum(r['reason'] == 'duplicate_endpoint' for r in excluded),
        unsupported_source_observations=sum(r['reason'] == 'unsupported_source_center' for r in excluded),
        other_support_exclusions=sum(r['reason'] == 'outside_common_curve_support' for r in excluded),
        mses_drag_cases=[j['pressure_case'] for j in jobs if j['mses_drag_available']],
        mses_pressure_cases=[j['pressure_case'] for j in jobs if j['mses_pressure_available']],
        targets=targets, largest_sampled_rises=[r for r in peaks if r['coordinate'] == 'table' and r['statistic'] == 'largest_rise'],
        exact_tested_geometry_established=False, numerical_convergence_established=False,
        qualification='Four fixed-parameter states; no causal intervention, fitted threshold, measured bubble endpoints or new MSES run.')
    require(summary['duplicated_endpoint_observations'] == 16 and summary['unsupported_source_observations'] == 12,
            'Expected physical tap exclusions changed')
    write(OUT/'summary.json', summary)
    write(OUT/'controls.json', controls)
    write(OUT/'initialization_histories.json', histories)
    diagnosis = []
    for target in targets:
        c = target['case']
        fore = next(r for r in layers if r['case'] == c and r['element'] == 'fore' and r['side'] == 'lower')
        aft = next(r for r in layers if r['case'] == c and r['element'] == 'aft' and r['side'] == 'upper')
        regions = [r for r in separations if r['case'] == c and r['element'] == 'aft' and r['side'] == 'upper']
        require(len(regions) == 1 and not regions[0]['reaches_TE'], 'Aft bubble summary needs explicit multiple-interval handling')
        selected_metrics = [r for r in metrics if r['case'] == c and r['surface'] in ('fore_lower', 'aft_upper')]
        coordinate_changes = [abs(r['rmse']-next(s['rmse'] for s in selected_metrics
            if s['surface'] == r['surface'] and s['region'] == r['region'] and s['method'] == r['method'] and s['coordinate'] == 'table'))
            for r in selected_metrics if r['coordinate'] == 'plotted']
        diagnosis.append(dict(case=c, target_cl=target['target_cl'], fore_lower_transition_x=fore['transition_x'],
            aft_upper_transition_x=aft['transition_x'], aft_upper_separation_x=regions[0]['start_x'],
            aft_upper_reattachment_x=regions[0]['end_x'],
            aft_upper_separated_x_extent=regions[0]['end_x']-regions[0]['start_x'],
            maximum_fore_lower_aft_upper_RMSE_coordinate_change=max(coordinate_changes)))
    contrasts = []
    for left, right in zip(targets, targets[1:]):
        row = dict(from_CL=left['target_cl'], to_CL=right['target_cl'],
            model_delta_CD_counts=(right['CD']-left['CD'])*1e4,
            experiment_delta_CD_counts=(right['experiment_CD']-left['experiment_CD'])*1e4,
            mses_delta_CD_counts=None if left['mses_CD'] is None or right['mses_CD'] is None
                else (right['mses_CD']-left['mses_CD'])*1e4)
        for element in ('fore', 'aft'):
            values = [next(r for r in forces if r['case'] == t['case'] and r['element'] == element) for t in (left, right)]
            row[element+'_delta_CD_counts'] = (values[1]['CD_wake']-values[0]['CD_wake'])*1e4
        contrasts.append(row)
    write(OUT/'diagnosis.json', dict(cases=diagnosis, adjacent_lift_contrasts=contrasts,
        qualification='Element contrasts decompose computed drag changes, not experimental drag error or a controlled wake effect.'))
    plot(jobs, states, taps, mses_pressure, targets, measured, theory, layers, separations, forces, intervals)
    for name in ('validation/report_s414_lower_bucket.py', 'validation/report_s414_reconstruction.py',
                 'validation/report_s702_matched_pressure.py'):
        inputs.add(ROOT/name)
    write(OUT/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=dict(inputs.files),
        analysis_runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})))
    print({k: v for k, v in summary.items() if k not in ('targets', 'largest_sampled_rises')})


def plot(jobs, states, taps, mses_pressure, targets, measured, theory, layers, separations, forces, intervals):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    def cp(ax, job, surface):
        c = job['pressure_case']
        element, side = surface.split('_')
        b = pressure_branch(states[job['name']]['surfaces'][element], side)
        ax.plot(b['x_over_cref'], b['Cp'], color=COLORS['multifoil'], label='MultiFoil N9')
        rows = [r for r in taps if r['case'] == c and r['surface'] == surface and not r['endpoint_duplicate']]
        ax.plot([r['table_x'] for r in rows], [float(r['Cp']) for r in rows], 'o', color=COLORS['experiment'],
                ms=4, mfc='white', label='Experiment (Table-II x)')
        rows = [r for r in mses_pressure if r['case'] == c and r['surface'] == surface]
        if job['mses_pressure_available']:
            ax.plot([float(r['x_over_c']) for r in rows], [float(r['Cp']) for r in rows],
                    color=COLORS['mses'], ls='--', label='Published MSES (CL 0.28)')
        ax.invert_yaxis()
        ax.grid(alpha=.2)
        ax.set(xlabel='x/c (total reference chord)', ylabel='Cp')

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), layout='constrained')
    for k, job in enumerate(jobs):
        for row, surface in enumerate(('fore_lower', 'aft_upper')):
            ax = axes[row, k]
            cp(ax, job, surface)
            ax.set_title(f"CL {job['target_cl']:.3f} · {surface.replace('_', ' ')}")
            ax.set_xlim((0, .81) if row == 0 else (.80, .98))
            ax.set_ylim((.55, -2.45) if row == 0 else (.2, -.9))
            layer = next(r for r in layers if r['case'] == job['pressure_case'] and r['element']+'_'+r['side'] == surface)
            tr = layer['transition_x']
            if tr is not None:
                ax.axvline(tr, color='#b45a00', ls=':', lw=1.3)
            note = 'Model transition: none before TE' if tr is None else f'Model transition x/c = {tr:.4f}'
            ax.text(.02, .04, note, transform=ax.transAxes, fontsize=8,
                    bbox=dict(facecolor='white', alpha=.8, edgecolor='none'))
            for r in separations:
                if r['case'] == job['pressure_case'] and r['element']+'_'+r['side'] == surface:
                    ax.axvspan(r['start_x'], r['end_x'], color=COLORS['multifoil'], alpha=.08)
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[1, 3].legend(fontsize=8, loc='upper left')
    fig.suptitle('S414 free-transition sequence · Re 1 million · M 0.10 · fixed nominal geometry\n'
                 'Shading: computed Cf < 0; dotted line: computed transition. MSES Cp available only at the final case.')
    save(fig, OUT, 'fore_transition_and_aft_recovery')

    fig, axes = plt.subplots(4, 4, figsize=(16, 13), layout='constrained')
    for k, job in enumerate(jobs):
        for row, surface in enumerate(SURFACES):
            cp(axes[row, k], job, surface)
            axes[row, k].set_title(f"CL {job['target_cl']:.3f} · {surface.replace('_', ' ')}")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle('All surface pressures · original source markers retained; duplicate endpoints omitted from plotting\n'
                 'Three unsupported tap centers per case appear for context and are excluded from scores.')
    save(fig, OUT, 'full_pressure_comparison')

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), layout='constrained')
    for ax, job in zip(axes, jobs):
        cp(ax, job, 'fore_lower')
        ax.set(xlim=(0, .06), ylim=(1.05, -2.2), title=f"Fore lower · CL {job['target_cl']:.3f}")
        tr = next(r['transition_x'] for r in layers if r['case'] == job['pressure_case'] and r['element'] == 'fore' and r['side'] == 'lower')
        if tr is not None and tr < .06:
            ax.axvline(tr, color='#b45a00', ls=':', label='Computed transition')
        ax.legend(fontsize=7, loc='lower right')
    fig.suptitle('Fore leading-edge pressure peak · same audited states and source taps\n'
                 'Dotted line: computed transition; no experimental transition coordinates are supplied.')
    save(fig, OUT, 'fore_leading_edge_pressure')

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    ax = axes[0, 0]
    for method, rows, label in [('experiment', measured, 'Experiment: appendix'), ('mses', theory, 'Published MSES: 2010 envelope')]:
        rows = [r for r in rows if float(r['CL']) <= .4]
        ax.plot([float(r['CL']) for r in rows], [float(r['CD'])*1e4 for r in rows], 'o-' if method == 'experiment' else '--',
                color=COLORS[method], ms=4, label=label)
    cls = [r['target_cl'] for r in targets]
    ax.plot(cls, [r['CD']*1e4 for r in targets], 's:', color=COLORS['multifoil'], label='MultiFoil: four states')
    ax.set(xlabel='CL', ylabel='Total CD (drag counts)', title='Matched-lift drag · lines connect samples')
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    ax.plot(cls, [r['error_CD_counts'] for r in targets], 's-', color=COLORS['multifoil'], label='MultiFoil − experiment')
    valid = [r for r in targets if r['mses_CD'] is not None]
    ax.plot([r['target_cl'] for r in valid], [r['mses_error_CD_counts'] for r in valid], 'D--', color=COLORS['mses'], label='MSES − experiment')
    ax.axhline(0, color='gray', lw=.8)
    ax.set(xlabel='CL', ylabel='Drag error (counts)', title='MSES drag missing at CL −0.072')
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    for element, color in [('fore', '#2364aa'), ('aft', '#b45a00')]:
        rows = [r for r in forces if r['element'] == element]
        ax.plot(cls, [r['CD_wake']*1e4 for r in rows], 'o-', color=color, label=element+' wake CD')
        ax.plot(cls, [r['CD_friction']*1e4 for r in rows], '--', color=color, alpha=.7, label=element+' friction diagnostic')
    ax.set(xlabel='CL', ylabel='Modeled contribution (counts)', title='Model decomposition; no measured element drag')
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    state = states[jobs[-1]['name']]
    for element in ('fore', 'aft'):
        xy = state['arrays'][element]['xy']
        ax.plot(xy[0], xy[1], label=element)
        wake = state['arrays'][element]['wake_xy']
        mask = wake[0] <= 1.05
        ax.plot(wake[0, mask], wake[1, mask], ':', color='gray')
    ax.set_aspect('equal', adjustable='box')
    ax.set(xlim=(-.03, 1.05), xlabel='x/c', ylabel='y/c', title='Frozen nominal reconstructed geometry')
    ax.text(.02, .03, 'Exact tested contours and grid convergence remain unqualified.', transform=ax.transAxes, fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle('S414 lower-bucket drag diagnostic · Ncrit 9 throughout · one drag count = 0.0001 CD')
    save(fig, OUT, 'geometry_and_drag')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for k, job in enumerate(jobs):
        ax = axes.flat[k]
        rows = [r for r in intervals if r['case'] == job['pressure_case'] and r['coordinate'] == 'table' and r['method'] == 'multifoil']
        x = [(r['left_x']+r['right_x'])/2 for r in rows]
        ax.plot(x, [r['experiment_delta_Cp'] for r in rows], 'o-', color=COLORS['experiment'], label='Experiment')
        ax.plot(x, [r['delta_Cp'] for r in rows], 's-', color=COLORS['multifoil'], label='MultiFoil')
        other = [r for r in intervals if r['case'] == job['pressure_case'] and r['coordinate'] == 'table' and r['method'] == 'mses']
        if other:
            ax.plot(x, [r['delta_Cp'] for r in other], 'D--', color=COLORS['mses'], label='Published MSES')
        ax.set(title=f"CL {job['target_cl']:.3f}", xlabel='Midpoint between successive taps, x/c', ylabel='Cp(right) − Cp(left)', ylim=(-.05, .39))
        ax.axhline(0, color='gray', lw=.7)
        ax.grid(alpha=.2)
    axes.flat[-1].legend(fontsize=8)
    fig.suptitle('Aft upper recovery · every supported adjacent interval in 0.80 ≤ x/c ≤ 0.98\n'
                 'Largest rise is a tap-sampling statistic, not a measured bubble endpoint.')
    save(fig, OUT, 'pressure_recovery_intervals')

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    palette = plt.cm.viridis(np.linspace(.1, .85, 4))
    for job, color in zip(jobs, palette):
        state = states[job['name']]
        label = f"CL {job['target_cl']:.3f}"
        fore = side_data(state['surfaces']['fore'], 'lower')
        aft = side_data(state['surfaces']['aft'], 'upper')
        axes[0, 0].plot(fore['x_over_cref'], fore['cf'], color=color, label=label)
        axes[0, 1].plot(aft['x_over_cref'], aft['cf'], color=color, label=label)
        wake = state['wakes']['fore']
        axes[1, 0].plot(wake['x_over_cref'], wake['delta_star'], color=color, label=label)
        axes[1, 1].plot(wake['x_over_cref'], wake['theta'], color=color, label=label)
    axes[0, 0].set(xlim=(0, .81), title='Fore lower computed skin friction', ylabel='Cf')
    axes[0, 1].set(xlim=(.80, .98), title='Aft upper computed separation and recovery', ylabel='Cf')
    axes[1, 0].set(xlim=(.79, 1.05), title='Fore-wake displacement thickness', ylabel='δ*/c')
    axes[1, 1].set(xlim=(.79, 1.05), title='Fore-wake momentum thickness', ylabel='θ/c')
    for ax in axes.flat:
        ax.set_xlabel('x/c (total chord)')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    for ax in axes[0]:
        ax.axhline(0, color='gray', lw=.8)
    fig.suptitle('Modeled boundary-layer and fore-wake states\n'
                 'Integral wake thicknesses are neither turbulence intensity nor measured wake width.')
    save(fig, OUT, 'boundary_layers_and_fore_wake')


if __name__ == '__main__':
    report()
