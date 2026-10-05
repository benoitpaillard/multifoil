#!/usr/bin/env python3
"""Report incomplete reconstruction coverage without inventing a missing force."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_s414_fore_reconstruction import ROOT, STUDY, OUT, LAUNCH, Inputs, load_states, require, sha, write, close, array_close
from run_s414_fore_reconstruction import geometry_elements
from report_s702_wake_extent import csv_write, audit_surface, side_data
from report_s702_matched_pressure import pressure_branch, save
from report_s414_reconstruction import interpolate
from report_s414_fore_nose import growth_parts, transition_root
from report_s414_wake_transfer import nearest_polyline
from viscous import new_model, mf

PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
SOURCE = ROOT/'validation/results/s414_trip_qualification_20261003'


def raw_speed(uk, param):
    """Stable positive inverse of the solver's Karman–Tsien velocity mapping."""
    lam = param.KTl if param.Minf > 0 else 0.
    return 2*uk/(1-lam+np.sqrt((1-lam)**2+4*lam*(uk/param.Vinf)**2))


def paired_drag(a, b):
    # This guard also applies to converged requests that have missed target lift:
    # callers must pass complete fixed-lift result objects, not diagnostic records.
    require(a['feasible'] and b['feasible'], 'Incomplete target cannot supply paired drag')
    return (b['CD']-a['CD'])*1e4


def lower_history(model, surface, layer):
    """Independent growth accounting from saved exports; no nonlinear solve."""
    chosen = np.flatnonzero(surface['side'] == 'lower')
    nodes = len(surface)-1-chosen
    values = surface[chosen]
    param = mf.build_param(model, 0)
    U = np.array([values['theta'], values['delta_star'], values['amplification_or_shear'],
                  raw_speed(values['ue_over_vinf'], param)])
    s = model.foil.s[nodes[0]]-model.foil.s[nodes]
    require(np.all(np.diff(s) > 0), 'Non-increasing flow distance')
    array_close(mf.get_cp(U[3], param)[0], values['Cp'], 'Pressure inversion mismatch', 3e-12)
    array_close(np.array([mf.get_Hk(u, param)[0] for u in U.T]), values['Hk'], 'Exported Hk mismatch', 3e-12)
    flags = values['turbulent'].astype(int)
    require(set(flags) <= {0, 1} and np.all(np.diff(flags) >= 0), 'Invalid transition flags')
    count = int(np.where(flags == 1)[0][0]) if flags.any() else len(values)
    ss, states = list(s[:count]), list(U[:, :count].T)
    transition = layer['transition']; root_error = None
    if transition is not None:
        require(0 < count < len(values), 'Transition must be bracketed')
        root, root_error = transition_root(U[:, count-1:count+1], s[count-1:count+1], param)
        fraction = (root-s[count-1])/(s[count]-s[count-1])
        xy = np.array([values['x_over_cref'], values['y_over_cref']]).T
        point = (1-fraction)*xy[count-1]+fraction*xy[count]
        array_close(point, np.array(transition['xy_over_cref']), 'Independent transition position', 2e-10)
        ut = (1-fraction)*U[:, count-1]+fraction*U[:, count]; ut[2] = param.ncrit
        ss.append(root); states.append(ut)
    g = np.array([growth_parts(u, param) for u in states])
    integrals = np.diff(ss)[:, None]*(g[:-1]+g[1:])*.5
    balance = float(np.max(abs(np.diff(np.array(states)[:, 2])-integrals[:, 0])))
    require(balance < 1e-7, 'Amplification balance fails')
    rows = [dict(x=float(v['x_over_cref']), y=float(v['y_over_cref']), s=float(s[j]),
        Cp=float(v['Cp']), Cf=float(v['cf']), Hk=float(v['Hk']), theta=float(v['theta']),
        delta_star=float(v['delta_star']), laminar=not bool(flags[j]),
        N=float(U[2, j]) if not flags[j] else None) for j, v in enumerate(values)]
    physical = [r for j, r in enumerate(rows) if r['laminar'] and nodes[j] <= np.argmin(model.foil.x[0])]
    xx = np.array([r['x'] for r in physical]); nn = np.array([r['N'] for r in physical])
    nose = [r for r in physical if r['x'] <= .06]
    stats = dict(maximum_growth_balance=balance, transition_root_residual=root_error,
        base_N=float(integrals[:, 1].sum()), regularization_N=float(integrals[:, 2].sum()),
        N_at_x0034=float(interpolate(xx, nn, [.034])[0]) if xx[0] <= .034 <= xx[-1] else None,
        maximum_laminar_nose_Hk=max(r['Hk'] for r in nose),
        fore_transition_x=None if transition is None else transition['xy_over_cref'][0])
    return rows, stats


def audit_request_exports(folder, index, record, diag, models, inputs, tolerance):
    """Export-level audit of a converged off-target request, not a selected target."""
    require(record['converged'] and record['admissible'], 'Cannot interpret a failed flow')
    require(record['residual'] <= tolerance and diag['status'] == 'converged', 'Request residual/status fails')
    close(record['residual'], max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
    for field, term in (('CL', 'CL'), ('CD', 'CD_wake'), ('CM', 'CM')):
        close(record[field], sum(f[term] for f in diag['element_forces'].values()), context='Request force sum')
    surfaces = {}
    for name, model in models.items():
        surface = inputs.array(folder/f'flow/surface_{index:03d}_{name}.csv')
        array_close(np.array([surface['x_over_cref'], surface['y_over_cref']])[:, ::-1],
                    model.foil.x, 'Request geometry differs', 2e-12)
        for side in ('lower', 'upper'): audit_surface(side_data(surface, side), diag['boundary_layers'][name][side])
        param = mf.build_param(model, 0)
        array_close(mf.get_cp(raw_speed(surface['ue_over_vinf'], param), param)[0], surface['Cp'], 'Request Cp mismatch', 3e-12)
        wake = inputs.array(folder/f'flow/wake_{index:03d}_{name}.csv')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid request wake')
        last = wake[-1]
        cd = 2*last['theta']*(last['ue_over_vinf']/param.Vinf)**((5+last['delta_star']/last['theta'])/2)
        close(cd, diag['element_forces'][name]['CD_wake'], context='Request terminal wake drag')
        surfaces[name] = surface
    return surfaces


def main():
    inputs = Inputs(); spec, states, controls = load_states(inputs)
    residual = inputs.json(OUT/'residual_recheck.json')
    require(residual['complete'] and residual['selected_states_audited'] == 3, 'Selected-state audit incomplete')
    for name, digest in residual['input_sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Residual audit input changed')
    prior = inputs.json(PRIOR/'verification.json')
    require(sha(PRIOR/'verification.json') == '3325127f24a971f345f9857529aa250ed6e7b05d10e588e93b2472720e03489f', 'Reference identity changed')
    for name in ('targets.csv', 'pressure_at_taps.csv'):
        require(sha(inputs.add(PRIOR/name)) == prior['output_sha256'][name], 'Reference changed: '+name)
    reference = {r['case']: r for r in inputs.csv(PRIOR/'targets.csv')}
    taps = inputs.csv(PRIOR/'pressure_at_taps.csv')
    targets, histories, growth, requests, pressures = [], [], [], [], []
    selected_profiles, request_profiles = {}, {}
    for job in spec['jobs']:
        name = job['name']; state = states[name]; r = state['result']; ex = reference[job['pressure_case']]
        models = {n: new_model(n, xy, spec['reynolds'], spec['mach'], spec['ncrit'], (.25, 0), 1e-5)
                  for n, xy in geometry_elements(job['scenario'], spec['panels'])}
        row = dict(job=name, case=job['pressure_case'], scenario=job['scenario'], target_cl=job['target_cl'],
            feasible=r['feasible'], CL=r['CL'], CD=r['CD'], alpha_deg=r['alpha_deg'],
            experiment_CD=float(ex['experiment_CD']), mses_CD=float(ex['mses_CD']),
            error_CD_counts=None if not r['feasible'] else (r['CD']-float(ex['experiment_CD']))*1e4,
            fore_transition_x=None, fore_CD=None, aft_CD=None, selected_path=r.get('selected_path'),
            requests=len(r['records']), failed_requests=sum(not q['admissible'] for q in r['records']),
            wall_seconds=r['wall_seconds'], error=r['error'])
        if r['feasible']:
            profile, stats = lower_history(models['fore'], state['surfaces']['fore'], r['boundary_layers']['fore']['lower'])
            selected_profiles[name] = profile
            growth.append(dict(job=name, evaluation=r['selected_flow_index'], matched_target=True, **stats))
            histories.extend(dict(job=name, evaluation=r['selected_flow_index'], matched_target=True, **v) for v in profile)
            row.update(fore_transition_x=stats['fore_transition_x'], fore_CD=r['element_forces']['fore']['CD_wake'],
                       aft_CD=r['element_forces']['aft']['CD_wake'])
            for t in taps:
                if t['case'] != job['pressure_case'] or t['method'] != 'multifoil' or t['included'] != 'True': continue
                element, side = t['surface'].split('_'); branch = pressure_branch(state['surfaces'][element], side)
                cp = float(interpolate(branch['x_over_cref'], branch['Cp'], [float(t['x_over_c'])])[0])
                pressures.append(dict(job=name, **{k:t[k] for k in ('case','surface','table_index','coordinate')},
                    x=float(t['x_over_c']), experiment_Cp=float(t['experiment_Cp']), model_Cp=cp,
                    error_Cp=cp-float(t['experiment_Cp'])))
        else:
            expected_paths = ['cold_local', 'warm_anchor', 'warm_entry']
            require([a['path'] for a in r['attempts']] == expected_paths and not any(a['feasible'] for a in r['attempts']), 'Failed path accounting differs')
            for i, (q, d) in enumerate(zip(r['records'], state['metadata']['diagnostics'])):
                item = dict(job=name, evaluation=i, path=q['initialization']['path'], phase=q['initialization']['phase'],
                    alpha_deg=q['alpha_deg'], CL=q['CL'], CD=q['CD'], residual=q['residual'], admissible=q['admissible'],
                    fore_transition_x=None, fore_CD=None, aft_CD=None, N_at_x0034=None,
                    regularization_N=None, error=q['failure_reason'], selected_target=False)
                if q['admissible']:
                    surfaces = audit_request_exports(STUDY/name, i, q, d, models, inputs, spec['residual_tolerance'])
                    profile, stats = lower_history(models['fore'], surfaces['fore'], d['boundary_layers']['fore']['lower'])
                    request_profiles[i] = profile
                    growth.append(dict(job=name, evaluation=i, matched_target=False, **stats))
                    histories.extend(dict(job=name, evaluation=i, matched_target=False, **v) for v in profile)
                    item.update({k: stats[k] for k in ('fore_transition_x', 'N_at_x0034', 'regularization_N')})
                    item.update(fore_CD=d['element_forces']['fore']['CD_wake'], aft_CD=d['element_forces']['aft']['CD_wake'])
                else:
                    require(all(q[k] is None for k in ('CL', 'CD', 'CM')), 'Failed request has a force')
                requests.append(item)
        targets.append(row)
    # The original nominal search contains exactly the same bracket-end angles.
    # Admit these as same-angle context, never as matched-lift target replacements.
    nominal_endpoint_profiles, endpoint_rows = {}, []
    state = states['cl0175_nominal']
    models = {n: new_model(n, xy, spec['reynolds'], spec['mach'], spec['ncrit'], (.25, 0), 1e-5)
              for n, xy in geometry_elements('nominal', spec['panels'])}
    for i in (0, 1):
        q = state['result']['records'][i]; d = state['metadata']['diagnostics'][i]
        surfaces = audit_request_exports(STUDY/'cl0175_nominal', i, q, d, models, inputs, spec['residual_tolerance'])
        profile, stats = lower_history(models['fore'], surfaces['fore'], d['boundary_layers']['fore']['lower'])
        nominal_endpoint_profiles[i] = profile
        growth.append(dict(job='cl0175_nominal', evaluation=i, matched_target=False, **stats))
        histories.extend(dict(job='cl0175_nominal', evaluation=i, matched_target=False, **v) for v in profile)
        cubic = next(r for r in requests if r['evaluation'] == i)
        require(q['alpha_deg'] == cubic['alpha_deg'] and q['initialization']['path'] == cubic['path'] == 'cold_local',
                'Endpoint angle or requested initialization differs')
        endpoint_rows.append(dict(evaluation=i, alpha_deg=q['alpha_deg'], nominal_CL=q['CL'], cubic_CL=cubic['CL'],
            nominal_CD=q['CD'], cubic_CD=cubic['CD'], delta_CD_counts=(cubic['CD']-q['CD'])*1e4,
            nominal_fore_transition_x=stats['fore_transition_x'], cubic_fore_transition_x=cubic['fore_transition_x'],
            nominal_nose_N=stats['N_at_x0034'], cubic_nose_N=cubic['N_at_x0034'],
            fore_delta_CD_counts=(cubic['fore_CD']-d['element_forces']['fore']['CD_wake'])*1e4,
            aft_delta_CD_counts=(cubic['aft_CD']-d['element_forces']['aft']['CD_wake'])*1e4,
            comparison_basis='same_angle_not_matched_lift'))
    a, b = states['cl0056_nominal'], states['cl0056_fore_source_cubic']
    require(np.array_equal(a['arrays']['aft']['xy'], b['arrays']['aft']['xy']), 'Aft geometry changed')
    branch0 = pressure_branch(a['surfaces']['fore'], 'lower'); branch1 = pressure_branch(b['surfaces']['fore'], 'lower')
    lo = max(branch0['x_over_cref'][0], branch1['x_over_cref'][0])
    xx = np.unique(np.r_[branch0['x_over_cref'], branch1['x_over_cref']]); xx = xx[(xx >= lo) & (xx <= .06)]
    cp0 = interpolate(branch0['x_over_cref'], branch0['Cp'], xx); cp1 = interpolate(branch1['x_over_cref'], branch1['Cp'], xx)
    comparison = dict(delta_CD_counts=paired_drag(a['result'], b['result']),
        fore_delta_CD_counts=(b['result']['element_forces']['fore']['CD_wake']-a['result']['element_forces']['fore']['CD_wake'])*1e4,
        aft_delta_CD_counts=(b['result']['element_forces']['aft']['CD_wake']-a['result']['element_forces']['aft']['CD_wake'])*1e4,
        maximum_nose_Cp_change=float(np.max(abs(cp1-cp0))),
        x_at_maximum_nose_Cp_change=float(xx[np.argmax(abs(cp1-cp0))]),
        nominal_minimum_nose_Cp=float(np.min(cp0)), cubic_minimum_nose_Cp=float(np.min(cp1)),
        same_selected_path=a['result']['selected_path'] == b['result']['selected_path'], aft_geometry_exact=True)
    xy0, xy1 = a['arrays']['fore']['xy'].T, b['arrays']['fore']['xy'].T
    comparison['nose_vertex_to_other_polyline_max_c'] = max(float(nearest_polyline(xy0[xy0[:,0] <= .06], xy1)[0].max()),
                                                         float(nearest_polyline(xy1[xy1[:,0] <= .06], xy0)[0].max()))
    comparison['nose_vertex_to_other_polyline_max_mm'] = comparison['nose_vertex_to_other_polyline_max_c']*457.2
    ref0, ref1 = reference['cl0056'], reference['cl0175']
    drops = dict(nominal=paired_drag(states['cl0056_nominal']['result'], states['cl0175_nominal']['result']),
        fore_source_cubic=None, experiment=(float(ref1['experiment_CD'])-float(ref0['experiment_CD']))*1e4,
        published_mses=(float(ref1['mses_CD'])-float(ref0['mses_CD']))*1e4)
    for name, rows in [('targets',targets),('fore_lower_stations',histories),('amplification_accounting',growth),
                       ('failed_target_requests',requests),('pressure_at_taps',pressures),
                       ('same_angle_endpoints',endpoint_rows)]: csv_write(OUT/(name+'.csv'), rows)
    summary = dict(analysis_complete=True, target_coverage_complete=False, required_targets=4, converged_targets=3,
        selected_states_freshly_audited=3, failed_target_converged_request_exports_audited=9,
        same_angle_nominal_request_exports_audited=2, same_angle_endpoints=endpoint_rows,
        requests=sum(r['requests'] for r in targets), failed_requests=sum(r['failed_requests'] for r in targets),
        exact_control_arrays=sum(len(r['array_comparisons']) for r in controls['jobs'].values()),
        controls=controls, targets=targets, low_lift_comparison=comparison, paired_drag_changes_counts=drops,
        maximum_fresh_residual=max(r['max_residual'] for r in residual['states']),
        maximum_growth_balance=max(r['maximum_growth_balance'] for r in growth),
        retrieved_files_verified=inputs.json(LAUNCH/'retrieval.json')['files_verified'],
        frozen_inputs_verified=len(spec['sha256']), new_nonlinear_flow_solves=0,
        exact_geometry_established=False, numerical_convergence_established=False, solution_uniqueness_established=False)
    write(OUT/'summary.json', summary)
    plots(targets, selected_profiles, request_profiles, requests, taps, endpoint_rows, nominal_endpoint_profiles)
    for name in ('validation/report_s414_fore_reconstruction.py','validation/report_s414_fore_nose.py',
                 'validation/report_s414_wake_transfer.py','validation/report_s414_reconstruction.py',
                 'validation/report_s702_matched_pressure.py','validation/report_s702_drag_mechanism.py',
                 'validation/report_design_initialization.py','validation/report_s702_wake_extent.py'):
        inputs.add(ROOT/name)
    write(OUT/'input_manifest.json', dict(files=inputs.files))
    print('Selected states: 3/4. Low-lift comparison:', comparison, flush=True)
    print('Paired drag changes:', drops, flush=True)


def plots(targets, selected, profiles, requests, taps, endpoints, endpoint_profiles):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    nominal = [r for r in targets if r['scenario'] == 'nominal']
    cubic = [r for r in targets if r['scenario'] != 'nominal' and r['feasible']]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout='constrained')
    for label, color, rows, key, marker in [('Experiment', '#222222', nominal, 'experiment_CD', 'o'),
        ('Published MSES (2010)', '#168a6e', nominal, 'mses_CD', '^'),
        ('MultiFoil nominal', '#2166ac', nominal, 'CD', 'o'),
        ('Fore cubic only', '#d6604d', cubic, 'CD', 's')]:
        axes[0].plot([r['target_cl'] for r in rows], [r[key]*1e4 for r in rows], marker=marker, color=color, label=label)
    for label, color, rows, marker in [('Nominal', '#2166ac', nominal, 'o'), ('Fore cubic only', '#d6604d', cubic, 's')]:
        axes[1].plot([r['target_cl'] for r in rows], [r['fore_transition_x'] for r in rows], marker=marker, color=color, label=label)
    for ax in axes:
        ax.set(xlabel='Target CL', xticks=[.056,.175], xlim=(.044,.187)); ax.grid(alpha=.2); ax.legend(fontsize=8)
        ax.text(.174, .87, 'Cubic CL 0.175 failed\nNo target value or paired drop', transform=ax.get_xaxis_transform(), ha='right', color='#a33a31', fontsize=9)
    axes[0].set(ylabel='CD [counts]', title='Low-lift drag error persists')
    axes[1].set(ylabel='Fore-lower transition x / total chord', title='Early transition remains near x/c = 0.022')
    fig.suptitle('Fore-only reconstruction: 3/4 targets converged; both baseline controls repeat exactly')
    save(fig, OUT, 'drag_and_transition')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), layout='constrained')
    for name, color, label in [('cl0056_nominal','#2166ac','Nominal'), ('cl0056_fore_source_cubic','#d6604d','Fore cubic only')]:
        rows = selected[name]; x = np.array([r['x'] for r in rows])
        geometric = np.arange(len(rows)) >= np.argmin(x)
        for ax, field in zip(axes, ('Cp','Hk','N')):
            y = np.array([np.nan if r[field] is None else r[field] for r in rows])
            mask = geometric & (x <= .06)
            if field in ('Hk','N'): mask &= np.array([r['laminar'] for r in rows])
            ax.plot(x[mask], y[mask], '.-', ms=3, color=color, label=label)
        last_lam = [r for r in rows if r['laminar']][-1]
        tr = next(r['fore_transition_x'] for r in targets if r['job'] == name)
        axes[2].plot([last_lam['x'], tr], [last_lam['N'], 9.], ':o', ms=4, mfc='white', color=color)
    chosen = [r for r in taps if r['case']=='cl0056' and r['surface']=='fore_lower' and r['coordinate']=='table'
              and r['included']=='True' and r['method']=='multifoil' and float(r['x_over_c']) <= .06]
    axes[0].plot([float(r['x_over_c']) for r in chosen], [float(r['experiment_Cp']) for r in chosen], 'ko', ms=5, label='Experiment taps')
    axes[0].invert_yaxis()
    for ax, label in zip(axes, ('Cp','Laminar Hk','Laminar amplification N')):
        ax.set(xlabel='x / total chord', ylabel=label, xlim=(0,.06)); ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[2].axhline(9, color='gray', ls='--'); axes[2].text(.035,9.1,'Ncrit = 9', fontsize=9)
    fig.suptitle('Matched CL 0.056: both reconstructions retain the nose event\nNo experimental transition coordinate or MSES nose-pressure curve is inferred')
    save(fig, OUT, 'low_lift_nose')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    palette = {0:'#9e5d00',1:'#168a6e',5:'#8b4bb0',6:'#d6604d'}
    chosen = [r for r in requests if r['evaluation'] in palette]
    for r in chosen:
        i = r['evaluation']; color = palette[i]; label = f'Request {i}: CL {r["CL"]:.6f}'
        for ax, field in zip((axes[0,0],axes[0,1],axes[1,0]), ('CL','CD','fore_transition_x')):
            ax.plot(r['alpha_deg'], r[field]*(1e4 if field=='CD' else 1), 'o', color=color, label=label)
        rows = [v for v in profiles[i] if v['laminar']]
        xx = [v['x'] for v in rows]+[r['fore_transition_x']]
        nn = [v['N'] for v in rows]+[9.]
        axes[1,1].plot(xx,nn,color=color,label=label)
    for ax, field in zip((axes[0,0],axes[0,1],axes[1,0]), ('nominal_CL','nominal_CD','nominal_fore_transition_x')):
        ax.plot([r['alpha_deg'] for r in endpoints], [r[field]*(1e4 if field=='nominal_CD' else 1) for r in endpoints],
                'kx', ms=8, mew=2, label='Nominal at the same angles')
    rows = [r for r in endpoint_profiles[0] if r['laminar']]
    axes[1,1].plot([r['x'] for r in rows]+[endpoints[0]['nominal_fore_transition_x']],
                   [r['N'] for r in rows]+[9.], 'k--', lw=1.5)
    axes[1,1].text(.11, .12, 'Dashed black: nominal at −2.56888°', transform=axes[1,1].transAxes, fontsize=9)
    axes[0,0].axhline(.175,color='gray',ls='--',label='Required CL = 0.175')
    axes[0,0].annotate('CL decreases as α increases\n(requests 5 → 6)', xy=(chosen[3]['alpha_deg'], chosen[3]['CL']),
        xytext=(-2.56,.1805), arrowprops=dict(arrowstyle='->', color='gray'), fontsize=8)
    for ax, label in zip((axes[0,0],axes[0,1],axes[1,0]), ('Achieved CL','CD [counts]','Fore-lower transition x/c')):
        ax.set(xlabel='Incidence [degrees]', ylabel=label, xlim=(-2.58,-2.46)); ax.grid(alpha=.2)
        ax.axvline(-2.5303895457,color='gray',ls=':',lw=1)
    axes[0,1].legend(fontsize=8)
    axes[1,1].set(xlabel='x / total chord', ylabel='Laminar amplification N', xlim=(0,.8))
    axes[1,1].axhline(9,color='gray',ls='--'); axes[1,1].grid(alpha=.2)
    fig.suptitle('Cubic CL 0.175 target not attained: transition-sensitive search interval\nConverged off-target requests only; dotted incidence: cold/entry failed, warm path converged')
    save(fig, OUT, 'failed_target_diagnostic')


if __name__ == '__main__': main()
