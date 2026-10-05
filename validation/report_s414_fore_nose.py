#!/usr/bin/env python3
"""Bounded geometry, pressure and amplification audit of archived S414 states."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import platform
import importlib.metadata
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_s414_lower_bucket import ROOT, STUDY, OUT as PRIOR, Inputs, load_states, sha, require, write
from report_s702_wake_extent import audit, csv_write, array_close
from report_s702_matched_pressure import pressure_branch, save
from report_s414_wake_transfer import nearest_polyline
from report_s414_reconstruction import interpolate
from s414_geometry import elements
from multifoil import read_dat
from viscous import CoupledSolver, new_model, mf

OUT = ROOT/'validation/results/s414_fore_nose_20261003'
QUAL = ROOT/'validation/results/s414_geometry_qualification_20261002'
OLD = ROOT/'validation/results/s414_reconstruction_20261002'
OLD_AUDIT = ROOT/'validation/results/s414_reconstruction_analysis_20261003'
PRIOR_HASH = '3325127f24a971f345f9857529aa250ed6e7b05d10e588e93b2472720e03489f'
CASES = ('cl0056_nominal', 'cl0175_nominal', 'cl028_control')
COLORS = ('#b45a00', '#2364aa', '#168a6e')


def laminar_indices(arrays):
    flags = arrays['turbulent'][arrays['lower']]
    require(set(np.unique(flags)) <= {0, 1}, 'Nonbinary turbulent flags')
    require(np.all(np.diff(flags.astype(int)) >= 0), 'Unexpected relaminarization')
    return arrays['lower'][flags == 0]  # Stored flags are integers, not booleans.


def growth_parts(U, param):
    total = float(mf.get_damp(U, param)[0])
    extra = float((1+np.tanh(5*(U[2]-param.ncrit)))*.001/U[0])
    return total, total-extra, extra


def transition_root(U, xi, param):
    """Independent bracketed scalar root of the implemented subcell equation."""
    left, right = xi
    def residual(s):
        q = (s-left)/(right-left)
        ut = (1-q)*U[:, 0]+q*U[:, 1]
        ut[2] = param.ncrit
        return param.ncrit-U[2, 0]-.5*(s-left)*(mf.get_damp(U[:, 0], param)[0]+mf.get_damp(ut, param)[0])
    require(residual(left)*residual(right) <= 0, 'Transition root not bracketed')
    value = brentq(residual, left, right, xtol=1e-14)
    return value, residual(value)


def rate_derivatives(U, param):
    analytic = mf.get_damp(U, param)[1]
    corrected = corrected_growth_derivative(U, param)
    rows = []
    for relative in (1e-5, 1e-6):
        for j in range(4):
            h = relative*max(abs(U[j]), 1e-6)
            a, b = U.copy(), U.copy()
            a[j] += h; b[j] -= h
            numerical = (mf.get_damp(a, param)[0]-mf.get_damp(b, param)[0])/(2*h)
            rows.append(dict(relative_step=relative, variable=('theta', 'delta_star', 'N', 'ue')[j],
                analytic=float(analytic[j]), finite_difference=float(numerical), diagnostic_corrected=float(corrected[j]),
                scaled_error=float(abs(numerical-analytic[j])/max(1, abs(analytic[j]))),
                corrected_scaled_error=float(abs(numerical-corrected[j])/max(1, abs(corrected[j])))))
    require(max(r['corrected_scaled_error'] for r in rows) < 2e-5,
            'Identified missing derivative term does not explain the discrepancy')
    return rows


def corrected_growth_derivative(U, param):
    """Diagnostic only: supply the missing derivative of 0.1 exp(-20/(Hk-1)).

    No function value, state or vendored implementation is changed.
    """
    derivative = mf.get_damp(U, param)[1].copy()
    hk, hk_u = mf.get_Hk(U, param)
    if hk < 1.05:
        hk, hk_u = 1.05, hk_u*0
    hmi = 1/(hk-1)
    ret = mf.get_Ret(U, param)[0]
    critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
    rn = (np.log10(ret)-(critical-.1))/.2
    if rn < 0: return derivative
    rf = 1 if rn >= 1 else 3*rn**2-2*rn**3
    da = .028*(hk-1)-.0345*np.exp(-(3.87*hmi-2.52)**2)
    missing_af_u = -np.exp(-20*hmi)*(-hmi**2*hk_u)
    return derivative+rf*da/U[0]*missing_af_u


def amplification(model, state, job):
    a = state['arrays']['fore']
    model.glob.U = a['U'].copy()
    model.glob.Nsys = a['U'].shape[1]
    model.vsol.turb = a['turbulent'].copy()
    model.vsol.Is = [a['lower'].copy(), a['upper'].copy(), np.arange(model.foil.N, model.glob.Nsys)]
    model.isol.xi = a['xi'].copy()
    param = mf.build_param(model, 0)
    ids, lam = a['lower'], laminar_indices(a)
    require(np.all(np.diff(a['xi'][ids]) > 0), 'Flow distance not increasing')
    surface = state['surfaces']['fore']
    cf = surface['cf'][::-1]
    stations, intervals, derivatives = [], [], []
    for k, j in enumerate(ids):
        u = a['U'][:, j]
        is_lam = a['turbulent'][j] == 0
        total, base, extra = growth_parts(u, param) if is_lam else (None, None, None)
        uk = mf.get_uk(u[3], param)[0]
        gradient = None
        if k:
            prev = ids[k-1]
            gradient = float((uk-mf.get_uk(a['U'][3, prev], param)[0])/(a['xi'][j]-a['xi'][prev]))
        stations.append(dict(case=job['pressure_case'], target_cl=job['target_cl'], flow_index=k, node=int(j),
            x=float(a['xy'][0, j]), y=float(a['xy'][1, j]), xi=float(a['xi'][j]),
            laminar=is_lam, N=float(u[2]) if is_lam else None,
            shear=float(u[2]) if not is_lam else None, theta=float(u[0]), delta_star=float(u[1]),
            Hk=float(mf.get_Hk(u, param)[0]), Re_theta=float(mf.get_Ret(u, param)[0]),
            uk=float(uk), Cp=float(mf.get_cp(u[3], param)[0]), Cf=float(cf[j]),
            interval_duk_dxi=gradient, dN_dxi=total, base_growth=base, regularization_growth=extra))
    pairs = [(int(j0), int(j1), None) for j0, j1 in zip(lam[:-1], lam[1:])]
    transition = None
    if len(lam) < len(ids):
        j0, j1 = ids[len(lam)-1:len(lam)+1]
        u, xi = a['U'][:, [j0, j1]], a['xi'][[j0, j1]]
        root, r = transition_root(u, xi, param)
        original = next(r for r in state['result']['transition_audit'] if r['element'] == 'fore' and r['side'] == 'lower')
        require(abs(root-original['transition_xi']) < 2e-11, 'Independent transition root differs')
        q = (root-xi[0])/(xi[1]-xi[0])
        ut = (1-q)*u[:, 0]+q*u[:, 1]; ut[2] = param.ncrit
        pairs.append((int(j0), int(j1), (root, ut, q)))
        transition = dict(xi=root, x=float((1-q)*a['xy'][0, j0]+q*a['xy'][0, j1]),
            root_residual=float(r), archived_xi=original['transition_xi'], fraction_of_cell=float(q),
            left_x=float(a['xy'][0, j0]), right_x=float(a['xy'][0, j1]),
            Hk=float(mf.get_Hk(ut, param)[0]), Re_theta=float(mf.get_Ret(ut, param)[0]))
    gq, gw = np.polynomial.legendre.leggauss(8)
    for j0, j1, subcell in pairs:
        u0, u1 = a['U'][:, j0].copy(), a['U'][:, j1].copy()
        s0, s1 = a['xi'][[j0, j1]]
        x0, x1 = a['xy'][0, [j0, j1]]
        if subcell is not None:
            s1, u1, q = subcell
            x1 = x0+q*(x1-x0)
        ds = s1-s0
        g0, b0, e0 = growth_parts(u0, param); g1, b1, e1 = growth_parts(u1, param)
        trapezoid = .5*(g0+g1)*ds
        gauss = .5*ds*sum(w*growth_parts((1-t)*u0+t*u1, param)[0]
                         for t, w in zip(.5*(gq+1), gw))
        balance = float(u1[2]-u0[2]-trapezoid)
        require(abs(balance) <= 1e-7, 'Laminar amplification balance failed')
        intervals.append(dict(case=job['pressure_case'], left_x=float(x0), right_x=float(x1),
            left_xi=float(s0), right_xi=float(s1), left_N=float(u0[2]), right_N=float(u1[2]),
            transition_cell=subcell is not None, delta_N=float(u1[2]-u0[2]), balance=balance,
            base_trapezoid=float(.5*(b0+b1)*ds), regularization_trapezoid=float(.5*(e0+e1)*ds),
            frozen_state_gauss8=float(gauss), gauss_minus_trapezoid=float(gauss-trapezoid)))
    active = [j for j in lam if growth_parts(a['U'][:, j], param)[1] > 1e-5]
    for j in sorted(set([active[0], active[len(active)//2], active[-1]])) if active else []:
        for row in rate_derivatives(a['U'][:, j], param):
            derivatives.append(dict(case=job['pressure_case'], node=int(j), x=float(a['xy'][0, j]), **row))
    saved_u = model.glob.U.copy(); saved_flags = model.vsol.turb.copy()
    try:
        last = mf.march_amplification(model, 0)
        difference = float(np.max(abs(model.glob.U[2, lam]-saved_u[2, lam])))
        require(last == len(lam)-1 and difference < 1e-6, 'Frozen-state amplification march differs')
    finally:
        model.glob.U[:] = saved_u; model.vsol.turb[:] = saved_flags
    require(np.array_equal(model.glob.U, a['U']) and np.array_equal(model.vsol.turb, a['turbulent']), 'Saved state mutated')
    geometric_lam = [r for r in stations if r['laminar'] and r['x'] >= a['xy'][0].min()
                     and r['node'] <= int(np.argmin(a['xy'][0]))]
    def n_at(x):
        xx = np.array([r['x'] for r in geometric_lam])
        if not xx[0] <= x <= xx[-1]: return None
        return float(interpolate(xx, [r['N'] for r in geometric_lam], [x])[0])
    summary = dict(case=job['pressure_case'], target_cl=job['target_cl'], laminar_stations=len(lam),
        transition=transition, first_base_growth_x=None if not active else float(a['xy'][0, active[0]]),
        maximum_laminar_nose_Hk=max(r['Hk'] for r in stations if r['laminar'] and r['x'] < .04),
        N_at_x0034=n_at(.034), N_at_x006=n_at(.06), maximum_laminar_N=float(a['U'][2, lam].max()),
        maximum_amplification_balance=max(abs(r['balance']) for r in intervals),
        total_regularization_N=sum(r['regularization_trapezoid'] for r in intervals),
        total_base_N=sum(r['base_trapezoid'] for r in intervals),
        frozen_state_gauss_minus_trapezoid=sum(r['gauss_minus_trapezoid'] for r in intervals),
        frozen_march_maximum_N_difference=difference, frozen_march_same_transition_cell=True,
        state_and_flags_restored_exactly=True,
        nose_laminar_intervals=sum(r['left_x'] >= 0 and r['right_x'] <= .04 for r in intervals))
    return stations, intervals, derivatives, summary


def geometry_audit(inputs, states):
    source = read_dat(inputs.add(QUAL/'report_test_fore.dat'))
    contours = {'source_polyline': source}
    for scenario in ('nominal', 'source_cubic'):
        raw = elements(scenario)[0][1]
        model = new_model('fore', raw, 1e6, .1, 9, (.25, 0.), 1e-5)
        contours[scenario+'_input'] = raw
        contours[scenario+'_solver'] = model.foil.x.T[::-1]
        if scenario == 'nominal':
            for state in states.values():
                array_close(model.foil.x, state['arrays']['fore']['xy'], 'Rebuilt fore geometry', 2e-12)
    nodes, turns, metrics = [], [], []
    for stage, xy in contours.items():
        for j, point in enumerate(xy):
            nodes.append(dict(stage=stage, node=j, x=float(point[0]), y=float(point[1])))
        nose = xy[xy[:, 0] <= .06]
        distance = nearest_polyline(nose, source)[0]
        lower = xy[int(np.argmin(xy[:, 0])):]
        segment = np.diff(lower, axis=0)
        lengths = np.linalg.norm(segment, axis=1)
        angles = np.unwrap(np.arctan2(segment[:, 1], segment[:, 0]))
        turning = np.diff(angles)
        for j, turn in enumerate(turning, 1):
            if lower[j, 0] <= .06:
                turns.append(dict(stage=stage, x=float(lower[j, 0]), y=float(lower[j, 1]),
                    left_spacing=float(lengths[j-1]), right_spacing=float(lengths[j]),
                    turning_degrees=float(np.rad2deg(turn)),
                    polygon_curvature=float(turn/(.5*(lengths[j-1]+lengths[j])))))
        metrics.append(dict(stage=stage, vertices=len(xy), nose_vertices=len(nose),
            minimum_x=float(xy[:, 0].min()), maximum_nose_vertex_distance_to_source_polyline=float(distance.max()),
            same_distance_mm=float(distance.max()*457.2),
            negative_lower_nose_turns=sum(r['turning_degrees'] < -.01 for r in turns if r['stage'] == stage),
            minimum_lower_nose_turn_degrees=min(r['turning_degrees'] for r in turns if r['stage'] == stage)))
    return contours, nodes, turns, metrics


def main():
    OUT.mkdir(exist_ok=True)
    inputs = Inputs()
    require(sha(inputs.add(PRIOR/'verification.json')) == PRIOR_HASH, 'Prior audit changed')
    for name, digest in inputs.json(PRIOR/'verification.json')['output_sha256'].items():
        require(sha(inputs.add(PRIOR/name)) == digest, 'Prior report changed: '+name)
    spec, all_states, controls = load_states(inputs)
    jobs = [next(j for j in spec['jobs'] if j['name'] == name) for name in CASES]
    states = {j['name']: all_states[j['name']] for j in jobs}
    contours, geometry_nodes, geometry_turns, geometry_metrics = geometry_audit(inputs, states)
    print('Frozen inputs and reconstructed nose geometry verified.', flush=True)
    # One unchanged linear inviscid operator; no wake building or nonlinear solve.
    solver = CoupledSolver(elements('nominal'), reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'])
    model = solver.models[0]
    array_close(model.foil.x, states[CASES[0]]['arrays']['fore']['xy'], 'Inviscid operator geometry', 2e-12)
    stations, intervals, derivatives, summaries, pressure = [], [], [], [], []
    taps = inputs.csv(PRIOR/'pressure_at_taps.csv')
    tap_rows = []
    for job in jobs:
        state = states[job['name']]
        alpha = state['result']['alpha_deg']
        gamma = model.isol.gamref @ np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
        cp_inv = mf.get_cp(gamma, model.param)[0][::-1]
        surface = state['surfaces']['fore']
        nose = int(np.argmin(surface['x_over_cref']))
        branch = pressure_branch(surface, 'lower')
        for k in range(nose, len(surface)):
            pressure.append(dict(case=job['pressure_case'], target_cl=job['target_cl'], x=float(surface['x_over_cref'][k]),
                y=float(surface['y_over_cref'][k]), inviscid_Cp=float(cp_inv[k]), viscous_Cp=float(surface['Cp'][k])))
        for tap in taps:
            if tap['case'] != job['pressure_case'] or tap['surface'] != 'fore_lower' or tap['method'] != 'multifoil' or tap['included'] != 'True': continue
            x = float(tap['x_over_c'])
            ci = float(interpolate(branch['x_over_cref'], cp_inv[nose:], [x])[0])
            tap_rows.append(dict(case=job['pressure_case'], coordinate=tap['coordinate'], table_index=int(tap['table_index']),
                x=x, experiment_Cp=float(tap['experiment_Cp']), viscous_Cp=float(tap['model_Cp']), inviscid_Cp=ci,
                viscous_error=float(tap['error_Cp']), inviscid_error=ci-float(tap['experiment_Cp'])))
        st, iv, de, su = amplification(model, state, job)
        stations.extend(st); intervals.extend(iv); derivatives.extend(de); summaries.append(su)
        print(job['pressure_case'], 'balances/root/march pass; original derivative criterion:',
              all(r['scaled_error'] < 2e-5 for r in de), flush=True)
    # Existing, already audited cubic-source control: no new sensitivity solve.
    old_manifest = inputs.json(OLD_AUDIT/'input_manifest.json')
    for name, digest in old_manifest['files'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Prior reconstruction input changed: '+name)
    old_spec = inputs.json(OLD/'study.json'); runtime = inputs.json(OLD/'runtime.json')
    cubic_job = next(j for j in old_spec['jobs'] if j['name'] == 'cl028_source_cubic')
    cubic = audit(OLD/cubic_job['name'], cubic_job, old_spec, runtime, inputs)
    cubic_branch = pressure_branch(cubic['surfaces']['fore'], 'lower')
    mses = [r for r in inputs.csv(QUAL/'pressure/mses_2010_envelope.csv')
            if r['case'] == 'cl028' and r['surface'] == 'fore_lower']
    cubic_rows = []
    for tap in [r for r in tap_rows if r['case'] == 'cl028']:
        value = float(interpolate(cubic_branch['x_over_cref'], cubic_branch['Cp'], [tap['x']])[0])
        published = float(interpolate([float(r['x_over_c']) for r in mses], [float(r['Cp']) for r in mses], [tap['x']])[0])
        cubic_rows.append(dict(tap, cubic_Cp=value, cubic_minus_nominal_Cp=value-tap['viscous_Cp'], mses_Cp=published))
    for name, rows in [('geometry_nodes', geometry_nodes), ('geometry_turning', geometry_turns), ('geometry_metrics', geometry_metrics),
        ('fore_lower_stations', stations), ('amplification_intervals', intervals), ('growth_derivative_checks', derivatives),
        ('fore_lower_pressure', pressure), ('pressure_at_taps', tap_rows), ('archived_cubic_control', cubic_rows)]:
        csv_write(OUT/(name+'.csv'), rows)
    summary = dict(complete=True, cases=summaries, geometry=geometry_metrics,
        unchanged_inviscid_operator_rebuilds=1, archived_incidence_inviscid_evaluations=3, new_nonlinear_flow_solves=0,
        amplification_balances_and_transition_roots_pass=True,
        original_growth_derivative_check_pass=max(r['scaled_error'] for r in derivatives) < 2e-5,
        maximum_original_growth_derivative_error=max(r['scaled_error'] for r in derivatives),
        diagnostic_corrected_derivative_check_pass=True,
        maximum_corrected_growth_derivative_error=max(r['corrected_scaled_error'] for r in derivatives),
        archived_cubic_control_CL=cubic['result']['CL'], archived_cubic_control_CD=cubic['result']['CD'],
        nominal_control_CD=states['cl028_control']['result']['CD'], ncrit=9,
        prior_exact_repeat=controls['passed'], exact_tested_geometry_established=False, numerical_convergence_established=False)
    write(OUT/'summary.json', summary)
    plot(jobs, contours, geometry_turns, stations, intervals, pressure, tap_rows, cubic_rows, summaries)
    for name in ('validation/report_s414_fore_nose.py', 'validation/S414_FORE_NOSE_PROTOCOL.md',
                 'validation/report_s414_wake_transfer.py', 'validation/report_s414_reconstruction.py'):
        inputs.add(ROOT/name)
    for name in ('source_manifest.json', 'author_source_check.json'):
        inputs.add(OUT/'sources'/name)
    write(OUT/'input_manifest.json', dict(files=inputs.files, prior_verification_sha256=PRIOR_HASH,
        analysis_runtime=dict(python=platform.python_version(), versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})))
    print('Report complete:', OUT, flush=True)


def plot(jobs, contours, turns, stations, intervals, pressure, taps, cubic, summaries):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    styles = {'source_polyline': ('#222222', 'o-'), 'nominal_input': ('#b45a00', ':'),
              'nominal_solver': ('#2364aa', '.-'), 'source_cubic_solver': ('#168a6e', '--')}
    for stage, (color, style) in styles.items():
        xy = contours[stage]
        for ax in axes[0]: ax.plot(xy[:, 0], xy[:, 1], style, color=color, ms=4, label=stage.replace('_', ' '))
        rows = [r for r in turns if r['stage'] == stage]
        axes[1, 0].plot([r['x'] for r in rows], [r['polygon_curvature'] for r in rows], style, color=color, ms=4)
        axes[1, 1].plot([r['x'] for r in rows], [r['right_spacing'] for r in rows], style, color=color, ms=4)
    axes[0, 0].set(xlim=(-.001, .04), ylim=(-.022, .04), title='Reported vertices and implemented nose')
    axes[0, 0].set_aspect('equal', adjustable='box')
    axes[0, 1].set(xlim=(-.0002, .006), ylim=(-.01, .01), title='Leading-edge detail')
    axes[0, 1].set_aspect('equal', adjustable='box')
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].set(xlim=(0, .04), title='Lower-surface polygon turning / arc length', ylabel='Discrete curvature × c')
    axes[1, 0].axhline(0, color='gray', lw=.8)
    axes[1, 1].set(xlim=(0, .06), title='Lower-surface station spacing', ylabel='Segment length / c', yscale='log')
    for ax in axes.flat: ax.grid(alpha=.2); ax.set_xlabel('x/c')
    fig.suptitle('S414 fore geometry pipeline · same 83 reported vertices\n'
                 'Polyline distance and discrete curvature are reconstruction diagnostics, not as-built measurements.')
    save(fig, OUT, 'geometry_and_sampling')

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    for k, job in enumerate(jobs):
        case = job['pressure_case']
        rows = [r for r in pressure if r['case'] == case]
        obs = [r for r in taps if r['case'] == case and r['coordinate'] == 'table']
        for ax in axes[:, k]:
            ax.plot([r['x'] for r in rows], [r['inviscid_Cp'] for r in rows], '--', color='#b45a00', label='Nominal inviscid, same α')
            ax.plot([r['x'] for r in rows], [r['viscous_Cp'] for r in rows], color='#2364aa', label='Saved nominal viscous')
            ax.plot([r['x'] for r in obs], [r['experiment_Cp'] for r in obs], 'o', ms=4, mfc='white', color='#222222', label='Experiment (Table-II x)')
            if case == 'cl028':
                cr = [r for r in cubic if r['coordinate'] == 'table']
                ax.plot([r['x'] for r in cr], [r['cubic_Cp'] for r in cr], 'D--', ms=3, color='#168a6e', label='Archived cubic-source control')
                ax.plot([r['x'] for r in cr], [r['mses_Cp'] for r in cr], ':', color='#8051a0', label='Published MSES, CL 0.28')
            ax.invert_yaxis(); ax.grid(alpha=.2); ax.set(xlabel='x/c', ylabel='Cp')
        axes[0, k].set(xlim=(0, .06), ylim=(1, -2.5), title=f"CL {job['target_cl']:.3f} · leading edge")
        axes[1, k].set(xlim=(0, .8), ylim=(.5, -2.5), title='Full fore lower surface')
    axes[0, 0].legend(fontsize=8); axes[0, 2].legend(fontsize=7)
    fig.suptitle('Pressure supplied to the boundary layer · unchanged nominal linear operator\n'
                 'Inviscid curves use archived viscous incidence, not matched inviscid lift; no new coupled-flow solve.')
    save(fig, OUT, 'pressure_input_diagnostic')

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    for job, color in zip(jobs, COLORS):
        rows = [r for r in stations if r['case'] == job['pressure_case']]
        lam = [r for r in rows if r['laminar']]
        label = f"CL {job['target_cl']:.3f}"
        for ax, field, chosen in [(axes[0, 0], 'Cp', rows), (axes[0, 1], 'Hk', lam),
            (axes[0, 2], 'Re_theta', lam), (axes[1, 0], 'dN_dxi', lam), (axes[1, 1], 'N', lam), (axes[1, 2], 'N', lam)]:
            ax.plot([r['x'] for r in chosen], [r[field] for r in chosen], '.-', ms=3, color=color, label=label)
        summary = next(s for s in summaries if s['case'] == job['pressure_case'])
        if summary['transition'] is not None:
            for ax in axes[1, 1:]: ax.plot(summary['transition']['x'], 9, '*', color=color, ms=11)
    for ax in axes.flat: ax.grid(alpha=.2); ax.set_xlabel('x/c'); ax.set_xlim(0, .06)
    axes[0, 0].invert_yaxis(); axes[0, 0].set(title='Computed pressure', ylabel='Cp')
    axes[0, 1].set(title='Laminar shape factor', ylabel='Hk')
    axes[0, 2].set(title='Laminar momentum-thickness Reynolds number', ylabel='Reθ', ylim=(0, 250))
    axes[1, 0].set(title='Laminar amplification rate', ylabel='dN/d(ξ/c)')
    axes[1, 1].set(title='Amplification through nose recovery', ylabel='N', ylim=(-.3, 9.7))
    axes[1, 2].set(title='Full fore-lower amplification history', ylabel='N', xlim=(0, .81), ylim=(-.3, 9.7))
    for ax in axes[1, 1:]: ax.axhline(9, color='gray', ls='--', lw=1)
    axes[0, 1].legend(fontsize=8)
    fig.suptitle('How the archived model reaches transition · Ncrit 9 unchanged\n'
                 'N and its rate stop at transition; downstream shear is never plotted as amplification. Stars mark subcell transition.')
    save(fig, OUT, 'amplification_mechanism')


if __name__ == '__main__': main()
