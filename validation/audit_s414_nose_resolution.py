#!/usr/bin/env python3
"""Audit the six-case local-resolution experiment without nonlinear solves.

Execution, evidence verification and scientific convergence remain separate.
Every converged request is reconstructed from saved state, including initialization
points; a failed solve never supplies a finite force or a selected comparison.
"""
import argparse
import csv
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from report_design_initialization import Inputs
from report_s702_drag_mechanism import require, close, audit_surface, side_data
from report_s702_wake_extent import array_close, thermo
from viscous import boundary_layer_summary, mf, te_info
from s414_nose_resolution import Solver, build_elements
from run_s414_transition_boundary import transition_audit, reached_by
import run_s414_nose_resolution as run

OUT = ROOT/'validation/results/s414_nose_resolution_analysis_20261005'
ELEMENTS = ('fore', 'aft')


def state_arrays(path, inputs):
    with np.load(inputs.add(path), allow_pickle=False) as saved:
        require(set(saved.files) == set(run.STATE_FIELDS), 'Saved state schema differs')
        result = {k: saved[k].copy() for k in saved.files}
    require(all(np.isfinite(v).all() for v in result.values()), 'Nonfinite saved state')
    return result


def audit_exports(folder, job, spec, runtime, inputs):
    result = inputs.json(folder/'result.json')
    require(result['complete'] and result['job'] == job['name'], 'Wrong/incomplete result')
    require(result['amplification_jacobian'] == job['amplification_jacobian'] and
            result['subdivision_factor'] == job['subdivision_factor'], 'Model/resolution differs')
    for name, digest in result['artifact_sha256'].items():
        require(run.sha(inputs.add(folder/name)) == digest, 'Worker artifact changed: '+name)
    if not result['records']:
        require(not result['feasible'] and result['execution_state'] == 'failed', 'Empty successful result')
        require(all(result[k] is None for k in ('CL', 'CD', 'CM', 'alpha_deg')), 'Failed job has force')
        return result, {}
    meta = inputs.json(folder/'flow/metadata.json'); polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(result['records']) == len(polar) == len(meta['diagnostics']) <= len(job['requests']), 'Record counts differ')
    require(len(result['records']) == len(job['requests']) or result['error'], 'Unexplained unevaluated requests')
    param = thermo(spec); states = {}; previous = None
    for index, (rec, p, diag, request) in enumerate(zip(result['records'], polar, meta['diagnostics'], job['requests'])):
        require(rec['index'] == index and rec['role'] == request['role'] and
                rec['requested_alpha_deg'] == request['alpha_deg'] == rec['alpha_deg'], 'Request history changed')
        require(rec['converged'] == bool(int(p['converged'])) == (diag['status'] == 'converged'), 'Convergence differs')
        for field in ('alpha_deg', 'CL', 'CD', 'CM', 'residual'):
            if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null field has finite export')
            else: close(rec[field], float(p[field]), context='record/export '+field)
        for field in ('mach', 'reynolds', 'ncrit'): close(diag[field], spec[field], context=field)
        require(diag['amplification_jacobian'] == job['amplification_jacobian'] and
                diag['ncrit_by_element'] == dict(fore=spec['ncrit'], aft=spec['ncrit']), 'Effective model differs')
        require(diag['nose_resolution'] == result['resolution_controls'], 'Resolution metadata differs')
        require(rec['iterations_final_step'] == diag['iterations'] == len(diag['history']) <= spec['iterations'], 'Iteration count differs')
        require(rec['solver_initialized_before'] == (previous is not None) and
                rec['previous_converged_alpha_deg'] == previous, 'Initialization state differs')
        require(rec['reached_by'] == reached_by(previous is not None, previous, request['alpha_deg'], rec['continuation']), 'Initialization description differs')
        require([(t['alpha_deg'], t['converged'], t['CL']) for t in rec['continuation']] ==
                [(t['alpha_deg'], t['converged'], t['CL']) for t in diag['continuation']], 'Continuation history differs')
        if not rec['converged']:
            require(all(rec[k] is None for k in ('CL', 'CD', 'CM')) and not rec['state_files'] and diag['failure_reason'], 'Failed evaluation has a force/state')
            require(not any((folder/f'flow/{kind}_{index:03d}_{e}.csv').exists()
                    for kind in ('surface', 'wake') for e in ELEMENTS), 'Failed evaluation has surface exports')
            continue
        previous = request['alpha_deg']
        require(rec['admissible'] and not diag['failure_reason'], 'Inadmissible converged record')
        require(rec['residual'] <= spec['residual_tolerance'], 'Recorded residual fails')
        close(rec['residual'], max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
        for key in ('element_forces', 'boundary_layers'): require(rec[key] == diag[key], 'Diagnostics differ')
        for field, term in (('CL', 'CL'), ('CD', 'CD_wake'), ('CM', 'CM')):
            close(rec[field], sum(f[term] for f in rec['element_forces'].values()), context='element sum '+field)
        arrays, surfaces = {}, {}
        for element in ELEMENTS:
            surface = inputs.array(folder/f'flow/surface_{index:03d}_{element}.csv')
            wake = inputs.array(folder/f'flow/wake_{index:03d}_{element}.csv')
            name = f'states/point_{index:02d}_{element}_state.npz'
            require(name in rec['state_files'], 'State file not recorded')
            a = state_arrays(folder/name, inputs); controls = result['resolution_controls'][element]; n = len(surface)
            require(n == controls['refined_nodes'] and len(wake) == controls['wake_nodes'] and
                    a['U'].shape == (4, n+len(wake)), 'State dimensions differ')
            for k, field in ((0, 'theta'), (1, 'delta_star'), (2, 'amplification_or_shear')):
                array_close(a['U'][k, :n][::-1], surface[field], 'Newton/surface '+field)
            for k, field in ((0, 'theta'), (1, 'delta_star')):
                array_close(a['U'][k, n:], wake[field], 'Newton/wake '+field)
            uk = mf.get_uk(a['U'][3], param)[0]
            array_close(uk[:n][::-1], surface['ue_over_vinf'], 'Surface corrected velocity')
            array_close(uk[n:], wake['ue_over_vinf'], 'Wake corrected velocity')
            hk = np.array([mf.get_Hk(a['U'][:, j], param)[0] for j in range(n+len(wake))])
            array_close(hk[:n][::-1], surface['Hk'], 'Surface Hk'); array_close(hk[n:], wake['Hk'], 'Wake Hk')
            require(np.array_equal(a['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence flags differ')
            labels = np.full(n, 'upper', dtype='<U5'); labels[a['lower']] = 'lower'
            require(np.array_equal(labels[::-1], surface['side']), 'Surface indices differ')
            require(set(a['lower']) | set(a['upper']) == set(range(n)) and not set(a['lower']) & set(a['upper']), 'Surface partition differs')
            for k, field in enumerate(('x_over_cref', 'y_over_cref')):
                array_close(a['xy'][k][::-1], surface[field], 'Surface geometry')
                array_close(a['wake_xy'][k], wake[field], 'Wake geometry')
            require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
            for side in ('upper', 'lower'): audit_surface(side_data(surface, side), rec['boundary_layers'][element][side])
            ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0)
            require(np.all(ds > 0), 'Nonpositive wake spacing')
            array_close(a['wake_s']-a['wake_s'][0], np.r_[0., np.cumsum(ds)], 'Wake arclength')
            array_close(np.linalg.norm(a['wake_t'], axis=0), np.ones(len(wake)), 'Wake tangent normalization')
            h = a['U'][1, n:]/a['U'][0, n:]
            sy = 2*a['U'][0, n:]*(uk[n:]/param.Vinf)**((5+h)/2)
            close(sy[-1], rec['element_forces'][element]['CD_wake'], context='terminal Squire–Young drag')
            arrays[element], surfaces[element] = a, surface
        states[index] = dict(record=rec, arrays=arrays, surfaces=surfaces)
    if result['feasible']:
        index = result['selected_flow_index']; rec = states[index]['record']
        require(index == len(job['requests'])-1 and rec['role'] == 'selected', 'Wrong selected request')
        for key in ('alpha_deg', 'CL', 'CD', 'CM', 'element_forces', 'boundary_layers', 'transition_audit'):
            require(result[key] == rec[key], 'Selected output differs: '+key)
        for element in ELEMENTS:
            final = state_arrays(folder/(element+'_state.npz'), inputs)
            require(all(np.array_equal(final[k], states[index]['arrays'][element][k]) for k in run.STATE_FIELDS), 'Final state differs from selected request')
    else:
        require(all(result[k] is None for k in ('CL', 'CD', 'CM', 'alpha_deg')) and
                'selected_flow_index' not in result, 'Failed target has selected force/state')
    return result, states


def recompute(job, spec, state, controls, remote_geometry=None):
    rec = state['record']; alpha = rec['alpha_deg']
    solver = Solver(build_elements(), subdivision_factor=job['subdivision_factor'], reynolds=spec['reynolds'],
        mach=spec['mach'], ncrit=spec['ncrit'], iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
        wake_length=spec['wake_length'], minimum_te_gap=spec['minimum_te_gap'], continuation_step=spec['continuation_step'],
        min_step=spec['min_step'], start_alpha=job['start_alpha'], amplification_jacobian=job['amplification_jacobian'])
    require(solver.resolution_controls == controls, 'Fresh geometry/resolution controls differ')
    if remote_geometry is not None:
        for model in solver.models:
            for label in ('x', 's', 't'):
                key = f"factor{job['subdivision_factor']}_{model.geom.name}_{label}"
                array_close(getattr(model.foil, label), remote_geometry[key], 'Actual remote preflight geometry '+key, 2e-12)
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = state['arrays'][model.geom.name]
        for field, actual, expected, tol in [('foil', model.foil.x, a['xy'], 2e-12),
                ('wake', model.wake.x, a['wake_xy'], 2e-11), ('wake tangent', model.wake.t, a['wake_t'], 2e-11),
                ('wake arclength', model.wake.s, a['wake_s'], 2e-11)]:
            array_close(actual, expected, 'Rebuilt '+field, tol)
        model.glob.U = a['U'].copy(); model.glob.Nsys = a['U'].shape[1]; model.vsol.turb = a['turbulent'].copy()
        model.isol.Istag = [int(a['lower'][0]), int(a['upper'][0])]
        model.isol.sgnue = np.ones(model.foil.N); model.isol.sgnue[a['lower']] = -1
        model.isol.sstag = float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model); mf.stagpoint_move(model)
        array_close(model.isol.xi, a['xi'], 'Restored station coordinates', 2e-12)
        for k, side in enumerate(('lower', 'upper')):
            require(np.array_equal(model.vsol.Is[k], a[side]), 'Restored branch indices differ')
        if abs(te_info(model.foil.x)[1]) < 1e-12: model.vsol.wgap = np.zeros(model.wake.N)
        else: mf.set_wake_gap(model)
    solver._build_source_coupling(); residual, _ = solver._system()
    blnorm = float(np.linalg.norm(residual[:3*solver.nsys])); coupling = float(np.max(abs(residual[3*solver.nsys:])))
    require(max(blnorm, coupling) <= spec['residual_tolerance'], 'Fresh Newton residual fails')
    roots = transition_audit(solver); old = rec['transition_audit']; shift = 0.
    require(len(roots) == len(old), 'Transition root count differs')
    for new, saved in zip(roots, old):
        require((new['element'], new['side']) == (saved['element'], saved['side']) and new['inside_cell'] and
                abs(new['amplification_residual']) <= 1e-7 and new['station_residual_max'] <= 1e-7, 'Fresh transition root fails')
        shift = max(shift, abs(new['transition_xi']-saved['transition_xi']))
    require(shift <= 1e-9, 'Fresh transition root moved')
    rows = []
    for model in solver.models:
        e = model.geom.name; surface = state['surfaces'][e]
        array_close(model.glob.U, state['arrays'][e]['U'], 'Audit changed state', 0.)
        mf.calc_force(model); mf.get_distributions(model)
        if model.oper.Ma > 0:
            for side in range(3):
                param = mf.build_param(model, side)
                for j in model.vsol.Is[side]:
                    rho, _ = mf.get_rho(model.glob.U[:, j], param); model.post.cf[j] *= rho/model.oper.rho
        array_close(model.post.cp[:model.foil.N][::-1], surface['Cp'], 'Fresh Cp', 2e-12)
        array_close(model.post.cf[:model.foil.N][::-1], surface['cf'], 'Fresh Cf', 2e-12)
        for name, field in (('CL', 'cl'), ('CM', 'cm'), ('CD_wake', 'cd'), ('CD_friction', 'cdf')):
            close(getattr(model.post, field), rec['element_forces'][e][name], context='Fresh '+e+' '+name)
        layers = boundary_layer_summary(model)
        for side in ('lower', 'upper'):
            actual, expected = layers[side], rec['boundary_layers'][e][side]
            require((actual['transition'] is None) == (expected['transition'] is None), 'Transition presence differs')
            if actual['transition']: close(actual['transition']['s_over_cref'], expected['transition']['s_over_cref'])
            for key in ('min_cf', 'max_Hk', 'separated_arclength_fraction'): close(actual[key], expected[key])
            require(len(actual['separation_regions']) == len(expected['separation_regions']), 'Separation count differs')
            for a, b in zip(actual['separation_regions'], expected['separation_regions']):
                require(a['reaches_trailing_edge'] == b['reaches_trailing_edge'], 'TE separation differs')
                for endpoint in ('start', 'end'): close(a[endpoint]['s_over_cref'], b[endpoint]['s_over_cref'])
        rows.append(dict(element=e, recomputed_CL=float(model.post.cl), recomputed_CD=float(model.post.cd),
            recomputed_CM=float(model.post.cm), maximum_Cp_difference=float(np.max(abs(model.post.cp[:model.foil.N][::-1]-surface['Cp'])))))
    return dict(job=job['name'], index=rec['index'], role=rec['role'], alpha_deg=alpha, BL_residual_l2=blnorm,
        coupling_residual_max=coupling, max_residual=max(blnorm, coupling), transition_cells=len(roots),
        transition_root_maximum_shift=shift, elements=rows, new_flow_solves=0,
        surface_s={m.geom.name: m.foil.s.tolist() for m in solver.models})


def audit(study, launch, out):
    inputs = Inputs(); spec = inputs.json(study/'study.json'); digest = run.sha(study/'study.json')
    launched = inputs.json(launch/'launch_receipt.json'); prepared = inputs.json(launch/'preparation.json')
    require(launched['study_sha256'] == prepared['study_sha256'] == digest, 'Launch/specification mismatch')
    for name, expected in spec['sha256'].items():
        require(run.sha(inputs.add(ROOT/name)) == expected and
                run.sha(inputs.add(study/'source_snapshot'/name)) == expected, 'Frozen input changed: '+name)
    hashes = inputs.json(launch/'retrieved_files_sha256.json'); retrieval = inputs.json(launch/'retrieval.json')
    require(retrieval['all_files_match'] and retrieval['files_verified'] == len(hashes), 'Incomplete retrieval')
    for name, expected in hashes.items(): require(run.sha(inputs.add(study/name)) == expected, 'Retrieved file changed: '+name)
    for name, expected in retrieval['launch_file_sha256'].items():
        require(run.sha(inputs.add(launch/('remote_'+name))) == expected, 'Remote launch/preflight artifact changed: '+name)
    preflight = inputs.json(launch/'remote_linear_preflight.json')
    require(preflight['passed'] and preflight['baseline_operator_exact'] and
            preflight['baseline_source_coupling_exact'] and preflight['nonlinear_solves'] == 0,
            'Remote operator/geometry preflight failed')
    require(run.sha(inputs.add(launch/'remote_linear_preflight.npz')) == preflight['geometry_sha256'],
            'Remote preflight geometry archive changed')
    for name, expected in preflight['source_sha256'].items(): require(spec['sha256'][name] == expected, 'Remote preflight source differs')
    with np.load(launch/'remote_linear_preflight.npz', allow_pickle=False) as saved:
        require(set(saved.files) == {f'factor{f}_{e}_{k}' for f in (1, 2, 4) for e in ELEMENTS for k in ('x', 's', 't')},
                'Remote preflight geometry coverage differs')
        remote_geometry = {k: saved[k].copy() for k in saved.files}
    require([r['factor'] for r in preflight['factors']] == [1, 2, 4] and
            all(r['residual'] <= 1e-10 and r['source_coupling']['passed']
                for r in preflight['factors']), 'Remote linear residual/coverage differs')
    runtime = inputs.json(study/'runtime.json'); status = inputs.json(study/'status.json')
    require(status['state'] in ('finished', 'execution_failed'), 'Study still running')
    require(runtime['python'] == spec['expected_python'] and runtime['versions'] == spec['expected_versions'] and
            set(runtime['thread_limits'].values()) == {'1'}, 'Runtime differs')
    pending, coverage, results = [], [], {}
    for job in spec['jobs']:
        folder = study/job['name']
        if not (folder/'result.json').exists():
            require(not folder.exists() or not any(folder.rglob('*.npz')),
                    'Interrupted attempt has unaudited saved states; preserve it for a separate recovery audit')
            coverage.append(dict(job=job['name'], execution='missing', selected_converged=False, evaluated_requests=0)); continue
        result, states = audit_exports(folder, job, spec, runtime, inputs); results[job['name']] = result
        if states:
            require(result['resolution_controls'] == next(r['resolution_controls'] for r in preflight['factors']
                    if r['factor'] == job['subdivision_factor']), 'Run/preflight resolution controls differ')
        coverage.append(dict(job=job['name'], execution=result['execution_state'], selected_converged=result['feasible'],
            evaluated_requests=len(result['records']), error=result['error']))
        pending.extend((job, state, result['resolution_controls']) for state in states.values())
    controls = {job['name']: run.compare_control(study/job['name'], job) for job in spec['jobs']
                if job['control'] and job['name'] in results}
    control_gate = len(controls) == 2 and all(c['passed'] for c in controls.values())
    if (study/'controls.json').exists():
        require(inputs.json(study/'controls.json') == dict(passed=control_gate, jobs=controls), 'Control gate changed')
    environment = dict(python=platform.python_version(), versions={k: importlib.metadata.version(k) for k in run.VERSIONS})
    receipt = dict(schema=1, complete=False, evidence_verified=False, input_sha256=inputs.files, study_sha256=digest,
        audit_runtime=environment, execution_state=status['state'], coverage=coverage, controls=controls,
        control_gate_passed=control_gate, saved_states=len(pending), states=[], new_flow_solves=0,
        criteria=dict(provenance=True, exports=True, fresh_saved_states=None, exact_controls=bool(control_gate)))
    out.mkdir(parents=True, exist_ok=True); path = out/'residual_recheck.json'
    if path.exists():
        old = run.read(path)
        require(old['input_sha256'] == receipt['input_sha256'] and old['audit_runtime'] == environment, 'Audit inputs changed; preserve prior receipt')
        receipt = old
    for job, state, controls in pending:
        if any(r['job'] == job['name'] and r['index'] == state['record']['index'] for r in receipt['states']): continue
        row = recompute(job, spec, state, controls, remote_geometry); receipt['states'].append(row); run.write(path, receipt)
        print(f"{len(receipt['states'])}/{len(pending)} {job['name']}[{row['index']}]: residual {row['max_residual']:.3g}", flush=True)
    receipt.update(complete=True, evidence_verified=True)
    receipt['criteria']['fresh_saved_states'] = len(receipt['states']) == len(pending)
    run.write(path, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=ROOT/run.STUDY)
    parser.add_argument('--launch', type=Path, default=ROOT/run.LAUNCH)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    receipt = audit(args.study.resolve(), args.launch.resolve(), args.out.resolve())
    print(json.dumps(dict(evidence_verified=receipt['evidence_verified'], criteria=receipt['criteria'],
        artifacts=[str(args.out/'residual_recheck.json')]), allow_nan=False))
    return 0


if __name__ == '__main__': raise SystemExit(main())
