#!/usr/bin/env python3
"""Independent audit of every saved V2 tranche-1 state; no nonlinear flow solve.

Each converged evaluation (controls, endpoint seeds and path points) gets a full
export/array consistency check and a fresh residual, force, Cp/Cf, station,
transition-root and wake-drag evaluation in a newly built solver. Failed
evaluations are checked to carry no force and no saved state.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, Inputs, sha, write
from report_s702_drag_mechanism import require, close, audit_surface, side_data
from report_s702_wake_extent import thermo, array_close
from viscous import CoupledSolver, boundary_layer_summary, mf, te_info
import run_s414_transition_boundary as deck

STUDY = ROOT/deck.STUDY
LAUNCH = ROOT/'validation/results/s414_transition_boundary_launch_20261004'
OUT = ROOT/'validation/results/s414_transition_boundary_analysis_20261004'
STUDY_HASH = 'b4d2be76dee5bff3d13e8dcf5b5fb7b2f8bfc12e424b9bb47830c5fbc986c965'
AUDIT_CODE = ['validation/audit_s414_transition_boundary.py', 'validation/report_s702_wake_extent.py',
              'validation/report_s702_drag_mechanism.py', 'validation/report_design_initialization.py']


def load(inputs):
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Unexpected study')
    launch = inputs.json(LAUNCH/'launch_receipt.json'); preparation = inputs.json(LAUNCH/'preparation.json')
    require(launch['study_sha256'] == preparation['study_sha256'] == STUDY_HASH, 'Launched study differs')
    spec = inputs.json(STUDY/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    receipt = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in receipt.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved output changed: '+name)
    retrieval = inputs.json(LAUNCH/'retrieval.json')
    require(retrieval['all_files_match'] and retrieval['files_verified'] == len(receipt), 'Retrieval coverage changed')
    runtime = inputs.json(STUDY/'runtime.json'); status = inputs.json(STUDY/'status.json')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == len(spec['jobs']),
            'Incomplete execution')
    require(runtime['python'] == spec['expected_python'] and runtime['versions'] == spec['expected_versions'],
            'Solver environment changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread limits changed')
    automatic = inputs.json(STUDY/'summary.json')
    # Pure comparison with the archive: never rewrite controls.json in the retrieved study.
    controls = deck.controls(STUDY, spec)
    require(controls['passed'] and controls == automatic['controls'] == inputs.json(STUDY/'controls.json'),
            'Exact control gate differs')
    for name in AUDIT_CODE: inputs.add(ROOT/name)
    return spec, runtime, automatic, controls, len(receipt)


def audit_exports(folder, job, spec, runtime, inputs):
    """Every evaluation of one job: records, polar, diagnostics and saved arrays must agree."""
    result = inputs.json(folder/'result.json'); meta = inputs.json(folder/'flow/metadata.json')
    polar = inputs.csv(folder/'flow/polar.csv'); param = thermo(spec); states = {}
    require(result['complete'] and result['job'] == job['name'] and result['amplification_jacobian'] == job['amplification_jacobian'],
            'Result does not belong to the declared job')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ from manifest')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(result['records']) == len(polar) == len(meta['diagnostics']) <= len(job['requests']), 'Evaluation counts differ')
    require(len(result['records']) == len(job['requests']) or result['error'] == 'Endpoint seed failed; path not started',
            'Unevaluated requests without a declared seed failure')
    previous = None
    for index, (rec, p, diag, request) in enumerate(zip(result['records'], polar, meta['diagnostics'], job['requests'])):
        require(rec['index'] == index and rec['role'] == request['role'] and rec['requested_alpha_deg'] == request['alpha_deg'],
                'Evaluation order differs from the declared requests')
        require(rec['converged'] == bool(int(p['converged'])) == (diag['status'] == 'converged'), 'Convergence differs')
        for field in ('CL', 'CD', 'CM', 'residual'):
            if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null evaluation has finite export')
            else: close(rec[field], float(p[field]), context='record/export '+field)
        close(float(p['alpha_deg']), request['alpha_deg'], context='exported incidence')
        for key in ('mach', 'reynolds', 'ncrit'): close(diag[key], spec[key], context=key)
        require(diag['ncrit_by_element'] == dict(fore=spec['ncrit'], aft=spec['ncrit']), 'Per-element Ncrit differs')
        require(diag['amplification_jacobian'] == job['amplification_jacobian'], 'Recorded derivative mode differs')
        require(rec['iterations_final_step'] == diag['iterations'] == len(diag['history']) <= spec['iterations'], 'Iteration count differs')
        require([(t['alpha_deg'], t['converged'], t['CL']) for t in rec['continuation']]
                == [(t['alpha_deg'], t['converged'], t['CL']) for t in diag['continuation']], 'Continuation trail differs')
        # The recorded start must agree with the solver's own trail; a warm start is never called cold.
        require(rec['solver_initialized_before'] == (previous is not None) and rec['previous_converged_alpha_deg'] == previous,
                'Recorded starting state differs from the evaluation sequence')
        require(rec['reached_by'] == deck.reached_by(previous is not None, previous, request['alpha_deg'], rec['continuation']),
                'Recorded start description differs')
        if not rec['converged']:
            require(all(rec[k] is None for k in ('CL', 'CD', 'CM')) and not rec['state_files'] and diag['failure_reason'],
                    'Failed evaluation has a force or a saved state')
            require(not any((folder/f'flow/{kind}_{index:03d}_{e}.csv').exists() for kind in ('surface', 'wake') for e in deck.ELEMENTS),
                    'Failed evaluation has surface exports')
            continue
        previous = request['alpha_deg']
        require(rec['admissible'] and not diag['failure_reason'], 'Converged record is inadmissible')
        residual = float(p['residual'])
        require(residual <= spec['residual_tolerance'], 'Tolerance failure')
        close(residual, max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
        for key in ('element_forces', 'boundary_layers'): require(rec[key] == diag[key], 'Recorded diagnostics differ')
        for field, term in (('CL', 'CL'), ('CD', 'CD_wake'), ('CM', 'CM')):
            close(rec[field], sum(f[term] for f in rec['element_forces'].values()), context='element sum '+field)
        for row in rec['transition_audit']:
            require(row['inside_cell'] and abs(row['amplification_residual']) <= 1e-7 and row['station_residual_max'] <= 1e-7,
                    'Subcell transition audit failed')
        arrays, surfaces = {}, {}
        for element, panels in zip(deck.ELEMENTS, spec['panels']):
            surface = inputs.array(folder/f'flow/surface_{index:03d}_{element}.csv')
            wake = inputs.array(folder/f'flow/wake_{index:03d}_{element}.csv')
            name = f'states/point_{index:02d}_{element}.npz'; require(name in rec['state_files'], 'State file not recorded')
            with np.load(inputs.add(folder/name), allow_pickle=False) as saved:
                a = {k: saved[k].copy() for k in saved.files}
            n = len(surface)
            require(n == panels+1 and a['U'].shape == (4, n+len(wake)), 'State dimensions differ')
            require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite Newton state')
            for k, field in ((0, 'theta'), (1, 'delta_star'), (2, 'amplification_or_shear')):
                array_close(a['U'][k, :n][::-1], surface[field], 'Newton/surface '+field)
            for k, field in ((0, 'theta'), (1, 'delta_star')):
                array_close(a['U'][k, n:], wake[field], 'Newton/wake '+field)
            uk = mf.get_uk(a['U'][3], param)[0]
            array_close(uk[:n][::-1], surface['ue_over_vinf'], 'Newton/surface corrected velocity')
            array_close(uk[n:], wake['ue_over_vinf'], 'Newton/wake corrected velocity')
            hk = np.array([mf.get_Hk(a['U'][:, j], param)[0] for j in range(n+len(wake))])
            array_close(hk[:n][::-1], surface['Hk'], 'Newton/surface Hk'); array_close(hk[n:], wake['Hk'], 'Newton/wake Hk')
            require(np.array_equal(a['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence mismatch')
            labels = np.full(n, 'upper', dtype='<U5'); labels[a['lower']] = 'lower'
            require(np.array_equal(labels[::-1], surface['side']), 'Surface indices differ')
            require(set(a['lower']) | set(a['upper']) == set(range(n)) and not set(a['lower']) & set(a['upper']), 'Surface partition differs')
            for k, field in enumerate(('x_over_cref', 'y_over_cref')):
                array_close(a['xy'][k][::-1], surface[field], 'Surface geometry differs')
                array_close(a['wake_xy'][k], wake[field], 'Wake geometry differs')
            require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
            for side in ('upper', 'lower'): audit_surface(side_data(surface, side), rec['boundary_layers'][element][side])
            ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0); arclength = np.r_[0., np.cumsum(ds)]
            array_close(a['wake_s']-a['wake_s'][0], arclength, 'Wake arclength differs')
            array_close(np.linalg.norm(a['wake_t'], axis=0), np.ones(len(wake)), 'Wake tangent normalization')
            info = rec['wake_discretization'][element]; require(info['nodes'] == len(wake), 'Wake node count differs')
            for k, value in (('arclength', arclength[-1]), ('first_spacing', ds[0]), ('last_spacing', ds[-1])):
                close(info[k], value, context=k)
            # Squire–Young wake drag rebuilt from the saved terminal wake station.
            h = a['U'][1, n:]/a['U'][0, n:]; sy = 2*a['U'][0, n:]*(uk[n:]/param.Vinf)**((5+h)/2)
            close(sy[-1], rec['element_forces'][element]['CD_wake'], context='terminal Squire–Young drag')
            arrays[element], surfaces[element] = a, surface
        states[index] = dict(record=rec, residual=residual, arrays=arrays, surfaces=surfaces)
    return result, states


def recompute(job, spec, state):
    """Fresh solver through the standard constructor; restore the saved state and re-evaluate everything."""
    rec = state['record']; alpha = rec['alpha_deg']
    solver = CoupledSolver(deck.geometry_elements(job['scenario'], spec['panels']), reynolds=spec['reynolds'],
        mach=spec['mach'], ncrit=spec['ncrit'], iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
        wake_length=job['wake_length'], amplification_jacobian=job['amplification_jacobian'])
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = state['arrays'][model.geom.name]
        array_close(model.foil.x, a['xy'], 'Rebuilt foil geometry', 2e-12)
        array_close(model.wake.x, a['wake_xy'], 'Rebuilt inviscid wake geometry', 2e-11)
        array_close(model.wake.t, a['wake_t'], 'Rebuilt wake tangents', 2e-11)
        array_close(model.wake.s, a['wake_s'], 'Rebuilt wake arclength', 2e-11)
        model.glob.U = a['U'].copy(); model.glob.Nsys = a['U'].shape[1]
        model.vsol.turb = a['turbulent'].copy()
        model.isol.Istag = [int(a['lower'][0]), int(a['upper'][0])]
        model.isol.sgnue = np.ones(model.foil.N); model.isol.sgnue[a['lower']] = -1
        model.isol.sstag = float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model); mf.stagpoint_move(model)
        array_close(model.isol.xi, a['xi'], 'Restored station coordinates', 2e-12)
        for k, label in enumerate(('lower', 'upper')):
            require(np.array_equal(model.vsol.Is[k], a[label]), 'Restored branch indices differ')
        if abs(te_info(model.foil.x)[1]) < 1e-12: model.vsol.wgap = np.zeros(model.wake.N)
        else: mf.set_wake_gap(model)
    solver._build_source_coupling()
    residual, _ = solver._system()
    blnorm = float(np.linalg.norm(residual[:3*solver.nsys])); coupling = float(np.max(abs(residual[3*solver.nsys:])))
    require(max(blnorm, coupling) <= spec['residual_tolerance'], 'Restored Newton residual fails tolerance')
    roots = deck.transition_audit(solver); old = rec['transition_audit']
    require(len(roots) == len(old), 'Transition cell count differs')
    root_difference = 0.
    for new, saved in zip(roots, old):
        require((new['element'], new['side']) == (saved['element'], saved['side']) and new['inside_cell'], 'Transition cell differs')
        require(abs(new['amplification_residual']) <= 1e-7 and new['station_residual_max'] <= 1e-7, 'Recomputed transition root fails')
        root_difference = max(root_difference, abs(new['transition_xi']-saved['transition_xi']))
    require(root_difference <= 1e-9, 'Recomputed transition root moved')
    rows = []
    for model in solver.models:
        e = model.geom.name; a = state['arrays'][e]; surface = state['surfaces'][e]
        array_close(model.glob.U, a['U'], 'Audit must not change saved solution', 0.)
        mf.calc_force(model); mf.get_distributions(model)
        if model.oper.Ma > 0:
            for side in range(3):
                param = mf.build_param(model, side)
                for j in model.vsol.Is[side]:
                    rho, _ = mf.get_rho(model.glob.U[:, j], param); model.post.cf[j] *= rho/model.oper.rho
        array_close(model.post.cp[:model.foil.N][::-1], surface['Cp'], 'Recomputed Cp', 2e-12)
        array_close(model.post.cf[:model.foil.N][::-1], surface['cf'], 'Recomputed Cf', 2e-12)
        forces = rec['element_forces'][e]
        for name, field in (('CL', 'cl'), ('CM', 'cm'), ('CD_wake', 'cd'), ('CD_friction', 'cdf')):
            close(getattr(model.post, field), forces[name], context='Recomputed '+e+' '+name)
        layers = boundary_layer_summary(model)
        for side in ('upper', 'lower'):
            expected = rec['boundary_layers'][e][side]; tr = layers[side]['transition']; before = expected['transition']
            require((tr is None) == (before is None), 'Restored transition presence differs')
            if tr: close(tr['s_over_cref'], before['s_over_cref'], context='Recomputed transition')
            for key in ('min_cf', 'max_Hk', 'separated_arclength_fraction'):
                close(layers[side][key], expected[key], context='Recomputed BL '+key)
            require(len(layers[side]['separation_regions']) == len(expected['separation_regions']), 'Separation interval count differs')
            for new, saved in zip(layers[side]['separation_regions'], expected['separation_regions']):
                require(new['reaches_trailing_edge'] == saved['reaches_trailing_edge'], 'TE separation flag differs')
                for endpoint in ('start', 'end'):
                    close(new[endpoint]['s_over_cref'], saved[endpoint]['s_over_cref'], context='Recomputed separation')
        rows.append(dict(element=e, recomputed_CL=float(model.post.cl), recomputed_CD=float(model.post.cd),
            maximum_Cp_difference=float(np.max(abs(model.post.cp[:model.foil.N][::-1]-surface['Cp'])))))
    return dict(job=job['name'], index=rec['index'], role=rec['role'], alpha_deg=alpha, mode=job['amplification_jacobian'],
        BL_residual_l2=blnorm, coupling_residual_max=coupling, recorded_residual=state['residual'],
        max_residual=max(blnorm, coupling), transition_cells=len(roots), transition_root_maximum_shift=root_difference,
        new_flow_solves=0, elements=rows)


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs()
    spec, runtime, automatic, controls, retrieved = load(inputs)
    audited = {}; failed = []
    for job in spec['jobs']:
        result, states = audit_exports(STUDY/job['name'], job, spec, runtime, inputs)
        audited[job['name']] = states
        failed += [dict(job=job['name'], index=r['index'], alpha_deg=r['requested_alpha_deg'], failure_reason=r['failure_reason'],
                        previous_converged_alpha_deg=r['previous_converged_alpha_deg']) for r in result['records'] if not r['converged']]
    environment = dict(python=platform.python_version(),
        versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})
    total = sum(len(s) for s in audited.values())
    receipt = dict(input_sha256=inputs.files, audit_runtime=environment, states=[], complete=False, new_flow_solves=0,
        retrieved_files_verified=retrieved, controls=controls, failed_evaluations=failed, saved_states=total,
        method='Full export/array audit of every evaluation and fresh coupled residual, force, Cp/Cf, station, '
               'transition-root and wake-drag evaluation of every saved state; no nonlinear solve')
    path = OUT/'residual_recheck.json'
    if path.exists():
        previous = json.loads(path.read_text())
        require(previous['input_sha256'] == receipt['input_sha256'] and previous['audit_runtime'] == environment,
                'Audit changed; preserve previous receipt')
        receipt = previous
    for job in spec['jobs']:
        for index, state in sorted(audited[job['name']].items()):
            if any(r['job'] == job['name'] and r['index'] == index for r in receipt['states']): continue
            row = recompute(job, spec, state); receipt['states'].append(row); write(path, receipt)
            print(f"{len(receipt['states'])}/{total} {job['name']}[{index}]: residual {row['max_residual']:.3g}", flush=True)
    receipt['complete'] = True; write(path, receipt)
    print('audit complete:', total, 'saved states;', len(failed), 'failed evaluation(s) without force')


if __name__ == '__main__':
    main()
