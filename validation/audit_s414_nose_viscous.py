#!/usr/bin/env python3
"""Independent audit of the V4 convex-nose viscous run; no nonlinear flow solve.

Fixed-lift targets get the unchanged full-export audit; fixed-incidence paths
get the unchanged V2 evaluation audit. Every saved state is then restored in a
newly built solver on the frozen geometry and its residual, forces, Cp/Cf,
stations, transition roots and wake drag are rebuilt.
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
from report_s702_drag_mechanism import require, close
from report_s702_wake_extent import audit as audit_fixed_lift, array_close
from audit_s414_transition_boundary import audit_exports as audit_path
from viscous import CoupledSolver, boundary_layer_summary, mf, te_info
import run_s414_transition_boundary as deck
import run_s414_nose_viscous as run

STUDY = ROOT/run.STUDY
LAUNCH = ROOT/'validation/results/s414_nose_viscous_launch_20261004'
OUT = ROOT/'validation/results/s414_nose_viscous_analysis_20261004'
STUDY_HASH = '5ffd7558f061ed7872a051d59e02f05d0a54827d2af30906acbdbf29274a7709'
AUDIT_CODE = ['validation/audit_s414_nose_viscous.py', 'validation/audit_s414_transition_boundary.py',
              'validation/report_s702_wake_extent.py', 'validation/report_s702_drag_mechanism.py',
              'validation/report_design_initialization.py']


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
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == len(spec['jobs']) == 6,
            'Incomplete execution')
    require(runtime['python'] == spec['expected_python'] and runtime['versions'] == spec['expected_versions'],
            'Solver environment changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread limits changed')
    automatic = inputs.json(STUDY/'summary.json')
    # Pure recomputation of the exact control gate; controls.json in the retrieved study is never rewritten.
    controls = run.controls(STUDY, spec)
    require(controls == automatic['controls'] == inputs.json(STUDY/'controls.json') and controls['passed'], 'Exact control gate differs')
    for name in AUDIT_CODE: inputs.add(ROOT/name)
    return spec, runtime, controls, len(receipt)


def recompute(job, spec, record, arrays, surfaces, label):
    """Fresh solver on the declared geometry; restore the saved state and re-evaluate everything."""
    alpha = record['alpha_deg']
    solver = CoupledSolver(run.geometry_elements(job['scenario'], spec['panels']), reynolds=spec['reynolds'], mach=spec['mach'],
        ncrit=spec['ncrit'], iterations=spec['iterations'], tolerance=spec['residual_tolerance'], wake_length=job['wake_length'],
        amplification_jacobian=job['amplification_jacobian'])
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = arrays[model.geom.name]
        array_close(model.foil.x, a['xy'], 'Rebuilt foil geometry', 2e-12)
        array_close(model.wake.x, a['wake_xy'], 'Rebuilt inviscid wake geometry', 2e-11)
        array_close(model.wake.t, a['wake_t'], 'Rebuilt wake tangents', 2e-11)
        array_close(model.wake.s, a['wake_s'], 'Rebuilt wake arclength', 2e-11)
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
    solver._build_source_coupling()
    residual, _ = solver._system()
    blnorm = float(np.linalg.norm(residual[:3*solver.nsys])); coupling = float(np.max(abs(residual[3*solver.nsys:])))
    require(max(blnorm, coupling) <= spec['residual_tolerance'], 'Restored Newton residual fails tolerance')
    roots = deck.transition_audit(solver); old = record['transition_audit']; shift = 0.
    require(len(roots) == len(old), 'Transition cell count differs')
    for new, saved in zip(roots, old):
        require((new['element'], new['side']) == (saved['element'], saved['side']) and new['inside_cell']
                and abs(new['amplification_residual']) <= 1e-7 and new['station_residual_max'] <= 1e-7, 'Recomputed transition root fails')
        shift = max(shift, abs(new['transition_xi']-saved['transition_xi']))
    require(shift <= 1e-9, 'Recomputed transition root moved')
    rows = []
    for model in solver.models:
        e = model.geom.name; surface = surfaces[e]
        array_close(model.glob.U, arrays[e]['U'], 'Audit must not change saved solution', 0.)
        mf.calc_force(model); mf.get_distributions(model)
        if model.oper.Ma > 0:
            for side in range(3):
                param = mf.build_param(model, side)
                for j in model.vsol.Is[side]:
                    rho, _ = mf.get_rho(model.glob.U[:, j], param); model.post.cf[j] *= rho/model.oper.rho
        array_close(model.post.cp[:model.foil.N][::-1], surface['Cp'], 'Recomputed Cp', 2e-12)
        array_close(model.post.cf[:model.foil.N][::-1], surface['cf'], 'Recomputed Cf', 2e-12)
        for name, field in (('CL', 'cl'), ('CM', 'cm'), ('CD_wake', 'cd'), ('CD_friction', 'cdf')):
            close(getattr(model.post, field), record['element_forces'][e][name], context='Recomputed '+e+' '+name)
        layers = boundary_layer_summary(model)
        for side in ('upper', 'lower'):
            expected = record['boundary_layers'][e][side]; tr = layers[side]['transition']; before = expected['transition']
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
    return dict(state=label, job=job['name'], scenario=job['scenario'], alpha_deg=alpha, BL_residual_l2=blnorm,
        coupling_residual_max=coupling, max_residual=max(blnorm, coupling), transition_cells=len(roots),
        transition_root_maximum_shift=shift, new_flow_solves=0, elements=rows)


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs()
    spec, runtime, controls, retrieved = load(inputs)
    pending, failed_targets, failed_points = [], [], []
    for job in spec['jobs']:
        folder = STUDY/job['name']
        if job['kind'] == 'fixed_lift':
            state = audit_fixed_lift(folder, job, spec, runtime, inputs); r = state['result']
            require(all(d['amplification_jacobian'] == job['amplification_jacobian'] for d in state['metadata']['diagnostics']), 'Recorded mode differs')
            if r['feasible']:
                pending.append((job['name'], job, r, state['arrays'], state['surfaces']))
            else:
                failed_targets.append(dict(job=job['name'], target_cl=job['target_cl'], error=r['error'], evaluations=len(r['records']),
                    attempts=[dict(path=a['path'], feasible=a['feasible'], error=a.get('error', '')) for a in r['attempts']]))
        else:
            result, states = audit_path(folder, job, spec, runtime, inputs)
            for index, state in sorted(states.items()):
                pending.append((f"{job['name']}[{index}]", job, state['record'], state['arrays'], state['surfaces']))
            failed_points += [dict(job=job['name'], index=r['index'], alpha_deg=r['requested_alpha_deg'], failure_reason=r['failure_reason'],
                                   previous_converged_alpha_deg=r['previous_converged_alpha_deg']) for r in result['records'] if not r['converged']]
    environment = dict(python=platform.python_version(),
        versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})
    receipt = dict(input_sha256=inputs.files, audit_runtime=environment, states=[], complete=False, new_flow_solves=0,
        retrieved_files_verified=retrieved, controls=controls, failed_targets=failed_targets, failed_path_points=failed_points,
        saved_states=len(pending),
        method='Unchanged fixed-lift and path export audits, then a fresh coupled residual, force, Cp/Cf, station, '
               'transition-root and wake-drag evaluation of every saved state on the frozen geometry; no nonlinear solve')
    path = OUT/'residual_recheck.json'
    if path.exists():
        previous = json.loads(path.read_text())
        require(previous['input_sha256'] == receipt['input_sha256'] and previous['audit_runtime'] == environment,
                'Audit changed; preserve previous receipt')
        receipt = previous
    for label, job, record, arrays, surfaces in pending:
        if any(r['state'] == label for r in receipt['states']): continue
        row = recompute(job, spec, record, arrays, surfaces, label); receipt['states'].append(row); write(path, receipt)
        print(f"{len(receipt['states'])}/{len(pending)} {label}: residual {row['max_residual']:.3g}", flush=True)
    receipt['complete'] = True; write(path, receipt)
    print('audit complete:', len(pending), 'saved states;', len(failed_targets), 'failed target(s);', len(failed_points), 'failed path point(s)')


if __name__ == '__main__':
    main()
