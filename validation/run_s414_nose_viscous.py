#!/usr/bin/env python3
"""V4 viscous check of the admitted convex-nose reconstruction; six restartable jobs.

One nominal control that must repeat an archived corrected state exactly, then
three fixed-lift targets and two fixed-incidence paths on the convex-nose
geometry. Protocol: validation/S414_NOSE_VISCOUS_PROTOCOL.md.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
from dataclasses import asdict
import fcntl
import importlib.metadata
import os
from pathlib import Path
import platform
import shutil
import sys
import time
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_s414_transition_boundary as deck
from run_s414_transition_boundary import read, sha, stamp, write, finite, prepare, VERSIONS, PYTHON, THREADS
from s414_geometry import elements as declared_elements
from multifoil import save_results, validate_geometry
from optimization import admissible
from fixed_lift_policy import InitializationPath, solve_with_initialization
from viscous import CoupledSolver

STUDY = 'validation/results/s414_nose_viscous_20261004'
GEOMETRY = 'validation/results/s414_nose_geometry_20261004'
FORE = GEOMETRY+'/fore_input_convex_nose.dat'
ARCHIVE = 'validation/results/s414_jacobian_20261003'
ARCHIVE_RECEIPT = 'validation/results/s414_jacobian_launch_20261003/retrieved_files_sha256.json'
LOWER = 'validation/results/s414_lower_bucket_20261003'
T2 = 'validation/results/s414_transition_boundary_t2_20261004'
PROTOCOL = 'validation/S414_NOSE_VISCOUS_PROTOCOL.md'
MODE = 'corrected_af'
STATE_FIELDS = ('U', 'xy', 'xi', 'turbulent', 'wake_xy', 'lower', 'upper', 'wake_s', 'wake_t')


def geometry_elements(scenario, panels):
    """Nominal, or the frozen convex-nose fore input with the unchanged nominal aft input."""
    nominal = declared_elements('nominal', panels)
    if scenario == 'nominal':
        return nominal
    if scenario != 'convex_nose' or list(panels) != [340, 240]:
        raise ValueError('Unknown scenario or panel count for the frozen fore input')
    fore = np.loadtxt(ROOT/FORE, skiprows=1)
    if fore.shape != nominal[0][1].shape or not np.array_equal(fore[[0, -1]], nominal[0][1][[0, -1]]):
        raise ValueError('Frozen fore input does not match the nominal trailing edge')
    result = [('fore', fore), nominal[1]]
    validate_geometry(result)
    return result


def paths_for(job, spec):
    """The archived ordered initialization policy (run_s702_wake_extent.paths_for)."""
    return [InitializationPath('cold_local', job['center']), InitializationPath('warm_anchor', spec['warm_start_alpha']),
            InitializationPath('warm_entry', spec['warm_start_alpha'], (job['center']+spec['entry_offset'],))]


def record(sol, label):
    """Same fields as the archived evaluation records (s702_benchmark.record)."""
    return dict(label=label, alpha_deg=sol.alpha, CL=finite(sol.cl), CD=finite(sol.cd), CM=finite(sol.cm),
        converged=bool(sol.converged), admissible=bool(admissible(sol)), residual=finite(sol.residual),
        failure_reason=sol.diagnostics.get('failure_reason', ''),
        local_supersonic=any(sol.diagnostics.get('local_supersonic_detected', {}).values()))


def specification():
    geometry = read(ROOT/GEOMETRY/'verification.json'); summary = read(ROOT/GEOMETRY/'summary.json')
    if not (geometry['complete'] and geometry['admission']['convex_nose'] and geometry['convex_nose_material']):
        raise ValueError('The convex-nose alternative is not admitted with a material pressure change; no viscous run')
    if sha(ROOT/FORE) != summary['frozen_inputs'][Path(FORE).name] or sha(ROOT/FORE) != geometry['output_sha256'][Path(FORE).name]:
        raise ValueError('Frozen convex-nose input changed')
    old = read(ROOT/ARCHIVE/'study.json'); receipt = read(ROOT/ARCHIVE_RECEIPT); lower = read(ROOT/LOWER/'study.json'); t2 = read(ROOT/T2/'study.json')
    reference = 'cl0056_corrected_af'; files = {PROTOCOL, FORE, GEOMETRY+'/verification.json', GEOMETRY+'/summary.json',
        ARCHIVE+'/study.json', ARCHIVE_RECEIPT, LOWER+'/study.json', T2+'/study.json', 'vendor/mfoil.py', 'vendor/__init__.py',
        deck.QUAL+'/report_test_fore.dat', deck.QUAL+'/report_test_aft.dat', deck.QUAL+'/verification.json'} | set(deck.TEST_DATA)
    for name in ('result.json', 'fore_state.npz', 'aft_state.npz'):
        if sha(ROOT/ARCHIVE/reference/name) != receipt[reference+'/'+name]:
            raise ValueError('Archived control changed: '+name)
        files.add(f'{ARCHIVE}/{reference}/{name}')
    archived = {j['name']: j for j in old['jobs']}; center = {j['target_cl']: j['center'] for j in lower['jobs']}
    if archived[reference]['amplification_jacobian'] != MODE or archived[reference]['center'] != center[.056]:
        raise ValueError('Archived control is not the corrected CL 0.056 job')
    common = dict(amplification_jacobian=MODE, wake_length=1.)
    lift = lambda name, scenario, target, **extra: dict(common, name=name, scenario=scenario, kind='fixed_lift',
        target_cl=target, center=center[target], **extra)
    angles = sorted(t2['nominal_new_incidences_deg']+t2['tranche1_angles_deg'])
    # The increasing path starts cold at its first incidence. A cold start at the upper end does not
    # converge for this geometry (preparation finding), so the decreasing path starts cold at the
    # archived warm-anchor incidence and reaches its first point by recorded continuation.
    path = lambda direction, order, start: dict(common, name='convex_nose_'+direction, scenario='convex_nose', kind='path',
        direction=direction, start_alpha=start, requests=[dict(alpha_deg=a, role='path') for a in order])
    jobs = [lift('nominal_cl0056_control', 'nominal', .056, control=True, reference_result=f'{ARCHIVE}/{reference}/result.json'),
            lift('convex_nose_cl0056', 'convex_nose', .056, control=False),
            lift('convex_nose_cl0175', 'convex_nose', .175, control=False),
            lift('convex_nose_cl0284', 'convex_nose', .284, control=False),
            dict(path('increasing', angles, angles[0]), control=False),
            dict(path('decreasing', angles[::-1], old['warm_start_alpha']), control=False)]
    for folder in (ROOT, ROOT/'validation'):
        files.update(str(p.relative_to(ROOT)) for p in folder.glob('*.py'))
    settings = {k: old[k] for k in ('reynolds', 'mach', 'ncrit', 'panels', 'reference_chord', 'model_chord_mm', 'aft_shift_c',
        'iterations', 'cl_tolerance', 'residual_tolerance', 'initial_half_width', 'alpha_bounds', 'max_evaluations', 'max_brackets',
        'warm_start_alpha', 'entry_offset', 'expected_versions', 'acceptance_control_CD_counts')}
    return dict(schema=1, experiment='S414 convex-nose reconstruction, V4 viscous check', **settings, expected_python=PYTHON,
        continuation_step=1., min_step=.25, minimum_te_gap=1e-5, path_angles_deg=angles, jobs=jobs, workers_limit=2,
        equivalence_screens=t2['equivalence_screens'], requested_path_points=16, fixed_lift_targets=3, control_targets=1,
        gate='The nominal CL 0.056 control repeats the archived corrected state with identical saved arrays before any convex-nose job',
        qualification='Reconstruction sensitivity on a figure-derived polygon; not the tested contour, not matched-geometry validation; audit required',
        sha256={name: sha(ROOT/name) for name in sorted(files)})


def save_states(folder, solver, prefix=''):
    info = {}; stream = np.array([np.cos(np.deg2rad(solver.last_alpha)), np.sin(np.deg2rad(solver.last_alpha))])
    for model in solver.models:
        np.savez_compressed(folder/f'{prefix}{model.geom.name}{"_state" if not prefix else ""}.npz', U=model.glob.U, xy=model.foil.x,
            xi=model.isol.xi, turbulent=model.vsol.turb, wake_xy=model.wake.x, lower=np.asarray(model.vsol.Is[0]),
            upper=np.asarray(model.vsol.Is[1]), wake_s=model.wake.s, wake_t=model.wake.t)
        spacing = np.linalg.norm(np.diff(model.wake.x, axis=1), axis=0)
        info[model.geom.name] = dict(nodes=int(model.wake.N), arclength=float(np.sum(spacing)), first_spacing=float(spacing[0]),
            last_spacing=float(spacing[-1]), streamwise_start=float(stream@model.wake.x[:, 0]), streamwise_end=float(stream@model.wake.x[:, -1]))
    return info


def make_solver(elements, job, spec, start_alpha):
    return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'], iterations=spec['iterations'],
        tolerance=spec['residual_tolerance'], wake_length=job['wake_length'], minimum_te_gap=spec['minimum_te_gap'],
        continuation_step=spec['continuation_step'], min_step=spec['min_step'], start_alpha=start_alpha,
        amplification_jacobian=job['amplification_jacobian'])


def lift_worker(folder, job, spec, elements):
    """The archived fixed-lift procedure, unchanged, through the standard solver with an explicit mode."""
    solutions = []
    result = dict(amplification_jacobian=job['amplification_jacobian'], scenario=job['scenario'], ncrit=spec['ncrit'],
        target_cl=job['target_cl'], wake_length=job['wake_length'], reused=False, feasible=False, CL=None, CD=None, CM=None,
        alpha_deg=None, error='', records=[], attempts=[], flow_directory=str((folder/'flow').relative_to(ROOT)))
    def checkpoint(sol, context):
        sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit']); solutions.append(sol)
        result['records'].append(dict(record(sol, context['phase']), initialization=context))
        save_results(folder/'flow', solutions, elements, dict(study=spec, job=job, target_cl=job['target_cl']), 'viscous')
        write(folder/'progress.json', result)
    outcome = solve_with_initialization(lambda alpha: make_solver(elements, job, spec, alpha), job['target_cl'], job['center'],
        spec['alpha_bounds'], paths_for(job, spec), half_width=spec['initial_half_width'], cl_tolerance=spec['cl_tolerance'],
        residual_tolerance=spec['residual_tolerance'], max_evaluations=spec['max_evaluations'], max_brackets=spec['max_brackets'],
        on_evaluation=checkpoint)
    result.update(attempts=outcome.attempts, error=outcome.error, initialization_paths=[asdict(p) for p in paths_for(job, spec)])
    if outcome.feasible:
        sol, solver = outcome.solution, outcome.solver
        if abs(solver.last_alpha-sol.alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
            raise ValueError('Returned solver does not match selected state')
        result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
            selected_flow_index=outcome.selected_flow_index, selected_path=outcome.selected_path,
            boundary_layers=sol.diagnostics['boundary_layers'], element_forces=sol.diagnostics['element_forces'],
            transition_audit=deck.transition_audit(solver), wake_discretization=save_states(folder, solver))
    return result


def path_worker(folder, job, spec, elements):
    """Fixed-incidence path with the V2 tranche-1 rules."""
    solver = make_solver(elements, job, spec, job['start_alpha']); solutions = []
    result = dict(job=job['name'], scenario=job['scenario'], amplification_jacobian=job['amplification_jacobian'],
        direction=job['direction'], ncrit=spec['ncrit'], records=[], error='', flow_directory=str((folder/'flow').relative_to(ROOT)))
    (folder/'states').mkdir(exist_ok=True)
    for index, request in enumerate(job['requests']):
        alpha = request['alpha_deg']; initialized, previous = bool(solver.initialized), solver.last_alpha
        sol = solver.solve(alpha); sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit'])
        solutions.append(sol); diag = sol.diagnostics
        trail = [dict(alpha_deg=t['alpha_deg'], converged=bool(t['converged']), CL=finite(t['CL'])) for t in diag['continuation']]
        entry = dict(index=index, role=request['role'], requested_alpha_deg=alpha, alpha_deg=sol.alpha, CL=finite(sol.cl), CD=finite(sol.cd),
            CM=finite(sol.cm), converged=bool(sol.converged), admissible=bool(admissible(sol)), residual=finite(sol.residual),
            failure_reason=diag.get('failure_reason', ''), iterations_final_step=diag['iterations'],
            solver_initialized_before=initialized, previous_converged_alpha_deg=previous,
            reached_by=deck.reached_by(initialized, previous, alpha, trail),
            approach=('cold' if not initialized else 'increasing' if alpha > previous else 'decreasing'),
            continuation=trail, wall_seconds=diag['wall_seconds'], state_files=[])
        if sol.converged:
            if abs(solver.last_alpha-alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
                raise ValueError('Solver state does not match the returned solution')
            entry.update(boundary_layers=diag['boundary_layers'], element_forces=diag['element_forces'],
                transition_audit=deck.transition_audit(solver),
                wake_discretization=save_states(folder/'states', solver, prefix=f'point_{index:02d}_'),
                state_files=[f'states/point_{index:02d}_{e}.npz' for e in deck.ELEMENTS])
        elif any(entry[k] is not None for k in ('CL', 'CD', 'CM')):
            raise ValueError('Failed evaluation carries a finite force')
        result['records'].append(entry)
        save_results(folder/'flow', solutions, elements, dict(study=spec, job=job), 'viscous'); write(folder/'progress.json', result)
        print(f"{index}: alpha={alpha!r} converged={sol.converged} CL={entry['CL']} CD={entry['CD']} {entry['reached_by']}", flush=True)
        if index == 0 and not sol.converged:
            result['error'] = 'Endpoint seed failed; path not started'; break
    return result


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter(); elements = geometry_elements(job['scenario'], spec['panels'])
        result = (lift_worker if job['kind'] == 'fixed_lift' else path_worker)(folder, job, spec, elements)
        result.update(kind=job['kind'], complete=True, wall_seconds=time.perf_counter()-begun)
        write(folder/'result.json', result)
        return result


def compare_control(folder, job, tolerance_counts):
    """Exact comparison with the archived corrected state; no tolerance on saved arrays in the pinned runtime."""
    reference = (ROOT/job['reference_result']).parent; actual, expected = read(folder/'result.json'), read(reference/'result.json')
    passed = bool(actual['feasible'] and expected['feasible']); delta = None; arrays = []
    if passed:
        delta = (actual['CD']-expected['CD'])*1e4
        passed = abs(delta) <= tolerance_counts and actual['selected_path'] == expected['selected_path'] and actual['alpha_deg'] == expected['alpha_deg']
        for element in deck.ELEMENTS:
            with np.load(folder/(element+'_state.npz')) as a, np.load(reference/(element+'_state.npz')) as b:
                if set(a.files) != set(b.files) or set(a.files) != set(STATE_FIELDS): raise ValueError('State fields changed')
                for field in a.files:
                    same = bool(np.array_equal(a[field], b[field])); passed = passed and same
                    arrays.append(dict(element=element, field=field, exactly_equal=same))
    return dict(passed=bool(passed), reference=job['reference_result'], CD_difference_counts=delta, array_comparisons=arrays)


def controls(out, spec):
    result = {j['name']: compare_control(out/j['name'], j, spec['acceptance_control_CD_counts']) for j in spec['jobs'] if j['control']}
    return dict(passed=all(r['passed'] for r in result.values()), jobs=result)


def summarize(out, spec):
    rows = []; points = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json'); meta = read(out/job['name']/'flow/metadata.json')
        if meta['case']['study'] != spec or meta['case']['job'] != job:
            raise ValueError('Exported study/job differs')
        if any(d['amplification_jacobian'] != MODE or d['ncrit_by_element'] != {'fore': 9., 'aft': 9.} for d in meta['diagnostics']):
            raise ValueError('Exported mode or Ncrit differs')
        if job['kind'] == 'fixed_lift':
            ok = r['feasible']
            if ok:
                selected = r['records'][r['selected_flow_index']]
                ok = bool(selected['admissible'] and selected['residual'] <= spec['residual_tolerance']
                    and abs(r['CL']-job['target_cl']) <= spec['cl_tolerance']
                    and abs(sum(v['CD_wake'] for v in r['element_forces'].values())-r['CD']) < 1e-12
                    and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7 and t['station_residual_max'] <= 1e-7
                            for t in r['transition_audit']))
            elif any(r[k] is not None for k in ('CL', 'CD', 'CM', 'alpha_deg')):
                raise ValueError('Failed target contains finite selected force')
            rows.append(dict(job=job['name'], kind='fixed_lift', scenario=job['scenario'], target_cl=job['target_cl'], feasible=r['feasible'],
                execution_check_pass=bool(ok) if r['feasible'] else True, CL=r['CL'], CD=r['CD'], alpha_deg=r['alpha_deg'],
                selected_path=r.get('selected_path'), evaluations=len(r['records']), error=r['error']))
        else:
            for rec in r['records']:
                ok = True
                if rec['converged']:
                    ok = bool(rec['admissible'] and rec['residual'] <= spec['residual_tolerance']
                        and abs(sum(v['CD_wake'] for v in rec['element_forces'].values())-rec['CD']) < 1e-12
                        and all((out/job['name']/n).exists() for n in rec['state_files']) and len(rec['state_files']) == 2
                        and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7 and t['station_residual_max'] <= 1e-7
                                for t in rec['transition_audit']))
                elif any(rec[k] is not None for k in ('CL', 'CD', 'CM')):
                    raise ValueError('Failed evaluation contains a finite force')
                points.append(dict(job=job['name'], direction=job['direction'], index=rec['index'], alpha_deg=rec['requested_alpha_deg'],
                    converged=rec['converged'], execution_check_pass=ok, CL=rec['CL'], CD=rec['CD'], reached_by=rec['reached_by'],
                    failure_reason=rec['failure_reason']))
            points += [dict(job=job['name'], direction=job['direction'], index=None, alpha_deg=q['alpha_deg'], converged=False,
                execution_check_pass=True, CL=None, CD=None, reached_by='not evaluated: endpoint seed failed', failure_reason=r['error'])
                for q in job['requests'][len(r['records']):]]
    control = controls(out, spec); targets = [r for r in rows if r['scenario'] != 'nominal']
    result = dict(required_jobs=len(spec['jobs']), completed_jobs=len(spec['jobs']),
        requested_path_points=len(points), converged_path_points=sum(p['converged'] for p in points),
        failed_path_points=sum(not p['converged'] for p in points),
        fixed_lift_targets=len(targets), converged_fixed_lift_targets=sum(r['feasible'] for r in targets),
        execution_checks_pass=bool(control['passed'] and all(r['execution_check_pass'] for r in rows+points)),
        controls=control, targets=rows, points=points,
        qualification='Execution only; failed targets and points are findings, not forces. Full audit and interpretation pending')
    write(out/'summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/STUDY)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 2:
        parser.error('Use a results directory and 1–2 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(f'Prepared {len(spec["jobs"])} convex-nose jobs; study sha256='+sha(out/'study.json')); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS}, thread_limits={k: os.environ[k] for k in THREADS})
        if runtime['versions'] != VERSIONS or runtime['python'] != PYTHON:
            raise ValueError('Use pinned Python 3.12.13 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        total = len(spec['jobs']); done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(), completed_jobs=done, total_jobs=total, **extra))
        try:
            for phase in ('nominal_control', 'convex_nose'):
                if phase == 'convex_nose':
                    control = controls(out, spec); write(out/'controls.json', control)
                    if not control['passed']:
                        status('execution_failed', error='Nominal control differs from the archived corrected state; convex-nose jobs not launched'); return 2
                pending = []
                for job in spec['jobs']:
                    if job['control'] != (phase == 'nominal_control'): continue
                    folder = out/job['name']
                    if (folder/'result.json').exists(): continue
                    if folder.exists():
                        dest = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                        dest.parent.mkdir(parents=True, exist_ok=True); shutil.move(folder, dest)
                    pending.append(job)
                status('running', phase=phase)
                with ProcessPoolExecutor(args.workers) as pool:
                    futures = {pool.submit(worker, out, j, spec): j for j in pending}
                    for future in as_completed(futures):
                        future.result(); done += 1
                        status('running', phase=phase, last_job=futures[future]['name'])
                        print(f'{done}/{total}: {futures[future]["name"]} done', flush=True)
            summary = summarize(out, spec); status('finished', execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}'); raise


if __name__ == '__main__':
    raise SystemExit(main())
