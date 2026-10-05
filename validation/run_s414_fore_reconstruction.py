#!/usr/bin/env python3
"""Four restartable fixed-lift jobs isolating fore coordinate reconstruction."""
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
from run_s414_reconstruction import ROOT, VERSIONS, prepare, read, sha, stamp, write, paths_for, record, transition_audit
from s414_geometry import elements as original_geometry
from multifoil import validate_geometry
from run_s414_lower_bucket import compare_control
from fixed_lift_policy import solve_with_initialization
from multifoil import save_results
from amplification_jacobian import AmplificationJacobianSolver, LEGACY

STUDY = 'validation/results/s414_fore_reconstruction_20261004'
OLD = 'validation/results/s414_lower_bucket_20261003'
SCENARIOS = ('nominal', 'fore_source_cubic')


def geometry_elements(scenario, panels):
    if scenario not in SCENARIOS:
        raise ValueError('Unknown fore reconstruction scenario')
    nominal = original_geometry('nominal', panels)
    if scenario == 'nominal':
        return nominal
    result = [original_geometry('source_cubic', panels)[0], nominal[1]]
    validate_geometry(result)
    return result


def specification():
    old = read(ROOT/OLD/'study.json')
    if sha(ROOT/OLD/'study.json') != '2c1d57eba8101dcbc8283f23032447bd18e656dc3f18812a8d0b4b17c82f8cec':
        raise ValueError('Archived comparison study changed')
    for name, digest in old['sha256'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('Historical input changed: '+name)
    receipt_name = 'validation/results/s414_lower_bucket_launch_20261003/retrieved_files_sha256.json'
    receipt = read(ROOT/receipt_name)
    jobs = []
    files = set(old['sha256']) | {OLD+'/study.json', receipt_name,
        'validation/S414_FORE_RECONSTRUCTION_PROTOCOL.md'}
    for scenario in SCENARIOS:
        for label in ('cl0056', 'cl0175'):
            prior_name = label+'_nominal'
            prior = next(j for j in old['jobs'] if j['name'] == prior_name)
            jobs.append(dict(prior, name=label+'_'+scenario, scenario=scenario,
                amplification_jacobian=LEGACY, control=scenario == 'nominal',
                reference_result=OLD+'/'+prior_name+'/result.json'))
            for suffix in ('result.json', 'fore_state.npz', 'aft_state.npz'):
                name = prior_name+'/'+suffix
                if sha(ROOT/OLD/name) != receipt[name]:
                    raise ValueError('Archived control changed: '+name)
                files.add(OLD+'/'+name)
    for folder in (ROOT, ROOT/'validation'):
        files.update(str(p.relative_to(ROOT)) for p in folder.glob('*.py'))
    spec = {k: v for k, v in old.items() if k not in ('jobs', 'sha256', 'reference_control', 'gate', 'qualification', 'workers_limit')}
    spec.update(experiment='S414 fore-only reconstruction at the lower-bucket transition event',
        jobs=jobs, workers_limit=2,
        gate='Both nominal jobs reproduce their own archived state exactly before fore reconstruction jobs',
        qualification='Conditional reconstruction sensitivity, not exact experimental geometry; audit required',
        sha256={name: sha(ROOT/name) for name in sorted(files)})
    return spec


def controls(out, spec):
    result = {j['name']: compare_control(out/j['name'], (ROOT/j['reference_result']).parent,
        spec['acceptance_control_CD_counts']) for j in spec['jobs'] if j['control']}
    return dict(passed=all(r['passed'] for r in result.values()), jobs=result)


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json'); ok = r['feasible']
        if ok:
            selected = r['records'][r['selected_flow_index']]
            ok = bool(selected['admissible'] and selected['residual'] <= spec['residual_tolerance']
                and abs(r['CL']-job['target_cl']) <= spec['cl_tolerance']
                and np.isfinite([r[k] for k in ('CL', 'CD', 'CM', 'alpha_deg')]).all()
                and abs(sum(v['CD_wake'] for v in r['element_forces'].values())-r['CD']) < 1e-12
                and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7
                    and t['station_residual_max'] <= 1e-7 for t in r['transition_audit']))
        elif any(r[k] is not None for k in ('CL', 'CD', 'CM', 'alpha_deg')):
            raise ValueError('Failed target contains finite selected force')
        metadata = out/job['name']/'flow/metadata.json'
        if metadata.exists():
            meta = read(metadata)
            if meta['case']['study'] != spec or meta['case']['job'] != job:
                raise ValueError('Exported study/job differs')
            if any(d['amplification_jacobian'] != job['amplification_jacobian']
                   or d['ncrit_by_element'] != {'fore': 9., 'aft': 9.} for d in meta['diagnostics']):
                raise ValueError('Exported mode or Ncrit differs')
        elif r['feasible']:
            raise ValueError('Missing flow export')
        rows.append(dict(job=job['name'], scenario=job['scenario'], mode=job['amplification_jacobian'], feasible=r['feasible'],
            execution_check_pass=ok, CL=r['CL'], CD=r['CD'], alpha_deg=r['alpha_deg'],
            selected_path=r.get('selected_path'), wall_seconds=r['wall_seconds'], error=r['error']))
    control = controls(out, spec)
    result = dict(required_jobs=4, converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=control['passed'] and all(r['execution_check_pass'] for r in rows),
        controls=control, jobs=rows, qualification='Execution only; full audit and interpretation pending')
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
            print('Prepared four fore reconstruction targets; study sha256='+sha(out/'study.json')); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != VERSIONS or runtime['python'] != '3.12.13':
            raise ValueError('Use pinned Python 3.12.13 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(),
                completed_jobs=done, total_jobs=4, **extra))
        try:
            for phase in ('nominal_controls', 'fore_reconstruction'):
                if phase == 'fore_reconstruction':
                    control = controls(out, spec); write(out/'controls.json', control)
                    if not control['passed']:
                        status('execution_failed', error='Nominal controls failed; reconstruction jobs not launched'); return 2
                pending = []
                for job in spec['jobs']:
                    if job['control'] != (phase == 'nominal_controls'): continue
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
                        result = future.result(); done += 1
                        status('running', phase=phase, last_job=futures[future]['name'], last_job_feasible=result['feasible'])
                        print(f'{done}/4: {futures[future]["name"]} feasible={result["feasible"]}', flush=True)
            summary = summarize(out, spec); status('finished', execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}'); raise


# The archived fixed-lift policy and solver settings are unchanged.
# Only the fore geometry input differs for the two intervention jobs.

def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = geometry_elements(job['scenario'], spec['panels'])
        def factory(alpha):
            return AmplificationJacobianSolver(elements, amplification_jacobian=job['amplification_jacobian'], reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
                iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
                wake_length=job['wake_length'], start_alpha=alpha)
        solutions = []
        result = dict(amplification_jacobian=job['amplification_jacobian'], scenario=job['scenario'], ncrit=spec['ncrit'], target_cl=job['target_cl'],
            wake_length=job['wake_length'], reused=False, feasible=False,
            CL=None, CD=None, CM=None, alpha_deg=None, error='', records=[], attempts=[],
            flow_directory=str((folder/'flow').relative_to(ROOT)))
        def checkpoint(sol, context):
            sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit'])
            solutions.append(sol)
            result['records'].append(dict(record(sol, context['phase']), initialization=context))
            save_results(folder/'flow', solutions, elements,
                         dict(study=spec, job=job, target_cl=job['target_cl']), 'viscous')
            write(folder/'progress.json', result)
        outcome = solve_with_initialization(factory, job['target_cl'], job['center'], spec['alpha_bounds'],
            paths_for(job, spec), half_width=spec['initial_half_width'], cl_tolerance=spec['cl_tolerance'],
            residual_tolerance=spec['residual_tolerance'], max_evaluations=spec['max_evaluations'],
            max_brackets=spec['max_brackets'], on_evaluation=checkpoint)
        result.update(attempts=outcome.attempts, error=outcome.error,
                      initialization_paths=[asdict(p) for p in paths_for(job, spec)])
        if outcome.feasible:
            sol, solver = outcome.solution, outcome.solver
            if abs(solver.last_alpha-sol.alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
                raise ValueError('Returned solver does not match selected state')
            result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                selected_flow_index=outcome.selected_flow_index, selected_path=outcome.selected_path,
                boundary_layers=sol.diagnostics['boundary_layers'], element_forces=sol.diagnostics['element_forces'],
                transition_audit=transition_audit(solver), wake_discretization={})
            stream = np.array([np.cos(np.deg2rad(sol.alpha)), np.sin(np.deg2rad(sol.alpha))])
            for model in solver.models:
                np.savez_compressed(folder/(model.geom.name+'_state.npz'), U=model.glob.U,
                    xy=model.foil.x, xi=model.isol.xi, turbulent=model.vsol.turb, wake_xy=model.wake.x,
                    lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]),
                    wake_s=model.wake.s, wake_t=model.wake.t)
                spacing = np.linalg.norm(np.diff(model.wake.x, axis=1), axis=0)
                result['wake_discretization'][model.geom.name] = dict(nodes=int(model.wake.N),
                    arclength=float(np.sum(spacing)), first_spacing=float(spacing[0]), last_spacing=float(spacing[-1]),
                    streamwise_start=float(stream@model.wake.x[:, 0]), streamwise_end=float(stream@model.wake.x[:, -1]))
        result['wall_seconds'] = time.perf_counter()-begun
        write(folder/'result.json', result)
        return result


if __name__ == '__main__':
    raise SystemExit(main())
