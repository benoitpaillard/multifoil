#!/usr/bin/env python3
"""Run/resume the predeclared seven-job common-initialization regression."""
import argparse
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import fcntl
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

from run_s702_repeatability import specification as baseline_specification, prepare, stamp, EXPECTED_VERSIONS
from s702_ncrit_sweep import ROOT, read, write, sha, tag, geometry, repanel, record, transition_audit
from report_s702_ncrit_sweep import csv_write
from multifoil import save_results
from viscous import CoupledSolver
from fixed_lift_policy import InitializationPath, solve_with_initialization
import numpy as np


def specification():
    # Keep the original source/calibration guards; the flow core is unchanged.
    original = baseline_specification()
    old = ROOT/'validation/results/s702_repeatability_20260915/study.json'
    if read(old) != original:
        raise ValueError('Original repeatability manifest differs from current baseline sources')
    refs = original['references']
    jobs = [dict(name=f"policy_d{r['flap']}_cl{tag(r['target_cl'])}",
                 flap=r['flap'], target_cl=r['target_cl'], center=r['local_seed_alpha'],
                 control=False, reference_result=r['result_path']) for r in refs]
    replay = next(j for j in original['jobs'] if j['path'] == 'warm_replay')
    ref = next(r for r in refs if r['flap'] == 10 and r['target_cl'] == .8)
    jobs.append(dict(name='replay_control_d10_cl0p8', flap=10, target_cl=.8,
                     center=ref['local_seed_alpha'], control=True,
                     reference_result=ref['result_path'], seed_alphas=replay['warm_seed_alphas']))
    sources = dict(original['sha256'])
    for name in ('fixed_lift_policy.py', 'test_fixed_lift_policy.py',
                 'validation/run_s702_initialization_policy.py',
                 'validation/S702_INITIALIZATION_POLICY_PROTOCOL.md',
                 'validation/results/s702_repeatability_20260915/study.json'):
        sources[name] = sha(ROOT/name)
    return dict(schema=1, experiment='S702 cold-then-warm initialization policy regression',
                jobs=jobs, references=refs, reynolds=700000., mach=.07, ncrit=6.25,
                panels=[260, 180], reference_chord=1., wake_length=1., iterations=100,
                cl_tolerance=2e-5, residual_tolerance=1e-7, initial_half_width=.05,
                alpha_bounds=[-6., 8.], max_evaluations=20, max_brackets=2,
                warm_start_alpha=-1., acceptance_CD_counts=.01,
                expected_versions=EXPECTED_VERSIONS,
                protocol='Declared cold local start then independent warm anchor; separate exact replay control',
                selection_protocol='First successful path; numerical lift only for incidence seeds',
                sha256=sources)


def paths_for(job, spec):
    if job['control']:
        return [InitializationPath('archived_replay', spec['warm_start_alpha'], tuple(job['seed_alphas']))]
    return [InitializationPath('cold_local', job['center']),
            InitializationPath('warm_anchor', spec['warm_start_alpha'])]


def worker(out, job, spec):
    folder = out/job['name']
    folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = [(name, repanel(xy, count)) for (name, xy), count in
                    zip(geometry(job['flap']), spec['panels'])]
        def factory(alpha):
            return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'],
                ncrit=spec['ncrit'], iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
                wake_length=spec['wake_length'], start_alpha=alpha)
        solutions = []
        result = dict(flap=job['flap'], ncrit=spec['ncrit'], target_cl=job['target_cl'],
                      reused=False, feasible=False, CL=None, CD=None, CM=None, alpha_deg=None,
                      error='', records=[], attempts=[],
                      flow_directory=str((folder/'flow').relative_to(ROOT)))
        def checkpoint(sol, context):
            sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit'])
            solutions.append(sol)
            result['records'].append(dict(record(sol, context['phase']), initialization=context))
            save_results(folder/'flow', solutions, elements,
                         dict(study=spec, job=job, target_cl=job['target_cl']), 'viscous')
            write(folder/'progress.json', result)
        outcome = solve_with_initialization(factory, job['target_cl'], job['center'],
            spec['alpha_bounds'], paths_for(job, spec), half_width=spec['initial_half_width'],
            cl_tolerance=spec['cl_tolerance'], residual_tolerance=spec['residual_tolerance'],
            max_evaluations=spec['max_evaluations'], max_brackets=spec['max_brackets'],
            on_evaluation=checkpoint)
        result.update(attempts=outcome.attempts, error=outcome.error,
                      initialization_paths=[asdict(p) for p in paths_for(job, spec)])
        if outcome.feasible:
            sol, solver = outcome.solution, outcome.solver
            if abs(solver.last_alpha-sol.alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
                raise ValueError('Returned solver does not match selected state')
            result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                selected_flow_index=outcome.selected_flow_index, selected_path=outcome.selected_path,
                boundary_layers=sol.diagnostics['boundary_layers'], element_forces=sol.diagnostics['element_forces'],
                transition_audit=transition_audit(solver))
            for model in solver.models:
                np.savez_compressed(folder/(model.geom.name+'_state.npz'), U=model.glob.U,
                    xy=model.foil.x, xi=model.isol.xi, turbulent=model.vsol.turb, wake_xy=model.wake.x,
                    lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]))
        result['wall_seconds'] = time.perf_counter()-begun
        write(folder/'result.json', result)
        return result


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        result = read(out/job['name']/'result.json')
        ref = read(ROOT/job['reference_result'])
        ok = result['feasible']
        if ok:
            selected = result['records'][result['selected_flow_index']]
            if (not selected['admissible'] or selected['residual'] > spec['residual_tolerance'] or
                    abs(result['CL']-job['target_cl']) > spec['cl_tolerance'] or
                    abs(sum(v['CD_wake'] for v in result['element_forces'].values())-result['CD']) > 1e-12):
                raise ValueError('Selected result violates regression criteria')
        delta = (result['CD']-ref['CD'])*1e4 if ok else None
        rows.append(dict(job=job['name'], control=job['control'], flap=job['flap'],
            target_cl=job['target_cl'], feasible=ok, selected_path=result.get('selected_path'),
            attempts=len(result['attempts']), requested_evaluations=len(result['records']),
            failed_requested_evaluations=sum(not r['admissible'] for r in result['records']),
            wall_seconds=result['wall_seconds'], difference_CD_counts=delta,
            force_regression_pass=ok and abs(delta) <= spec['acceptance_CD_counts']))
    csv_write(out/'comparison.csv', rows)
    result = dict(completed_utc=stamp(), required_primary=6,
        converged_primary=sum(r['feasible'] for r in rows if not r['control']),
        required_control=1, converged_control=sum(r['feasible'] for r in rows if r['control']),
        force_regression_pass=all(r['force_regression_pass'] for r in rows),
        requested_evaluations=sum(r['requested_evaluations'] for r in rows),
        failed_requested_evaluations=sum(r['failed_requested_evaluations'] for r in rows),
        qualification='Implementation regression only. Audit selected full states before adopting policy in optimization.')
    write(out/'summary.json', result)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 7:
        ap.error('Use a results directory and 1–7 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(json.dumps(dict(state='prepared', jobs=len(spec['jobs']), study_sha256=sha(out/'study.json'))))
            return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(),
            versions={n: importlib.metadata.version(n) for n in EXPECTED_VERSIONS},
            executable=sys.executable,
            thread_limits={n: os.environ[n] for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != EXPECTED_VERSIONS or sys.version_info[:2] != (3, 12):
            raise ValueError('Use Python 3.12 and the pinned numerical dependencies')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        pending = []
        for job in spec['jobs']:
            folder = out/job['name']
            if (folder/'result.json').exists():
                continue
            if folder.exists():
                archived = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                archived.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(folder, archived)
            pending.append(job)
        done = len(spec['jobs'])-len(pending)
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), completed_jobs=done,
                                          total_jobs=len(spec['jobs']), **extra))
        status('running')
        try:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(worker, out, job, spec): job for job in pending}
                for future in as_completed(futures):
                    result = future.result(); done += 1
                    status('running', last_completed_job=futures[future]['name'],
                           last_job_success=result['feasible'])
            summary = summarize(out, spec)
            status('finished', force_regression_pass=summary['force_regression_pass'])
            return 0 if summary['force_regression_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}')
            raise


if __name__ == '__main__':
    raise SystemExit(main())
