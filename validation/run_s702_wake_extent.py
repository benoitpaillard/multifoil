#!/usr/bin/env python3
"""Run/resume five frozen S702 wake-extent jobs or monitor their remote progress."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
from dataclasses import asdict
import fcntl
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from run_design_initialization import ROOT, VERSIONS, prepare, read, sha, stamp
from s702_ncrit_sweep import write, tag, geometry, repanel, record, transition_audit
from multifoil import save_results
from viscous import CoupledSolver
from fixed_lift_policy import InitializationPath, solve_with_initialization
import numpy as np

OLD = 'validation/results/s702_initialization_policy_20260921'
AUDIT = 'validation/results/s702_initialization_policy_analysis'
STUDY = 'validation/results/s702_wake_extent_20260922'
REMOTE = '/home/jony/multifoil-runs/s702-wake-extent-20260922'
CODE = ['multifoil.py', 'viscous.py', 'optimization.py', 'paneling.py', 'vendor/mfoil.py',
        'fixed_lift_policy.py', 'shape_optimization.py', 'validation/run_design_initialization.py',
        'validation/s702_ncrit_sweep.py', 'validation/s702_geometry.py',
        'validation/s702_benchmark.py', 'validation/s702_transition_study.py']


def specification():
    if sha(ROOT/OLD/'study.json') != '5388fef33ad7833c66b9d8db43da4c6304a54a4d50a4aae38f1da2a2d27219cb':
        raise ValueError('Unexpected archived 1c study')
    previous = read(ROOT/OLD/'study.json')
    if not read(ROOT/AUDIT/'summary.json')['full_state_audit_pass']:
        raise ValueError('Archived selected-state audit did not pass')
    audit = read(ROOT/AUDIT/'input_manifest.json')['files']
    for name in CODE:
        if name in previous['sha256'] and sha(ROOT/name) != previous['sha256'][name]:
            raise ValueError('Verified numerical source changed: '+name)
    files = CODE+['validation/run_s702_wake_extent.py', 'validation/S702_WAKE_EXTENT_PROTOCOL.md',
        'test_s702_wake_extent.py', 'test_fixed_lift_policy.py', 'test_optimization.py', 'test_paneling.py',
        OLD+'/study.json', OLD+'/runtime.json', AUDIT+'/summary.json', AUDIT+'/input_manifest.json']
    for name in ('fore.dat', 'aft_d5.dat', 'aft_d10.dat', 'source.json'):
        name = 'validation/data/s702/'+name
        if sha(ROOT/name) != previous['sha256'][name]:
            raise ValueError('Archived geometry source changed: '+name)
        files.append(name)
    jobs = []
    for flap in (5, 10):
        for cl in (.8, 1.2):
            prior = next(j for j in previous['jobs'] if j['flap'] == flap and j['target_cl'] == cl and not j['control'])
            folder = ROOT/OLD/prior['name']
            ref = read(folder/'result.json')
            if not ref['feasible'] or abs(ref['CL']-cl) > 2e-5:
                raise ValueError('Missing feasible 1c reference')
            # Preserve all audited selected inputs, including full saved Newton states.
            for name in audit:
                if Path(name).is_relative_to(Path(OLD)/prior['name']):
                    if sha(ROOT/name) != audit[name]['sha256']:
                        raise ValueError('Audited reference changed: '+name)
                    files.append(name)
            jobs.append(dict(name=f'd{flap}_cl{tag(cl)}_wake2', flap=flap, target_cl=cl,
                center=prior['center'], wake_length=2., control=False,
                reference_result=str((folder/'result.json').relative_to(ROOT))))
    control = dict(next(j for j in jobs if j['flap'] == 10 and j['target_cl'] == .8))
    control.update(name='d10_cl0p8_wake1_control', wake_length=1., control=True)
    jobs.append(control)
    return dict(schema=1, experiment='S702 operational wake extent: four 2c cases and one fresh 1c control',
        jobs=jobs, reynolds=700000., mach=.07, ncrit=6.25, panels=[260, 180], reference_chord=1.,
        iterations=100, cl_tolerance=2e-5, residual_tolerance=1e-7, initial_half_width=.05,
        alpha_bounds=[-6., 8.], max_evaluations=20, max_brackets=2, warm_start_alpha=-1.,
        entry_offset=-.4, acceptance_control_CD_counts=.01, design_screening_counts=.5,
        expected_versions=VERSIONS, archived_wake_length=1.,
        qualification='Extent and derived wake-node discretization sensitivity; not independent aerodynamic validation.',
        selection_protocol='First successful cold-local, warm-anchor, then explicit-entry path; no measured drag selection.',
        sha256={name: sha(ROOT/name) for name in sorted(set(files))})


def paths_for(job, spec):
    return [InitializationPath('cold_local', job['center']),
            InitializationPath('warm_anchor', spec['warm_start_alpha']),
            InitializationPath('warm_entry', spec['warm_start_alpha'], (job['center']+spec['entry_offset'],))]


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = [(name, repanel(xy, count)) for (name, xy), count in zip(geometry(job['flap']), spec['panels'])]
        def factory(alpha):
            return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
                iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
                wake_length=job['wake_length'], start_alpha=alpha)
        solutions = []
        result = dict(flap=job['flap'], ncrit=spec['ncrit'], target_cl=job['target_cl'],
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


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        result, ref = read(out/job['name']/'result.json'), read(ROOT/job['reference_result'])
        ok = result['feasible']
        if ok:
            selected = result['records'][result['selected_flow_index']]
            if (not selected['admissible'] or selected['residual'] is None or
                selected['residual'] > spec['residual_tolerance'] or
                not np.isfinite([result[k] for k in ('CL', 'CD', 'CM', 'alpha_deg')]).all() or
                abs(result['CL']-job['target_cl']) > spec['cl_tolerance'] or
                abs(sum(v['CD_wake'] for v in result['element_forces'].values())-result['CD']) > 1e-12):
                raise ValueError('Selected state violates admissibility/force criteria')
        elif any(result[k] is not None for k in ('CL', 'CD', 'CM', 'alpha_deg')):
            raise ValueError('Failed target has selected finite force/incidence')
        delta = (result['CD']-ref['CD'])*1e4 if ok else None
        rows.append(dict(job=job['name'], control=job['control'], flap=job['flap'], target_cl=job['target_cl'],
            feasible=ok, selected_path=result.get('selected_path'), difference_from_1c_CD_counts=delta,
            control_pass=bool(ok and abs(delta) <= spec['acceptance_control_CD_counts']) if job['control'] else None,
            evaluations=len(result['records']), failed_evaluations=sum(not r['admissible'] for r in result['records']),
            wall_seconds=result['wall_seconds']))
    complete = all(r['feasible'] for r in rows)
    controls = [r for r in rows if r['control']]
    control_pass = len(controls) == 1 and all(r['control_pass'] for r in controls)
    summary = dict(completed_utc=stamp(), required_jobs=len(rows), feasible_jobs=sum(r['feasible'] for r in rows),
        coverage_complete=complete, control_pass=control_pass, execution_checks_pass=complete and control_pass,
        jobs=rows, qualification='Full selected-state audit and physical interpretation remain required.')
    write(out/'summary.json', summary)
    return summary


PROBE = r'''
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); out=root/sys.argv[2]
def read(p): return json.loads(p.read_text()) if p.exists() else {}
status=read(out/'status.json'); spec=read(out/'study.json')
try:
    cmd=Path('/proc/'+str(status.get('pid'))+'/cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive='run_s702_wake_extent.py' in cmd and sys.argv[2] in cmd
except OSError: alive=False
jobs=[]
for j in spec.get('jobs',[]):
    f=out/j['name']; r=read(f/'result.json'); p=r or read(f/'progress.json'); records=p.get('records',[])
    jobs.append(dict(job=j['name'], state=('CONVERGED' if r['feasible'] else 'FAILED') if r else ('RUNNING' if f.exists() else 'PENDING'),
        evaluations=len(records), failed=sum(not a['admissible'] for a in records),
        path=p.get('selected_path') or (records[-1]['initialization']['path'] if records else ''), error=p.get('error','')))
print(json.dumps(dict(status=status,alive=alive,jobs=jobs,summary=read(out/'summary.json'))))
'''


def remote_status(watch=False, as_json=False):
    try:
        while True:
            try:
                proc = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
                    '-o', 'StrictHostKeyChecking=yes', 'jony4-anywhere',
                    'python3 - '+shlex.quote(REMOTE)+' '+shlex.quote(STUDY)],
                    input=PROBE, text=True, capture_output=True, timeout=25)
                if proc.returncode: raise RuntimeError(proc.stderr.strip())
                data = json.loads(proc.stdout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
                print('Status connection failed: '+str(exc), flush=True)
                if not watch: return 3
                time.sleep(10); continue
            status = data['status']; state = status.get('state', 'not started')
            if as_json: print(json.dumps(data, indent=2), flush=True)
            else:
                print(f"\nS702 wake extent / jony4: {state}; {status.get('completed_jobs', 0)}/5 completed")
                for row in data['jobs']:
                    print(f"{row['job']:27s} {row['state']:10s} solves={row['evaluations']:2d} failed={row['failed']:2d} {row['path']}")
                    if row['error']: print('  '+row['error'])
                if data['summary']:
                    s = data['summary']
                    print(f"Feasible: {s['feasible_jobs']}/5; fresh 1c control: {'PASS' if s['control_pass'] else 'FAIL'}")
                    print('Selected-state audit and physical interpretation remain required.')
                if status.get('error'): print(status['error'])
                if state not in ('finished', 'execution_failed') and not data['alive']:
                    print('Runner inactive; inspect launch.log before resuming.')
                print('Log: '+REMOTE+'/launch.log', flush=True)
            if state == 'finished': return 0 if data['summary'].get('execution_checks_pass') else 2
            if state == 'execution_failed' or not data['alive']: return 1
            if not watch: return 0
            time.sleep(10)
    except KeyboardInterrupt:
        print('\nMonitor stopped; remote job is unaffected.')
        return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/STUDY)
    ap.add_argument('--workers', type=int, default=5)
    ap.add_argument('--prepare-only', action='store_true')
    ap.add_argument('--remote-status', action='store_true')
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()
    if args.remote_status: return remote_status(args.watch, args.json)
    if args.watch or args.json: ap.error('--watch and --json require --remote-status')
    out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 5:
        ap.error('Use a results directory and 1–5 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(json.dumps(dict(state='prepared', jobs=len(spec['jobs']), study_sha256=sha(out/'study.json'))))
            return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != VERSIONS or sys.version_info[:2] != (3, 12):
            raise ValueError('Use pinned Python 3.12 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(),
                completed_jobs=done, total_jobs=len(spec['jobs']), **extra))
        status('running')
        try:
            pending = []
            for job in spec['jobs']:
                folder = out/job['name']
                if (folder/'result.json').exists(): continue
                if folder.exists():
                    target = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                    target.parent.mkdir(parents=True, exist_ok=True); shutil.move(folder, target)
                pending.append(job)
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(worker, out, j, spec):j for j in pending}
                for future in as_completed(futures):
                    r = future.result(); done += 1
                    status('running', last_job=futures[future]['name'], last_job_feasible=r['feasible'])
                    print(f"{done}/{len(spec['jobs'])} {futures[future]['name']} feasible={r['feasible']}", flush=True)
            summary = summarize(out, spec)
            status('finished', execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}')
            raise


if __name__ == '__main__':
    raise SystemExit(main())
