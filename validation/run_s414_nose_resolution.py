#!/usr/bin/env python3
"""Frozen six-case nose-spacing diagnostic; replay incidence histories, never retarget lift."""
import argparse
import ast
import contextlib
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
for key in THREADS:
    os.environ[key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
STUDY = 'validation/results/s414_nose_resolution_20261005'
LAUNCH = 'validation/results/s414_nose_resolution_launch_20261005'
PROTOCOL = 'validation/S414_NOSE_RESOLUTION_PROTOCOL.md'
ARCHIVE = 'validation/results/s414_jacobian_20261003'
RECEIPT = 'validation/results/s414_jacobian_launch_20261003/retrieved_files_sha256.json'
QUAL = 'validation/results/s414_geometry_qualification_20261002'
MODE = 'corrected_af'
VERSIONS = {'numpy': '2.3.5', 'scipy': '1.17.0', 'matplotlib': '3.10.8'}
PYTHON = '3.12.13'
STATE_FIELDS = ('U', 'xy', 'xi', 'turbulent', 'wake_xy', 'lower', 'upper', 'wake_s', 'wake_t')
SOURCE_FILES = ('multifoil.py', 'viscous.py', 'paneling.py', 'amplification_jacobian.py',
    'diagnostic_discretization.py', 'optimization.py', 'validation/s414_geometry.py',
    'validation/run_s414_transition_boundary.py', 'validation/s414_nose_resolution.py',
    'validation/run_s414_nose_resolution.py', 'validation/prepare_s414_nose_resolution_payload.py',
    'validation/s414_nose_resolution_remote.py', 'test_s414_nose_resolution.py',
    'test_s414_nose_resolution_runner.py', 'validation/check_s414_nose_resolution.py',
    'validation/run_s414_nose_resolution_launch.py')


def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def stamp(): return datetime.now(timezone.utc).isoformat()
def finite(x): return float(x) if x is not None and math.isfinite(x) else None


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n'); temp.replace(path)


def archived_requests(result):
    """Freeze the successful archived path without adaptively reproducing its lift search."""
    index = result['selected_flow_index']; rows = result['records'][:index+1]
    if not result['feasible'] or result['selected_path'] != 'cold_local' or index != 4 or len(rows) != 5:
        raise ValueError('Unexpected archived selected path or request budget')
    if any(not r['converged'] or r['initialization']['path'] != 'cold_local' for r in rows):
        raise ValueError('Archived replay contains a failed or different initialization path')
    if rows[-1]['alpha_deg'] != result['alpha_deg']:
        raise ValueError('Selected archive is not the final replay request')
    path = result['initialization_paths'][0]
    if path['name'] != 'cold_local' or path['seed_alphas']:
        raise ValueError('Unexpected archived cold initialization')
    return path['start_alpha'], [dict(alpha_deg=r['alpha_deg'], role='selected' if i == index else 'initialization',
        archived_index=i, archived_initialization=r['initialization']) for i, r in enumerate(rows)]


def local_dependencies(seeds):
    """Pin the local import closure, including imports inside audit functions; never import historical runners."""
    pending = list(seeds); found = set()
    while pending:
        name = pending.pop()
        if name in found: continue
        found.add(name)
        if not name.endswith('.py'): continue
        tree = ast.parse((ROOT/name).read_text(), filename=name)
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import): modules = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules = [node.module]+[node.module+'.'+a.name for a in node.names]
            for module in modules:
                relative = Path(*module.split('.'))
                for base in (ROOT, ROOT/'validation'):
                    for candidate in (base/relative.with_suffix('.py'), base/relative/'__init__.py'):
                        if candidate.is_file(): pending.append(str(candidate.relative_to(ROOT)))
    return found


def specification():
    old = read(ROOT/ARCHIVE/'study.json'); receipt = read(ROOT/RECEIPT)
    if sha(ROOT/ARCHIVE/'study.json') != receipt['study.json']:
        raise ValueError('Archived study manifest changed')
    expected = dict(reynolds=1e6, mach=.1, ncrit=9., panels=[340, 240], iterations=100, residual_tolerance=1e-7)
    if any(old[k] != value for k, value in expected.items()) or old['expected_versions'] != VERSIONS:
        raise ValueError('Archived baseline differs from the admitted settings')
    geometry_receipt = read(ROOT/QUAL/'verification.json')
    for name in ('report_test_fore.dat', 'report_test_aft.dat'):
        if sha(ROOT/QUAL/name) != geometry_receipt['output_sha256'][name]:
            raise ValueError('Qualified geometry input changed: '+name)
    files = {PROTOCOL, ARCHIVE+'/study.json', RECEIPT, QUAL+'/report_test_fore.dat',
             QUAL+'/report_test_aft.dat', QUAL+'/verification.json', 'vendor/mfoil.py', 'vendor/__init__.py'}
    jobs = []
    for pressure_case, reference in [('cl0056', 'cl0056_corrected_af'), ('cl0284', 'cl028_corrected_af')]:
        archived_job = next(j for j in old['jobs'] if j['name'] == reference)
        if archived_job['amplification_jacobian'] != MODE or archived_job['scenario'] != 'nominal':
            raise ValueError('Archive physics or geometry differs')
        for name in ('result.json', 'fore_state.npz', 'aft_state.npz'):
            relative = f'{ARCHIVE}/{reference}/{name}'
            if sha(ROOT/relative) != receipt[reference+'/'+name]:
                raise ValueError('Archived corrected control changed: '+relative)
            files.add(relative)
        result = read(ROOT/ARCHIVE/reference/'result.json'); start, requests = archived_requests(result)
        for factor in (1, 2, 4):
            jobs.append(dict(name=f'{pressure_case}_factor{factor}', pressure_case=pressure_case,
                subdivision_factor=factor, control=factor == 1, start_alpha=start, requests=requests,
                alpha_deg=result['alpha_deg'], archived_target_cl=archived_job['target_cl'],
                reference_result=f'{ARCHIVE}/{reference}/result.json', amplification_jacobian=MODE))
    sources = list(SOURCE_FILES)
    for name in ('validation/audit_s414_nose_resolution.py', 'validation/report_s414_nose_resolution.py',
                 'test_s414_nose_resolution_analysis.py'):
        if (ROOT/name).exists(): sources.append(name)
    files.update(local_dependencies(sources))
    settings = {k: old[k] for k in ('reynolds', 'mach', 'ncrit', 'panels', 'reference_chord', 'model_chord_mm',
                                   'iterations', 'residual_tolerance')}
    return dict(schema=1, experiment='S414 local coupled fore-lower nose resolution', **settings,
        wake_length=1., minimum_te_gap=1e-5, continuation_step=1., min_step=.25, jobs=jobs,
        expected_python=PYTHON, expected_versions=VERSIONS, workers_limit=2,
        requested_selected_cases=6, requested_evaluations=30, attempts_per_job=1,
        case_timeout_seconds=1800, study_timeout_seconds=5400,
        gate='Both factor1 controls reproduce archived forces and all18 saved state arrays exactly before refinements',
        qualification='Fixed-incidence resolution diagnostic on the nominal reconstructed polygon; not matched-lift validation',
        sha256={n: sha(ROOT/n) for n in sorted(files)})


def check_inputs(spec):
    for name, digest in spec['sha256'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('Frozen source/input changed: '+name)


def prepare(out, spec):
    out.mkdir(parents=True, exist_ok=True)
    if (out/'study.json').exists():
        if read(out/'study.json') != spec:
            raise ValueError('Study differs; use a new revision/directory')
    elif any(p.name != 'run.lock' for p in out.iterdir()):
        raise ValueError('Nonempty output without study manifest')
    else:
        write(out/'study.json', spec)
    check_inputs(spec)
    for name, digest in spec['sha256'].items():
        destination = out/'source_snapshot'/name
        if not destination.exists():
            destination.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/name, destination)
        if sha(destination) != digest:
            raise ValueError('Frozen snapshot differs: '+name)


def save_states(folder, solver, prefix=''):
    import numpy as np
    saved = []
    for model in solver.models:
        name = f'{prefix}{model.geom.name}_state.npz'
        np.savez_compressed(folder/name, U=model.glob.U, xy=model.foil.x, xi=model.isol.xi,
            turbulent=model.vsol.turb, wake_xy=model.wake.x, lower=np.asarray(model.vsol.Is[0]),
            upper=np.asarray(model.vsol.Is[1]), wake_s=model.wake.s, wake_t=model.wake.t)
        saved.append(name)
    return saved


def worker(out, job, spec):
    from s414_nose_resolution import Solver, build_elements
    from run_s414_transition_boundary import transition_audit, reached_by
    from multifoil import save_results
    from optimization import admissible
    folder = out/job['name']; folder.mkdir(exist_ok=True); (folder/'states').mkdir(exist_ok=True)
    started = time.monotonic()
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        result = dict(job=job['name'], complete=False, execution_state='running', feasible=False, records=[],
            CL=None, CD=None, CM=None, alpha_deg=None, error='', amplification_jacobian=MODE,
            subdivision_factor=job['subdivision_factor'], qualification=spec['qualification'])
        try:
            solver = Solver(build_elements(), subdivision_factor=job['subdivision_factor'],
                reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'], iterations=spec['iterations'],
                tolerance=spec['residual_tolerance'], wake_length=spec['wake_length'], minimum_te_gap=spec['minimum_te_gap'],
                continuation_step=spec['continuation_step'], min_step=spec['min_step'], start_alpha=job['start_alpha'],
                amplification_jacobian=MODE)
            result['resolution_controls'] = solver.resolution_controls
            result['surface_node_sha256'] = {m.geom.name: hashlib.sha256(m.foil.x.tobytes()).hexdigest() for m in solver.models}
            write(folder/'progress.json', result); solutions = []
            for index, request in enumerate(job['requests']):
                initialized, previous = bool(solver.initialized), solver.last_alpha
                sol = solver.solve(request['alpha_deg']); diag = sol.diagnostics
                diag['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit']); solutions.append(sol)
                row = dict(index=index, role=request['role'], requested_alpha_deg=request['alpha_deg'], alpha_deg=sol.alpha,
                    CL=finite(sol.cl), CD=finite(sol.cd), CM=finite(sol.cm), converged=bool(sol.converged),
                    admissible=bool(admissible(sol)), residual=finite(sol.residual),
                    failure_reason=diag.get('failure_reason', ''), iterations_final_step=diag['iterations'],
                    solver_initialized_before=initialized, previous_converged_alpha_deg=previous,
                    reached_by=reached_by(initialized, previous, sol.alpha, diag['continuation']),
                    continuation=[dict(alpha_deg=t['alpha_deg'], converged=bool(t['converged']), CL=finite(t['CL']))
                                  for t in diag['continuation']], wall_seconds=diag['wall_seconds'], state_files=[])
                if sol.converged:
                    if solver.last_alpha != sol.alpha or solver.last_solution.cl != sol.cl:
                        raise ValueError('Solver state does not match returned point')
                    row.update(boundary_layers=diag['boundary_layers'], element_forces=diag['element_forces'],
                        transition_audit=transition_audit(solver),
                        state_files=['states/'+n for n in save_states(folder/'states', solver, f'point_{index:02d}_')])
                elif any(row[k] is not None for k in ('CL', 'CD', 'CM')):
                    raise ValueError('Failed evaluation carries finite force')
                result['records'].append(row)
                save_results(folder/'flow', solutions, solver.elements, dict(study=spec, job=job), 'viscous')
                write(folder/'progress.json', result)
                print(f"request {index+1}/5 alpha={sol.alpha!r} converged={sol.converged} CL={row['CL']} CD={row['CD']}", flush=True)
                if not row['converged'] or not row['admissible'] or row['residual'] > spec['residual_tolerance']:
                    result['error'] = 'Recorded incidence path failed; remaining requests not evaluated'; break
                if request['role'] == 'selected':
                    result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                        boundary_layers=diag['boundary_layers'], element_forces=diag['element_forces'],
                        transition_audit=row['transition_audit'], selected_flow_index=index,
                        state_files=save_states(folder, solver))
            result.update(complete=True, execution_state='complete')
        except Exception as exc:
            result.update(complete=True, execution_state='failed', error=f'{type(exc).__name__}: {exc}')
            import traceback
            traceback.print_exc()
        result['wall_seconds'] = time.monotonic()-started
        result['artifact_sha256'] = {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*'))
            if p.is_file() and p.name not in ('result.json', 'progress.json', 'run.log', 'process.log')}
        write(folder/'result.json', result)
    return result


def compare_control(folder, job, root=ROOT):
    import numpy as np
    actual = read(folder/'result.json'); reference = root/job['reference_result']; expected = read(reference)
    passed = actual['feasible'] and expected['feasible']; arrays = []
    forces = {k: actual.get(k) == expected[k] for k in ('alpha_deg', 'CL', 'CD', 'CM')}
    passed = passed and all(forces.values())
    if actual['feasible']:
        for element in ('fore', 'aft'):
            with np.load(folder/(element+'_state.npz')) as a, np.load(reference.parent/(element+'_state.npz')) as b:
                if set(a.files) != set(STATE_FIELDS) or set(b.files) != set(STATE_FIELDS):
                    raise ValueError('Saved state schema differs')
                for field in STATE_FIELDS:
                    same = bool(np.array_equal(a[field], b[field])); passed = passed and same
                    arrays.append(dict(element=element, field=field, exactly_equal=same))
    return dict(passed=bool(passed), exact_forces=forces, array_comparisons=arrays, reference=str(job['reference_result']))


def verify_result(folder, spec):
    result = read(folder/'result.json')
    if not result.get('complete'):
        raise ValueError('Incomplete result cannot be resumed as completed')
    for name, digest in result.get('artifact_sha256', {}).items():
        if sha(folder/name) != digest:
            raise ValueError('Completed artifact changed: '+str(folder/name))
    return result


def run_phase(out, spec, jobs, workers, deadline, status):
    """One attempt per job. On restart, completed work is reused; interrupted work is never silently retried."""
    pending = []; active = {}
    for job in jobs:
        folder = out/job['name']
        if (folder/'result.json').exists():
            verify_result(folder, spec); continue
        if folder.exists():
            raise ValueError('Interrupted attempt preserved; explicit new revision required: '+job['name'])
        pending.append(job)
    try:
        while pending or active:
            if time.monotonic() > deadline:
                raise TimeoutError('Frozen study wall-clock limit reached')
            while pending and len(active) < workers:
                job = pending.pop(0); folder = out/job['name']; folder.mkdir()
                write(folder/'attempt.json', dict(job=job['name'], attempt=1, started_utc=stamp(), study_sha256=sha(out/'study.json')))
                log = (folder/'process.log').open('ab', buffering=0)
                process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--out', str(out), '--worker', job['name']],
                    cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                write(folder/'attempt.json', dict(job=job['name'], attempt=1, started_utc=stamp(),
                    study_sha256=sha(out/'study.json'), pid=process.pid, timeout_seconds=spec['case_timeout_seconds']))
                log.close(); active[job['name']] = (process, time.monotonic(), job)
            for name, (process, begun, job) in list(active.items()):
                if time.monotonic()-begun > spec['case_timeout_seconds']:
                    raise TimeoutError('Frozen case wall-clock limit reached: '+name)
                code = process.poll()
                if code is not None:
                    if code or not (out/name/'result.json').exists():
                        raise RuntimeError(f'Worker exited without usable result: {name}, exit={code}')
                    verify_result(out/name, spec); del active[name]
            status('running', active_jobs=list(active), pending_jobs=[j['name'] for j in pending])
            if active:
                time.sleep(2)
    finally:
        for process, _, _ in active.values():
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
        for process, _, _ in active.values():
            try: process.wait(timeout=5)
            except subprocess.TimeoutExpired: os.killpg(process.pid, signal.SIGKILL); process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/STUDY)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--worker')
    args = parser.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 2:
        parser.error('Use a results directory and at most two workers')
    if args.worker:
        spec = read(out/'study.json')
        # Each detached worker retains its own timeout even if the supervising process dies.
        signal.setitimer(signal.ITIMER_REAL, spec['case_timeout_seconds'])
        check_inputs(spec)
        job = next(j for j in spec['jobs'] if j['name'] == args.worker)
        worker(out, job, spec); return 0
    def interrupted(signum, frame):
        raise TimeoutError(f'Supervisor interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(json.dumps(dict(study_sha256=sha(out/'study.json'), cases=6, requested_evaluations=30))); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS}, thread_limits={k: os.environ[k] for k in THREADS})
        if runtime['python'] != PYTHON or runtime['versions'] != VERSIONS:
            raise ValueError('Pinned numerical environment differs')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed')
        write(out/'runtime.json', runtime)
        budget_path = out/'budget.json'
        if budget_path.exists():
            budget = read(budget_path)
            if budget['study_sha256'] != sha(out/'study.json') or budget['limit_seconds'] != spec['study_timeout_seconds']:
                raise ValueError('Persisted compute budget differs')
        else:
            budget = dict(started_utc=stamp(), deadline_unix=time.time()+spec['study_timeout_seconds'],
                limit_seconds=spec['study_timeout_seconds'], study_sha256=sha(out/'study.json'))
            write(budget_path, budget)
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(), total_jobs=6,
                completed_jobs=sum((out/j['name']/'result.json').exists() for j in spec['jobs']), **extra))
        try:
            deadline = time.monotonic()+max(0., budget['deadline_unix']-time.time())
            controls = [j for j in spec['jobs'] if j['control']]
            run_phase(out, spec, controls, args.workers, deadline, status)
            gates = {j['name']: compare_control(out/j['name'], j) for j in controls}
            gate = dict(passed=all(v['passed'] for v in gates.values()), jobs=gates)
            write(out/'controls.json', gate)
            if not gate['passed']:
                status('execution_failed', error='Exact archived baseline gate failed; refinements not launched'); return 2
            run_phase(out, spec, [j for j in spec['jobs'] if not j['control']], args.workers, deadline, status)
            results = [verify_result(out/j['name'], spec) for j in spec['jobs']]
            execution_ok = all(r['execution_state'] == 'complete' for r in results)
            summary = dict(completed_jobs=len(results), required_jobs=6, converged_selected_cases=sum(r['feasible'] for r in results),
                requested_selected_cases=6, requested_evaluations=30, evaluated_requests=sum(len(r['records']) for r in results),
                execution_checks_pass=execution_ok, controls_passed=True, evidence_audit='pending', physical_interpretation='pending',
                qualification=spec['qualification'])
            write(out/'summary.json', summary)
            status('finished' if execution_ok else 'execution_failed', execution_checks_pass=execution_ok)
            return 0 if execution_ok else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}'); raise


if __name__ == '__main__': raise SystemExit(main())
