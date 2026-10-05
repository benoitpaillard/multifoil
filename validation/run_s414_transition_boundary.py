#!/usr/bin/env python3
"""V2 tranche 1: fixed-incidence paths through the S414 transition-boundary interval.

Six restartable jobs: two legacy endpoint controls replaying the archived
initialization, then four corrected_af paths (two fore reconstructions, both
incidence directions) on the same five angles. No fixed-lift search is used.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import csv
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
import sys
import time
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s414_geometry import elements as original_geometry, QUAL
from multifoil import save_results, validate_geometry
from optimization import admissible
from viscous import CoupledSolver, mf

STUDY = 'validation/results/s414_transition_boundary_20261004'
ARCHIVE = 'validation/results/s414_fore_reconstruction_20261004'
RECEIPT = 'validation/results/s414_fore_reconstruction_launch_20261004/retrieved_files_sha256.json'
PROTOCOL = 'validation/S414_TRANSITION_BOUNDARY_PROTOCOL.md'
INTEGRATION = 'validation/results/s414_default_integration_20261004/verification.json'
ARCHIVED_JOBS = {'nominal': 'cl0175_nominal', 'fore_source_cubic': 'cl0175_fore_source_cubic'}
ELEMENTS = ('fore', 'aft')
VERSIONS = {'numpy': '2.3.5', 'scipy': '1.17.0', 'matplotlib': '3.10.8'}
PYTHON = '3.12.13'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
# Inputs needed by the frozen preflight tests of this payload.
TEST_DATA = ['validation/results/s414_lower_bucket_20261003/cl0056_nominal/fore_state.npz',
             'examples/s702_incompressible_design.json', 'validation/data/s702/fore.dat',
             'validation/data/s702/aft_d5.dat', 'validation/data/nlr7301/NLR7301MainECARP.dat',
             'validation/data/nlr7301/NLR7301FlapECARP.dat', 'validation/data/nlr7301/LICENSE']


def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def stamp(): return datetime.now(timezone.utc).isoformat()


def write(path, value):
    path = Path(path); temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n'); temporary.replace(path)


def finite(value):
    return float(value) if value is not None and math.isfinite(value) else None


def geometry_elements(scenario, panels):
    """Nominal, or the archived fore-only source-cubic reconstruction; aft input unchanged."""
    if scenario not in ARCHIVED_JOBS:
        raise ValueError('Unknown fore reconstruction scenario')
    nominal = original_geometry('nominal', panels)
    if scenario == 'nominal':
        return nominal
    result = [original_geometry('source_cubic', panels)[0], nominal[1]]
    validate_geometry(result)
    return result


def control_files(index):
    return [f'flow/{kind}_{index:03d}_{element}.csv' for kind in ('surface', 'wake') for element in ELEMENTS]


def archived_endpoints():
    """Full-precision bracket ends and start incidence, read from verified archived exports."""
    receipt = read(ROOT/RECEIPT); found = {}; files = {RECEIPT}
    for scenario, job in ARCHIVED_JOBS.items():
        names = ['result.json', 'flow/polar.csv', 'flow/metadata.json']+control_files(0)+control_files(1)
        for name in names:
            if sha(ROOT/ARCHIVE/job/name) != receipt[job+'/'+name]:
                raise ValueError('Archived endpoint changed: '+job+'/'+name)
            files.add(ARCHIVE+'/'+job+'/'+name)
        result = read(ROOT/ARCHIVE/job/'result.json'); records = result['records'][:2]
        for i, rec in enumerate(records):
            if rec['initialization'] != dict(index=i, path='cold_local', bracket=0, phase='fixed_lift') or not rec['converged']:
                raise ValueError('Archived endpoint is not a converged recorded cold-local evaluation')
        path = result['initialization_paths'][0]
        if path['name'] != 'cold_local' or path['seed_alphas']:
            raise ValueError('Archived initialization differs')
        found[scenario] = (records[0]['alpha_deg'], records[1]['alpha_deg'], path['start_alpha'])
    if len(set(found.values())) != 1:
        raise ValueError('Archived endpoint incidences differ between geometries')
    return found['nominal'], sorted(files)


def specification():
    (lo, hi, start), archive_files = archived_endpoints()
    if not lo < start < hi:
        raise ValueError('Archived start incidence lies outside its bracket')
    integration = read(ROOT/INTEGRATION)
    if not (integration['complete'] and integration['gate_passed']):
        raise ValueError('V1 integration gate has not passed; no nonlinear run is admitted')
    old = read(ROOT/ARCHIVE/'study.json')
    # Endpoints are the archived values; interior angles divide the bracket equally.
    angles = [lo]+[lo+(hi-lo)*k/4 for k in (1, 2, 3)]+[hi]
    jobs = []
    for scenario in ARCHIVED_JOBS:
        jobs.append(dict(name=scenario+'_legacy_endpoints', scenario=scenario, amplification_jacobian='legacy',
            control=True, direction='recorded', start_alpha=start, wake_length=1.,
            archive=ARCHIVE+'/'+ARCHIVED_JOBS[scenario],
            requests=[dict(alpha_deg=a, role='control_endpoint') for a in (lo, hi)]))
    for scenario in ARCHIVED_JOBS:
        common = dict(scenario=scenario, amplification_jacobian='corrected_af', control=False, start_alpha=start, wake_length=1.)
        jobs.append(dict(common, name=scenario+'_corrected_increasing', direction='increasing',
            requests=[dict(alpha_deg=a, role='path') for a in angles]))
        # Recorded initialization of the upper endpoint passes through the lower one.
        jobs.append(dict(common, name=scenario+'_corrected_decreasing', direction='decreasing',
            requests=[dict(alpha_deg=lo, role='endpoint_seed')]+[dict(alpha_deg=a, role='path') for a in angles[::-1]]))
    files = set(archive_files) | {PROTOCOL, INTEGRATION, ARCHIVE+'/study.json', 'vendor/mfoil.py', 'vendor/__init__.py',
        QUAL+'/report_test_fore.dat', QUAL+'/report_test_aft.dat', QUAL+'/verification.json'} | set(TEST_DATA)
    evidence = read(ROOT/QUAL/'verification.json')['output_sha256']
    for name in ('report_test_fore.dat', 'report_test_aft.dat'):
        if sha(ROOT/QUAL/name) != evidence[name]:
            raise ValueError('Qualified geometry changed: '+name)
    for folder in (ROOT, ROOT/'validation'):
        files.update(str(p.relative_to(ROOT)) for p in folder.glob('*.py'))
    settings = {k: old[k] for k in ('reynolds', 'mach', 'ncrit', 'panels', 'reference_chord', 'model_chord_mm',
                                    'aft_shift_c', 'iterations', 'residual_tolerance')}
    return dict(schema=1, experiment='S414 fixed-incidence transition-boundary paths, V2 tranche 1',
        **settings, continuation_step=1., min_step=.25, minimum_te_gap=1e-5,
        bracket=dict(lower_alpha_deg=lo, upper_alpha_deg=hi, archived_start_alpha_deg=start),
        angles_deg=angles, jobs=jobs, workers_limit=2, expected_versions=VERSIONS, expected_python=PYTHON,
        requested_path_points=20, requested_control_points=4, endpoint_seed_evaluations=2,
        controls_require_identical_exports=True,
        equivalence_screens=dict(CL=2e-5, CD_counts=.01, Cp=1e-4, transition_x_over_c=1e-4),
        gate='Both legacy endpoint controls reproduce the archived exports exactly before any corrected path starts',
        qualification='Numerical path diagnostic at fixed incidence; not matched-lift drag validation, '
                      'not physical hysteresis, not exact tested geometry; audit required',
        sha256={name: sha(ROOT/name) for name in sorted(files)})


def prepare(out, spec):
    out.mkdir(parents=True, exist_ok=True); manifest = out/'study.json'
    if manifest.exists():
        if read(manifest) != spec:
            raise ValueError('Code, inputs or settings changed; use a new study directory')
    else:
        if any(p.name != 'run.lock' for p in out.iterdir()):
            raise ValueError('Nonempty output directory without a study manifest')
        write(manifest, spec)
    for name, digest in spec['sha256'].items():
        saved = out/'source_snapshot'/name
        if not saved.exists():
            saved.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/name, saved)
        if sha(saved) != digest:
            raise ValueError('Source snapshot changed: '+name)


def transition_audit(solver):
    """Subcell e^N equation at each transition cell; same check as s702_transition_study."""
    rows = []
    for model in solver.models:
        for side, label in enumerate(('lower', 'upper')):
            ids = np.asarray(model.vsol.Is[side], dtype=int); flags = model.vsol.turb[ids]
            for i in np.flatnonzero(np.diff(flags.astype(int)) == 1)+1:
                js = ids[[i-1, i]]; x = model.isol.xi[js]; u = model.glob.U[:, js]
                param = mf.build_param(model, side); param.turb = False
                r, _, _ = mf.residual_transition(model, param, x, u, np.zeros(2))
                xt = model.vsol.xt
                ut = u[:, 0]+(xt-x[0])/(x[1]-x[0])*(u[:, 1]-u[:, 0]); ut[2] = param.ncrit
                d1, _ = mf.get_damp(u[:, 0], param); dt, _ = mf.get_damp(ut, param)
                rows.append(dict(element=model.geom.name, side=label, bracket=x.tolist(), transition_xi=float(xt),
                    inside_cell=bool(x[0] <= xt <= x[1]),
                    amplification_residual=float(param.ncrit-u[2, 0]-.5*(xt-x[0])*(d1+dt)),
                    station_residual_max=float(np.max(abs(r)))))
    return rows


def reached_by(initialized, previous, alpha, trail):
    """State the actual start of each requested point; a warm start is never called cold."""
    if initialized:
        return f'continuation from previous converged incidence {previous!r}'
    if not trail:
        return 'not reached'
    if trail[0]['alpha_deg'] == alpha:
        return 'cold solve at requested incidence'
    if trail[0]['converged']:
        return f"continuation from cold solve at start incidence {trail[0]['alpha_deg']!r}"
    return f"cold solve at requested incidence after failed cold solve at start incidence {trail[0]['alpha_deg']!r}"


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = geometry_elements(job['scenario'], spec['panels'])
        solver = CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
            iterations=spec['iterations'], tolerance=spec['residual_tolerance'], wake_length=job['wake_length'],
            minimum_te_gap=spec['minimum_te_gap'], continuation_step=spec['continuation_step'],
            min_step=spec['min_step'], start_alpha=job['start_alpha'],
            amplification_jacobian=job['amplification_jacobian'])
        if solver.amplification_jacobian != job['amplification_jacobian']:
            raise ValueError('Solver mode differs from declared job')
        solutions = []
        result = dict(job=job['name'], scenario=job['scenario'], amplification_jacobian=job['amplification_jacobian'],
            direction=job['direction'], ncrit=spec['ncrit'], complete=False, records=[], error='',
            flow_directory=str((folder/'flow').relative_to(ROOT)))
        (folder/'states').mkdir(exist_ok=True)
        for index, request in enumerate(job['requests']):
            alpha = request['alpha_deg']
            initialized, previous = bool(solver.initialized), solver.last_alpha
            sol = solver.solve(alpha)
            sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit'])
            solutions.append(sol); diag = sol.diagnostics
            entry = dict(index=index, role=request['role'], requested_alpha_deg=alpha, alpha_deg=sol.alpha,
                CL=finite(sol.cl), CD=finite(sol.cd), CM=finite(sol.cm), converged=bool(sol.converged),
                admissible=bool(admissible(sol)), residual=finite(sol.residual),
                failure_reason=diag.get('failure_reason', ''), iterations_final_step=diag['iterations'],
                solver_initialized_before=initialized, previous_converged_alpha_deg=previous,
                reached_by=reached_by(initialized, previous, alpha, diag['continuation']),
                approach=('cold' if not initialized else 'increasing' if alpha > previous else 'decreasing'),
                # Strict JSON: a failed internal step has no lift, never NaN.
                continuation=[dict(alpha_deg=t['alpha_deg'], converged=bool(t['converged']), CL=finite(t['CL']))
                              for t in diag['continuation']],
                wall_seconds=diag['wall_seconds'], state_files=[])
            if sol.converged:
                if abs(solver.last_alpha-alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
                    raise ValueError('Solver state does not match the returned solution')
                entry.update(boundary_layers=diag['boundary_layers'], element_forces=diag['element_forces'],
                    transition_audit=transition_audit(solver), wake_discretization={})
                stream = np.array([np.cos(np.deg2rad(sol.alpha)), np.sin(np.deg2rad(sol.alpha))])
                for model in solver.models:
                    name = f'states/point_{index:02d}_{model.geom.name}.npz'
                    np.savez_compressed(folder/name, U=model.glob.U, xy=model.foil.x, xi=model.isol.xi,
                        turbulent=model.vsol.turb, wake_xy=model.wake.x, lower=np.asarray(model.vsol.Is[0]),
                        upper=np.asarray(model.vsol.Is[1]), wake_s=model.wake.s, wake_t=model.wake.t)
                    entry['state_files'].append(name)
                    spacing = np.linalg.norm(np.diff(model.wake.x, axis=1), axis=0)
                    entry['wake_discretization'][model.geom.name] = dict(nodes=int(model.wake.N),
                        arclength=float(np.sum(spacing)), first_spacing=float(spacing[0]), last_spacing=float(spacing[-1]),
                        streamwise_start=float(stream@model.wake.x[:, 0]), streamwise_end=float(stream@model.wake.x[:, -1]))
            elif any(entry[k] is not None for k in ('CL', 'CD', 'CM')):
                raise ValueError('Failed evaluation carries a finite force')
            result['records'].append(entry)
            save_results(folder/'flow', solutions, elements, dict(study=spec, job=job), 'viscous')
            write(folder/'progress.json', result)
            print(f"{index}: alpha={alpha!r} converged={sol.converged} CL={entry['CL']} CD={entry['CD']} "
                  f"{entry['reached_by']}", flush=True)
            first_path = request['role'] == 'path' and not any(r['role'] == 'path' for r in result['records'][:-1])
            if (request['role'] == 'endpoint_seed' or first_path) and not sol.converged:
                # No substitute seed: a path without its own endpoint state stays explicitly unevaluated,
                # so later points cannot be approached from the wrong side or by unrelated cold starts.
                result['error'] = 'Endpoint seed failed; path not started'
                break
        result.update(complete=True, wall_seconds=time.perf_counter()-begun)
        write(folder/'result.json', result)
        return result


def table(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def compare_control(folder, job):
    """Exact comparison with the archived endpoint exports; no tolerance in the pinned runtime."""
    archive = ROOT/job['archive']; result = read(folder/'result.json'); rows = []; passed = True
    new_meta, old_meta = read(folder/'flow/metadata.json'), read(archive/'flow/metadata.json')
    new_polar, old_polar = table(folder/'flow/polar.csv'), table(archive/'flow/polar.csv')
    for index, rec in enumerate(result['records']):
        same = bool(rec['converged'] and new_polar[index] == old_polar[index])
        rows.append(dict(evaluation=index, item='polar row', identical=same)); passed = passed and same
        for key in ('iterations', 'continuation', 'history', 'boundary_layers', 'element_forces'):
            same = new_meta['diagnostics'][index].get(key) == old_meta['diagnostics'][index].get(key)
            if key == 'continuation':
                same = ([(t['alpha_deg'], t['converged'], t['CL'] if t['converged'] else None) for t in new_meta['diagnostics'][index][key]]
                        == [(t['alpha_deg'], t['converged'], t['CL'] if t['converged'] else None) for t in old_meta['diagnostics'][index][key]])
            rows.append(dict(evaluation=index, item='diagnostic '+key, identical=bool(same))); passed = passed and same
        for name in control_files(index):
            same = (folder/name).exists() and sha(folder/name) == sha(archive/name)
            rows.append(dict(evaluation=index, item=name, identical=bool(same))); passed = passed and same
    return dict(passed=bool(passed and len(result['records']) == 2), reference=job['archive'], comparisons=rows)


def controls(out, spec):
    result = {j['name']: compare_control(out/j['name'], j) for j in spec['jobs'] if j['control']}
    return dict(passed=all(r['passed'] for r in result.values()), jobs=result)


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json')
        meta = read(out/job['name']/'flow/metadata.json')
        if meta['case']['study'] != spec or meta['case']['job'] != job:
            raise ValueError('Exported study/job differs')
        for rec, diag in zip(r['records'], meta['diagnostics']):
            if (diag['amplification_jacobian'] != job['amplification_jacobian']
                    or diag['ncrit_by_element'] != {'fore': 9., 'aft': 9.}):
                raise ValueError('Exported mode or Ncrit differs')
            ok = True
            if rec['converged']:
                ok = bool(rec['admissible'] and rec['residual'] <= spec['residual_tolerance']
                    and abs(rec['alpha_deg']-rec['requested_alpha_deg']) <= 1e-12
                    and abs(sum(v['CD_wake'] for v in rec['element_forces'].values())-rec['CD']) < 1e-12
                    and all((out/job['name']/name).exists() for name in rec['state_files']) and len(rec['state_files']) == 2
                    and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7
                            and t['station_residual_max'] <= 1e-7 for t in rec['transition_audit']))
            elif any(rec[k] is not None for k in ('CL', 'CD', 'CM')):
                raise ValueError('Failed evaluation contains a finite force')
            rows.append(dict(job=job['name'], scenario=job['scenario'], mode=job['amplification_jacobian'],
                direction=job['direction'], index=rec['index'], role=rec['role'], alpha_deg=rec['requested_alpha_deg'],
                converged=rec['converged'], execution_check_pass=ok, CL=rec['CL'], CD=rec['CD'],
                reached_by=rec['reached_by'], iterations_final_step=rec['iterations_final_step'],
                failure_reason=rec['failure_reason']))
        requested = len(job['requests']); evaluated = len(r['records'])
        for request in job['requests'][evaluated:]:
            rows.append(dict(job=job['name'], scenario=job['scenario'], mode=job['amplification_jacobian'],
                direction=job['direction'], index=None, role=request['role'], alpha_deg=request['alpha_deg'],
                converged=False, execution_check_pass=True, CL=None, CD=None, reached_by='not evaluated: endpoint seed failed',
                iterations_final_step=None, failure_reason=r['error']))
    control = controls(out, spec); path = [r for r in rows if r['role'] == 'path']
    result = dict(required_jobs=len(spec['jobs']), completed_jobs=len(spec['jobs']),
        requested_path_points=len(path), converged_path_points=sum(r['converged'] for r in path),
        failed_path_points=sum(not r['converged'] for r in path),
        control_points=sum(r['role'] == 'control_endpoint' for r in rows),
        endpoint_seed_evaluations=sum(r['role'] == 'endpoint_seed' for r in rows),
        execution_checks_pass=bool(control['passed'] and all(r['execution_check_pass'] for r in rows)),
        controls=control, points=rows,
        qualification='Execution only; nonconverged points are findings, not forces. Full audit and interpretation pending')
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
            print('Prepared six transition-boundary jobs; study sha256='+sha(out/'study.json')); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k: os.environ[k] for k in THREADS})
        if runtime['versions'] != VERSIONS or runtime['python'] != PYTHON:
            raise ValueError('Use pinned Python 3.12.13 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        total = len(spec['jobs']); done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(),
                completed_jobs=done, total_jobs=total, **extra))
        try:
            for phase in ('legacy_endpoint_controls', 'corrected_paths'):
                if phase == 'corrected_paths':
                    control = controls(out, spec); write(out/'controls.json', control)
                    if not control['passed']:
                        status('execution_failed', error='Legacy endpoint controls differ from archive; corrected paths not launched'); return 2
                pending = []
                for job in spec['jobs']:
                    if job['control'] != (phase == 'legacy_endpoint_controls'): continue
                    folder = out/job['name']
                    if (folder/'result.json').exists(): continue
                    if folder.exists():
                        # A continuation path is not a portable checkpoint: restart the whole job.
                        dest = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                        dest.parent.mkdir(parents=True, exist_ok=True); shutil.move(folder, dest)
                    pending.append(job)
                status('running', phase=phase)
                with ProcessPoolExecutor(args.workers) as pool:
                    futures = {pool.submit(worker, out, j, spec): j for j in pending}
                    for future in as_completed(futures):
                        result = future.result(); done += 1
                        status('running', phase=phase, last_job=futures[future]['name'],
                               last_job_converged_points=sum(r['converged'] for r in result['records']),
                               last_job_points=len(result['records']))
                        print(f'{done}/{total}: {futures[future]["name"]} '
                              f'{sum(r["converged"] for r in result["records"])}/{len(result["records"])} converged', flush=True)
            summary = summarize(out, spec); status('finished', execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}'); raise


if __name__ == '__main__':
    raise SystemExit(main())
