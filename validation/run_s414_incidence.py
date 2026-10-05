#!/usr/bin/env python3
"""Two exact corrected controls, followed by two reported-incidence diagnostics."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
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
from run_s414_jacobian import ROOT, VERSIONS, prepare, read, sha, stamp, write, geometry_elements, record, transition_audit
from run_s414_jacobian import worker as fixed_lift_worker, controls
from amplification_jacobian import AmplificationJacobianSolver
from multifoil import save_results
from optimization import admissible

STUDY = 'validation/results/s414_incidence_20261003'
OLD = 'validation/results/s414_jacobian_20261003'
ANALYSIS = 'validation/results/s414_jacobian_analysis_20261003'
PREFLIGHT = 'validation/results/s414_incidence_preflight_20261003'


def specification():
    if sha(ROOT/OLD/'study.json') != '9eaafe63d282d2c89a3cfc24fcd726f26894d573bf42cc63911da9b0a0194e19':
        raise ValueError('Prior derivative study changed')
    old = read(ROOT/OLD/'study.json')
    report = read(ROOT/ANALYSIS/'summary.json'); receipt = read(ROOT/ANALYSIS/'verification.json')
    if not report['complete'] or not report['full_state_audit_pass'] or not receipt['complete']:
        raise ValueError('Prior scientific audit incomplete')
    files = {OLD+'/study.json', ANALYSIS+'/verification.json'}
    for name, digest in receipt['output_sha256'].items():
        path = ANALYSIS+'/'+name
        if sha(ROOT/path) != digest: raise ValueError('Audited analysis changed: '+path)
        files.add(path)
    for manifest in (ANALYSIS+'/input_manifest.json', PREFLIGHT+'/input_manifest.json'):
        files.add(manifest)
        for name, digest in read(ROOT/manifest)['files'].items():
            if sha(ROOT/name) != digest: raise ValueError('Audited input changed: '+name)
            files.add(name)
    for name in ('summary.json', 'pressure_at_taps.csv', 'pressure_metrics.csv'):
        files.add(PREFLIGHT+'/'+name)
    if not read(ROOT/PREFLIGHT/'summary.json')['complete']: raise ValueError('Incidence preflight incomplete')
    jobs = []
    for control in (True, False):
        for case in ('cl0056', 'cl028'):
            prior = next(j for j in old['jobs'] if j['name'] == case+'_corrected_af')
            job = dict(prior, name=case+('_control' if control else '_measured_alpha'), control=control,
                comparison_basis='fixed_lift' if control else 'fixed_alpha',
                reference_result=OLD+'/'+prior['name']+'/result.json')
            if not control:
                job['reference_experimental_cl'] = job.pop('target_cl')
                job['target_alpha'] = prior['experiment_alpha']; job['center'] = job['target_alpha']
                job['initialization_policy'] = 'single local cold solve at measured caption incidence'
            jobs.append(job)
    files.add('validation/S414_INCIDENCE_PROTOCOL.md')
    for folder in (ROOT, ROOT/'validation'):
        files.update(str(p.relative_to(ROOT)) for p in folder.glob('*.py'))
    spec = {k: v for k, v in old.items() if k not in ('jobs', 'sha256', 'gate', 'qualification', 'experiment', 'jacobian_check_sha256')}
    spec.update(experiment='S414 measured-incidence pressure and transition diagnostic', jobs=jobs,
        gate='Both corrected fixed-lift controls repeat exactly before measured-incidence jobs',
        qualification='Fixed-angle cases are not matched-lift drag validations; no angle fit or geometry change',
        sha256={name: sha(ROOT/name) for name in sorted(files)})
    return spec


def accepted_fixed_angle(sol, job, spec):
    return bool(admissible(sol) and np.isfinite(sol.residual) and sol.residual <= spec['residual_tolerance']
        and np.isfinite(sol.alpha) and abs(sol.alpha-job['target_alpha']) <= 1e-12)


def worker(out, job, spec):
    if job['comparison_basis'] == 'fixed_lift':
        return fixed_lift_worker(out, job, spec)
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        started = time.perf_counter(); elements = geometry_elements(job['scenario'], spec['panels'])
        result = dict(comparison_basis='fixed_alpha', target_alpha=job['target_alpha'],
            reference_experimental_cl=job['reference_experimental_cl'], scenario=job['scenario'],
            amplification_jacobian=job['amplification_jacobian'], ncrit=spec['ncrit'], wake_length=job['wake_length'],
            feasible=False, CL=None, CD=None, CM=None, alpha_deg=None, error='', records=[], attempts=[],
            flow_directory=str((folder/'flow').relative_to(ROOT)), reused=False)
        write(folder/'progress.json', result)
        solver = AmplificationJacobianSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
            amplification_jacobian=job['amplification_jacobian'], iterations=spec['iterations'],
            tolerance=spec['residual_tolerance'], wake_length=job['wake_length'], start_alpha=job['target_alpha'])
        sol = solver.solve(job['target_alpha'])
        sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit'])
        context = dict(index=0, path='cold_measured_incidence', phase='fixed_alpha')
        result['records'].append(dict(record(sol, 'fixed_alpha'), initialization=context))
        save_results(folder/'flow', [sol], elements, dict(study=spec, job=job,
            comparison_basis='fixed_alpha', target_alpha=job['target_alpha']), 'viscous')
        result['feasible'] = accepted_fixed_angle(sol, job, spec)
        result['attempts'] = [dict(path='cold_measured_incidence', start_alpha=job['target_alpha'],
            requested_alpha=job['target_alpha'], feasible=result['feasible'])]
        if result['feasible']:
            result.update(CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                selected_flow_index=0, selected_path='cold_measured_incidence',
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
        else:
            result['error'] = sol.diagnostics.get('failure_reason') or 'Fixed-angle admissibility or residual check failed'
        result['wall_seconds'] = time.perf_counter()-started
        write(folder/'result.json', result)
        return result


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json'); ok = r['feasible']
        if ok:
            selected = r['records'][r['selected_flow_index']]
            ok = bool(selected['admissible'] and selected['residual'] <= spec['residual_tolerance']
                and np.isfinite([r[k] for k in ('CL', 'CD', 'CM', 'alpha_deg')]).all()
                and abs(sum(v['CD_wake'] for v in r['element_forces'].values())-r['CD']) < 1e-12
                and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7
                    and t['station_residual_max'] <= 1e-7 for t in r['transition_audit']))
            if job['comparison_basis'] == 'fixed_lift': ok = ok and abs(r['CL']-job['target_cl']) <= spec['cl_tolerance']
            else: ok = ok and abs(r['alpha_deg']-job['target_alpha']) <= 1e-12 and r['comparison_basis'] == 'fixed_alpha'
        elif any(r[k] is not None for k in ('CL', 'CD', 'CM', 'alpha_deg')):
            raise ValueError('Failed target contains finite selected force')
        metadata = out/job['name']/'flow/metadata.json'
        if metadata.exists():
            meta = read(metadata)
            if meta['case']['study'] != spec or meta['case']['job'] != job:
                raise ValueError('Exported study/job differs')
            if any(d['amplification_jacobian'] != 'corrected_af' or d['ncrit_by_element'] != {'fore': 9., 'aft': 9.} for d in meta['diagnostics']):
                raise ValueError('Exported mode or Ncrit differs')
        elif r['feasible']: raise ValueError('Missing flow export')
        rows.append(dict(job=job['name'], comparison_basis=job['comparison_basis'], feasible=r['feasible'],
            execution_check_pass=bool(ok), CL=r['CL'], CD=r['CD'], alpha_deg=r['alpha_deg'],
            selected_path=r.get('selected_path'), wall_seconds=r['wall_seconds'], error=r['error']))
    check = controls(out, spec)
    result = dict(required_jobs=4, converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=check['passed'] and all(r['execution_check_pass'] for r in rows),
        controls=check, jobs=rows, qualification='Mixed comparison bases; scientific audit pending')
    write(out/'summary.json', result)
    return result


# Restart, locking, immutable runtime and gated orchestration below are adapted
# without numerical changes from the completed derivative comparison.

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
            print('Prepared four incidence targets; study sha256='+sha(out/'study.json')); return 0
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
            for phase in ('corrected_controls', 'measured_incidence'):
                if phase == 'measured_incidence':
                    control = controls(out, spec); write(out/'controls.json', control)
                    if not control['passed']:
                        status('execution_failed', error='Corrected controls failed; measured-incidence jobs not launched'); return 2
                pending = []
                for job in spec['jobs']:
                    if job['control'] != (phase == 'corrected_controls'): continue
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

if __name__ == '__main__':
    raise SystemExit(main())
