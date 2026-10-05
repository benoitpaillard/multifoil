#!/usr/bin/env python3
"""Nine guarded S414 reconstruction jobs at two measured pressure lifts."""
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

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from run_design_initialization import ROOT, VERSIONS, prepare, read, sha, stamp
from s702_ncrit_sweep import write, record, transition_audit
from s414_geometry import elements as geometry_elements, SCENARIOS, SHIFT, QUAL
from run_s702_wake_extent import CODE, paths_for
from fixed_lift_policy import solve_with_initialization
from multifoil import save_results
from viscous import CoupledSolver
import numpy as np

STUDY='validation/results/s414_reconstruction_20261002'
DRAG='validation/results/s414_drag_reference_20261002'
PARENT='validation/results/s702_matched_pressure_20261002'


def specification():
    pins={PARENT+'/study.json':'786752736ff78e82005be0fb22009bbdba50015db57d595ab33bc51126bb26c4',
        QUAL+'/verification.json':'a7e4257b8aebd5a5cb3f78096c59a937fa19fe82fe7b60039a8627c8eaa8eb0a',
        DRAG+'/source.json':'9aea45e09958ed077976f6fc1c3663c923606cbb716f92708bd2b25bc593ff5c'}
    for name,digest in pins.items():
        if sha(ROOT/name)!=digest:raise ValueError('Qualified source changed: '+name)
    previous=read(ROOT/PARENT/'study.json')
    for name in CODE:
        if name in previous['sha256'] and sha(ROOT/name)!=previous['sha256'][name]:
            raise ValueError('Archived numerical code changed: '+name)
    files=set(CODE)|set(pins)|{
        'validation/run_s702_wake_extent.py','validation/run_s414_reconstruction.py',
        'validation/s414_geometry.py','validation/S414_RECONSTRUCTION_PROTOCOL.md',
        'validation/extract_s414_drag.py','validation/qualify_s414_geometry.py','validation/extract_s702.py',
        'test_s414_reconstruction.py','test_fixed_lift_policy.py','test_optimization.py','test_paneling.py'}
    evidence=read(ROOT/QUAL/'verification.json')['output_sha256']
    names=['report_test_fore.dat','report_test_aft.dat','geometry_metrics.json','source_manifest.json',
        'pressure/experiment.csv','pressure/mses_2010_envelope.csv','pressure/source.json',
        'pressure/tap_coordinate_checks.csv','pressure_tap_geometry_support.csv','FINDINGS.md']
    for name in names:
        path=QUAL+'/'+name
        if sha(ROOT/path)!=evidence[name]:raise ValueError('Qualified geometry/pressure changed: '+name)
        files.add(path)
    source=read(ROOT/DRAG/'source.json')
    for name,digest in source['input_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Drag input changed: '+name)
        files.add(name)
    for name,digest in source['files_sha256'].items():
        path=DRAG+'/'+name
        if sha(ROOT/path)!=digest:raise ValueError('Drag extraction changed: '+name)
        files.add(path)
    jobs=[]
    for scenario in SCENARIOS:
        for target in source['targets']:
            jobs.append(dict(name=target['case']+'_'+scenario,pressure_case=target['case'],
                scenario=scenario,target_cl=target['CL'],experiment_alpha=target['alpha_deg'],
                center=target['mses_seed_alpha'],wake_length=1.,control=scenario=='nominal'))
    repeat=dict(jobs[0],name='cl028_nominal_repeat',repeat_of=jobs[0]['name'])
    jobs.insert(2,repeat)
    return dict(schema=1,experiment='S414 bounded reconstruction sensitivity at measured pressure lifts',
        jobs=jobs,reynolds=1e6,mach=.1,ncrit=9.,panels=[340,240],reference_chord=1.,
        model_chord_mm=457.2,aft_shift_c=SHIFT,iterations=100,cl_tolerance=2e-5,
        residual_tolerance=1e-7,initial_half_width=.05,alpha_bounds=[-6.,8.],
        max_evaluations=20,max_brackets=2,warm_start_alpha=-1.,entry_offset=-.4,
        workers_limit=4,expected_versions=VERSIONS,acceptance_control_CD_counts=.01,
        controls_require_exact_saved_arrays=True,
        qualification='Exploratory nominal reconstruction, not exact tested geometry or a probabilistic uncertainty bound; no Ncrit fit.',
        sha256={name:sha(ROOT/name) for name in sorted(files)})


def check_controls(out,spec):
    repeat=next(j for j in spec['jobs'] if 'repeat_of' in j)
    folder=out/repeat['name'];reference=out/repeat['repeat_of']
    actual,expected=read(folder/'result.json'),read(reference/'result.json')
    ok=bool(actual['feasible'] and expected['feasible']);delta=None;arrays=[]
    if ok:
        delta=(actual['CD']-expected['CD'])*1e4
        ok=abs(delta)<=spec['acceptance_control_CD_counts'] and actual['selected_path']==expected['selected_path']
        for element in ('fore','aft'):
            with np.load(folder/(element+'_state.npz')) as a,np.load(reference/(element+'_state.npz')) as b:
                if set(a.files)!=set(b.files):raise ValueError('Repeat state fields changed')
                for field in a.files:
                    same=np.array_equal(a[field],b[field]);ok=ok and same
                    arrays.append(dict(element=element,field=field,exactly_equal=same))
    return dict(passed=bool(ok),repeat=repeat['name'],reference=repeat['repeat_of'],
        CD_difference_counts=delta,array_comparisons=arrays)


def summarize(out,spec):
    rows=[]
    for job in spec['jobs']:
        result=read(out/job['name']/'result.json');ok=result['feasible']
        if ok:
            selected=result['records'][result['selected_flow_index']]
            ok=bool(selected['admissible'] and selected['residual'] is not None
                and selected['residual']<=spec['residual_tolerance']
                and np.isfinite([result[k] for k in ('CL','CD','CM','alpha_deg')]).all()
                and abs(result['CL']-job['target_cl'])<=spec['cl_tolerance']
                and abs(sum(v['CD_wake'] for v in result['element_forces'].values())-result['CD'])<1e-12)
            ok=ok and all(t['inside_cell'] and abs(t['amplification_residual'])<=1e-7
                and t['station_residual_max']<=1e-7 for t in result['transition_audit'])
        elif any(result[k] is not None for k in ('CL','CD','CM','alpha_deg')):
            raise ValueError('Failed target contains selected finite forces')
        metadata=out/job['name']/'flow/metadata.json'
        if metadata.exists():
            meta=read(metadata)
            if meta['case']['study']!=spec or meta['case']['job']!=job:raise ValueError('Exported scenario changed')
            if any(d['ncrit_by_element']!={'fore':9.,'aft':9.} for d in meta['diagnostics']):
                raise ValueError('Exported Ncrit changed')
        elif result['feasible']:raise ValueError('Missing selected-state exports')
        rows.append(dict(job=job['name'],scenario=job['scenario'],feasible=result['feasible'],
            execution_check_pass=bool(ok),CL=result['CL'],CD=result['CD'],alpha_deg=result['alpha_deg'],
            selected_path=result.get('selected_path'),wall_seconds=result['wall_seconds'],error=result['error']))
    controls=check_controls(out,spec)
    summary=dict(required_jobs=len(rows),converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=controls['passed'] and all(r['execution_check_pass'] for r in rows),
        controls=controls,jobs=rows,qualification='Execution gates only; full state/residual audit and physical interpretation pending.')
    write(out/'summary.json',summary);return summary


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = geometry_elements(job['scenario'], spec['panels'])
        def factory(alpha):
            return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
                iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
                wake_length=job['wake_length'], start_alpha=alpha)
        solutions = []
        result = dict(scenario=job['scenario'], ncrit=spec['ncrit'], target_cl=job['target_cl'],
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


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/STUDY)
    ap.add_argument('--workers',type=int,default=4);ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args();out=args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1<=args.workers<=4:
        ap.error('Use a results directory and 1–4 workers')
    out.mkdir(parents=True,exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        spec=specification();prepare(out,spec)
        if args.prepare_only:
            print('Prepared nine S414 reconstruction targets; study sha256='+sha(out/'study.json'));return 0
        runtime=dict(python=platform.python_version(),hostname=platform.node(),executable=sys.executable,
            versions={k:importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions']!=VERSIONS or runtime['python']!='3.12.13':
            raise ValueError('Use pinned Python 3.12.13 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json')!=runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json',runtime)
        done=sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state,**extra):
            write(out/'status.json',dict(state=state,updated_utc=stamp(),pid=os.getpid(),
                completed_jobs=done,total_jobs=len(spec['jobs']),**extra))
        status('running')
        try:
            # Exact repeat gates the sensitivity cases, including on resume.
            for phase in ('baselines','sensitivities'):
                if phase=='sensitivities':
                    checks=check_controls(out,spec)
                    write(out/'controls.json',checks)
                    if not checks['passed']:
                        status('execution_failed',error='Exact nominal repeat failed; sensitivities were not launched')
                        return 2
                pending=[]
                for j in spec['jobs']:
                    if bool(j['control']) != (phase=='baselines'):continue
                    folder=out/j['name']
                    if (folder/'result.json').exists():continue
                    if folder.exists():
                        dest=out/'interrupted_attempts'/stamp().replace(':','')/j['name']
                        dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(folder,dest)
                    pending.append(j)
                status('running',phase=phase)
                with ProcessPoolExecutor(args.workers) as pool:
                    futures={pool.submit(worker,out,j,spec):j for j in pending}
                    for future in as_completed(futures):
                        r=future.result();done+=1
                        status('running',phase=phase,last_job=futures[future]['name'],last_job_feasible=r['feasible'])
                        print(f"{done}/{len(spec['jobs'])}: {futures[future]['name']} feasible={r['feasible']}",flush=True)
            summary=summarize(out,spec);status('finished',execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed',error=f'{type(exc).__name__}: {exc}');raise


if __name__=='__main__':raise SystemExit(main())
