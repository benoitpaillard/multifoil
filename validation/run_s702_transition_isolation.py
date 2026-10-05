#!/usr/bin/env python3
"""Ten guarded S702 jobs separating fore and aft transition sensitivity."""
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
from s702_ncrit_sweep import write, geometry, repanel, record, transition_audit
from run_s702_wake_extent import CODE, paths_for
from fixed_lift_policy import solve_with_initialization
from multifoil import save_results
from viscous import CoupledSolver
import numpy as np

STUDY='validation/results/s702_transition_isolation_20261002'
PARENT='validation/results/s702_matched_pressure_20261002'
ANALYSIS='validation/results/s702_matched_pressure_analysis_20261002'
GEO='validation/results/s702_geometry_reconciliation_20261002'
QUAL='validation/results/s702_pressure_qualification_20261002'


def specification():
    pins={PARENT+'/study.json':'786752736ff78e82005be0fb22009bbdba50015db57d595ab33bc51126bb26c4',
        ANALYSIS+'/verification.json':'f6cbfc3276bdd9a683cf209c169dac108d684cf77590b1c32615b3969adb6c7a',
        ANALYSIS+'/input_manifest.json':'f7abaa8c2fa8bccc418677bf604c81ef4808c9f034db63a7ba083044e9db2b18'}
    for name,digest in pins.items():
        if sha(ROOT/name)!=digest:raise ValueError('Audited baseline changed: '+name)
    old=read(ROOT/PARENT/'study.json');audited=read(ROOT/ANALYSIS/'input_manifest.json')['files']
    for name,digest in old['sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Original code/input changed: '+name)
    verification=read(ROOT/ANALYSIS/'verification.json')
    for name in ('summary.json','residual_recheck.json'):
        if sha(ROOT/ANALYSIS/name)!=verification['output_sha256'][name]:raise ValueError('Audit evidence changed')
    summary=read(ROOT/ANALYSIS/'summary.json')
    if not (summary['full_state_audit_pass'] and summary['full_residual_recomputation_pass'] and summary['converged_targets']==8):
        raise ValueError('All archived states must pass full audit')
    files=set(old['sha256'])|set(pins)|{PARENT+'/runtime.json',ANALYSIS+'/summary.json',ANALYSIS+'/residual_recheck.json',
        'validation/run_s702_transition_isolation.py','validation/S702_TRANSITION_ISOLATION_PROTOCOL.md',
        'test_s702_transition_isolation.py','validation/qualify_s702_geometry.py'}
    source=read(ROOT/GEO/'source_manifest.json')
    for name,digest in source['input_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Geometry source changed: '+name)
        files.add(name)
    for name,digest in source['files_sha256'].items():
        path=GEO+'/'+name
        if sha(ROOT/path)!=digest:raise ValueError('Geometry qualification changed: '+path)
        files.add(path)
    files.update({GEO+'/source_manifest.json',GEO+'/FINDINGS.md',GEO+'/public_access_manifest.json'})
    for access in read(ROOT/GEO/'public_access_manifest.json')['sources']:
        if access['download_status']=='complete':
            path=GEO+'/'+access['file']
            if sha(ROOT/path)!=access['sha256']:raise ValueError('Public source snapshot changed')
            files.add(path)
    references=[]
    for j in old['jobs']:
        folder=ROOT/PARENT/j['name'];r=read(folder/'result.json')
        if not r['feasible']:raise ValueError('Missing archived endpoint')
        # Preserve every consumed export and selected state for later full audit.
        selected=r['selected_flow_index']
        needed=['result.json','fore_state.npz','aft_state.npz','flow/metadata.json','flow/polar.csv',
            'flow/geometry_fore.dat','flow/geometry_aft.dat']
        needed += [f'flow/{kind}_{selected:03d}_{e}.csv' for kind in ('surface','wake') for e in ('fore','aft')]
        for name in needed:
            path=str((folder/name).relative_to(ROOT))
            if path not in audited or sha(ROOT/path)!=audited[path]:raise ValueError('Audited reference changed: '+path)
            files.add(path)
        references.append(dict(job=j['name'],pressure_case=j['pressure_case'],ncrit=j['ncrit'],result=PARENT+'/'+j['name']+'/result.json'))
    jobs=[]
    for case in read(ROOT/QUAL/'source.json')['cases']:
        endpoints=[j for j in old['jobs'] if j['pressure_case']==case['name']]
        center=sum(j['center'] for j in endpoints)/2
        for label,nvalues in [('fore_only',dict(fore=6.25,aft=9.)),('aft_only',dict(fore=9.,aft=6.25))]:
            jobs.append(dict(name=case['name']+'_'+label,pressure_case=case['name'],flap=case['flap'],
                target_cl=case['target_cl'],experiment_alpha=case['alpha'],center=center,
                ncrit_by_element=nvalues,wake_length=1.,control=False,
                center_source_jobs=[j['name'] for j in endpoints]))
    controls=[]
    for name in ('d5_a2p03_n9','d5_a6p11_n6p25'):
        j=next(v for v in old['jobs'] if v['name']==name)
        controls.append(dict(name=name+'_control',pressure_case=j['pressure_case'],flap=j['flap'],
            target_cl=j['target_cl'],experiment_alpha=j['experiment_alpha'],center=j['center'],
            ncrit_by_element=dict(fore=j['ncrit'],aft=j['ncrit']),wake_length=1.,control=True,
            reference_result=PARENT+'/'+name+'/result.json'))
    settings={k:old[k] for k in ('reynolds','mach','panels','reference_chord','iterations','cl_tolerance',
        'residual_tolerance','initial_half_width','alpha_bounds','max_evaluations','max_brackets',
        'warm_start_alpha','entry_offset','workers_limit','expected_versions')}
    return dict(schema=1,experiment='S702 2x2 element-specific transition sensitivity at measured pressure lifts',
        jobs=controls+jobs,references=references,**settings,acceptance_control_CD_counts=.01,
        controls_require_exact_saved_arrays=True,
        qualification='Counterfactual transition sensitivity on fixed reconstructed geometry; no measured Ncrit, wake-turbulence model or new calibration.',
        sha256={name:sha(ROOT/name) for name in sorted(files)})


def make_solver(elements,job,spec,alpha):
    values=job['ncrit_by_element']
    if set(values)!={'fore','aft'} or any(not np.isfinite(v) or v<=0 for v in values.values()):
        raise ValueError('Invalid declared per-element thresholds')
    solver=CoupledSolver(elements,reynolds=spec['reynolds'],mach=spec['mach'],ncrit=values['fore'],
        iterations=spec['iterations'],tolerance=spec['residual_tolerance'],wake_length=job['wake_length'],start_alpha=alpha)
    if {m.geom.name for m in solver.models}!=set(values):raise ValueError('Solver element names differ')
    for model in solver.models:model.param.ncrit=values[model.geom.name]
    return solver


def check_controls(out,spec):
    rows=[]
    for j in spec['jobs']:
        if not j['control']:continue
        folder=out/j['name'];r=read(folder/'result.json');refpath=ROOT/j['reference_result'];ref=read(refpath)
        differences=[];ok=bool(r['feasible'] and ref['feasible']);delta=None
        if ok:
            delta=(r['CD']-ref['CD'])*1e4
            ok=abs(delta)<=spec['acceptance_control_CD_counts'] and r['selected_path']==ref['selected_path']
            for element in ('fore','aft'):
                with np.load(folder/(element+'_state.npz')) as actual,np.load(refpath.parent/(element+'_state.npz')) as expected:
                    if set(actual.files)!=set(expected.files):raise ValueError('Control state fields changed')
                    for field in actual.files:
                        a,b=actual[field],expected[field];same=np.array_equal(a,b)
                        differences.append(dict(element=element,field=field,exactly_equal=same))
                        ok=ok and same
        rows.append(dict(job=j['name'],passed=bool(ok),CD_difference_counts=delta,array_comparisons=differences))
    return dict(passed=len(rows)==2 and all(r['passed'] for r in rows),controls=rows)


def summarize(out,spec):
    rows=[]
    for j in spec['jobs']:
        r=read(out/j['name']/'result.json');ok=r['feasible']
        if ok:
            selected=r['records'][r['selected_flow_index']]
            ok=bool(selected['admissible'] and selected['residual']<=spec['residual_tolerance']
                and abs(r['CL']-j['target_cl'])<=spec['cl_tolerance']
                and np.isfinite([r[k] for k in ('CL','CD','CM','alpha_deg')]).all()
                and r['selected_solver_ncrit_by_element']==j['ncrit_by_element']
                and abs(sum(v['CD_wake'] for v in r['element_forces'].values())-r['CD'])<1e-12)
            ok=ok and all(t['inside_cell'] and abs(t['amplification_residual'])<=1e-7 and t['station_residual_max']<=1e-7 for t in r['transition_audit'])
        elif any(r[k] is not None for k in ('CL','CD','CM','alpha_deg')):
            raise ValueError('Failed target contains finite forces')
        meta=read(out/j['name']/'flow/metadata.json')
        if meta['case']['study']!=spec or meta['case']['job']!=j:raise ValueError('Exported scenario changed')
        for d in meta['diagnostics']:
            if d['ncrit_by_element']!=j['ncrit_by_element'] or d['ncrit']!=j['ncrit_by_element']['fore']:
                raise ValueError('Incorrect per-element Ncrit provenance')
        rows.append(dict(job=j['name'],feasible=r['feasible'],execution_check_pass=bool(ok),
            CL=r['CL'],CD=r['CD'],alpha_deg=r['alpha_deg'],selected_path=r.get('selected_path'),
            ncrit_by_element=j['ncrit_by_element'],wall_seconds=r['wall_seconds'],error=r['error']))
    controls=check_controls(out,spec)
    result=dict(required_jobs=len(rows),converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=controls['passed'] and all(r['execution_check_pass'] for r in rows),
        controls=controls,jobs=rows,
        qualification='Execution and fresh controls only; full selected-state residual audit and pressure/BL factorial interpretation pending.')
    write(out/'summary.json',result);return result


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements = [(name, repanel(xy, count)) for (name, xy), count in zip(geometry(job['flap']), spec['panels'])]
        def factory(alpha):
            return make_solver(elements, job, spec, alpha)
        solutions = []
        result = dict(flap=job['flap'], ncrit=job['ncrit_by_element']['fore'],
            ncrit_by_element=job['ncrit_by_element'], target_cl=job['target_cl'],
            wake_length=job['wake_length'], reused=False, feasible=False,
            CL=None, CD=None, CM=None, alpha_deg=None, error='', records=[], attempts=[],
            flow_directory=str((folder/'flow').relative_to(ROOT)))
        def checkpoint(sol, context):
            sol.diagnostics['ncrit_by_element'] = dict(job['ncrit_by_element'])
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
            actual = {m.geom.name: m.param.ncrit for m in solver.models}
            if actual != job['ncrit_by_element']:
                raise ValueError('Selected solver thresholds differ from declared scenario')
            result['selected_solver_ncrit_by_element'] = actual
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
            print('Prepared 10 transition-isolation targets; study sha256='+sha(out/'study.json'));return 0
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
            # Controls gate the mixed cases, including on resume.
            for phase in ('controls','mixed'):
                if phase=='mixed':
                    checks=check_controls(out,spec)
                    write(out/'controls.json',checks)
                    if not checks['passed']:
                        status('execution_failed',error='Fresh controls failed; mixed states were not launched')
                        return 2
                pending=[]
                for j in spec['jobs']:
                    if bool(j['control']) != (phase=='controls'):continue
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
