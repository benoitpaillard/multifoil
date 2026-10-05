#!/usr/bin/env python3
"""Run/resume the separated S702 discretization study and monitor jony4."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from s702_ncrit_sweep import write, tag
import run_s702_panel_convergence as previous
import run_s702_wake_extent as flow

STUDY = 'validation/results/s702_discretization_20261001'
REMOTE = '/home/jony/multifoil-runs/s702-discretization-20261001'
PARENT = 'validation/results/s702_panel_convergence_20260922'
AUDIT = 'validation/results/s702_panel_convergence_analysis'
STAGES = ('fine', 'contour', 'surface', 'wake_count', 'wake_spacing')

from dataclasses import asdict
import contextlib
import numpy as np
from multifoil import save_results, repanel
from viscous import CoupledSolver, new_model
from diagnostic_discretization import DiscretizedSolver
from s702_ncrit_sweep import geometry, record, transition_audit
from fixed_lift_policy import solve_with_initialization


def specification():
    prior = previous.specification()
    if sha(ROOT/PARENT/'study.json') != 'b522cead604443e3470eb7cb38df9efdd1863ea15a6bf8cc6440cbf6eeb43833' or read(ROOT/PARENT/'study.json') != prior:
        raise ValueError('Frozen parent sources/settings changed')
    audit = read(ROOT/AUDIT/'summary.json')
    if not audit['full_state_audit_pass'] or not audit['control_saved_arrays_exactly_equal']:
        raise ValueError('Parent state audit/control failed')
    pinned = read(ROOT/AUDIT/'input_manifest.json')['files']
    files = list(prior['sha256']) + [PARENT+'/study.json', PARENT+'/runtime.json',
        AUDIT+'/summary.json', AUDIT+'/input_manifest.json', AUDIT+'/FINDINGS.md',
        'diagnostic_discretization.py', 'validation/run_s702_discretization.py',
        'validation/S702_DISCRETIZATION_PROTOCOL.md', 'test_diagnostic_discretization.py',
        'test_s702_discretization.py']
    jobs = []
    for flap in (5, 10):
        old = next(j for j in prior['jobs'] if j['flap'] == flap and j['target_cl'] == .8 and not j['control'])
        fine = old['reference_result']
        finer = str(Path(PARENT)/old['name']/'result.json')
        for ref in (fine, finer):
            folder = Path(ref).parent
            for name, digest in pinned.items():
                if Path(name).is_relative_to(folder):
                    if sha(ROOT/name) != digest:
                        raise ValueError('Audited reference changed: '+name)
                    files.append(name)
        for stage in STAGES:
            jobs.append(dict(name=f'd{flap}_{stage}', flap=flap, target_cl=.8, center=old['center'],
                stage=stage, control=False, source_panels=[340,240] if stage == 'fine' else [440,320],
                panels=[340,240] if stage in ('fine','contour') else [440,320], wake_length=1.,
                wake_count_recipe='finer' if stage in ('wake_count','wake_spacing') else 'fine',
                wake_spacing_recipe='finer' if stage == 'wake_spacing' else 'fine',
                reference_result=fine, finer_reference_result=finer))
    control = dict(next(j for j in jobs if j['name']=='d10_fine'))
    control.update(name='d10_legacy_control', stage='legacy_control', control=True)
    jobs.append(control)
    inherited = {k:v for k,v in prior.items() if k not in
        ('last_automatic_grid_extension','per_job_panel_override','grid_sequence','archived_wake_length')}
    return dict(inherited, experiment='S702 contour/surface/wake factor separation', jobs=jobs,
        parent_study=PARENT, independent_discretization_controls=True,
        stages=list(STAGES), panels=[340,240], workers_limit=4,
        qualification='Conditional numerical factor contributions; not experimental validation or a new grid extension.',
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def numerical_setup(job, spec):
    raw = geometry(job['flap'])
    recipes = {}
    for label, counts in (('fine',[340,240]), ('finer',[440,320])):
        models = [new_model(name, repanel(xy,count), spec['reynolds'], spec['mach'],
            spec['ncrit'], np.array([.25,0.]), 1e-5) for (name,xy),count in zip(raw,counts)]
        tips = [.5*(m.foil.x[:,0]+m.foil.x[:,-1]) for m in models]
        recipes[label] = dict(nodes={}, first_spacing={})
        for model, tip in zip(models,tips):
            recipes[label]['nodes'][model.geom.name] = max(16,int(np.ceil(model.foil.N/10+10+20*np.log(
                max(1.,(1.+max(np.linalg.norm(tip-other) for other in tips))/model.local_chord)))))
            recipes[label]['first_spacing'][model.geom.name] = float(.5*((model.foil.s[1]-model.foil.s[0])+
                (model.foil.s[-1]-model.foil.s[-2])))
    elements = [(name,repanel(xy,count)) for (name,xy),count in zip(raw,job['source_panels'])]
    settings = dict(surface_nodes={name:count+1 for (name,_),count in zip(raw,job['panels'])},
        wake_nodes=recipes[job['wake_count_recipe']]['nodes'],
        wake_first_spacing=recipes[job['wake_spacing_recipe']]['first_spacing'])
    return elements, settings, recipes


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a',buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter()
        elements, settings, recipes = numerical_setup(job,spec)
        for name,xy in elements:
            np.save(folder/('source_'+name+'.npy'),xy)
        write(folder/'numerical_setup.json',dict(requested=settings,recipes=recipes,
            applied={} if job['control'] else settings, production_solver=job['control'],
            source_sha256={name:sha(folder/('source_'+name+'.npy')) for name,_ in elements}))
        def factory(alpha):
            kwargs = dict(reynolds=spec['reynolds'],mach=spec['mach'],ncrit=spec['ncrit'],
                iterations=spec['iterations'],tolerance=spec['residual_tolerance'],
                wake_length=job['wake_length'],start_alpha=alpha)
            return CoupledSolver(elements,**kwargs) if job['control'] else DiscretizedSolver(elements,**settings,**kwargs)
        solutions = []
        effective = dict(spec,panels=job['panels'],source_panels=job['source_panels'])
        result = dict(flap=job['flap'],ncrit=spec['ncrit'],target_cl=job['target_cl'],
            wake_length=job['wake_length'],stage=job['stage'],reused=False,feasible=False,
            CL=None,CD=None,CM=None,alpha_deg=None,error='',records=[],attempts=[],
            flow_directory=str((folder/'flow').relative_to(ROOT)))
        def checkpoint(sol,context):
            sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'],aft=spec['ncrit'])
            solutions.append(sol)
            result['records'].append(dict(record(sol,context['phase']),initialization=context))
            save_results(folder/'flow',solutions,elements,dict(study=effective,job=job,target_cl=job['target_cl']),'viscous')
            write(folder/'progress.json',result)
        outcome = solve_with_initialization(factory,job['target_cl'],job['center'],spec['alpha_bounds'],
            flow.paths_for(job,spec),half_width=spec['initial_half_width'],cl_tolerance=spec['cl_tolerance'],
            residual_tolerance=spec['residual_tolerance'],max_evaluations=spec['max_evaluations'],
            max_brackets=spec['max_brackets'],on_evaluation=checkpoint)
        result.update(attempts=outcome.attempts,error=outcome.error,
            initialization_paths=[asdict(p) for p in flow.paths_for(job,spec)])
        if outcome.feasible:
            sol,solver = outcome.solution,outcome.solver
            if solver.last_alpha != sol.alpha or solver.last_solution.cl != sol.cl:
                raise ValueError('Returned solver differs from selected state')
            result.update(feasible=True,CL=sol.cl,CD=sol.cd,CM=sol.cm,alpha_deg=sol.alpha,
                selected_flow_index=outcome.selected_flow_index,selected_path=outcome.selected_path,
                boundary_layers=sol.diagnostics['boundary_layers'],element_forces=sol.diagnostics['element_forces'],
                transition_audit=transition_audit(solver),wake_discretization={})
            stream = np.array([np.cos(np.deg2rad(sol.alpha)),np.sin(np.deg2rad(sol.alpha))])
            for model in solver.models:
                np.savez_compressed(folder/(model.geom.name+'_state.npz'),U=model.glob.U,xy=model.foil.x,
                    xi=model.isol.xi,turbulent=model.vsol.turb,wake_xy=model.wake.x,
                    lower=np.asarray(model.vsol.Is[0]),upper=np.asarray(model.vsol.Is[1]),
                    wake_s=model.wake.s,wake_t=model.wake.t,foil_s=model.foil.s,foil_t=model.foil.t,
                    source_xy=model.geom.xpoint)
                spacing = np.linalg.norm(np.diff(model.wake.x,axis=1),axis=0)
                result['wake_discretization'][model.geom.name] = dict(nodes=int(model.wake.N),
                    arclength=float(sum(spacing)),first_spacing=float(spacing[0]),last_spacing=float(spacing[-1]),
                    streamwise_start=float(stream@model.wake.x[:,0]),streamwise_end=float(stream@model.wake.x[:,-1]))
        result['wall_seconds'] = time.perf_counter()-begun
        write(folder/'result.json',result)
        return result


def summarize(out,spec):
    summary = flow.summarize(out,spec)
    results = {j['name']:read(out/j['name']/'result.json') for j in spec['jobs']}
    endpoints = []
    for job in spec['jobs']:
        if job['stage'] not in ('fine','wake_spacing'): continue
        r = results[job['name']]
        ref = read(ROOT/(job['reference_result'] if job['stage']=='fine' else job['finer_reference_result']))
        delta = (r['CD']-ref['CD'])*1e4 if r['feasible'] else None
        endpoints.append(dict(job=job['name'],difference_CD_counts=delta,
            passed=delta is not None and abs(delta)<=spec['acceptance_control_CD_counts']))
    pairs = []
    for stage in STAGES:
        a,b = results['d5_'+stage],results['d10_'+stage]
        pairs.append(dict(stage=stage,complete=a['feasible'] and b['feasible'],
            flap_delta_CD_counts=(b['CD']-a['CD'])*1e4 if a['feasible'] and b['feasible'] else None))
    exact = results['d10_fine']['feasible'] and results['d10_legacy_control']['feasible']
    if exact:
        for name in ('fore','aft'):
            with np.load(out/'d10_fine'/(name+'_state.npz')) as a, np.load(out/'d10_legacy_control'/(name+'_state.npz')) as b:
                exact = exact and set(a.files)==set(b.files) and all(np.array_equal(a[k],b[k]) for k in a.files)
    for row in summary['jobs']:
        row['difference_from_fine_CD_counts'] = row.pop('difference_from_1c_CD_counts')
    increments = []
    for a,b in zip(pairs[:-1],pairs[1:]):
        complete = a['complete'] and b['complete']
        increments.append(dict(from_stage=a['stage'],to_stage=b['stage'],complete=complete,
            paired_change_CD_counts=b['flap_delta_CD_counts']-a['flap_delta_CD_counts'] if complete else None,
            individual_change_CD_counts={str(f):(results[f'd{f}_'+b['stage']]['CD']-
                results[f'd{f}_'+a['stage']]['CD'])*1e4 for f in (5,10)} if complete else None))
    summary.update(endpoint_checks=endpoints,paired_stages=pairs,control_saved_arrays_exactly_equal=bool(exact),
        conditional_increments=increments,
        execution_checks_pass=bool(summary['execution_checks_pass'] and exact and all(e['passed'] for e in endpoints)),
        qualification='Provisional numerical comparison; full state audit, factor increments and interpretation remain required.')
    write(out/'summary.json',summary)
    return summary


PROBE = r'''
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); out=root/sys.argv[2]
def read(p): return json.loads(p.read_text()) if p.exists() else {}
status=read(out/'status.json'); spec=read(out/'study.json')
try:
    cmd=Path('/proc/'+str(status.get('pid'))+'/cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive='run_s702_discretization.py' in cmd and sys.argv[2] in cmd
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
                print(f"\nS702 separated discretization / jony4: {state}; {status.get('completed_jobs', 0)}/11 completed")
                for row in data['jobs']:
                    print(f"{row['job']:27s} {row['state']:10s} solves={row['evaluations']:2d} failed={row['failed']:2d} {row['path']}")
                    if row['error']: print('  '+row['error'])
                if data['summary']:
                    s = data['summary']
                    print(f"Feasible: {s['feasible_jobs']}/11; original-solver control: {'PASS' if s['control_pass'] else 'FAIL'}")
                    print(f"Endpoint reproduction: {sum(e['passed'] for e in s['endpoint_checks'])}/4; "
                          f"control arrays exact: {s['control_saved_arrays_exactly_equal']}")
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
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--prepare-only', action='store_true')
    ap.add_argument('--remote-status', action='store_true')
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()
    if args.remote_status: return remote_status(args.watch, args.json)
    if args.watch or args.json: ap.error('--watch and --json require --remote-status')
    out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 4:
        ap.error('Use a results directory and 1–4 workers')
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
