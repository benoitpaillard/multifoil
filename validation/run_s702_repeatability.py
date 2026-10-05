#!/usr/bin/env python3
"""Run/resume the predeclared S702 initialization matrix with unchanged workers."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import traceback

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

from s702_ncrit_sweep import ROOT, read, write, sha, tag, local_seed, run_group
from report_s702_ncrit_sweep import csv_write, flow_directory

RESULTS = ROOT/'validation/results'
EXPECTED_VERSIONS = {'numpy': '2.3.5', 'scipy': '1.17.0', 'matplotlib': '3.10.8'}
CODE = ['multifoil.py', 'viscous.py', 'optimization.py', 'paneling.py', 'vendor/mfoil.py',
        'validation/s702_ncrit_sweep.py', 'validation/s702_geometry.py',
        'validation/s702_benchmark.py', 'validation/s702_transition_study.py',
        'validation/report_s702_ncrit_sweep.py', 'validation/s702_sweep_metrics.py',
        'validation/run_s702_repeatability.py']


def stamp():
    return datetime.now(timezone.utc).isoformat()


def specification():
    import csv
    analysis = RESULTS/'s702_ncrit_analysis'
    frozen = read(analysis/'frozen_selection.json')
    if frozen['selected_ncrit'] != 6.25 or sha(analysis/'points.csv') != frozen['calibration_points_sha256']:
        raise ValueError('Frozen S702 calibration changed')
    archive = {f['path']: f['sha256'] for f in read(ROOT/'handover/archive_manifest.json')['files']}
    inputs = set(CODE + ['validation/S702_REPEATABILITY_PROTOCOL.md',
        'validation/data/s702/source.json', 'validation/data/s702/fore.dat',
        'validation/data/s702/aft_d5.dat', 'validation/data/s702/aft_d10.dat',
        'validation/results/s702_ncrit_analysis/holdout_after_retries.csv',
        'validation/results/s702_ncrit_analysis/frozen_selection.json'])
    references, centers = [], {}
    with (analysis/'holdout_after_retries.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    for flap in (5, 10):
        for cl in (.5, .8, 1.2):
            row = next(r for r in rows if float(r['ncrit']) == 6.25 and int(r['flap']) == flap and float(r['target_cl']) == cl)
            path = ROOT/row['result_path']
            result = read(path)
            if not result['feasible'] or abs(result['CL']-cl) > 2e-5:
                raise ValueError('Archived selected reference is not feasible')
            prior_spec = read(path.parents[2]/'study.json')
            if (prior_spec['reynolds'], prior_spec['mach'], prior_spec['panels']) != (700000., .07, [260,180]):
                raise ValueError('Reference physical protocol changed')
            for name in ('multifoil.py','viscous.py','optimization.py','paneling.py','vendor/mfoil.py'):
                if sha(ROOT/name) != prior_spec['sha256'][name]:
                    raise ValueError('Numerical core differs from baseline: '+name)
            flow = flow_directory(result)
            inputs.update([row['result_path'], str((path.parents[2]/'study.json').relative_to(ROOT)),
                           str((flow/'metadata.json').relative_to(ROOT)), str((flow/'polar.csv').relative_to(ROOT))])
            index = result['selected_flow_index']
            for element in ('fore','aft'):
                inputs.add(str((flow/f'geometry_{element}.dat').relative_to(ROOT)))
                for kind in ('surface','wake'):
                    inputs.add(str((flow/f'{kind}_{index:03d}_{element}.csv').relative_to(ROOT)))
            center = result['alpha_deg']
            if row['local_retry_of']:
                inputs.add(row['local_retry_of'])
                original = read(ROOT/row['local_retry_of'])
                center = local_seed(original)
            centers[(flap, cl)] = center
            references.append(dict(flap=flap, target_cl=cl, result_path=row['result_path'],
                original_failure=row['local_retry_of'], CL=result['CL'], CD=result['CD'], CM=result['CM'],
                alpha_deg=result['alpha_deg'], boundary_layers=result['boundary_layers'],
                element_forces=result['element_forces'], local_seed_alpha=center))
    jobs = []
    for flap in (5,10):
        for cl in (.5,.8,1.2):
            jobs.append(dict(name=f'cold_d{flap}_cl{tag(cl)}', path='cold', flap=flap, targets=[cl]))
        jobs.append(dict(name=f'increasing_d{flap}', path='increasing', flap=flap, targets=[.5,.8,1.2]))
        jobs.append(dict(name=f'decreasing_d{flap}', path='decreasing', flap=flap, targets=[1.2,.8,.5]))
    failed = read(RESULTS/'s702_ncrit_holdout_selected/n6p25_d+10/cl0p8/result.json')
    warm_prefix = [r['alpha_deg'] for r in failed['records'][:2] if r['admissible']]
    if len(warm_prefix) != 2:
        raise ValueError('Expected the two original admissible warm replay states')
    jobs.append(dict(name='warm_replay_d10_cl0p8', path='warm_replay', flap=10, targets=[.8], warm_seed_alphas=warm_prefix))
    for job in jobs:
        job.update(ncrit=6.25, seed_alphas={str(cl): centers[(job['flap'],cl)] for cl in job['targets']})
    for name in inputs:
        if name in archive and sha(ROOT/name) != archive[name]:
            raise ValueError('Archived execution dependency changed: '+name)
    return dict(schema=1, experiment='S702 existing-policy initialization repeatability',
        jobs=jobs, references=references, ncrit_values=[6.25], reynolds=700000., mach=.07,
        panels=[260,180], reference_chord=1., wake_length=1., cl_tolerance=2e-5,
        residual_tolerance=1e-7, iterations=100, reuse_legacy=False,
        initial_half_width=.05, start_at_seed=True,
        protocol='Existing three-attempt worker; local cold / increasing target lift / decreasing target lift / archived warm replay',
        selection_protocol='First successful numerical attempt; no measured drag enters path selection',
        restart_policy='Restart incomplete whole jobs from their original initialization; archive all partial outputs',
        expected_versions=EXPECTED_VERSIONS, python_major_minor='3.12',
        nominal_fixed_lift_targets=sum(len(j['targets']) for j in jobs),
        sha256={name:sha(ROOT/name) for name in sorted(inputs)})


def prepare(out, spec):
    out.mkdir(parents=True, exist_ok=True)
    manifest = out/'study.json'
    if manifest.exists():
        if read(manifest) != spec:
            raise ValueError('Study or sources changed; choose a fresh output directory')
    else:
        if any(p.name != 'run.lock' for p in out.iterdir()):
            raise ValueError('Nonempty directory without study manifest')
        write(manifest, spec)
    for name, expected in spec['sha256'].items():
        if name.endswith(('.py','.md')):
            target = out/'source_snapshot'/name
            target.parent.mkdir(parents=True,exist_ok=True)
            if not target.exists():
                shutil.copyfile(ROOT/name,target)
            if sha(target) != expected:
                raise ValueError('Source snapshot changed: '+name)


def worker(out, job, spec):
    group_spec = dict(spec, targets=job['targets'], initialization_path=job['path'])
    if job['path'] == 'warm_replay':
        group_spec.update(start_at_seed=False, warm_seed_alphas=job['warm_seed_alphas'])
    # The archived worker returns dict(**job, targets=actual_results).
    worker_job = {k:v for k,v in job.items() if k != 'targets'}
    return run_group(out, worker_job, group_spec)


def compare(out, spec):
    references = {(r['flap'],r['target_cl']):r for r in spec['references']}
    rows, layers = [], []
    for job in spec['jobs']:
        group = read(out/job['name']/'result.json')
        for result in group['targets']:
            ref = references[job['flap'], result['target_cl']]
            ok = result['feasible']
            row = dict(job=job['name'], path=job['path'], flap=job['flap'], target_cl=result['target_cl'],
                feasible=ok, first_attempt_feasible=result['attempts'][0]['feasible'],
                used_retry=len(result['attempts']) > 1, attempts=len(result['attempts']),
                requested_evaluations=len(result['records']),
                failed_requested_evaluations=sum(not r['admissible'] for r in result['records']),
                wall_seconds=result['wall_seconds'], error=result['error'],
                result_path=str((out/job['name']/f"cl{tag(result['target_cl'])}"/'result.json').relative_to(ROOT)),
                reference_result=ref['result_path'])
            for field in ('CL','CD','CM','alpha_deg'):
                row[field] = result[field] if ok else None
                row['reference_'+field] = ref[field]
                row['difference_'+field] = result[field]-ref[field] if ok else None
            row['difference_CD_counts'] = row['difference_CD']*1e4 if ok else None
            row['residual'] = result['records'][result['selected_flow_index']]['residual'] if ok else None
            if ok:
                if row['residual'] > spec['residual_tolerance'] or abs(result['CL']-result['target_cl']) > spec['cl_tolerance']:
                    raise ValueError('Selected target violates tolerance')
                if abs(sum(v['CD_wake'] for v in result['element_forces'].values())-result['CD']) > 1e-12:
                    raise ValueError('Element wake terms differ from total CD')
                for element in ('fore','aft'):
                    for side in ('upper','lower'):
                        now, old = result['boundary_layers'][element][side], ref['boundary_layers'][element][side]
                        layer = dict(job=job['name'],flap=job['flap'],target_cl=result['target_cl'],element=element,side=side)
                        for prefix, bl in [('',now),('reference_',old)]:
                            layer[prefix+'transition_x_local'] = bl['transition']['x_over_local_chord'] if bl['transition'] else None
                            for key in ('separated_arclength_fraction','min_cf','max_Hk'):
                                layer[prefix+key]=bl[key]
                        layers.append(layer)
            rows.append(row)
    csv_write(out/'comparison.csv', rows)
    csv_write(out/'boundary_layer_comparison.csv', layers)
    spans, flap_deltas = [], []
    for flap in (5,10):
        for cl in (.5,.8,1.2):
            group=[r for r in rows if r['flap']==flap and r['target_cl']==cl and r['path']!='warm_replay']
            complete=all(r['feasible'] for r in group)
            spans.append(dict(flap=flap,target_cl=cl,required=3,converged=sum(r['feasible'] for r in group),
                complete=complete, all_first_attempts=all(r['first_attempt_feasible'] for r in group),
                CD_span_counts=(max(r['CD'] for r in group)-min(r['CD'] for r in group))*1e4 if complete else None))
    for path in ('cold','increasing','decreasing'):
        for cl in (.5,.8,1.2):
            pair=[next(r for r in rows if r['path']==path and r['flap']==f and r['target_cl']==cl) for f in (5,10)]
            complete=all(r['feasible'] for r in pair)
            baseline=(references[10,cl]['CD']-references[5,cl]['CD'])*1e4
            delta=(pair[1]['CD']-pair[0]['CD'])*1e4 if complete else None
            flap_deltas.append(dict(path=path,target_cl=cl,complete=complete,delta_10_minus_5_counts=delta,
                archived_delta_counts=baseline,difference_counts=delta-baseline if complete else None))
    csv_write(out/'path_spans.csv',spans)
    csv_write(out/'paired_flap_deltas.csv',flap_deltas)
    result=dict(completed_utc=stamp(),required_targets=len(rows),converged_targets=sum(r['feasible'] for r in rows),
        first_attempt_successes=sum(r['first_attempt_feasible'] for r in rows),path_spans=spans,
        requested_evaluations=sum(r['requested_evaluations'] for r in rows),
        failed_requested_evaluations=sum(r['failed_requested_evaluations'] for r in rows),
        summed_target_wall_seconds=sum(r['wall_seconds'] for r in rows),
        qualification='Numerical path dependence only; separate first-attempt and fallback success; full state/plot audit still required')
    write(out/'summary.json',result)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args(); out=args.out.resolve()
    if not out.is_relative_to(RESULTS) or not 1 <= args.workers <= 8:
        ap.error('Use an output directory under validation/results and 1–8 workers')
    out.mkdir(parents=True,exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        spec=specification(); prepare(out,spec)
        if args.prepare_only:
            print(json.dumps(dict(status='prepared',jobs=len(spec['jobs']),targets=spec['nominal_fixed_lift_targets'],out=str(out))))
            return 0
        runtime=dict(python=platform.python_version(), hostname=platform.node(),
            versions={n:importlib.metadata.version(n) for n in EXPECTED_VERSIONS},
            thread_limits={n:os.environ[n] for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != EXPECTED_VERSIONS or sys.version_info[:2] != (3,12):
            raise ValueError('Use Python 3.12 with requirements-handover.txt; numerical dependencies must match baseline')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed since launch; use a new study directory')
        write(out/'runtime.json',runtime)
        interruptions=[]
        pending=[]
        for job in spec['jobs']:
            folder=out/job['name']
            if (folder/'result.json').exists():
                continue
            if folder.exists():
                archived=out/'interrupted_attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')/job['name']
                archived.parent.mkdir(parents=True,exist_ok=True)
                shutil.move(folder,archived)
                interruptions.append(str(archived.relative_to(ROOT)))
            pending.append(job)
        invocation=dict(started_utc=stamp(),workers=args.workers,pending_jobs=len(pending),archived_interruptions=interruptions)
        write(out/f"invocation_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json",invocation)
        done=len(spec['jobs'])-len(pending)
        def status(state, **extra):
            write(out/'status.json',dict(state=state,updated_utc=stamp(),completed_jobs=done,total_jobs=len(spec['jobs']),**extra))
        status('running',pid=os.getpid())
        print(json.dumps(dict(event='started',**invocation)),flush=True)
        try:
            with ProcessPoolExecutor(args.workers) as pool:
                futures={pool.submit(worker,out,job,spec):job for job in pending}
                for future in as_completed(futures):
                    job=futures[future]; group=future.result(); done+=1
                    status('running',last_completed_job=job['name'],last_job_success=group['complete_success'])
                    print(json.dumps(dict(event='job_completed',job=job['name'],complete_success=group['complete_success'],wall_seconds=group['wall_seconds'])),flush=True)
            summary=compare(out,spec)
            status('finished',converged_targets=summary['converged_targets'],required_targets=summary['required_targets'],followup='Review comparison.csv, path_spans.csv and selected states before changing the policy')
            print(json.dumps(summary),flush=True)
            return 0 if summary['converged_targets']==summary['required_targets'] else 2
        except Exception:
            status('execution_error',traceback=traceback.format_exc())
            raise


if __name__=='__main__':
    raise SystemExit(main())
