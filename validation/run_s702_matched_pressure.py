#!/usr/bin/env python3
"""Eight restartable pressure-matched targets, using the frozen existing worker."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import fcntl
import importlib.metadata
import os
from pathlib import Path
import platform
import shutil
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from run_design_initialization import ROOT, VERSIONS, prepare, read, sha, stamp
from s702_ncrit_sweep import write, tag
from run_s702_wake_extent import worker as frozen_worker, CODE
import numpy as np

STUDY='validation/results/s702_matched_pressure_20261002'
QUAL='validation/results/s702_pressure_qualification_20261002'
PARENT='validation/results/s702_wake_extent_20260922/study.json'
SEEDS='validation/results/s702_drag_mechanism/selected_states.csv'


def csv_rows(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))


def seed(rows, flap, ncrit, target):
    values=[r for r in rows if float(r['reynolds'])==1e6 and float(r['mach'])==.1
            and int(r['flap'])==flap and float(r['ncrit'])==ncrit]
    near=sorted(values,key=lambda r:abs(float(r['CL'])-target))[:2]
    if len(near)!=2:raise ValueError('Missing archived numerical lift seeds')
    a,b=sorted(near,key=lambda r:float(r['CL']))
    cla,clb=float(a['CL']),float(b['CL']);aa,ab=float(a['alpha_deg']),float(b['alpha_deg'])
    slope=(clb-cla)/(ab-aa)
    if not (.02<=slope<=.3 and cla-.1<=target<=clb+.1):
        raise ValueError('Numerical seed exceeds declared slope/extrapolation limits')
    return aa+(target-cla)/slope,[a['state_id'],b['state_id']]


def specification():
    pins={PARENT:'d79a1add8303d9f5b5de97d18396f7360a5dada3e44f7e88593ab3693c4d6f96',
          SEEDS:'bbad31dd9f5657a234e702506c1fdd83737b1bc87e86d579f4374ee2d27da0bd',
          'validation/run_s702_wake_extent.py':'1b0746093013b6c451ccb9e9536065152105f31ae4efe2f5c8ea7d53cdf05d1a',
          'validation/run_design_initialization.py':'952f88d657e3aad8b66c9352b1f569c6dfd2eb768767621ca3428c0aa3eec839'}
    for name,digest in pins.items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen reference changed: '+name)
    previous=read(ROOT/PARENT)
    for name in CODE:
        if name in previous['sha256'] and sha(ROOT/name)!=previous['sha256'][name]:
            raise ValueError('Verified numerical code changed: '+name)
    source=read(ROOT/QUAL/'source.json')
    files=list(pins)+CODE+[QUAL+'/source.json','validation/run_s702_matched_pressure.py',
        'validation/S702_MATCHED_PRESSURE_PROTOCOL.md','test_s702_matched_pressure.py',
        'test_fixed_lift_policy.py','test_optimization.py','test_paneling.py',
        'validation/extract_s702.py','validation/extract_s702_flap_pressure.py']
    for prefix,key in (('', 'input_sha256'),(QUAL+'/', 'files_sha256')):
        for name,digest in source[key].items():
            path=prefix+name
            if sha(ROOT/path)!=digest:raise ValueError('Qualified pressure source changed: '+path)
            files.append(path)
    for name in ('fore.dat','aft_d5.dat','aft_d10.dat','source.json'):
        name='validation/data/s702/'+name
        if sha(ROOT/name)!=previous['sha256'][name]:raise ValueError('Geometry source changed: '+name)
        files.append(name)
    rows=csv_rows(ROOT/SEEDS);jobs=[]
    for c in source['cases']:
        for ncrit in (9.,6.25):
            center,ids=seed(rows,c['flap'],ncrit,c['target_cl'])
            jobs.append(dict(name=c['name']+'_n'+tag(ncrit),flap=c['flap'],ncrit=ncrit,
                target_cl=c['target_cl'],center=center,seed_state_ids=ids,wake_length=1.,
                pressure_case=c['name'],experiment_alpha=c['alpha']))
    return dict(schema=1,experiment='S702 pressure at four experimental snapshot lifts',
        reynolds=1e6,mach=.1,panels=[340,240],reference_chord=1.,jobs=jobs,
        iterations=100,cl_tolerance=2e-5,residual_tolerance=1e-7,initial_half_width=.05,
        alpha_bounds=[-6.,8.],max_evaluations=20,max_brackets=2,warm_start_alpha=-1.,
        entry_offset=-.4,workers_limit=4,expected_versions=VERSIONS,
        qualification='Calibrated-airfoil physical observable check; no new calibration, convergence or independent-geometry validation claim.',
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def worker(out,job,spec):
    # The unchanged archived worker checkpoints requested/failed states and the
    # complete selected Newton state. Only the declared common Ncrit varies.
    return frozen_worker(out,job,dict(spec,ncrit=job['ncrit']))


def summarize(out,spec):
    rows=[]
    for j in spec['jobs']:
        r=read(out/j['name']/'result.json');ok=r['feasible']
        if ok:
            selected=r['records'][r['selected_flow_index']]
            ok=bool(selected['admissible'] and selected['residual']<=spec['residual_tolerance']
                and abs(r['CL']-j['target_cl'])<=spec['cl_tolerance']
                and np.isfinite([r[k] for k in ('CL','CD','CM','alpha_deg')]).all())
        else:
            if any(r[k] is not None for k in ('CL','CD','CM','alpha_deg')):
                raise ValueError('Failed target contains finite forces')
        rows.append(dict(job=j['name'],feasible=r['feasible'],execution_check_pass=ok,
            CL=r['CL'],CD=r['CD'],alpha_deg=r['alpha_deg'],selected_path=r.get('selected_path'),
            wall_seconds=r['wall_seconds'],error=r['error']))
    result=dict(required_jobs=len(rows),converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=all(r['execution_check_pass'] for r in rows),jobs=rows,
        qualification='Provisional completion checks only; full selected-state audit and pressure interpretation pending.')
    write(out/'summary.json',result);return result


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
            print('Prepared 8 pressure targets; study sha256='+sha(out/'study.json'));return 0
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
            pending=[]
            for j in spec['jobs']:
                folder=out/j['name']
                if (folder/'result.json').exists():continue
                if folder.exists():
                    dest=out/'interrupted_attempts'/stamp().replace(':','')/j['name']
                    dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(folder,dest)
                pending.append(j)
            with ProcessPoolExecutor(args.workers) as pool:
                futures={pool.submit(worker,out,j,spec):j for j in pending}
                for future in as_completed(futures):
                    r=future.result();done+=1
                    status('running',last_job=futures[future]['name'],last_job_feasible=r['feasible'])
                    print(f"{done}/{len(spec['jobs'])}: {futures[future]['name']} feasible={r['feasible']}",flush=True)
            summary=summarize(out,spec);status('finished',execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed',error=f'{type(exc).__name__}: {exc}');raise


if __name__=='__main__':raise SystemExit(main())
