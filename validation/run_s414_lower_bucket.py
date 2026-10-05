#!/usr/bin/env python3
"""Four guarded S414 free-transition cases, with an archived-state control gate."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
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
import numpy as np
from run_s414_reconstruction import worker,ROOT,VERSIONS,prepare,read,sha,stamp,write
from extract_s414_drag import ascending_branch,interpolate

STUDY='validation/results/s414_lower_bucket_20261003'
OLD='validation/results/s414_reconstruction_20261002'
AUDIT='validation/results/s414_reconstruction_analysis_20261003'
SOURCE='validation/results/s414_trip_qualification_20261003'
CONTROL='cl028_control'
REFERENCE=OLD+'/cl028_nominal'
PINS={OLD+'/study.json':'56e91de0f2c19c16c2d732c084e3141cd19c765d38bcd29db3d074ae5685d9f9',
      AUDIT+'/verification.json':'f3e1bd12e648eeacd95f55ea63deb16d7eb88f3f5c70007fba4943c279f65453',
      SOURCE+'/verification.json':'afd1cc00c3a83326760d7097b854a53e38f54a638fb809e1deaa42a08f731aec'}


def csv_rows(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))


def specification():
    files=set(PINS)
    for name,digest in PINS.items():
        if sha(ROOT/name)!=digest:raise ValueError('Audited reference changed: '+name)
    old=read(ROOT/OLD/'study.json');receipt=read(ROOT/AUDIT/'verification.json')
    for name,digest in old['sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen numerical/source input changed: '+name)
        files.add(name)
    for name in ('input_manifest.json','residual_recheck.json','summary.json'):
        path=AUDIT+'/'+name
        if sha(ROOT/path)!=receipt['output_sha256'][name]:raise ValueError('State audit changed: '+path)
        files.add(path)
    archived=read(ROOT/AUDIT/'input_manifest.json')['files']
    for name,digest in archived.items():
        if name.startswith(REFERENCE+'/') or name==OLD+'/runtime.json':
            if sha(ROOT/name)!=digest:raise ValueError('Archived control changed: '+name)
            files.add(name)
    source_receipt=read(ROOT/SOURCE/'verification.json')
    for name in ('free_boundary_pressure.csv','experimental_tables.csv','published_mses_drag.csv',
                 'evidence_qualification.json','summary.json','source_manifest.json',
                 'pressure_tap_checks.csv','pressure_recovery_intervals.csv','repeated_pressure_crosscheck.csv'):
        path=SOURCE+'/'+name
        if sha(ROOT/path)!=source_receipt['output_sha256'][name]:raise ValueError('Qualified measurement changed: '+path)
        files.add(path)
    source=read(ROOT/SOURCE/'source_manifest.json')
    for name,digest in source['input_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Measurement extraction input changed: '+name)
        files.add(name)
    prior=next(j for j in old['jobs'] if j['name']=='cl028_nominal')
    control=dict(prior,name=CONTROL,reference_result=REFERENCE+'/result.json',
        pressure_figure='7(d)',center_method='Exact archived nominal job center',
        mses_drag_available=True,mses_pressure_available=True)
    jobs=[control]
    lift=[{k:float(v) if k in ('CL','alpha_deg') else v for k,v in r.items()}
          for r in csv_rows(ROOT/'validation/results/s414_drag_reference_20261002/mses_lift.csv')]
    branch=ascending_branch(lift)
    pressure=csv_rows(ROOT/SOURCE/'free_boundary_pressure.csv')
    for case,figure,cl,alpha in [('clm0072','7(a)',-.072,-4.09),('cl0056','7(b)',.056,-3.07),('cl0175','7(c)',.175,-2.06)]:
        rows=[r for r in pressure if r['figure']==figure]
        if len(rows)!=94 or any(float(r['CL'])!=cl or float(r['alpha_deg'])!=alpha for r in rows):
            raise ValueError('Pressure case association changed')
        available=branch[0]['CL']<=cl<=branch[-1]['CL']
        center=interpolate(branch,cl,y='alpha_deg') if available else alpha
        jobs.append(dict(name=case+'_nominal',pressure_case=case,pressure_figure=figure,
            scenario='nominal',target_cl=cl,experiment_alpha=alpha,center=center,wake_length=1.,control=False,
            center_method='Published MSES increasing-lift interpolation' if available else 'Measured incidence seed only; no MSES extrapolation',
            mses_drag_available=cl>=0,mses_pressure_available=False))
    spec={k:v for k,v in old.items() if k not in ('jobs','sha256','experiment','qualification','workers_limit')}
    spec.update(experiment='S414 free-transition lower-bucket pressure/transition sequence',jobs=jobs,
        workers_limit=3,reference_control=REFERENCE,archived_study_sha256=PINS[OLD+'/study.json'],
        qualification='Bounded nominal reconstruction; no exact geometry admission, numerical-convergence claim or closure fit.',
        gate='Exact archived CL 0.284 state must repeat before all three new targets')
    files.update(['validation/run_s414_lower_bucket.py','validation/S414_LOWER_BUCKET_PROTOCOL.md',
                  'test_s414_lower_bucket.py'])
    spec['sha256']={name:sha(ROOT/name) for name in sorted(files)}
    return spec


def compare_control(folder,reference,tolerance_counts):
    actual,expected=read(folder/'result.json'),read(reference/'result.json')
    passed=bool(actual['feasible'] and expected['feasible']);delta=None;arrays=[]
    if passed:
        delta=(actual['CD']-expected['CD'])*1e4
        passed=abs(delta)<=tolerance_counts and actual['selected_path']==expected['selected_path']
        for element in ('fore','aft'):
            with np.load(folder/(element+'_state.npz')) as a,np.load(reference/(element+'_state.npz')) as b:
                if set(a.files)!=set(b.files):raise ValueError('Archived state fields changed')
                for field in a.files:
                    same=np.array_equal(a[field],b[field]);passed=passed and same
                    arrays.append(dict(element=element,field=field,exactly_equal=same))
    return dict(passed=bool(passed),control=CONTROL,reference=str(reference.relative_to(ROOT)) if reference.is_relative_to(ROOT) else str(reference),
                CD_difference_counts=delta,array_comparisons=arrays)


def check_control(out,spec):
    result=compare_control(out/CONTROL,ROOT/spec['reference_control'],spec['acceptance_control_CD_counts'])
    write(out/'controls.json',result)
    if not result['passed']:raise ValueError('Archived control mismatch; new pressure cases are not launched')
    return result


def summarize(out,spec):
    rows=[]
    for job in spec['jobs']:
        result=read(out/job['name']/'result.json');valid=bool(result['feasible'])
        if valid:
            selected=result['records'][result['selected_flow_index']]
            valid=bool(selected['admissible'] and selected['residual'] is not None and selected['residual']<=spec['residual_tolerance']
                and np.isfinite([result[k] for k in ('CL','CD','CM','alpha_deg')]).all()
                and abs(result['CL']-job['target_cl'])<=spec['cl_tolerance']
                and abs(sum(f['CD_wake'] for f in result['element_forces'].values())-result['CD'])<1e-12
                and all(t['inside_cell'] and abs(t['amplification_residual'])<=1e-7 and t['station_residual_max']<=1e-7
                        for t in result['transition_audit']))
        elif any(result[k] is not None for k in ('CL','CD','CM','alpha_deg')):
            raise ValueError('Failed target contains selected finite forces')
        metadata=out/job['name']/'flow/metadata.json'
        if metadata.exists():
            meta=read(metadata)
            if meta['case']['study']!=spec or meta['case']['job']!=job:raise ValueError('Exported case differs')
            if any(d['ncrit_by_element']!={'fore':9.,'aft':9.} for d in meta['diagnostics']):raise ValueError('Ncrit changed')
        elif result['feasible']:raise ValueError('Missing flow exports')
        rows.append(dict(job=job['name'],feasible=result['feasible'],execution_check_pass=bool(valid),
            CL=result['CL'],CD=result['CD'],alpha_deg=result['alpha_deg'],selected_path=result.get('selected_path'),
            wall_seconds=result['wall_seconds'],error=result['error']))
    controls=check_control(out,spec)
    summary=dict(required_jobs=4,converged_jobs=sum(r['feasible'] for r in rows),
        execution_checks_pass=controls['passed'] and all(r['execution_check_pass'] for r in rows),
        controls=controls,jobs=rows,qualification='Execution gates only; full saved-state audit and scientific comparison pending.')
    write(out/'summary.json',summary);return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/STUDY);ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args();out=args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1<=args.workers<=3:ap.error('Use a results directory and 1–3 workers')
    out.mkdir(parents=True,exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);spec=specification();prepare(out,spec)
        if args.prepare_only:print('Prepared four S414 lower-bucket jobs; study sha256='+sha(out/'study.json'));return 0
        runtime=dict(python=platform.python_version(),hostname=platform.node(),executable=sys.executable,
            versions={k:importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions']!=VERSIONS or runtime['python']!='3.12.13':raise ValueError('Use pinned Python 3.12.13 environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json')!=runtime:raise ValueError('Runtime changed; use a new directory')
        write(out/'runtime.json',runtime)
        done=sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state,**extra):write(out/'status.json',dict(state=state,updated_utc=stamp(),pid=os.getpid(),completed_jobs=done,total_jobs=4,**extra))
        def pending(job):
            folder=out/job['name']
            if (folder/'result.json').exists():return False
            if folder.exists():
                dest=out/'interrupted_attempts'/stamp().replace(':','')/job['name'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(folder,dest)
            return True
        status('running',phase='archived_control')
        try:
            control=next(j for j in spec['jobs'] if j['control'])
            if pending(control):
                result=worker(out,control,spec);done+=1
                status('running',phase='archived_control',last_job=control['name'],last_job_feasible=result['feasible'])
                print(f'{done}/4: {control["name"]} feasible={result["feasible"]}',flush=True)
            check_control(out,spec)
            jobs=[j for j in spec['jobs'] if not j['control'] and pending(j)]
            status('running',phase='lower_bucket',archived_control_pass=True)
            with ProcessPoolExecutor(args.workers) as pool:
                futures={pool.submit(worker,out,j,spec):j for j in jobs}
                for future in as_completed(futures):
                    result=future.result();done+=1
                    status('running',phase='lower_bucket',last_job=futures[future]['name'],last_job_feasible=result['feasible'])
                    print(f'{done}/4: {futures[future]["name"]} feasible={result["feasible"]}',flush=True)
            summary=summarize(out,spec);status('finished',execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed',error=f'{type(exc).__name__}: {exc}');raise


if __name__=='__main__':raise SystemExit(main())
