#!/usr/bin/env python3
"""Run/resume five matched S702 reserved-condition panel checks on jony4."""
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
import run_s702_wake_extent as baseline

STUDY = 'validation/results/s702_holdout_panels_20260922'
REMOTE = '/home/jony/multifoil-runs/s702-holdout-panels-20260922'
WAKE_STUDY = 'validation/results/s702_wake_extent_20260922/study.json'
AUDIT = 'validation/results/s702_wake_extent_analysis'


def specification():
    old = baseline.specification()
    if sha(ROOT/WAKE_STUDY) != 'd79a1add8303d9f5b5de97d18396f7360a5dada3e44f7e88593ab3693c4d6f96' or read(ROOT/WAKE_STUDY) != old:
        raise ValueError('Frozen wake-study inputs changed; do not bypass its source guards')
    audit = read(ROOT/AUDIT/'summary.json')
    if not audit['full_state_audit_pass'] or not audit['control_shared_arrays_exactly_equal']:
        raise ValueError('Wake-study audit/control did not pass')
    jobs = []
    for job in old['jobs']:
        j = dict(job, wake_length=1., panels=[260,180] if job['control'] else [340,240])
        j['name'] = f"d{j['flap']}_cl{tag(j['target_cl'])}_"+('base_control' if j['control'] else 'fine')
        jobs.append(j)
    files = list(old['sha256']) + [WAKE_STUDY, AUDIT+'/summary.json', AUDIT+'/surface_diagnostics.csv',
        AUDIT+'/input_manifest.json', AUDIT+'/FINDINGS.md', AUDIT+'/source_manifest.json', AUDIT+'/xfoil_doc.txt',
        'validation/report_s702_wake_extent.py', 'validation/run_s702_holdout_panels.py',
        'validation/S702_HOLDOUT_PANELS_PROTOCOL.md', 'test_s702_holdout_panels.py']
    return dict(old, experiment='S702 reserved-condition panel sensitivity: four fine cases plus one base control',
        jobs=jobs, panels=[340,240], base_panels=[260,180], wake_length=1.,
        per_job_panel_override=True,
        qualification='Surface-panel setting and its derived wake discretization sensitivity; no independent validation.',
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def effective_spec(job, spec):
    expected = spec['base_panels'] if job['control'] else spec['panels']
    if job['panels'] != expected or job['wake_length'] != spec['wake_length']:
        raise ValueError('Job settings differ from declared panel protocol')
    return dict(spec, panels=list(job['panels']))


def worker(out, job, spec):
    # The unchanged worker exports this explicit job-specific configuration.
    return baseline.worker(out, job, effective_spec(job,spec))


def summarize(out, spec):
    summary = baseline.summarize(out,spec)
    for row in summary['jobs']:
        row['difference_from_base_CD_counts'] = row.pop('difference_from_1c_CD_counts')
    summary['qualification'] = 'Full selected-state audit, local surface resolution and paired grid interpretation remain required.'
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
    alive='run_s702_holdout_panels.py' in cmd and sys.argv[2] in cmd
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
                print(f"\nS702 holdout panels / jony4: {state}; {status.get('completed_jobs', 0)}/5 completed")
                for row in data['jobs']:
                    print(f"{row['job']:27s} {row['state']:10s} solves={row['evaluations']:2d} failed={row['failed']:2d} {row['path']}")
                    if row['error']: print('  '+row['error'])
                if data['summary']:
                    s = data['summary']
                    print(f"Feasible: {s['feasible_jobs']}/5; fresh base control: {'PASS' if s['control_pass'] else 'FAIL'}")
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
