#!/usr/bin/env python3
"""Run/resume the final bounded S702 third-grid check and monitor jony4."""
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
import run_s702_holdout_panels as previous
import run_s702_wake_extent as flow

STUDY = 'validation/results/s702_panel_convergence_20260922'
REMOTE = '/home/jony/multifoil-runs/s702-panel-convergence-20260922'
FINE = 'validation/results/s702_holdout_panels_20260922'
AUDIT = 'validation/results/s702_holdout_panels_analysis'


def specification():
    prior = previous.specification()
    if sha(ROOT/FINE/'study.json') != 'cd861525d01a02d788bbaddabaacb0c81f9896e153ab2229743190b775fea6f1' or read(ROOT/FINE/'study.json') != prior:
        raise ValueError('Frozen fine-grid inputs changed; do not bypass source guards')
    audit = read(ROOT/AUDIT/'summary.json')
    if not audit['full_state_audit_pass'] or not audit['control_shared_arrays_exactly_equal']:
        raise ValueError('Fine-grid audit/control did not pass')
    pinned = read(ROOT/AUDIT/'input_manifest.json')['files']
    files = list(prior['sha256']) + [FINE+'/study.json', FINE+'/runtime.json', AUDIT+'/summary.json',
        AUDIT+'/input_manifest.json', AUDIT+'/FINDINGS.md', AUDIT+'/targets.csv', AUDIT+'/flap_deltas.csv',
        AUDIT+'/geometry_comparison.csv', AUDIT+'/surface_diagnostics.csv',
        'validation/run_s702_panel_convergence.py', 'validation/S702_PANEL_CONVERGENCE_PROTOCOL.md',
        'validation/report_s702_holdout_panels.py', 'test_s702_panel_convergence.py']
    jobs = []
    for job in prior['jobs']:
        if job['control']: continue
        folder = Path(FINE)/job['name']
        refpath = str(folder/'result.json')
        ref = read(ROOT/refpath)
        if not ref['feasible'] or abs(ref['CL']-job['target_cl']) > prior['cl_tolerance']:
            raise ValueError('Missing feasible fine-grid reference')
        for name,digest in pinned.items():
            if Path(name).is_relative_to(folder):
                if sha(ROOT/name) != digest: raise ValueError('Audited fine-grid input changed: '+name)
                files.append(name)
        jobs.append(dict(job,name=f"d{job['flap']}_cl{tag(job['target_cl'])}_finer",panels=[440,320],
            reference_result=refpath,coarse_reference_result=job['reference_result']))
    control = dict(next(j for j in jobs if j['flap']==10 and j['target_cl']==.8))
    control.update(name='d10_cl0p8_fine_control',panels=[340,240],control=True)
    jobs.append(control)
    return dict(prior,experiment='S702 third-grid check: four 440/320 cases plus one 340/240 control',
        jobs=jobs,panels=[440,320],base_panels=[340,240],grid_sequence=[[260,180],[340,240],[440,320]],
        last_automatic_grid_extension=True,
        qualification='Three-grid operational sensitivity; no formal error bound or independent aerodynamic validation.',
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def worker(out,job,spec):
    return flow.worker(out,job,previous.effective_spec(job,spec))


def summarize(out,spec):
    summary = flow.summarize(out,spec)
    for row in summary['jobs']:
        row['difference_from_fine_CD_counts'] = row.pop('difference_from_1c_CD_counts')
    summary['qualification'] = 'Full selected-state audit and three-grid contraction decision remain required.'
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
    alive='run_s702_panel_convergence.py' in cmd and sys.argv[2] in cmd
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
                print(f"\nS702 third grid / jony4: {state}; {status.get('completed_jobs', 0)}/5 completed")
                for row in data['jobs']:
                    print(f"{row['job']:27s} {row['state']:10s} solves={row['evaluations']:2d} failed={row['failed']:2d} {row['path']}")
                    if row['error']: print('  '+row['error'])
                if data['summary']:
                    s = data['summary']
                    print(f"Feasible: {s['feasible_jobs']}/5; fresh 340/240 control: {'PASS' if s['control_pass'] else 'FAIL'}")
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
