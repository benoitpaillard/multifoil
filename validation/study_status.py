#!/usr/bin/env python3
"""Read jony4's current design integration progress without controlling the job."""
import argparse
import json
import shlex
import subprocess
import time

REMOTE = '/home/jony/multifoil-runs/design-initialization-20260921'
STUDY = 'validation/results/design_initialization_20260921'
PROBE = r'''
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); out=root/sys.argv[2]
def read(p): return json.loads(p.read_text()) if p.exists() else {}
status=read(out/'status.json'); spec=read(out/'study.json')
pid=status.get('pid')
try:
    command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive='run_design_initialization.py' in command and sys.argv[2] in command
except OSError: alive=False
jobs=[]
for j in spec.get('jobs',[]):
    f=out/j['name']; result=read(f/'result.json'); progress=result or read(f/'progress.json')
    records=progress.get('records',[])
    jobs.append(dict(job=j['name'], state=('PASS' if result['feasible'] else 'FAILED') if result else ('RUNNING' if f.exists() else 'PENDING'),
        evaluations=len(records), failed=sum(not r['admissible'] for r in records),
        path=progress.get('selected_path') or (records[-1]['initialization']['path'] if records else ''),
        error=progress.get('failure_reason','')))
print(json.dumps(dict(status=status,alive=alive,jobs=jobs,summary=read(out/'summary.json'),remote_root=str(root))))
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--host', default='jony4-anywhere')
    ap.add_argument('--remote-root', default=REMOTE)
    ap.add_argument('--study', default=STUDY)
    args = ap.parse_args()
    try:
        while True:
            try:
                probe = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8',
                    '-o', 'StrictHostKeyChecking=yes', args.host,
                    'python3 - '+shlex.quote(args.remote_root)+' '+shlex.quote(args.study)],
                    input=PROBE, text=True, capture_output=True, timeout=20)
                if probe.returncode: raise RuntimeError(probe.stderr.strip())
                data = json.loads(probe.stdout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
                print('Status connection failed: '+str(exc), flush=True)
                if not args.watch: return 3
                time.sleep(10); continue
            status = data['status']; terminal = status.get('state') in ('finished', 'execution_failed')
            if args.json:
                print(json.dumps(data, indent=2), flush=True)
            else:
                print(f"\nMultiFoil / jony4: {status.get('state', 'not started')} — "
                      f"{status.get('completed_jobs', 0)}/{status.get('total_jobs', 12)} completed", flush=True)
                for j in data['jobs']:
                    print(f"{j['job']:24s} {j['state']:8s} solves={j['evaluations']:2d} failed={j['failed']:2d} {j['path']}")
                    if j['error']: print('  '+j['error'])
                if data['summary']:
                    s = data['summary']
                    print(f"Feasible: {s['converged_jobs']}/{s['required_jobs']}; pairs passed: {s['passed_pairs']}/{s['required_pairs']}")
                    print('Selected-state review remains required.')
                if status.get('error'): print(status['error'])
                if not data['alive'] and not terminal: print('Runner inactive; inspect launch.log before resuming.')
                print('Log: '+args.remote_root+'/launch.log', flush=True)
            if terminal:
                return (0 if status.get('regression_pass') else 2) if status['state'] == 'finished' else 1
            if not data['alive']: return 1
            if not args.watch: return 0
            time.sleep(10)
    except KeyboardInterrupt:
        print('\nMonitor stopped; remote job is unaffected.')
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
