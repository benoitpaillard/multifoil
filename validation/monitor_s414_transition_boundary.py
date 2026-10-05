#!/usr/bin/env python3
"""Read-only monitor for the S414 transition-boundary paths (V2 tranche 1).

This independent monitor does not import or modify the frozen numerical runner.
"""
import argparse
import json
import subprocess
import time

REMOTE = '/home/jony/multifoil-runs/s414-transition-boundary-20261004'
PROBE = r'''
import json
from pathlib import Path
root=Path('/home/jony/multifoil-runs/s414-transition-boundary-20261004')
out=root/'validation/results/s414_transition_boundary_20261004'
def read(p):
    try: return json.loads(p.read_text()) if p.exists() else {}
    except ValueError: return {}
status=read(out/'status.json'); spec=read(out/'study.json')
pid=status.get('pid',read(root/'launch_receipt.json').get('pid'))
try:
    proc=Path('/proc')/str(pid)
    cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive=(proc/'cwd').resolve()==root and ('run_s414_transition_boundary.py' in cmd or 'run_s414_transition_boundary.sh' in cmd)
except OSError: alive=False
jobs=[]
for j in spec.get('jobs',[]):
    f=out/j['name']; result=read(f/'result.json'); progress=result or read(f/'progress.json')
    records=progress.get('records',[])
    jobs.append(dict(job=j['name'],state='DONE' if result else (('RUNNING' if alive else 'INTERRUPTED') if f.exists() else 'PENDING'),
        requested=len(j['requests']),evaluated=len(records),converged=sum(r['converged'] for r in records),
        failed=sum(not r['converged'] for r in records),error=progress.get('error','')))
preflight=root/'preflight.log'; tail=preflight.read_text().strip().splitlines()[-1] if preflight.exists() and preflight.read_text().strip() else ''
print(json.dumps(dict(status=status,alive=alive,pid=pid,jobs=jobs,controls=read(out/'controls.json').get('passed'),
    summary=read(out/'summary.json'),preflight_last_line=tail)))
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    try:
        while True:
            try:
                response = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
                    '-o', 'ConnectTimeout=15', 'jony4-anywhere', 'python3 -'], input=PROBE,
                    text=True, capture_output=True, timeout=25, check=True)
                data = json.loads(response.stdout)
            except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
                print('Status connection failed: '+str(exc), flush=True)
                if not args.watch: return 3
                time.sleep(10); continue
            state = data['status'].get('state', 'preflight' if data['alive'] else 'inactive')
            done = sum(j['state'] == 'DONE' for j in data['jobs'])
            if args.json: print(json.dumps(data, indent=2), flush=True)
            else:
                print(f"\nS414 transition boundary / jony4: {state}; {done}/{len(data['jobs'])} jobs completed; process alive={data['alive']}")
                if state == 'preflight': print('Preflight: '+data['preflight_last_line'])
                for job in data['jobs']:
                    print(f"{job['job']:40s} {job['state']:11s} evaluated={job['evaluated']}/{job['requested']} "
                          f"converged={job['converged']} failed={job['failed']}")
                    if job['error']: print('  '+job['error'])
                if data['controls'] is not None:
                    print('Legacy endpoint controls identical to archive: '+('PASS' if data['controls'] else 'FAIL'))
                if data['summary']:
                    s = data['summary']
                    print(f"Path points converged/failed/requested: {s['converged_path_points']}/{s['failed_path_points']}/{s['requested_path_points']}")
                    print('Execution gates: '+('PASS' if s['execution_checks_pass'] else 'FAIL'))
                    print('Nonconverged points are findings, not forces. Full saved-state audit and interpretation still required.')
                if data['status'].get('error'): print(data['status']['error'])
                print('Log: '+REMOTE+'/launch.log', flush=True)
            if state == 'finished': return 0 if data['summary'].get('execution_checks_pass') else 2
            if state == 'execution_failed' or not data['alive']: return 1
            if not args.watch: return 0
            time.sleep(10)
    except KeyboardInterrupt:
        print('\nMonitor stopped; remote jobs are unaffected.')
        return 0


if __name__ == '__main__': raise SystemExit(main())
