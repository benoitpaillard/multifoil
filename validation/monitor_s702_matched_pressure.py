#!/usr/bin/env python3
"""Read-only monitor for the S702 matched pressure study.

This independent monitor does not import or modify the frozen numerical runner.
"""
import argparse
import json
import subprocess
import time

REMOTE = '/home/jony/multifoil-runs/s702-matched-pressure-20261002'
PROBE = r'''
import json
from pathlib import Path
root=Path('/home/jony/multifoil-runs/s702-matched-pressure-20261002')
out=root/'validation/results/s702_matched_pressure_20261002'
def read(p): return json.loads(p.read_text()) if p.exists() else {}
status=read(out/'status.json'); spec=read(out/'study.json')
pid=status.get('pid',read(root/'launch_receipt.json').get('pid'))
try:
    proc=Path('/proc')/str(pid)
    cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive=(proc/'cwd').resolve()==root and ('run_s702_matched_pressure.py' in cmd or 'run_matched_pressure.sh' in cmd)
except OSError: alive=False
jobs=[]
for j in spec.get('jobs',[]):
    f=out/j['name']; result=read(f/'result.json'); progress=result or read(f/'progress.json')
    records=progress.get('records',[])
    jobs.append(dict(job=j['name'],state=('CONVERGED' if result['feasible'] else 'FAILED') if result else ('RUNNING' if f.exists() and alive else 'PENDING'),
        evaluations=len(records),failed=sum(not r['admissible'] for r in records),error=progress.get('error','')))
print(json.dumps(dict(status=status,alive=alive,pid=pid,jobs=jobs,summary=read(out/'summary.json'))))
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch',action='store_true')
    parser.add_argument('--json',action='store_true')
    args=parser.parse_args()
    try:
        while True:
            try:
                response=subprocess.run(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',
                    '-o','ConnectTimeout=15','jony4-anywhere','python3 -'],input=PROBE,
                    text=True,capture_output=True,timeout=25,check=True)
                data=json.loads(response.stdout)
            except (OSError,subprocess.SubprocessError,json.JSONDecodeError) as exc:
                print('Status connection failed: '+str(exc),flush=True)
                if not args.watch:return 3
                time.sleep(10);continue
            state=data['status'].get('state','preflight' if data['alive'] else 'inactive')
            done=sum(j['state'] in ('CONVERGED','FAILED') for j in data['jobs'])
            if args.json:print(json.dumps(data,indent=2),flush=True)
            else:
                print(f"\nS702 matched pressure / jony4: {state}; {done}/{len(data['jobs'])} completed; process alive={data['alive']}")
                for job in data['jobs']:
                    print(f"{job['job']:25s} {job['state']:10s} evaluations={job['evaluations']:2d} failed={job['failed']:2d}")
                    if job['error']:print('  '+job['error'])
                if data['summary']:
                    print('Execution gates: '+('PASS' if data['summary']['execution_checks_pass'] else 'FAIL'))
                    print('Full selected-state audit and scientific interpretation still required.')
                if data['status'].get('error'):print(data['status']['error'])
                print('Log: '+REMOTE+'/launch.log',flush=True)
            if state=='finished':return 0 if data['summary'].get('execution_checks_pass') else 2
            if state=='execution_failed' or not data['alive']:return 1
            if not args.watch:return 0
            time.sleep(10)
    except KeyboardInterrupt:
        print('\nMonitor stopped; remote jobs are unaffected.')
        return 0


if __name__=='__main__':raise SystemExit(main())
