#!/usr/bin/env python3
"""Read live jony4 progress. --watch refreshes until completion; Ctrl-C only stops watching."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REMOTE = '/home/jony/multifoil-runs/s702-policy-20260921'
TERMINAL = {'finished', 'finished_with_regression_failures', 'execution_failed'}
PROBE = r'''
from pathlib import Path
from datetime import datetime,timezone
import json,sys
r=Path(sys.argv[1]); rt=r/'handover/policy_runtime'; study=r/'validation/results/s702_initialization_policy_20260921'
def read(p):
    return json.loads(p.read_text()) if p.exists() else None
launch=read(rt/'launch_status.json') or {}
pid=launch.get('launcher_pid')
try:
    command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\x00', b' ').decode()
    alive='resume_s702_policy.py' in command or 'launch_s702_policy.sh' in command
except (OSError,TypeError):
    alive=False
spec=read(study/'study.json') or {'jobs':[]}
jobs=[]
for j in spec['jobs']:
    folder=study/j['name']; result=read(folder/'result.json'); progress=result or read(folder/'progress.json') or {}
    records=progress.get('records',[])
    status=('CONVERGED' if result['feasible'] else 'FAILED') if result else ('RUNNING' if folder.exists() else 'PENDING')
    init=records[-1].get('initialization',{}) if records else {}
    jobs.append({'job':j['name'],'status':status,'evaluations':len(records),
                 'failed_evaluations':sum(not x['admissible'] for x in records),
                 'path':progress.get('selected_path') or init.get('path',''),
                 'error':progress.get('error','')})
print(json.dumps({'observed_utc':datetime.now(timezone.utc).isoformat(),'launch':launch,
    'launcher_alive':alive,'study':read(study/'status.json'),'summary':read(study/'summary.json'),
    'known_baseline_failure':(rt/'baseline_classification.json').exists(),'jobs':jobs,'remote_root':str(r)}))
'''


def fetch(host='jony4-anywhere'):
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8',
        '-o', 'StrictHostKeyChecking=yes', host, 'python3 - '+shlex.quote(REMOTE)],
        input=PROBE, text=True, capture_output=True, timeout=20)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or 'Remote status probe failed')
    return json.loads(result.stdout)


def render(data):
    phase = data['launch'].get('phase', 'not launched')
    print('MultiFoil — jony4 — '+phase)
    print('Observed: '+data['observed_utc'])
    status = data['study'] or {}
    print(f"Completed: {status.get('completed_jobs', 0)}/{status.get('total_jobs', 7)} jobs")
    if data['known_baseline_failure']:
        print('Known baseline issue: the 100-chord element-order test remains failed.')
    if phase not in TERMINAL and not data['launcher_alive']:
        print('ATTENTION: launcher is no longer active; inspect the log before resuming.')
    if data['launch'].get('error'):
        print('Error: '+data['launch']['error'])
    print(f"{'Job':32s} {'State':10s} {'Solves':>6s} {'Failed':>6s}  Initialization")
    for j in data['jobs']:
        print(f"{j['job']:32s} {j['status']:10s} {j['evaluations']:6d} {j['failed_evaluations']:6d}  {j['path']}")
    summary = data['summary']
    if summary:
        print(f"Primary converged: {summary['converged_primary']}/6; replay control: {summary['converged_control']}/1")
        print('Force regression: '+('PASS' if summary['force_regression_pass'] else 'FAIL'))
        print('Full selected-state review is still required before optimizer integration.')
    elif phase not in TERMINAL:
        print('Counts are requested solves; internal continuation steps are additional.')
    print('Log: '+REMOTE+'/handover/policy_runtime/continuation.log')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('--interval', type=float, default=10.)
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()
    if not 2 <= args.interval <= 300:
        ap.error('--interval must be between 2 and 300 seconds')
    try:
        while True:
            try:
                data = fetch()
            except (RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
                print(f'Status unavailable: {exc}', file=sys.stderr)
                if not args.watch:
                    return 3
                time.sleep(args.interval)
                continue
            if args.watch and sys.stdout.isatty() and not args.json:
                print('\033[2J\033[H', end='')
            if args.json:
                print(json.dumps(data, indent=2), flush=True)
            else:
                render(data); sys.stdout.flush()
            phase = data['launch'].get('phase')
            if phase in TERMINAL:
                return data['launch'].get('exit_code') or 0
            if not data['launcher_alive']:
                return 1
            if not args.watch:
                return 0
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print('\nStopped watching; the remote simulation continues.')
        return 130


if __name__ == '__main__':
    raise SystemExit(main())
