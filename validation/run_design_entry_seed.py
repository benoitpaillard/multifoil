#!/usr/bin/env python3
"""Run the entry diagnostic or a fresh combined-policy regression; inspect remote progress."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import copy
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
from shape_optimization import CODE, ShapeFamily, check_config, flow_job, write_json
from multifoil import read_dat
import numpy as np

OLD = 'validation/results/design_initialization_20260921'
STUDY = 'validation/results/design_entry_seed_20260921'
REMOTE = '/home/jony/multifoil-runs/design-entry-seed-20260921'
COMBINED_STUDY = 'validation/results/design_combined_policy_20260921'
COMBINED_REMOTE = '/home/jony/multifoil-runs/design-combined-policy-20260921'
ENTRY_AUDIT = 'validation/results/design_entry_seed_analysis/summary.json'


def specification(combined=False):
    if sha(ROOT/OLD/'study.json') != 'f95ec243e739b69d883062ae9ecb37133b0e71910d62e352d6cdafd103a201f6':
        raise ValueError('Unexpected original integration study')
    original = read(ROOT/OLD/'study.json')
    for name in CODE:
        if sha(ROOT/name) != original['sha256'][name]:
            raise ValueError('Verified evaluator or flow core changed: '+name)
    jobs = []
    references = []
    for flap in ((5, 10) if combined else (10,)):
        for label in ('baseline', 'minus', 'plus'):
            prior = next(j for j in original['jobs'] if j['name'] == f'd{flap}_{label}_legacy')
            path = f'{OLD}/{prior["name"]}/result.json'
            ref = read(ROOT/path)
            if not ref['feasible']:
                raise ValueError('Missing feasible legacy reference')
            jobs.append(dict(name=f'd{flap}_{label}_'+('combined' if combined else 'entry'),
                flap=flap, amplitude=prior['amplitude'], center=prior['center'],
                seed_alpha=prior['center']-original['legacy_half_width'], reference_result=path))
            references.append(path)
    if combined:
        if not read(ROOT/ENTRY_AUDIT)['full_state_audit_pass']:
            raise ValueError('Entry diagnostic requires a passed full-state audit')
        repeat = copy.deepcopy(next(j for j in jobs if j['name']=='d10_baseline_combined'))
        repeat['name'] = 'd10_baseline_repeat_combined'; repeat['repeat_of'] = 'd10_baseline_combined'
        jobs.append(repeat)
    files = CODE+['validation/run_design_initialization.py', 'validation/run_design_entry_seed.py',
        'validation/DESIGN_ENTRY_SEED_PROTOCOL.md', 'test_shape_initialization.py', 'test_shape_optimization.py',
        'test_fixed_lift_policy.py', 'test_optimization.py', 'test_paneling.py',
        'examples/s702_incompressible_design.json', 'validation/data/s702/fore.dat',
        'validation/data/s702/aft_d10.dat', OLD+'/study.json']+references
    if combined:
        files += ['validation/data/s702/aft_d5.dat', ENTRY_AUDIT,
                  'validation/DESIGN_COMBINED_POLICY_PROTOCOL.md', 'test_design_combined_policy.py']
    return dict(schema=1, experiment=('Combined cold/warm/entry policy on six contours plus baseline repeat' if combined
                                     else 'Explicit lower-entry seed after Mach-zero baseline initialization failure'),
        combined_policy=combined,
        jobs=jobs, config=original['config'], target_cl=.8, warm_start_alpha=-1., alpha_domain=[-6., 8.],
        half_width=.05, max_brackets=2, max_evaluations=20, acceptance_CD_counts=.01,
        expected_versions=VERSIONS, original_study_sha256=sha(ROOT/OLD/'study.json'),
        qualification=('Fresh numerical regression of the complete ordered policy; independent aerodynamic ranking remains unvalidated.' if combined
                       else 'Diagnostic of one additional path only; does not replace the failed twelve-job matrix.'),
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def configuration(job, spec):
    cfg = copy.deepcopy(spec['config'])
    cfg['elements'][1]['file'] = f"../validation/data/s702/aft_d{job['flap']}.dat"
    paths = [dict(name='warm_entry', start_alpha=spec['warm_start_alpha'], seed_alphas=[job['seed_alpha']])]
    if spec.get('combined_policy'):
        paths = [dict(name='cold_local', start_alpha=job['center']),
                 dict(name='warm_anchor', start_alpha=spec['warm_start_alpha'])]+paths
    cfg['points'] = [dict(name='cl08', cl=spec['target_cl'], weight=1., start_alpha=spec['warm_start_alpha'],
        alpha_bounds=spec['alpha_domain'], initialization=dict(center=job['center'], half_width=spec['half_width'],
            max_brackets=spec['max_brackets'], max_evaluations=spec['max_evaluations'],
            paths=paths))]
    check_config(cfg)
    return cfg


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        cfg = configuration(job, spec)
        raw = [(e['name'], read_dat(ROOT/'examples'/e['file'])) for e in cfg['elements']]
        raw = ShapeFamily(raw, cfg['modes']).deform([job['amplitude']])
        write_json(folder/'configuration.json', cfg)
        for name, xy in raw:
            np.savetxt(folder/f'input_{name}.dat', xy, header=name, comments='', fmt='%.16g')
        return flow_job((folder, raw, cfg, cfg['panels'], cfg['scenarios'][0], cfg['points'][0]))


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json'); ref = read(ROOT/job['reference_result'])
        good = r['feasible']
        delta = (r['CD']-ref['CD'])*1e4 if good else None
        passed = good and abs(delta) <= spec['acceptance_CD_counts'] and abs(r['CL']-spec['target_cl']) <= 2e-5 and r['residual'] <= 1e-7
        rows.append(dict(job=job['name'], feasible=good, difference_CD_counts=delta, passed=bool(passed),
            selected_path=r.get('selected_path'),
            evaluations=r['flow_evaluations'], failed_evaluations=sum(not a['admissible'] for a in r['records'])))
    repeats = []
    for job in spec['jobs']:
        if 'repeat_of' not in job: continue
        a, b = read(out/job['name']/'result.json'), read(out/job['repeat_of']/'result.json')
        good = a['feasible'] and b['feasible']
        delta = (a['CD']-b['CD'])*1e4 if good else None
        repeats.append(dict(job=job['name'], repeat_of=job['repeat_of'],
                            difference_CD_counts=delta, passed=bool(good and abs(delta) <= spec['acceptance_CD_counts'])))
    result = dict(completed_utc=stamp(), required_jobs=len(rows), feasible_jobs=sum(r['feasible'] for r in rows),
        passed_jobs=sum(r['passed'] for r in rows), regression_pass=all(r['passed'] for r in rows) and all(r['passed'] for r in repeats),
        jobs=rows, repeats=repeats,
        qualification='Review all selected states, entry histories and repeatability before accepting the combined policy.' if spec.get('combined_policy')
                      else 'Review selected states, entry convergence and branch differences before adding this path to the common policy.')
    write_json(out/'summary.json', result)
    return result


PROBE = r'''
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); out=root/sys.argv[2]
def read(p): return json.loads(p.read_text()) if p.exists() else {}
status=read(out/'status.json'); spec=read(out/'study.json')
try:
    cmd=Path('/proc/'+str(status.get('pid'))+'/cmdline').read_bytes().replace(b'\0',b' ').decode()
    alive='run_design_entry_seed.py' in cmd and sys.argv[2] in cmd
except OSError: alive=False
jobs=[]
for job in spec.get('jobs',[]):
    folder=out/job['name']; r=read(folder/'result.json'); progress=r or read(folder/'progress.json')
    jobs.append(dict(job=job['name'], state=('PASS' if r['feasible'] else 'FAILED') if r else ('RUNNING' if folder.exists() else 'PENDING'),
        evaluations=len(progress.get('records',[])), error=progress.get('failure_reason','')))
print(json.dumps(dict(status=status,alive=alive,jobs=jobs,summary=read(out/'summary.json'),experiment=spec.get('experiment','Study'))))
'''


def remote_status(watch=False, as_json=False, combined=False):
    remote, study = (COMBINED_REMOTE, COMBINED_STUDY) if combined else (REMOTE, STUDY)
    try:
        while True:
            try:
                proc = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
                    '-o', 'StrictHostKeyChecking=yes', 'jony4-anywhere',
                    'python3 - '+shlex.quote(remote)+' '+shlex.quote(study)],
                    input=PROBE, text=True, capture_output=True, timeout=25)
                if proc.returncode: raise RuntimeError(proc.stderr.strip())
                data = json.loads(proc.stdout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
                print('Connection failed: '+str(exc), flush=True)
                if not watch: return 3
                time.sleep(10); continue
            state = data['status'].get('state', 'not started')
            if as_json: print(json.dumps(data, indent=2), flush=True)
            else:
                print(f"\n{data['experiment']}: {state}; {data['status'].get('completed_jobs', 0)}/{data['status'].get('total_jobs', len(data['jobs']))} completed")
                for row in data['jobs']:
                    print(f"{row['job']:24s} {row['state']:8s} solves={row['evaluations']} {row['error']}")
                if data['summary']:
                    print(f"Force/tolerance checks: {data['summary']['passed_jobs']}/{data['summary']['required_jobs']}. Selected-state review still required.")
                    print('Overall numerical check: '+('PASS' if data['summary']['regression_pass'] else 'FAIL'))
                    for repeat in data['summary'].get('repeats', []):
                        print(f"Repeat {repeat['job']}: difference={repeat['difference_CD_counts']} CD counts; "
                              +('PASS' if repeat['passed'] else 'FAIL'))
                if state not in ('finished', 'execution_failed') and not data['alive']:
                    print('Runner inactive; inspect launch.log before resuming.')
                if data['status'].get('error'): print(data['status']['error'])
                print('Log: '+remote+'/launch.log', flush=True)
            if state == 'finished': return 0 if data['summary'].get('regression_pass') else 2
            if state == 'execution_failed' or not data['alive']: return 1
            if not watch: return 0
            time.sleep(10)
    except KeyboardInterrupt:
        print('\nMonitor stopped; remote job is unaffected.')
        return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--prepare-only', action='store_true')
    ap.add_argument('--remote-status', action='store_true')
    ap.add_argument('--watch', action='store_true')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--combined', action='store_true', help='Run or inspect the fresh seven-job combined-policy regression')
    args = ap.parse_args()
    if args.remote_status: return remote_status(args.watch, args.json, args.combined)
    if args.watch or args.json: ap.error('--watch and --json require --remote-status')
    out = (args.out or ROOT/(COMBINED_STUDY if args.combined else STUDY)).resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= (6 if args.combined else 3):
        ap.error('Use a results directory and 1–3 workers (up to 6 for --combined)')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(args.combined); prepare(out, spec)
        if args.prepare_only:
            print(json.dumps(dict(state='prepared', jobs=len(spec['jobs']), study_sha256=sha(out/'study.json'))))
            return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k:importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != VERSIONS or sys.version_info[:2] != (3, 12):
            raise ValueError('Use pinned Python 3.12 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a fresh study directory')
        write_json(out/'runtime.json', runtime)
        done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write_json(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(),
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
            status('finished', regression_pass=summary['regression_pass'])
            return 0 if summary['regression_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}')
            raise


if __name__ == '__main__':
    raise SystemExit(main())
