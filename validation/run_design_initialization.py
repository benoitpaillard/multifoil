#!/usr/bin/env python3
"""Restartable Mach-zero integration check through the actual design flow evaluator."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import copy
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import read_dat
from shape_optimization import CODE, ShapeFamily, check_config, flow_job, write_json

OLD = 'validation/results/s702_initialization_policy_20260921/study.json'
AUDIT = 'validation/results/s702_initialization_policy_analysis/summary.json'
VERSIONS = {'numpy': '2.3.5', 'scipy': '1.17.0', 'matplotlib': '3.10.8'}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def specification():
    if sha(ROOT/OLD) != '5388fef33ad7833c66b9d8db43da4c6304a54a4d50a4aae38f1da2a2d27219cb':
        raise ValueError('Unexpected prior policy study')
    previous = read(ROOT/OLD)
    if not read(ROOT/AUDIT)['full_state_audit_pass']:
        raise ValueError('Previous full-state audit did not pass')
    for name in ('multifoil.py', 'viscous.py', 'paneling.py', 'optimization.py', 'vendor/mfoil.py', 'fixed_lift_policy.py'):
        if sha(ROOT/name) != previous['sha256'][name]:
            raise ValueError('Flow core or verified helper changed: '+name)
    cfg = read(ROOT/'examples/s702_incompressible_design.json')
    cfg.update(reynolds=700000., scenarios=[dict(name='ncrit6p25', ncrit=6.25)])
    jobs = []
    for flap in (5, 10):
        center = next(j['center'] for j in previous['jobs'] if j['flap'] == flap and j['target_cl'] == .8)
        for label, amplitude in (('minus', -.001), ('baseline', 0.), ('plus', .001)):
            for mode in ('legacy', 'policy'):
                jobs.append(dict(name=f'd{flap}_{label}_{mode}', flap=flap, amplitude=amplitude,
                                 mode=mode, center=center, pair=f'd{flap}_{label}'))
    files = CODE + ['validation/run_design_initialization.py', 'validation/DESIGN_INITIALIZATION_PROTOCOL.md',
        'validation/study_status.py', 'test_shape_initialization.py', 'test_shape_optimization.py',
        'test_fixed_lift_policy.py', 'test_optimization.py', 'test_paneling.py',
        'examples/s702_incompressible_design.json', OLD, AUDIT,
        'validation/data/s702/fore.dat', 'validation/data/s702/aft_d5.dat', 'validation/data/s702/aft_d10.dat']
    return dict(schema=1, experiment='Mach-zero design evaluator initialization integration', config=cfg,
        jobs=jobs, target_cl=.8, legacy_half_width=.4, policy_half_width=.05,
        alpha_domain=[-6., 8.], warm_start_alpha=-1., acceptance_CD_counts=.01,
        expected_versions=VERSIONS, qualification='Implementation comparison, not independent aerodynamic validation or an optimization campaign.',
        sha256={f: sha(ROOT/f) for f in sorted(set(files))})


def prepare(out, spec):
    out.mkdir(parents=True, exist_ok=True)
    manifest = out/'study.json'
    if manifest.exists():
        if read(manifest) != spec:
            raise ValueError('Code, inputs or settings changed; use a new study directory')
    else:
        if any(p.name != 'run.lock' for p in out.iterdir()):
            raise ValueError('Nonempty output directory without a study manifest')
        write_json(manifest, spec)
    for name, digest in spec['sha256'].items():
        saved = out/'source_snapshot'/name
        if not saved.exists():
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, saved)
        if sha(saved) != digest:
            raise ValueError('Source snapshot changed: '+name)


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        cfg = copy.deepcopy(spec['config'])
        cfg['elements'][0]['file'] = '../validation/data/s702/fore.dat'
        cfg['elements'][1]['file'] = f"../validation/data/s702/aft_d{job['flap']}.dat"
        center = job['center']
        point = dict(name='cl08', cl=spec['target_cl'], weight=1., start_alpha=spec['warm_start_alpha'],
                     alpha_bounds=[center-spec['legacy_half_width'], center+spec['legacy_half_width']])
        if job['mode'] == 'policy':
            point.update(alpha_bounds=spec['alpha_domain'], initialization=dict(center=center,
                half_width=spec['policy_half_width'], max_brackets=2, max_evaluations=20,
                paths=[dict(name='cold_local', start_alpha=center),
                       dict(name='warm_anchor', start_alpha=spec['warm_start_alpha'])]))
        cfg['points'] = [point]
        check_config(cfg)
        raw = [(e['name'], read_dat(ROOT/'examples'/e['file'])) for e in cfg['elements']]
        raw = ShapeFamily(raw, cfg['modes']).deform([job['amplitude']])
        for name, xy in raw:
            np.savetxt(folder/f'input_{name}.dat', xy, header=name, comments='', fmt='%.16g')
        write_json(folder/'configuration.json', cfg)
        return flow_job((folder, raw, cfg, cfg['panels'], cfg['scenarios'][0], point))


def summarize(out, spec):
    results = {j['name']: read(out/j['name']/'result.json') for j in spec['jobs']}
    pairs = []
    for name in sorted({j['pair'] for j in spec['jobs']}):
        old, new = results[name+'_legacy'], results[name+'_policy']
        feasible = old['feasible'] and new['feasible']
        difference = (new['CD']-old['CD'])*1e4 if feasible else None
        valid = feasible and abs(difference) <= spec['acceptance_CD_counts']
        for r in (old, new):
            if r['feasible']:
                valid = valid and abs(r['CL']-spec['target_cl']) <= 2e-5 and r['residual'] <= 1e-7
        pairs.append(dict(pair=name, complete=feasible, difference_CD_counts=difference,
                          passed=bool(valid), policy_path=new.get('selected_path')))
    summary = dict(completed_utc=stamp(), required_jobs=len(spec['jobs']),
        converged_jobs=sum(r['feasible'] for r in results.values()), required_pairs=len(pairs),
        passed_pairs=sum(r['passed'] for r in pairs), pairs=pairs,
        regression_pass=all(r['passed'] for r in pairs),
        requested_evaluations=sum(r['flow_evaluations'] for r in results.values()),
        qualification='Force/tolerance checks only. Audit selected surfaces, transition and saved policy states before broader design use.')
    write_json(out/'summary.json', summary)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 6:
        ap.error('Use a results directory and 1–6 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(json.dumps(dict(state='prepared', jobs=len(spec['jobs']), study_sha256=sha(out/'study.json'))))
            return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != VERSIONS or sys.version_info[:2] != (3, 12):
            raise ValueError('Use pinned Python 3.12 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
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
                    archived = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                    archived.parent.mkdir(parents=True, exist_ok=True); shutil.move(folder, archived)
                pending.append(job)
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(worker, out, j, spec): j for j in pending}
                for future in as_completed(futures):
                    r = future.result(); done += 1
                    status('running', last_job=futures[future]['name'], last_job_feasible=r['feasible'])
                    print(f"{done}/{len(spec['jobs'])}: {futures[future]['name']} feasible={r['feasible']}", flush=True)
            summary = summarize(out, spec)
            status('finished', regression_pass=summary['regression_pass'])
            return 0 if summary['regression_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}')
            raise


if __name__ == '__main__':
    raise SystemExit(main())
