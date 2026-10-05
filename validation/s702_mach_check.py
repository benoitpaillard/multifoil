#!/usr/bin/env python3
"""Isolate the effect of Mach 0.1 -> 0 at fixed lift on the S702 ranking.

The Mach-zero points are a sensitivity experiment, not matched-condition
validation against the published Mach-0.1 MSES curves.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from multifoil import repanel, save_results
from viscous import CoupledSolver
from optimization import solve_fixed_lift, LiftTargetError
from s702_geometry import geometry, FILES
from s702_benchmark import record


def run_job(out, flap, mach, spec):
    folder = Path(out)/f'flap{flap}_mach{mach:g}'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/'result.json'
    if path.exists():
        return json.loads(path.read_text())
    elements = [(n, repanel(x, p)) for (n, x), p in zip(geometry(flap), spec['panels'])]
    bounds = (1., 2.) if flap == 5 else (-1., 0.)
    result = dict(flap=flap, mach=mach, target_cl=.8, alpha_bounds=bounds,
                  feasible=False, CL=None, CD=None, CM=None, alpha_deg=None, error='', records=[])
    solutions = []

    def checkpoint(sol):
        solutions.append(sol)
        result['records'].append(record(sol, 'fixed_lift'))
        save_results(folder/'flow', solutions, elements, spec, 'viscous')
        (folder/'progress.json').write_text(json.dumps(result, indent=2)+'\n')

    try:
        solver = CoupledSolver(elements, reynolds=1e6, mach=mach, ncrit=9.,
                               iterations=100, tolerance=1e-7, start_alpha=-1.)
        sol, _ = solve_fixed_lift(solver, .8, bounds, 2e-5, 20, checkpoint)
        result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                      boundary_layers=sol.diagnostics['boundary_layers'],
                      element_forces=sol.diagnostics['element_forces'])
    except (LiftTargetError, ValueError, FloatingPointError) as exc:
        result['error'] = f'{type(exc).__name__}: {exc}'
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_mach_check')
    ap.add_argument('--workers', type=int, default=2)
    args = ap.parse_args()
    if not 1 <= args.workers <= 4:
        ap.error('workers must be between 1 and 4')
    code = ['multifoil.py', 'viscous.py', 'paneling.py', 'optimization.py', 'vendor/mfoil.py',
            'validation/s702_geometry.py', 'validation/s702_benchmark.py',
            'validation/s702_mach_check.py']
    data = ['validation/data/s702/'+f for f in ('fore.dat', FILES[5], FILES[10], 'source.json')]
    spec = dict(schema=1, reynolds=1e6, ncrit=9., panels=[260, 180], start_alpha=-1.,
                target_cl=.8, cl_tolerance=2e-5, iterations=100, residual_tolerance=1e-7,
                wake_length=1., jobs=[(f, m) for f in (5, 10) for m in (.1, 0.)],
                reference_chord=1., moment_reference=[.25, 0.],
                reference_status='S702 vector geometry; M=0 is a sensitivity calculation',
                sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in code+data})
    spec = json.loads(json.dumps(spec))
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != spec:
            raise ValueError('Study/code changed; choose a new output directory')
    elif any(args.out.iterdir()):
        raise ValueError('Nonempty directory lacks a manifest')
    else:
        manifest.write_text(json.dumps(spec, indent=2)+'\n')
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(run_job, str(args.out), f, m, spec) for f, m in spec['jobs']]
        for future in as_completed(futures):
            r = future.result()
            print({k: r[k] for k in ('flap', 'mach', 'feasible', 'CD', 'error')}, flush=True)
    records = [json.loads((args.out/f'flap{f}_mach{m:g}/result.json').read_text()) for f, m in spec['jobs']]
    (args.out/'summary.json').write_text(json.dumps(records, indent=2)+'\n')
    return 0 if all(r['feasible'] for r in records) else 2


if __name__ == '__main__':
    raise SystemExit(main())
