#!/usr/bin/env python3
"""Run/resume NLR7301 fixed-lift contour and placement sensitivity checks.

This is numerical verification of design screening, not independent validation
of geometry rankings. Original MSES is not run. Outputs include every requested
flow evaluation, exact modified input coordinates and the solver's surfaces.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import read_dat, repanel, validate_geometry, save_results
from viscous import CoupledSolver
from optimization import solve_fixed_lift, camber_bump, minimum_clearance, LiftTargetError
from validate import verify_data, write_csv

GEOMETRIES = {
    'baseline': (0., 0.),
    'camber_minus': (-.001, 0.), 'camber_plus': (.001, 0.),
    'flap_down': (0., -.002), 'flap_up': (0., .002),
    'camber_minus_half': (-.0005, 0.), 'camber_plus_half': (.0005, 0.)}
SETTINGS = {
    'coarse': dict(panels=[180, 120], wake_length=1., ncrit=9, start_alpha=4.),
    'base': dict(panels=[260, 180], wake_length=1., ncrit=9, start_alpha=4.),
    'fine': dict(panels=[340, 240], wake_length=1., ncrit=9, start_alpha=4.),
    'long_wake': dict(panels=[260, 180], wake_length=2., ncrit=9, start_alpha=4.),
    'ncrit6': dict(panels=[260, 180], wake_length=1., ncrit=6, start_alpha=4.),
    'start0': dict(panels=[260, 180], wake_length=1., ncrit=9, start_alpha=0.)}


def specification():
    jobs = [(g, s) for s in ('coarse', 'base', 'fine', 'long_wake')
            for g in list(GEOMETRIES)[:5]]
    jobs += [(g, 'base') for g in list(GEOMETRIES)[5:]]
    jobs += [(g, 'ncrit6') for g in ('baseline', 'camber_minus', 'camber_plus')]
    jobs += [('baseline', 'start0')]
    files = ['multifoil.py', 'viscous.py', 'vendor/mfoil.py', 'optimization.py',
             'validation/optimization_pilot.py',
             'validation/data/nlr7301/NLR7301MainECARP.dat',
             'validation/data/nlr7301/NLR7301FlapECARP.dat']
    return dict(schema=1, target_cl=2.484, reynolds=2510000, mach=.185,
                reference_chord=1., moment_reference=[.25, 0.], alpha_bounds=[4., 7.],
                cl_tolerance=2e-5, iterations=100, residual_tolerance=1e-7,
                minimum_clearance=.015, geometries=GEOMETRIES, settings=SETTINGS,
                jobs=jobs, sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files})


def one_job(out_string, geom, setting, spec):
    out = Path(out_string)/f'{geom}__{setting}'
    out.mkdir(parents=True, exist_ok=True)
    result_path = out/'result.json'
    if result_path.exists():
        return json.loads(result_path.read_text())
    start = time.perf_counter()
    params = dict(spec['settings'][setting])
    bump, dy = spec['geometries'][geom]
    data = ROOT/'validation/data/nlr7301'
    main = camber_bump(read_dat(data/'NLR7301MainECARP.dat'), bump)
    flap = read_dat(data/'NLR7301FlapECARP.dat') + np.array([0., dy])
    raw = [('main', main), ('flap', flap)]
    validate_geometry(raw)
    clearance = minimum_clearance(main, flap)
    for name, xy in raw:
        np.savetxt(out/f'input_{name}.dat', xy, header=name, comments='', fmt='%.16g')
    record = dict(geometry=geom, setting=setting, camber_amplitude=bump, flap_dy=dy,
                  clearance=clearance, feasible=False, alpha_deg=None, CL=None, CD=None, CM=None,
                  lift_error=None, failure_reason='', flow_evaluations=0)
    evaluations = []
    elements = [(name, repanel(xy, n)) for (name, xy), n in zip(raw, params['panels'])]
    cfg = dict(spec, active_geometry=geom, active_setting=setting, solver_settings=params)

    def checkpoint(solution):
        evaluations.append(solution)
        save_results(out/'flow', evaluations, elements, cfg, 'viscous')

    try:
        validate_geometry(elements)
        if clearance < spec['minimum_clearance']:
            raise LiftTargetError('Input clearance below the declared design constraint')
        solver = CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'],
            ncrit=params['ncrit'], iterations=spec['iterations'], tolerance=spec['residual_tolerance'],
            wake_length=params['wake_length'], start_alpha=params['start_alpha'])
        solution, _ = solve_fixed_lift(solver, spec['target_cl'], spec['alpha_bounds'],
                                       spec['cl_tolerance'], on_evaluation=checkpoint)
        solved = [(s['name'], s['xy']) for s in solution.surfaces]
        validate_geometry(solved)
        solved_clearance = minimum_clearance(solved[0][1], solved[1][1])
        if solved_clearance < spec['minimum_clearance']:
            raise LiftTargetError('Repaneled clearance below the declared design constraint')
        record.update(feasible=True, alpha_deg=solution.alpha, CL=solution.cl, CD=solution.cd,
                      CM=solution.cm, lift_error=solution.cl-spec['target_cl'],
                      solved_clearance=solved_clearance,
                      max_separated_arclength_fraction=max(
                          side['separated_arclength_fraction']
                          for bl in solution.diagnostics['boundary_layers'].values() for side in bl.values()))
    except (LiftTargetError, ValueError, FloatingPointError) as exc:
        record['failure_reason'] = f'{type(exc).__name__}: {exc}'
    record['flow_evaluations'] = len(evaluations)
    record['wall_seconds'] = time.perf_counter()-start
    result_path.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    return record


def summarize(out, spec):
    records = []
    for geom, setting in spec['jobs']:
        path = out/f'{geom}__{setting}'/'result.json'
        if path.exists():
            records.append(json.loads(path.read_text()))
    by_key = {(r['geometry'], r['setting']): r for r in records}
    for r in records:
        base = by_key.get(('baseline', r['setting']))
        r['delta_cd_counts'] = ((r['CD']-base['CD'])*1e4
                               if r['feasible'] and base and base['feasible'] else None)
    if records:
        fields = list(dict.fromkeys(k for r in records for k in r))
        write_csv(out/'results.csv', records, fields)
    directions = []
    for setting in SETTINGS:
        for mode, minus, plus, step in [('camber', 'camber_minus', 'camber_plus', .001),
                                      ('flap_y', 'flap_down', 'flap_up', .002),
                                      ('camber_half', 'camber_minus_half', 'camber_plus_half', .0005)]:
            a, b = by_key.get((minus, setting)), by_key.get((plus, setting))
            if a and b:
                directions.append(dict(setting=setting, mode=mode, step=step,
                    dCD_dp=((b['CD']-a['CD'])/(2*step) if a['feasible'] and b['feasible'] else None)))
    if directions:
        write_csv(out/'directional_differences.csv', directions, list(directions[0]))
    comparison = []
    for geom in list(GEOMETRIES)[1:5]:
        values = [by_key.get((geom, s)) for s in ('coarse', 'base', 'fine', 'long_wake')]
        complete = all(r and r['delta_cd_counts'] is not None for r in values)
        delta = [r['delta_cd_counts'] for r in values] if complete else []
        comparison.append(dict(geometry=geom, complete=complete,
            min_delta_cd_counts=min(delta) if delta else None,
            max_delta_cd_counts=max(delta) if delta else None,
            sign_stable=(all(v < 0 for v in delta) or all(v > 0 for v in delta)) if delta else None,
            improves_all_settings=all(v < 0 for v in delta) if delta else None))
    summary = dict(requested=len(spec['jobs']), completed=len(records),
        feasible=sum(r['feasible'] for r in records), geometry_comparison=comparison,
        independent_geometry_ranking_validated=False,
        note='Resolution/wake spread is an observed sensitivity range, not a statistical error bar.')
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'validation/results/optimization_pilot')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('workers must be between 1 and 8')
    verify_data()
    spec = specification()
    # Normalize tuples before comparing with JSON readback.
    spec = json.loads(json.dumps(spec))
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest = out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != spec:
            raise ValueError('Study or code changed: choose a new output directory')
    elif any(out.iterdir()):
        raise ValueError('Nonempty output directory lacks a study manifest')
    else:
        manifest.write_text(json.dumps(spec, indent=2)+'\n')
    if not args.summarize_only:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(one_job, str(out), g, s, spec) for g, s in spec['jobs']]
            for future in as_completed(futures):
                r = future.result()
                print(f"{r['geometry']} / {r['setting']}: feasible={r['feasible']} "
                      f"CL={r['CL']} CD={r['CD']} {r['failure_reason']}", flush=True)
                summarize(out, spec)
    summary = summarize(out, spec)
    print(json.dumps(summary, indent=2))
    return 0 if summary['feasible'] == summary['requested'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
