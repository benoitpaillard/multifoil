#!/usr/bin/env python3
"""Reproducible S702 polar or fixed-lift trials on published vector geometry.

Downward and upward polar branches start independently at the specified center.
Fixed-lift brackets come from adjacent admissible points in a preserved polar.
Failure is archived, never converted to stall or a finite drag prediction.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import repanel, save_results
from viscous import CoupledSolver
from optimization import admissible, solve_fixed_lift, LiftTargetError
from s702_geometry import geometry, FILES

SETTINGS = {'coarse': [180, 120], 'base': [260, 180], 'fine': [340, 240]}
CODE = ['multifoil.py', 'viscous.py', 'paneling.py', 'vendor/mfoil.py',
        'optimization.py', 'validation/s702_geometry.py', 'validation/s702_benchmark.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value):
    return float(value) if value is not None and math.isfinite(value) else None


def record(sol, label):
    return dict(label=label, alpha_deg=sol.alpha, CL=finite(sol.cl), CD=finite(sol.cd),
        CM=finite(sol.cm), converged=bool(sol.converged), admissible=bool(admissible(sol)),
        residual=finite(sol.residual), failure_reason=sol.diagnostics.get('failure_reason', ''),
        local_supersonic=any(sol.diagnostics.get('local_supersonic_detected', {}).values()))


def specification(args):
    inputs = CODE + ['validation/data/s702/'+f for f in ['fore.dat', *FILES.values(), 'source.json']]
    spec = dict(schema_version=1, reference='NASA/CR-20220014964, S702, Figures 1 and 23-26',
        geometry_status='vector-figure reconstruction; not original author coordinates',
        validation_pass_claim=False, mode=args.mode, reynolds=args.reynolds, mach=args.mach,
        ncrit=9., reference_chord=1., moment_reference=[.25, 0.],
        tolerance=1e-7, iterations=100, wake_length=1.,
        center_alpha=args.center, angles=args.angles, targets=args.targets,
        cl_tolerance=2e-5, fixed_lift_budget=20,
        settings={s: SETTINGS[s] for s in args.settings},
        jobs=[dict(deflection=d, setting=s) for s in args.settings for d in args.deflections],
        sha256={f: digest(ROOT/f) for f in inputs})
    if args.mode == 'fixed':
        if args.seed_polar is None:
            raise ValueError('Fixed-lift mode requires --seed-polar')
        seed = args.seed_polar.resolve()
        prior = json.loads((seed/'study.json').read_text())
        keys = ['reynolds', 'mach', 'ncrit', 'reference_chord', 'moment_reference',
                'tolerance', 'iterations', 'wake_length']
        if prior['mode'] != 'polar' or any(prior[k] != spec[k] for k in keys):
            raise ValueError('Seed polar conditions differ')
        for f in inputs:
            if prior['sha256'].get(f) != spec['sha256'][f]:
                raise ValueError('Seed polar source differs: '+f)
        if any(prior['settings'].get(s) != SETTINGS[s] for s in args.settings):
            raise ValueError('Seed panel settings differ')
        # Embed seed records: the manifest does not depend on machine-specific paths.
        spec['seed_study_sha256'] = digest(seed/'study.json')
        spec['seed_results'] = {}
        for j in spec['jobs']:
            name = f"d{j['deflection']:+d}__{j['setting']}"
            path = seed/name/'result.json'
            spec['seed_results'][name] = dict(sha256=digest(path), data=json.loads(path.read_text()))
    return spec


def run_job(root_string, job, spec):
    d, setting = job['deflection'], job['setting']
    name = f'd{d:+d}__{setting}'
    out = Path(root_string)/name
    out.mkdir(parents=True, exist_ok=True)
    result_path = out/'result.json'
    if result_path.exists():
        return json.loads(result_path.read_text())
    started = time.perf_counter()
    result = dict(**job, mode=spec['mode'], records=[], targets=[], error='')
    try:
        raw = geometry(d)
        elements = [(n, repanel(x, count)) for (n, x), count in zip(raw, spec['settings'][setting])]
        for n, x in raw:
            np.savetxt(out/f'input_{n}.dat', x, fmt='%.16g', header=n, comments='')
        def solver():
            return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'],
                ncrit=spec['ncrit'], iterations=spec['iterations'], tolerance=spec['tolerance'],
                wake_length=spec['wake_length'], start_alpha=spec['center_alpha'])
        def collector(label):
            solutions = []
            def checkpoint(sol):
                solutions.append(sol)
                result['records'].append(record(sol, label))
                save_results(out/label, solutions, elements, dict(study=spec, job=job), 'viscous')
                result_path.with_suffix('.partial.json').write_text(json.dumps(result, indent=2)+'\n')
            return checkpoint
        if spec['mode'] == 'polar':
            center = spec['center_alpha']
            branches = [('down', sorted([a for a in spec['angles'] if a <= center], reverse=True)),
                        ('up', sorted([a for a in spec['angles'] if a > center]))]
            for label, angles in branches:
                if not angles:
                    continue
                flow = solver()
                checkpoint = collector(label)
                for alpha in angles:
                    checkpoint(flow.solve(alpha))
        else:
            seed = spec['seed_results'][name]['data']['records']
            # Keep the requested incidence order but require adjacent valid points;
            # never bracket across a failed or inadmissible polar point.
            ordered = sorted(seed, key=lambda r: r['alpha_deg'])
            for target in spec['targets']:
                label = f'cl{target:.4f}'
                entry = dict(target_cl=target, feasible=False, CL=None, CD=None,
                             CM=None, alpha_deg=None, reason='')
                try:
                    brackets = [(a, b) for a, b in zip(ordered, ordered[1:])
                        if a['admissible'] and b['admissible'] and a['CL'] <= target <= b['CL']
                        and a['CL'] < b['CL']]
                    if len(brackets) != 1:
                        raise LiftTargetError('No unique adjacent increasing-lift bracket in seed polar')
                    left, right = brackets[0]
                    bounds = [left['alpha_deg'], right['alpha_deg']]
                    entry['bounds'] = bounds
                    answer, _ = solve_fixed_lift(solver(), target, bounds, spec['cl_tolerance'],
                        max_evaluations=spec['fixed_lift_budget'], on_evaluation=collector(label))
                    entry.update(feasible=True, CL=answer.cl, CD=answer.cd,
                                 CM=answer.cm, alpha_deg=answer.alpha)
                except (ValueError, RuntimeError, FloatingPointError) as exc:
                    entry['reason'] = f'{type(exc).__name__}: {exc}'
                result['targets'].append(entry)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        result['error'] = f'{type(exc).__name__}: {exc}'
    result.update(flow_evaluations=len(result['records']),
                  admissible_evaluations=sum(r['admissible'] for r in result['records']),
                  wall_seconds=time.perf_counter()-started)
    result['complete_success'] = (not result['error'] and bool(result['records']) and
        (all(t['feasible'] for t in result['targets']) if spec['mode'] == 'fixed' else
         all(r['admissible'] for r in result['records'])))
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    result_path.with_suffix('.partial.json').unlink(missing_ok=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=['polar', 'fixed'], default='polar')
    p.add_argument('--settings', default='base')
    p.add_argument('--deflections', default='-2,0,5,10')
    p.add_argument('--angles', default='-4,-3,-2,-1,0,1,2,3,4,6,8,10,12')
    p.add_argument('--targets', default='.3,.6,1.0')
    p.add_argument('--center', type=float, default=-1.)
    p.add_argument('--reynolds', type=float, default=1e6)
    p.add_argument('--mach', type=float, default=.1)
    p.add_argument('--seed-polar', type=Path)
    p.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_polar')
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    try:
        args.settings = args.settings.split(',')
        args.deflections = [int(v) for v in args.deflections.split(',')]
        args.angles = [float(v) for v in args.angles.split(',')]
        args.targets = [float(v) for v in args.targets.split(',')]
        if not 1 <= args.workers <= 4 or not set(args.settings) <= SETTINGS.keys() or not set(args.deflections) <= FILES.keys():
            raise ValueError('Use 1..4 workers, known settings and published deflections')
        for values in [args.settings, args.deflections, args.angles, args.targets]:
            if not values or len(values) != len(set(values)):
                raise ValueError('Selections must be nonempty and contain no duplicates')
        if not np.isfinite(args.angles+args.targets+[args.center, args.reynolds, args.mach]).all():
            raise ValueError('Inputs must be finite')
        if args.reynolds <= 0 or not 0 <= args.mach < .7:
            raise ValueError('Require Re > 0 and 0 <= Mach < .7')
        spec = specification(args)
    except (ValueError, OSError, KeyError) as exc:
        p.error(str(exc))
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out/'study.json'
    if manifest.exists() and json.loads(manifest.read_text()) != spec:
        p.error('Inputs or source hashes changed; select a new output directory')
    manifest.write_text(json.dumps(spec, indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run_job, str(args.out), j, spec) for j in spec['jobs']]):
            r = future.result()
            print(json.dumps({k: r[k] for k in ['deflection', 'setting', 'flow_evaluations',
                'admissible_evaluations', 'complete_success', 'wall_seconds', 'error']}), flush=True)
    records = [json.loads((args.out/f"d{j['deflection']:+d}__{j['setting']}"/'result.json').read_text())
               for j in spec['jobs']]
    summary = dict(jobs=len(records), successful_jobs=sum(r['complete_success'] for r in records),
        requested_evaluations=sum(r['flow_evaluations'] for r in records),
        admissible_evaluations=sum(r['admissible_evaluations'] for r in records),
        fixed_lift_targets=sum(len(r['targets']) for r in records),
        successful_targets=sum(t['feasible'] for r in records for t in r['targets']),
        validation_pass_claim=False)
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return 0 if all(r['complete_success'] for r in records) else 2


if __name__ == '__main__':
    raise SystemExit(main())
