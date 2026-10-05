#!/usr/bin/env python3
"""Controlled S702 Ncrit sweep; measured drag never enters the flow solver.

One solver per (Ncrit, flap) follows increasing lift targets. Existing N=6/9
results can be reused at exactly matching conditions. Their different incidence
paths remain explicit. All failed evaluations and bracket retries are saved.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import repanel, save_results
from optimization import solve_fixed_lift, admissible
from viscous import CoupledSolver
from s702_geometry import geometry, FILES
from s702_benchmark import record
from s702_transition_study import transition_audit


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tag(value):
    return f'{value:g}'.replace('.', 'p')


def seed_alpha(flap, cl, ncrit):
    """Interpolate numerical incidences only; no experimental force fitting."""
    nine = read(ROOT/f'validation/results/s702_fixed/d{flap:+d}__base/result.json')['targets']
    xs = np.array([r['target_cl'] for r in nine])
    a9 = np.array([r['alpha_deg'] for r in nine])
    a6 = np.array([read(ROOT/f'validation/results/s702_transition_diagnostic/d{flap:+d}_cl{c:g}_both/result.json')['alpha_deg'] for c in xs])
    center = np.interp(cl, xs, a9+(ncrit-9.)*(a9-a6)/3.)
    return float(center)


def local_seed(result):
    """Estimate incidence from two admissible model lift values, never drag."""
    points = {}
    for r in result['records']:
        if r['admissible'] and np.isfinite([r['alpha_deg'], r['CL']]).all():
            points[r['alpha_deg']] = r['CL']
    nearest = sorted(points, key=lambda a: abs(points[a]-result['target_cl']))[:2]
    if len(nearest) < 2:
        raise ValueError('Local retry requires two distinct admissible incidences')
    a, b = sorted(nearest)
    slope = (points[b]-points[a])/(b-a)
    if not .02 <= slope <= .3:
        raise ValueError('Computed lift slope is unsuitable for a local subcritical retry')
    center = a+(result['target_cl']-points[a])/slope
    if not a-1. <= center <= b+1.:
        raise ValueError('Local retry would extrapolate incidence by more than one degree')
    return float(center)


def legacy(flap, cl, ncrit, spec):
    if (not spec['reuse_legacy'] or spec['reynolds'] != 1e6 or
        spec['mach'] != .1 or spec['panels'] != [260, 180] or
        cl not in (.5, .8, 1.2) or ncrit not in (6., 9.)):
        return None
    if ncrit == 9.:
        path = ROOT/f'validation/results/s702_fixed/d{flap:+d}__base/result.json'
        raw = read(path)
        result = next(r for r in raw['targets'] if r['target_cl'] == cl).copy()
        result['records'] = [r for r in raw['records'] if r['label'] == f'cl{cl:.4f}']
    else:
        path = ROOT/f'validation/results/s702_transition_diagnostic/d{flap:+d}_cl{cl:g}_both/result.json'
        result = read(path)
    if not result['feasible'] or abs(result['CL']-cl) > spec['cl_tolerance']:
        raise ValueError('Legacy target is not admissible')
    if not all(r['admissible'] and r['residual'] <= spec['residual_tolerance'] for r in result['records']):
        raise ValueError('Legacy flow does not meet the requested convergence criteria')
    result.update(flap=flap, ncrit=ncrit, reused=True,
        source_path=str(path.relative_to(ROOT)), source_sha256=sha(path),
        flow_directory=str((path.parent/(f'cl{cl:.4f}' if ncrit == 9. else 'flow')).relative_to(ROOT)))
    return result


def run_group(out, job, spec):
    folder = Path(out)/job['name']
    folder.mkdir(exist_ok=True)
    if (folder/'result.json').exists():
        return read(folder/'result.json')
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        return _run_group(folder, job, spec)


def _run_group(folder, job, spec):
    start = time.perf_counter()
    elements = [(name, repanel(xy, count)) for (name, xy), count in
                zip(geometry(job['flap']), spec['panels'])]

    def fresh(start_alpha=-1.):
        return CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'],
            ncrit=job['ncrit'], iterations=100, tolerance=spec['residual_tolerance'],
            wake_length=1., start_alpha=start_alpha)

    solver = None
    targets = []
    for cl in spec['targets']:
        target_dir = folder/f'cl{tag(cl)}'
        target_dir.mkdir(exist_ok=True)
        result_path = target_dir/'result.json'
        if result_path.exists():
            targets.append(read(result_path))
            solver = None  # A resumed run does not pretend to restore a flow state.
            continue
        result = legacy(job['flap'], cl, job['ncrit'], spec)
        if result is not None:
            write(result_path, result)
            targets.append(result)
            continue
        result = dict(flap=job['flap'], ncrit=job['ncrit'], target_cl=cl,
            reused=False, feasible=False, CL=None, CD=None, CM=None, alpha_deg=None,
            error='', records=[], attempts=[],
            flow_directory=str((target_dir/'flow').relative_to(ROOT)))
        solutions = []
        selected_states = {}
        begun = time.perf_counter()
        phase = 'fixed_lift'

        def checkpoint(sol):
            sol.diagnostics['ncrit_by_element'] = dict(fore=job['ncrit'], aft=job['ncrit'])
            solutions.append(sol)
            result['records'].append(record(sol, phase))
            if admissible(sol) and abs(sol.cl-cl) <= spec['cl_tolerance']:
                selected_states[sol.alpha] = (copy.deepcopy(solver.models), copy.deepcopy(sol))
            save_results(target_dir/'flow', solutions, elements,
                         dict(study=spec, job=job, target_cl=cl), 'viscous')
            write(target_dir/'progress.json', result)

        center = job['seed_alphas'][str(cl)]
        widths = [spec.get('initial_half_width', .3), .75, .15]
        if solver is None:
            solver = fresh(center) if spec.get('start_at_seed', False) else fresh()
        for attempt, width in enumerate(widths):
            phase = 'fixed_lift'
            if attempt and 'not bracketed' in result['error']:
                width = max(.75, result['attempts'][-1]['bounds'][1]-result['attempts'][-1]['bounds'][0])
            if attempt and 'not bracketed' not in result['error']:
                # Cold restart near the already predicted target, avoiding an
                # inadmissible remote bracket endpoint. Keep the failed attempt.
                width = .15 if attempt == 1 else .05
                solver = fresh(center)
            bounds = [center-width, center+width]
            trial = dict(bounds=bounds, first_record=len(solutions),
                         initialized=solver.initialized, start_alpha=solver.start_alpha)
            result['attempts'].append(trial)
            try:
                if attempt == 0 and spec.get('warm_seed_alphas'):
                    phase = 'warm_seed'
                    for alpha in spec['warm_seed_alphas']:
                        seed_solution = solver.solve(alpha)
                        checkpoint(seed_solution)
                        if not admissible(seed_solution):
                            raise RuntimeError('Recorded initialization path could not be replayed')
                    phase = 'fixed_lift'
                sol, _ = solve_fixed_lift(solver, cl, bounds,
                    spec['cl_tolerance'], 20, checkpoint)
                # The fixed-lift helper may return the first endpoint after it
                # also evaluates the second. Restore the *selected* exact state.
                models, selected = selected_states[sol.alpha]
                solver.models, solver.last_solution = models, selected
                solver.last_alpha, solver.initialized = sol.alpha, True
                result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm,
                    alpha_deg=sol.alpha, error='',
                    selected_flow_index=next(i for i in range(len(solutions)-1, -1, -1)
                                             if solutions[i].alpha == sol.alpha),
                    boundary_layers=sol.diagnostics['boundary_layers'],
                    element_forces=sol.diagnostics['element_forces'],
                    transition_audit=transition_audit(solver))
                for model in solver.models:
                    np.savez_compressed(target_dir/(model.geom.name+'_state.npz'),
                        U=model.glob.U, xy=model.foil.x, xi=model.isol.xi,
                        turbulent=model.vsol.turb, wake_xy=model.wake.x,
                        lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]))
                trial.update(feasible=True, last_record=len(solutions)-1)
                break
            except (ValueError, RuntimeError, FloatingPointError, AssertionError) as exc:
                result['error'] = f'{type(exc).__name__}: {exc}'
                trial.update(feasible=False, error=result['error'], last_record=len(solutions)-1)
                write(target_dir/'progress.json', result)
        if not result['feasible']:
            solver = None
        result['wall_seconds'] = time.perf_counter()-begun
        write(result_path, result)
        targets.append(result)
        print({k: result[k] for k in ('ncrit', 'flap', 'target_cl', 'feasible', 'CD', 'error')}, flush=True)
    group = dict(**job, targets=targets, complete_success=all(r['feasible'] for r in targets),
                 wall_seconds=time.perf_counter()-start)
    write(folder/'result.json', group)
    return group


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_ncrit_coarse')
    parser.add_argument('--ncrit', default='3,4,5,6,7,8,9,10,11,12')
    parser.add_argument('--flaps', default='-2,0,5,10')
    parser.add_argument('--targets', default='.5,.8,1.2')
    parser.add_argument('--panels', default='260,180')
    parser.add_argument('--reynolds', type=float, default=1e6)
    parser.add_argument('--mach', type=float, default=.1)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--no-reuse', action='store_true')
    parser.add_argument('--seed-from', type=Path,
                        help='Separate retry using model lift records from one failed target')
    parser.add_argument('--warm-seed', action='store_true',
                        help='Replay up to two initial admissible source evaluations before the local bracket')
    args = parser.parse_args()
    args.out = args.out.resolve()
    ns = sorted(set(float(v) for v in args.ncrit.split(',')))
    flaps = sorted(set(int(v) for v in args.flaps.split(',')))
    targets = sorted(set(float(v) for v in args.targets.split(',')))
    panels = [int(v) for v in args.panels.split(',')]
    if args.warm_seed and not args.seed_from:
        parser.error('--warm-seed requires --seed-from')
    if (not ns or not np.isfinite(ns+targets).all() or min(ns) <= 0 or max(ns) > 20 or
        not set(flaps) <= FILES.keys() or not targets or min(targets) < .5 or max(targets) > 1.2 or
        len(panels) != 2 or min(panels) < 50 or sum(panels) > 1200 or
        not 1 <= args.workers <= 8 or not 0 <= args.mach < .7 or not 1e5 <= args.reynolds <= 1e7):
        parser.error('Invalid conditions or sweep selections')
    for name in ('s702_fixed', 's702_transition_diagnostic'):
        for path, expected in read(ROOT/f'validation/results/{name}/study.json')['sha256'].items():
            if sha(ROOT/path) != expected:
                raise ValueError('Archived seed/reuse source changed: '+path)
    jobs = [dict(name=f'n{tag(n)}_d{f:+d}', ncrit=n, flap=f,
                 seed_alphas={str(cl): seed_alpha(f, cl, n) for cl in targets})
            for n in ns for f in flaps]
    retry_source = None
    if args.seed_from:
        retry_source = args.seed_from.resolve()
        if not retry_source.is_relative_to(ROOT) or len(jobs) != 1 or len(targets) != 1:
            parser.error('A local retry requires one Ncrit, flap and lift target inside this project')
        prior = read(retry_source)
        prior_spec = read(retry_source.parents[2]/'study.json')
        if prior['feasible'] or (prior['ncrit'], prior['flap'], prior['target_cl']) != (ns[0],flaps[0],targets[0]):
            parser.error('Local retry must refer to the identical failed target')
        if (prior_spec['reynolds'],prior_spec['mach'],prior_spec['panels']) != (args.reynolds,args.mach,panels):
            parser.error('Local retry cannot change physical conditions or panel counts')
        jobs[0]['seed_alphas'][str(targets[0])] = local_seed(prior)
    inputs = ['multifoil.py', 'viscous.py', 'paneling.py', 'optimization.py', 'vendor/mfoil.py',
        'validation/s702_geometry.py', 'validation/s702_benchmark.py',
        'validation/s702_transition_study.py', 'validation/s702_ncrit_sweep.py',
        'validation/data/s702/source.json', 'validation/data/s702/fore.dat']
    inputs += ['validation/data/s702/'+FILES[f] for f in flaps]
    inputs += [str(p.relative_to(ROOT)) for name in ('s702_fixed', 's702_transition_diagnostic')
               for p in sorted((ROOT/'validation/results'/name).glob('**/result.json'))]
    if retry_source:
        inputs += [str(retry_source.relative_to(ROOT)), str((retry_source.parents[2]/'study.json').relative_to(ROOT))]
    spec = dict(schema=1, ncrit_values=ns, flaps=flaps, targets=targets, jobs=jobs,
        reynolds=args.reynolds, mach=args.mach, panels=panels, wake_length=1.,
        cl_tolerance=2e-5, residual_tolerance=1e-7, iterations=100,
        reuse_legacy=not args.no_reuse,
        protocol='Common Ncrit, ascending CL continuation; local brackets with recorded retries',
        selection_protocol='Equal-point MAPE; complete coverage required; leave-one-flap-out selection',
        qualification='Calibration and transfer study; no contour-gradient validation claim',
        sha256={p: sha(ROOT/p) for p in inputs})
    if retry_source:
        spec.update(local_retry_of=str(retry_source.relative_to(ROOT)), start_at_seed=True,
                    initial_half_width=.05, reuse_legacy=False,
                    protocol='Separate numerical retry: cold start at incidence extrapolated from computed lift; original failed target retained')
        if args.warm_seed:
            prefix = []
            for r in prior['records'][:2]:
                if not r['admissible']:
                    break
                prefix.append(r['alpha_deg'])
            if not prefix:
                raise ValueError('No admissible initial source path to replay')
            spec.update(start_at_seed=False, warm_seed_alphas=prefix,
                        protocol='Separate numerical retry: replay original admissible initialization, then use a local bracket from computed lift')
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out/'study.json'
    if manifest.exists():
        if read(manifest) != spec:
            raise ValueError('Study or source changed; select a new output directory')
    elif any(args.out.iterdir()):
        raise ValueError('Nonempty output directory has no manifest')
    else:
        write(manifest, spec)
    for relative, expected in spec['sha256'].items():
        if relative.endswith('.py'):
            snapshot = args.out/'source_snapshot'/relative
            snapshot.parent.mkdir(parents=True, exist_ok=True)
            if not snapshot.exists():
                snapshot.write_bytes((ROOT/relative).read_bytes())
            if sha(snapshot) != expected:
                raise ValueError('Execution source snapshot differs from manifest')
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(run_group, args.out, job, spec) for job in jobs]
        for future in as_completed(futures):
            group = future.result()
            print(dict(name=group['name'], success=group['complete_success'],
                       wall_seconds=group['wall_seconds']), flush=True)
    groups = [read(args.out/job['name']/'result.json') for job in jobs]
    write(args.out/'summary.json', groups)
    return 0 if all(g['complete_success'] for g in groups) else 2


if __name__ == '__main__':
    raise SystemExit(main())
