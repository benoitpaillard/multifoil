#!/usr/bin/env python3
"""Constrained, derivative-free incompressible contour optimization.

Coordinates must already share one reference chord. A C2 camber deformation
preserves thickness at fixed chordwise position and fixes the ends. Each shape
and operating point gets a fresh viscous solver. No data-fitted drag correction
or shape adjoint is used. See validation/INCOMPRESSIBLE_DESIGN.md.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from multifoil import read_dat, repanel, validate_geometry, save_results, signed_area
from optimization import admissible, solve_fixed_lift, LiftTargetError, minimum_clearance
from fixed_lift_policy import InitializationPath, solve_with_initialization
from viscous import CoupledSolver

ROOT = Path(__file__).resolve().parent
CODE = ['shape_optimization.py', 'fixed_lift_policy.py', 'optimization.py', 'multifoil.py', 'viscous.py',
        'amplification_jacobian.py', 'paneling.py', 'vendor/mfoil.py']


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def camber_mode(x, start, end):
    """Unit-peak compact polynomial; value, slope and curvature vanish at ends."""
    if not 0 <= start < end <= 1:
        raise ValueError('Camber support must satisfy 0 <= start < end <= 1')
    t = np.clip((np.asarray(x)-start)/(end-start), 0., 1.)
    return 64*t**3*(1-t)**3


class ShapeFamily:
    """Apply all modes in frames frozen from the undeformed coordinates.

    Frames are used only for deformation. Coordinates, element incidence,
    assembly chord and coefficient normalization are never rescaled.
    """
    def __init__(self, elements, modes):
        validate_geometry(elements)
        self.elements = [(name, np.array(xy, copy=True)) for name, xy in elements]
        self.modes = modes
        self.frames, self.basis = {}, []
        names = [n for n, _ in elements]
        if len(set(names)) != len(names):
            raise ValueError('Element names must be unique')
        for name, xy in self.elements:
            te = (xy[0]+xy[-1])/2
            le = xy[np.argmax(np.linalg.norm(xy-te, axis=1))]
            chord = np.linalg.norm(te-le)
            axis = (te-le)/chord
            self.frames[name] = (le.copy(), axis, np.array([-axis[1], axis[0]]), chord)
        for m in modes:
            if m['element'] not in names:
                raise ValueError('Mode element does not exist')
            xy = dict(self.elements)[m['element']]
            origin, axis, normal, chord = self.frames[m['element']]
            x = (xy-origin)@axis/chord
            self.basis.append(chord*camber_mode(x, *m['support'])[:, None]*normal)
        self.areas = {n: abs(signed_area(x)) for n, x in self.elements}

    def deform(self, parameters):
        p = np.asarray(parameters, dtype=float)
        if p.shape != (len(self.modes),) or not np.isfinite(p).all():
            raise ValueError('One finite amplitude is required per shape mode')
        changed = {n: x.copy() for n, x in self.elements}
        for value, mode, basis in zip(p, self.modes, self.basis):
            if not mode['bounds'][0] <= value <= mode['bounds'][1]:
                raise ValueError('Shape parameter outside declared bounds')
            changed[mode['element']] += value*basis
        result = [(n, changed[n]) for n, _ in self.elements]
        validate_geometry(result)
        return result


def geometry_metrics(elements, baseline_areas):
    validate_geometry(elements)
    distances = [minimum_clearance(a, b) for i, (_, a) in enumerate(elements)
                 for _, b in elements[i+1:]]
    return dict(minimum_clearance=min(distances) if distances else None,
                relative_area_change={n: abs(signed_area(x))/baseline_areas[n]-1
                                      for n, x in elements})


def initialization_settings(point):
    """Optional declared local-search policy; no state or seed shared across shapes."""
    if 'initialization' not in point:
        return None
    policy = point['initialization']
    allowed = {'center', 'half_width', 'max_evaluations', 'max_brackets', 'paths'}
    if not isinstance(policy, dict) or set(policy)-allowed or not {'center', 'paths'} <= set(policy):
        raise ValueError('Initialization requires center and paths, with no unknown options')
    options = dict(center=policy['center'], half_width=policy.get('half_width', .05),
                   max_evaluations=policy.get('max_evaluations', 20),
                   max_brackets=policy.get('max_brackets', 2))
    lo, hi = point['alpha_bounds']
    if (not np.isfinite([options['center'], options['half_width']]).all() or
            not lo <= options['center'] <= hi or options['half_width'] <= 0 or
            type(options['max_evaluations']) is not int or options['max_evaluations'] < 3 or
            type(options['max_brackets']) is not int or not 1 <= options['max_brackets'] <= 2 or
            not isinstance(policy['paths'], list) or not policy['paths']):
        raise ValueError('Invalid initialization bounds, paths or budgets')
    paths = []
    for entry in policy['paths']:
        if (not isinstance(entry, dict) or set(entry)-{'name', 'start_alpha', 'seed_alphas'} or
                not {'name', 'start_alpha'} <= set(entry) or
                not isinstance(entry['name'], str) or not re.fullmatch(r'[A-Za-z0-9_-]+', entry['name']) or
                not isinstance(entry.get('seed_alphas', []), list)):
            raise ValueError('Invalid initialization path')
        path = InitializationPath(entry['name'], entry['start_alpha'], tuple(entry.get('seed_alphas', [])))
        if not all(np.isfinite(a) and lo <= a <= hi for a in (path.start_alpha, *path.seed_alphas)):
            raise ValueError('Initialization start or seed exceeds alpha_bounds')
        paths.append(path)
    if len({p.name for p in paths}) != len(paths):
        raise ValueError('Initialization path names must be unique')
    return dict(paths=paths, **options)


def check_config(cfg):
    if cfg.get('amplification_jacobian', 'corrected_af') not in ('corrected_af', 'legacy'):
        raise ValueError('Unknown amplification Jacobian')
    if cfg.get('mach') != 0:
        raise ValueError('This design workflow requires an explicit mach: 0')
    if cfg.get('reference_chord') != 1:
        raise ValueError('Coordinates must use a shared reference chord of 1')
    if not np.isfinite(cfg['reynolds']) or cfg['reynolds'] <= 0:
        raise ValueError('Reynolds number must be positive')
    if not cfg['modes'] or not cfg['points'] or not cfg['scenarios']:
        raise ValueError('Modes, lift points and transition scenarios are required')
    for category in ('elements', 'points', 'scenarios', 'modes'):
        names = [x['name'] for x in cfg[category]]
        if len(set(names)) != len(names) or any(not re.fullmatch(r'[A-Za-z0-9_-]+', n) for n in names):
            raise ValueError('Names must be unique safe filenames within each category')
    for m in cfg['modes']:
        lo, hi = m['bounds']
        if not np.isfinite([lo, hi, m['step']]).all() or not lo < 0 < hi or m['step'] <= 0:
            raise ValueError('Mode bounds must contain zero; step must be positive')
        camber_mode(np.array([0.]), *m['support'])
    for p in cfg['points']:
        lo, hi = p['alpha_bounds']
        if not np.isfinite([lo, hi, p['cl'], p['weight'], p['start_alpha']]).all() or lo >= hi or p['weight'] <= 0:
            raise ValueError('Invalid fixed-lift point')
        initialization_settings(p)
    for s in cfg['scenarios']:
        if not np.isfinite(s['ncrit']) or s['ncrit'] <= 0:
            raise ValueError('Ncrit must be positive')
    for key in ('panels', 'verification_panels'):
        if len(cfg[key]) != len(cfg['elements']) or any(type(n) is not int or not 40 <= n <= 1000 for n in cfg[key]):
            raise ValueError('One valid panel count per element is required')
        if sum(n+1 for n in cfg[key]) > 1200:
            raise ValueError('Panel total exceeds the viscous solver limit')
    if any(f <= b for f, b in zip(cfg['verification_panels'], cfg['panels'])):
        raise ValueError('Verification must use more panels on every element')
    if type(cfg['search']['max_designs']) is not int or cfg['search']['max_designs'] < 1+2*len(cfg['modes']):
        raise ValueError('Budget must cover the initial complete coordinate poll')
    if not 0 < cfg['search']['minimum_step_fraction'] <= 1:
        raise ValueError('minimum_step_fraction must be in (0,1]')
    for key in ('acceptance_counts', 'max_point_drag_increase_counts', 'max_moment_change',
                'max_separated_fraction', 'minimum_clearance', 'max_relative_area_change'):
        value = cfg['constraints'][key]
        if not np.isfinite(value) or value < 0:
            raise ValueError('Constraint values must be finite and nonnegative')
    if not 0 < cfg['constraints']['max_separated_fraction'] < 1:
        raise ValueError('Separation fraction must be between zero and one')
    if cfg['constraints']['acceptance_counts'] <= 0:
        raise ValueError('Use a positive minimum drag improvement')
    if not np.isfinite(cfg['moment_reference']).all() or len(cfg['moment_reference']) != 2:
        raise ValueError('moment_reference needs two finite coordinates')


def flow_job(payload):
    """Fresh flow state for each geometry, transition scenario and lift point."""
    folder, raw, cfg, panels, scenario, point = payload
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/'result.json'
    if path.exists():
        return json.loads(path.read_text())
    raw = [(n, np.asarray(x)) for n, x in raw]
    result = dict(scenario=scenario['name'], point=point['name'], ncrit=scenario['ncrit'],
                  target_cl=point['cl'], feasible=False, alpha_deg=None, CL=None, CD=None,
                  CM=None, failure_reason='', flow_evaluations=0, records=[])
    solutions = []
    started = time.perf_counter()
    try:
        elements = [(n, repanel(x, p)) for (n, x), p in zip(raw, panels)]
        validate_geometry(elements)
        if geometry_metrics(elements, {n: abs(signed_area(x)) for n, x in raw})['minimum_clearance'] is not None:
            if min(minimum_clearance(a, b) for i, (_, a) in enumerate(elements) for _, b in elements[i+1:]) < cfg['constraints']['minimum_clearance']:
                raise LiftTargetError('Repaneled gap violates clearance constraint')
        def factory(start_alpha):
            return CoupledSolver(elements, reynolds=cfg['reynolds'], mach=0.,
                amplification_jacobian=cfg.get('amplification_jacobian', 'corrected_af'),
                ncrit=scenario['ncrit'], iterations=100, tolerance=1e-7, wake_length=1.,
                moment_reference=cfg['moment_reference'], start_alpha=start_alpha)

        def checkpoint(sol, context=None):
            solutions.append(sol)
            values = {k: getattr(sol, attr) for k, attr in
                      (('alpha_deg', 'alpha'), ('CL', 'cl'), ('CD', 'cd'), ('CM', 'cm'), ('residual', 'residual'))}
            result['records'].append(dict(
                **{k: float(v) if v is not None and np.isfinite(v) else None for k, v in values.items()},
                admissible=bool(admissible(sol)), initialization=context or dict(
                    index=len(solutions)-1, path='legacy', bracket=0, phase='fixed_lift')))
            save_results(folder/'flow', solutions, elements,
                         dict(config=cfg, panels=panels, scenario=scenario, point=point,
                              initialization_history=[r['initialization'] for r in result['records']]), 'viscous')
            write_json(folder/'progress.json', result)

        policy = initialization_settings(point)
        if policy is None:
            solver = factory(point['start_alpha'])
            sol, _ = solve_fixed_lift(solver, point['cl'], point['alpha_bounds'],
                                     cl_tolerance=2e-5, max_evaluations=20, on_evaluation=checkpoint)
            result.update(selected_path='legacy', selected_flow_index=next(i for i, s in enumerate(solutions) if s is sol))
        else:
            outcome = solve_with_initialization(factory, point['cl'], alpha_bounds=point['alpha_bounds'],
                cl_tolerance=2e-5, residual_tolerance=1e-7, on_evaluation=checkpoint, **policy)
            result['initialization_attempts'] = outcome.attempts
            if not outcome.feasible:
                raise LiftTargetError(outcome.error)
            sol, solver = outcome.solution, outcome.solver
            result.update(selected_path=outcome.selected_path, selected_flow_index=outcome.selected_flow_index)
            # Only this path guarantees restoration of the complete selected solver.
            for model in solver.models:
                np.savez_compressed(folder/(model.geom.name+'_state.npz'), U=model.glob.U,
                    xy=model.foil.x, xi=model.isol.xi, turbulent=model.vsol.turb, wake_xy=model.wake.x,
                    lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]))
        solved = [(s['name'], s['xy']) for s in sol.surfaces]
        measures = geometry_metrics(solved, {n: abs(signed_area(x)) for n, x in raw})
        if measures['minimum_clearance'] is not None and measures['minimum_clearance'] < cfg['constraints']['minimum_clearance']:
            raise LiftTargetError('Actual solver surface violates clearance constraint')
        if max(abs(v) for v in measures['relative_area_change'].values()) > cfg['constraints']['max_relative_area_change']:
            raise LiftTargetError('Solver surface area drift exceeds geometry tolerance')
        separation = max(s['separated_arclength_fraction']
                         for bl in sol.diagnostics['boundary_layers'].values() for s in bl.values())
        result.update(CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha,
                      residual=sol.residual, max_separated_fraction=separation,
                      actual_geometry=measures, boundary_layers=sol.diagnostics['boundary_layers'],
                      element_forces=sol.diagnostics['element_forces'])
        if separation > cfg['constraints']['max_separated_fraction']:
            raise LiftTargetError('Separated surface fraction exceeds design restriction')
        result['feasible'] = True
    except (LiftTargetError, ValueError, FloatingPointError) as exc:
        result['failure_reason'] = f'{type(exc).__name__}: {exc}'
    result.update(flow_evaluations=len(solutions),
                  admissible_evaluations=sum(admissible(s) for s in solutions),
                  wall_seconds=time.perf_counter()-started)
    write_json(path, result)
    return result


def paired_score(candidate, baseline, cfg):
    """Worst scenario's weighted drag CHANGE; every point keeps its own CL.

    Paired differences prevent a high-drag transition scenario from hiding a
    deterioration in another one. Infeasible points never become an objective.
    """
    expected = {(s['name'], p['name']) for s in cfg['scenarios'] for p in cfg['points']}
    current = {(r['scenario'], r['point']): r for r in candidate}
    base = {(r['scenario'], r['point']): r for r in baseline}
    if (len(current) != len(candidate) or len(base) != len(baseline) or
            set(current) != expected or set(base) != expected):
        return dict(feasible=False, score=None, reason='Missing or duplicate operating points')
    if not all(r['feasible'] and all(r[k] is not None and math.isfinite(r[k]) for k in ('CL', 'CD', 'CM'))
               for r in candidate+baseline):
        return dict(feasible=False, score=None, reason='An operating point is infeasible')
    weights = {p['name']: p['weight'] for p in cfg['points']}
    limits = cfg['constraints']
    changes, means = [], {}
    for s in cfg['scenarios']:
        delta = 0.
        for point in cfg['points']:
            key = s['name'], point['name']
            a, b = current[key], base[key]
            if abs(a['CL']-point['cl']) > 2e-5 or abs(b['CL']-point['cl']) > 2e-5:
                return dict(feasible=False, score=None, reason='Lift is not matched')
            dc = a['CD']-b['CD']
            changes.append(dict(scenario=key[0], point=key[1], delta_cd_counts=dc*1e4,
                                delta_cm=a['CM']-b['CM']))
            if dc*1e4 > limits['max_point_drag_increase_counts']:
                return dict(feasible=False, score=None, reason='Off-design drag restriction exceeded', changes=changes)
            if abs(a['CM']-b['CM']) > limits['max_moment_change']:
                return dict(feasible=False, score=None, reason='Moment change restriction exceeded', changes=changes)
            delta += weights[key[1]]*dc
        means[s['name']] = delta/sum(weights.values())
    return dict(feasible=True, score=max(means.values()), reason='',
                scenario_delta_cd_counts={s: v*1e4 for s, v in means.items()}, changes=changes)


def pattern_search(evaluate, modes, max_designs, min_fraction, acceptance):
    """Bounded coordinate polling with cached points and an explicit budget.

    A complete poll is needed before claiming step convergence. Budget exhaustion
    is reported separately, including when it cuts a poll short.
    """
    initial = tuple(0. for _ in modes)
    cache = {}

    def get(p):
        p = tuple(round(float(x), 14) for x in p)
        if p not in cache:
            cache[p] = evaluate(p, len(cache))
        return cache[p]

    best = get(initial)
    if not best['feasible']:
        return dict(status='baseline_infeasible', parameters=list(initial), best=best,
                    history=list(cache.values()), complete_poll=False)
    p = initial
    steps = np.array([m['step'] for m in modes])
    fraction, complete = 1., True
    status = 'step_limit'
    while fraction >= min_fraction-1e-12:
        center = p
        polled = []
        complete = True
        for j, mode in enumerate(modes):
            for sign in (-1., 1.):
                q = list(center)
                q[j] = float(np.clip(q[j]+sign*steps[j]*fraction, *mode['bounds']))
                key = tuple(round(x, 14) for x in q)
                if key == center:
                    continue
                if key not in cache and len(cache) >= max_designs:
                    complete = False
                    break
                r = get(key)
                if r['feasible']:
                    polled.append((key, r))
            if not complete:
                break
        eligible = [(q, r) for q, r in polled if r['score'] < best['score']-acceptance]
        if eligible:
            p, best = min(eligible, key=lambda item: item[1]['score'])
        elif complete:
            fraction *= .5
        if not complete:
            status = 'budget_exhausted'
            break
        if len(cache) >= max_designs and fraction >= min_fraction-1e-12:
            status = 'budget_exhausted'
            break
    return dict(status=status, parameters=list(p), best=best, history=list(cache.values()),
                complete_poll=complete, final_step_fraction=fraction)


def run(config_path, out, workers):
    cfg = json.loads(config_path.read_text())
    cfg.setdefault('amplification_jacobian', 'corrected_af')
    check_config(cfg)
    inputs = [(e['name'], (config_path.parent/e['file']).resolve()) for e in cfg['elements']]
    family = ShapeFamily([(n, read_dat(p)) for n, p in inputs], cfg['modes'])
    spec = dict(schema=1, config=cfg, code_sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in CODE},
                input_sha256={n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in inputs},
                derivative_method='none; bounded coordinate polling', mach=0.,
                validation_claim=False)
    out.mkdir(parents=True, exist_ok=True)
    path = out/'study.json'
    if path.exists():
        if json.loads(path.read_text()) != spec:
            raise ValueError('Code, input or configuration changed; use a new output directory')
    elif any(out.iterdir()):
        raise ValueError('Nonempty output directory has no study manifest')
    else:
        write_json(path, spec)
    baseline = []

    def evaluate(p, number, panels=None, label=None):
        folder = out/(label or f'design_{number:03d}')
        folder.mkdir(parents=True, exist_ok=True)
        result_path = folder/'result.json'
        if result_path.exists():
            return json.loads(result_path.read_text())
        r = dict(parameters=list(p), folder=folder.name, feasible=False, score=None, reason='', flows=[])
        try:
            raw = family.deform(p)
            measures = geometry_metrics(raw, family.areas)
            r['geometry'] = measures
            for n, x in raw:
                np.savetxt(folder/f'{n}.dat', x, header=n, comments='', fmt='%.16g')
            if measures['minimum_clearance'] is not None and measures['minimum_clearance'] < cfg['constraints']['minimum_clearance']:
                raise ValueError('Deformed input violates minimum clearance')
            if max(abs(v) for v in measures['relative_area_change'].values()) > cfg['constraints']['max_relative_area_change']:
                raise ValueError('Input area drift exceeds geometry tolerance')
            tasks = [(str(folder/f"{s['name']}__{point['name']}"), raw, cfg, panels or cfg['panels'], s, point)
                     for s in cfg['scenarios'] for point in cfg['points']]
            with ProcessPoolExecutor(workers) as pool:
                r['flows'] = list(pool.map(flow_job, tasks))
            r.update(paired_score(r['flows'], baseline or r['flows'], cfg))
        except ValueError as exc:
            r['reason'] = str(exc)
        write_json(result_path, r)
        print(f"{folder.name}: p={list(p)}, feasible={r['feasible']}, "
              f"worst mean delta CD={None if r['score'] is None else r['score']*1e4} counts; {r['reason']}", flush=True)
        return r

    def objective(p, number):
        r = evaluate(p, number)
        if number == 0:
            baseline.extend(r['flows'])
        return r

    search = pattern_search(objective, cfg['modes'], cfg['search']['max_designs'],
                            cfg['search']['minimum_step_fraction'], cfg['constraints']['acceptance_counts']*1e-4)
    write_json(out/'search.json', search)
    selected = search['parameters']
    verification = dict(status='no_shape_change_selected', accepted_for_screening=False)
    if any(abs(v) > 1e-14 for v in selected):
        # Pair the fine candidate with the fine baseline, never with base-grid CD.
        base_fine = evaluate([0.]*len(selected), 0, cfg['verification_panels'], 'verification_baseline')
        prior_baseline = baseline[:]
        baseline[:] = base_fine['flows']
        candidate_fine = evaluate(selected, 0, cfg['verification_panels'], 'verification_candidate')
        fine = paired_score(candidate_fine['flows'], base_fine['flows'], cfg)
        # base_fine was evaluated with the base baseline in scope; only its raw
        # flow records matter here. Recompute all comparisons on the fine grid.
        base_fine.update(paired_score(base_fine['flows'], base_fine['flows'], cfg))
        write_json(out/'verification_baseline/result.json', base_fine)
        base_means = search['best'].get('scenario_delta_cd_counts', {})
        fine_means = fine.get('scenario_delta_cd_counts', {})
        spread = {s: abs(fine_means[s]-base_means[s]) for s in base_means if s in fine_means}
        margin = cfg['constraints']['acceptance_counts']
        passed = bool(fine['feasible'] and fine['score']*1e4 < -margin and
                      len(spread) == len(cfg['scenarios']) and
                      all(-fine_means[s] > spread[s]+margin and -base_means[s] > spread[s]+margin for s in spread))
        verification = dict(status='passed_numerical_screen' if passed else 'not_verified',
                            accepted_for_screening=passed, fine_comparison=fine,
                            differential_grid_spread_counts=spread,
                            note='Observed two-grid sensitivity is not a physical error bound.')
        baseline[:] = prior_baseline
    write_json(out/'verification.json', verification)
    delivered = selected if verification['accepted_for_screening'] else [0.]*len(selected)
    final_dir = out/'selected_geometry'
    final_dir.mkdir(exist_ok=True)
    for n, x in family.deform(delivered):
        np.savetxt(final_dir/f'{n}.dat', x, header=n, comments='', fmt='%.16g')
    summary = dict(search_status=search['status'], evaluated_designs=len(search['history']),
                   candidate_parameters=selected, delivered_parameters=delivered,
                   numerical_screen_passed=verification['accepted_for_screening'],
                   independently_validated_shape_improvement=False,
                   geometry_folder='selected_geometry', verification=verification)
    write_json(out/'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 2 if search['status'] == 'baseline_infeasible' else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('config', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=2)
    args = ap.parse_args()
    if not 1 <= args.workers <= 4:
        ap.error('workers must be between 1 and 4')
    return run(args.config.resolve(), args.out.resolve(), args.workers)


if __name__ == '__main__':
    raise SystemExit(main())
