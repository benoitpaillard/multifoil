#!/usr/bin/env python3
"""Identify S702 drag sensitivity to transition, without fitting drag offsets.

The N=6 scenario is a diagnostic disturbance assumption, not a measurement of
this tunnel. The same assumption is transferred to all four configurations and
three lift targets. Fore-only and aft-only changes isolate each element's role.
Every requested flow evaluation and the final transition residual are retained.
"""
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
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import repanel, save_results
from optimization import solve_fixed_lift
from viscous import CoupledSolver
from vendor import mfoil as mf
from s702_geometry import geometry, FILES
from s702_benchmark import record


def transition_audit(solver):
    """Check the actual subcell e^N equation after the global solve converged."""
    rows = []
    for model in solver.models:
        for side, label in enumerate(('lower', 'upper')):
            ids = np.asarray(model.vsol.Is[side], dtype=int)
            flags = model.vsol.turb[ids]
            transitions = np.flatnonzero(np.diff(flags.astype(int)) == 1)+1
            for i in transitions:
                js = ids[[i-1, i]]
                x = model.isol.xi[js]
                u = model.glob.U[:, js]
                param = mf.build_param(model, side)
                param.turb = False
                r, _, _ = mf.residual_transition(model, param, x, u, np.zeros(2))
                xt = model.vsol.xt
                ut = u[:, 0]+(xt-x[0])/(x[1]-x[0])*(u[:, 1]-u[:, 0])
                ut[2] = param.ncrit
                d1, _ = mf.get_damp(u[:, 0], param)
                dt, _ = mf.get_damp(ut, param)
                ramp = param.ncrit-u[2, 0]-.5*(xt-x[0])*(d1+dt)
                rows.append(dict(element=model.geom.name, side=label,
                    bracket=x.tolist(), transition_xi=float(xt),
                    inside_cell=bool(x[0] <= xt <= x[1]),
                    amplification_residual=float(ramp),
                    station_residual_max=float(np.max(abs(r)))))
    return rows


def run_job(out, job, spec):
    folder = Path(out)/job['name']
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/'result.json'
    if path.exists():
        return json.loads(path.read_text())
    start = time.perf_counter()
    result = dict(**job, feasible=False, CL=None, CD=None, CM=None,
                  alpha_deg=None, error='', records=[])
    elements = [(name, repanel(x, p)) for (name, x), p in
                zip(geometry(job['flap']), spec['panels'])]
    solutions = []

    def checkpoint(sol):
        # The legacy solver's scalar diagnostic records the fore-element Ncrit.
        # Preserve the complete scenario explicitly on every flow record.
        sol.diagnostics['ncrit_by_element'] = job['ncrit_by_element']
        solutions.append(sol)
        result['records'].append(record(sol, 'fixed_lift'))
        save_results(folder/'flow', solutions, elements,
                     dict(study=spec, job=job), 'viscous')
        (folder/'progress.json').write_text(json.dumps(result, indent=2)+'\n')

    try:
        solver = CoupledSolver(elements, reynolds=spec['reynolds'],
            mach=spec['mach'], ncrit=9., iterations=100, tolerance=1e-7,
            wake_length=spec['wake_length'], start_alpha=-1.)
        for model in solver.models:
            model.param.ncrit = job['ncrit_by_element'][model.geom.name]
        sol, _ = solve_fixed_lift(solver, job['target_cl'], job['bounds'],
                                 2e-5, 20, checkpoint)
        result.update(feasible=True, CL=sol.cl, CD=sol.cd, CM=sol.cm,
            alpha_deg=sol.alpha, boundary_layers=sol.diagnostics['boundary_layers'],
            element_forces=sol.diagnostics['element_forces'],
            transition_audit=transition_audit(solver))
        for model in solver.models:
            np.savez_compressed(folder/(model.geom.name+'_state.npz'),
                U=model.glob.U, xy=model.foil.x, xi=model.isol.xi,
                turbulent=model.vsol.turb, wake_xy=model.wake.x,
                lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]))
    except (ValueError, RuntimeError, FloatingPointError, AssertionError) as exc:
        result['error'] = f'{type(exc).__name__}: {exc}'
        result['feasible'] = False
    result['wall_seconds'] = time.perf_counter()-start
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_transition_diagnostic')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--ncrit', type=float, default=6.)
    ap.add_argument('--panels', default='260,180')
    ap.add_argument('--flaps', default='-2,0,5,10')
    ap.add_argument('--targets', default='.5,.8,1.2')
    ap.add_argument('--isolate', action='store_true')
    ap.add_argument('--mach', type=float, default=.1)
    args = ap.parse_args()
    panels = [int(v) for v in args.panels.split(',')]
    flaps = [int(v) for v in args.flaps.split(',')]
    targets = [float(v) for v in args.targets.split(',')]
    if (len(panels) != 2 or min(panels) < 50 or sum(panels) > 1200 or
        not 1 <= args.workers <= 4 or not set(flaps) <= FILES.keys() or
        not set(targets) <= {.5, .8, 1.2} or not np.isfinite(args.ncrit) or
        args.ncrit <= 0 or not 0 <= args.mach < .7):
        ap.error('Invalid conditions or selections')
    jobs = []
    seed_paths = []
    for flap in flaps:
        seed_path = ROOT/f'validation/results/s702_fixed/d{flap:+d}__base/result.json'
        seed_paths.append(str(seed_path.relative_to(ROOT)))
        seed = json.loads(seed_path.read_text())['targets']
        for target in targets:
            prior = next(r for r in seed if r['target_cl'] == target)
            # Previous bracketing points, not previous solutions, initialize this
            # changed-physics calculation. Expand one degree for sensitivity.
            bounds = [prior['bounds'][0]-1., prior['bounds'][1]+1.]
            scenarios = [('both', dict(fore=args.ncrit, aft=args.ncrit))]
            if args.isolate and flap in (5, 10) and target == .8:
                scenarios += [('fore', dict(fore=args.ncrit, aft=9.)),
                              ('aft', dict(fore=9., aft=args.ncrit))]
            for label, values in scenarios:
                jobs.append(dict(name=f'd{flap:+d}_cl{target:g}_{label}',
                    flap=flap, target_cl=target, scenario=label,
                    ncrit_by_element=values, bounds=bounds))
    inputs = ['multifoil.py', 'viscous.py', 'paneling.py', 'vendor/mfoil.py',
              'optimization.py', 'validation/s702_geometry.py',
              'validation/s702_benchmark.py', 'validation/s702_transition_study.py',
              'validation/data/s702/source.json', 'validation/data/s702/fore.dat']
    inputs += ['validation/data/s702/'+FILES[f] for f in flaps]+seed_paths
    spec = dict(schema=1, reynolds=1e6, mach=args.mach, ncrit=args.ncrit,
        panels=panels, wake_length=1., cl_tolerance=2e-5,
        residual_tolerance=1e-7, iterations=100, start_alpha=-1.,
        validation_claim=False, jobs=jobs,
        interpretation='Transition sensitivity; Ncrit is not a measured tunnel value',
        sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != spec:
            raise ValueError('Study/source changed; select a new output directory')
    elif any(args.out.iterdir()):
        raise ValueError('Nonempty output directory lacks a manifest')
    else:
        manifest.write_text(json.dumps(spec, indent=2)+'\n')
    with ProcessPoolExecutor(args.workers) as pool:
        for fut in as_completed([pool.submit(run_job, str(args.out), j, spec) for j in jobs]):
            r = fut.result()
            print({k: r[k] for k in ('name', 'feasible', 'CD', 'error', 'wall_seconds')}, flush=True)
    rows = [json.loads((args.out/j['name']/'result.json').read_text()) for j in jobs]
    (args.out/'summary.json').write_text(json.dumps(rows, indent=2)+'\n')
    return 0 if all(r['feasible'] for r in rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
