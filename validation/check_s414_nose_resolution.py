#!/usr/bin/env python3
"""Small prelaunch linear-operator check; no boundary-layer/Newton solves."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s414_nose_resolution import Solver, build_elements
from viscous import CoupledSolver, mf


def source_check(solver):
    """Build the real source coupling; independently check all block indices."""
    alpha = -3.459015324304847
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha
        model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        mf.stagpoint_find(model)
        mf.identify_surfaces(model)
    solver._build_source_coupling()
    fs, ss, sc = 0, 0, 0
    blocks = []
    mass = 1e-4*(1.+np.sin(np.arange(solver.nsys)*.37))
    strengths, signs = [], []
    for model, foil, system, source in zip(solver.models, solver.fs, solver.ss, solver.sources):
        n, nw = model.foil.N, model.wake.N
        expected = [(fs, fs+n), (ss, ss+n+nw), (sc, sc+n+nw-2)]
        if [(s.start, s.stop) for s in (foil, system, source)] != expected:
            raise ValueError('Source/system/surface block indices differ')
        local = mass[ss:ss+n+nw]
        strengths.extend(np.diff(model.isol.sgnue*local[:n])/np.diff(model.foil.s))
        strengths.extend(np.diff(local[n:])/np.diff(model.wake.s))
        signs.extend(model.isol.sgnue)
        signs.extend(np.ones(nw))
        if not np.array_equal(model.vsol.ue_sigma, solver.ue_sigma[system, source]):
            raise ValueError('Local source operator is not its declared global block')
        blocks.append(dict(element=model.geom.name, surface=expected[0], system=expected[1], source=expected[2]))
        fs, ss, sc = expected[0][1], expected[1][1], expected[2][1]
    if (fs, ss, sc) != (solver.nfoil, solver.nsys, solver.nsource):
        raise ValueError('Source/system blocks do not cover the whole operator')
    if solver.ue_sigma.shape != (ss, sc) or solver.D.shape != (ss, ss):
        raise ValueError('Source coupling dimensions differ')
    if not np.isfinite(solver.ue_sigma).all() or not np.isfinite(solver.D).all():
        raise ValueError('Nonfinite source coupling')
    direct = np.asarray(signs)*(solver.ue_sigma@np.asarray(strengths))
    error = float(np.max(abs(solver.D@mass-direct)))
    if error > 1e-10*max(1., float(np.max(abs(direct)))):
        raise ValueError('Global mass/source mapping fails the independent action check')
    return dict(passed=True, alpha_deg=alpha, surface_nodes=fs, system_nodes=ss,
                sources=sc, blocks=blocks, independent_mass_action_error=error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    geometry_path = args.output.with_suffix('.npz')
    if args.output.exists() or geometry_path.exists():
        raise ValueError('Preserve existing check receipt; choose a new output')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    geometry = {}
    for factor in (1, 2, 4):
        solver = Solver(build_elements(), subdivision_factor=factor, reynolds=1e6, mach=.1, ncrit=9.)
        # Streamfunction constants are the last element-count unknowns.
        from scipy.linalg import lu_solve
        full = lu_solve(solver.lu, solver.rhs)
        residual = float(np.max(abs(solver.aic@full-solver.rhs)))
        if not np.isfinite(full).all() or residual > 1e-10:
            raise ValueError('Refined linear operator fails its residual criterion')
        if factor == 1:
            native = CoupledSolver(build_elements(), reynolds=1e6, mach=.1, ncrit=9.)
            for label in ('aic', 'rhs', 'gref'):
                if not np.array_equal(getattr(solver, label), getattr(native, label)):
                    raise ValueError('Factor-1 operator differs from the standard solver: '+label)
        source = source_check(solver)
        if factor == 1:
            source_check(native)
            for label in ('ue_sigma', 'D'):
                if not np.array_equal(getattr(solver, label), getattr(native, label)):
                    raise ValueError('Factor-1 source coupling differs from the standard solver: '+label)
        rows.append(dict(factor=factor, matrix_size=solver.aic.shape[0], residual=residual,
                         resolution_controls=solver.resolution_controls, source_coupling=source))
        for model in solver.models:
            foil = model.foil
            prefix = f'factor{factor}_{model.geom.name}'
            for label in ('x', 's', 't'):
                geometry[prefix+'_'+label] = getattr(foil, label)
        print(f'factor {factor}: {solver.nfoil} surface nodes, linear residual {residual:.3g}', flush=True)
    np.savez_compressed(geometry_path, **geometry)
    value = dict(passed=True, factors=rows, nonlinear_solves=0,
        geometry_archive=geometry_path.name, geometry_sha256=hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        baseline_operator_exact=True, baseline_source_coupling_exact=True, source_sha256={str(p.relative_to(ROOT)):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [Path(__file__), ROOT/'validation/s414_nose_resolution.py', ROOT/'viscous.py',
         ROOT/'diagnostic_discretization.py', ROOT/'paneling.py', ROOT/'vendor/mfoil.py']})
    args.output.write_text(json.dumps(value, indent=2)+'\n')


if __name__ == '__main__':
    main()
