#!/usr/bin/env python3
"""Single-element shape-change verification against the bundled native mfoil.

Both implementations use the identical actual surface coordinates and closures.
Native mfoil uses its own viscous interaction/Newton solve and lift constraint.
This isolates implementation effects; it is not independent physical validation.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from multifoil import naca4, repanel, save_results
from viscous import CoupledSolver, mf
from optimization import solve_fixed_lift
from shape_optimization import ShapeFamily, write_json, CODE


def velocity_basis(model):
    """Correct the native Python 2xN/Nx2 layout bug in get_ueinvref only."""
    air = model.isol.sgnue[:, None]*model.isol.gamref
    if model.oper.viscous and model.wake.N:
        wake = model.isol.uewiref.copy()
        wake[0] = air[-1]
        return np.vstack((air, wake))
    return air


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/incompressible_shape_check')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    spec = dict(mach=0., reynolds=1e6, ncrit=9., target_cl=.5, panels=200,
                amplitudes=[-.001, 0., .001], support=[.05, .95],
                native_reference='bundled Python mfoil, independent Newton/lift constraint; identical surfaces',
                native_adapters=['sparse Jacobian reallocation', 'get_ueinvref Nx2 array layout, checked against get_ueinv',
                                 'scalar lift residual packed as a length-one array for concatenate'],
                sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest()
                        for f in CODE+['validation/incompressible_shape_check.py']})
    path = args.out/'study.json'
    if path.exists():
        if json.loads(path.read_text()) != spec:
            raise ValueError('Study/code changed: choose a new output folder')
    elif any(args.out.iterdir()):
        raise ValueError('Nonempty output folder lacks study manifest')
    else:
        write_json(path, spec)
    family = ShapeFamily([('foil', naca4('0012', 401))],
                         [dict(element='foil', support=spec['support'], bounds=[-.002, .002])])
    records = []
    for amplitude in spec['amplitudes']:
        folder = args.out/f'camber_{amplitude:+.4f}'
        folder.mkdir(exist_ok=True)
        record_path = folder/'result.json'
        if record_path.exists():
            records.append(json.loads(record_path.read_text()))
            continue
        elements = [(n, repanel(x, 200)) for n, x in family.deform([amplitude])]
        solver = CoupledSolver(elements, reynolds=1e6, mach=0., ncrit=9.,
                               iterations=100, tolerance=1e-8, start_alpha=3.)
        solutions = []
        def checkpoint(sol):
            solutions.append(sol)
            save_results(folder/'multifoil', solutions, elements, spec, 'viscous')
        sol, _ = solve_fixed_lift(solver, .5, (3., 6.), 2e-6, 20, checkpoint)
        model = copy.deepcopy(solver.models[0])
        exact_surface = model.foil.x.copy()
        model.oper.givencl = True
        model.oper.cltgt = .5
        model.oper.initbl = True
        model.param.rtol = 1e-8
        model.param.niglob = 100
        original = mf.build_glob_sys
        original_basis = mf.get_ueinvref
        original_lift = mf.clalpha_residual
        # Check this layout adapter against the native velocity and its
        # incidence derivative before using the native coupled lift trim.
        old_alpha = model.oper.alpha
        h = 1e-4
        model.oper.alpha = old_alpha+h
        plus = mf.get_ueinv(model)
        model.oper.alpha = old_alpha-h
        minus = mf.get_ueinv(model)
        model.oper.alpha = old_alpha
        basis = velocity_basis(model)
        a = np.deg2rad(old_alpha)
        np.testing.assert_allclose(basis@np.array([np.cos(a), np.sin(a)]), mf.get_ueinv(model), atol=1e-13)
        derivative = basis@np.array([-np.sin(a), np.cos(a)])*np.pi/180
        difference = float(np.max(abs(derivative-(plus-minus)/(2*h))))
        if difference > 1e-9:
            raise RuntimeError('Velocity-basis derivative check failed')
        def sparse_compatibility(m):
            # Same sparse reallocation compatibility fix as test_viscous.py.
            m.glob.realloc = True
            original(m)
        def lift_layout(m):
            residual, incidence, state = original_lift(m)
            return np.atleast_1d(residual), incidence, state
        try:
            mf.build_glob_sys = sparse_compatibility
            mf.get_ueinvref = velocity_basis
            mf.clalpha_residual = lift_layout
            mf.solve_viscous(model)
        finally:
            mf.build_glob_sys = original
            mf.get_ueinvref = original_basis
            mf.clalpha_residual = original_lift
        if not model.glob.conv or abs(model.post.cl-.5) > 2e-6:
            raise RuntimeError('Native mfoil did not converge to the prescribed lift')
        np.testing.assert_array_equal(model.foil.x, exact_surface)
        native = dict(CL=float(model.post.cl), CD=float(model.post.cd),
                      CM=float(model.post.cm), alpha_deg=float(model.oper.alpha))
        r = dict(amplitude=amplitude, multifoil=dict(CL=sol.cl, CD=sol.cd, CM=sol.cm, alpha_deg=sol.alpha),
                 native_mfoil=native, flow_evaluations=len(solutions),
                 velocity_basis_derivative_max_error=difference,
                 delta_cd_counts=(sol.cd-model.post.cd)*1e4)
        mf.get_distributions(model)
        np.savetxt(folder/'native_surface.csv', np.column_stack((model.foil.x.T, model.post.cp[:model.foil.N],
                   model.post.cf[:model.foil.N], model.glob.U[:, :model.foil.N].T)),
                   delimiter=',', header='x,y,cp,cf,theta,delta_star,amplification_or_shear,ue', comments='')
        write_json(record_path, r)
        records.append(r)
        print(r, flush=True)
    derivatives = {method: (records[2][method]['CD']-records[0][method]['CD'])/.002
                   for method in ('multifoil', 'native_mfoil')}
    summary = dict(records=records, central_dcd_damplitude=derivatives,
                   derivative_difference=derivatives['multifoil']-derivatives['native_mfoil'],
                   same_derivative_sign=derivatives['multifoil']*derivatives['native_mfoil']>0,
                   independent_aerodynamic_validation=False)
    write_json(args.out/'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
