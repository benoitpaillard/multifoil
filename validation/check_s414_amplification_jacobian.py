#!/usr/bin/env python3
"""Check the assembled coupled Jacobian at two unchanged archived S414 states."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import numpy as np
from audit_s414_lower_bucket import Inputs, load_states, ROOT, sha, require, write
from report_s414_wake_transfer import restore
from report_s702_wake_extent import csv_write
from amplification_jacobian import configure, LEGACY, CORRECTED, mf

OUT = ROOT/'validation/results/s414_jacobian_check_20261003'
PRIOR = ROOT/'validation/results/s414_fore_nose_20261003'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(); out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    require(not (out/'summary.json').exists(), 'Preserve completed check; choose a new output directory')
    inputs = Inputs()
    receipt = inputs.json(PRIOR/'verification.json')
    require(receipt['complete'] and not receipt['original_growth_derivative_check_pass'], 'Prior finding changed')
    for name, digest in receipt['output_sha256'].items():
        require(sha(inputs.add(PRIOR/name)) == digest, 'Prior diagnostic changed: '+name)
    for name, digest in inputs.json(PRIOR/'input_manifest.json')['files'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Prior frozen input changed: '+name)
    spec, states, _ = load_states(inputs)
    rows, cases = [], []
    for name in ('cl0056_nominal', 'cl028_control'):
        job = next(j for j in spec['jobs'] if j['name'] == name)
        solver = restore(job, spec, states[name]); fore = solver.models[0]
        saved = [m.glob.U.copy() for m in solver.models]; flags = [m.vsol.turb.copy() for m in solver.models]
        configure(solver.models, LEGACY); r0, j0 = solver._system()
        configure(solver.models, CORRECTED); r1, j1 = solver._system()
        residual_change = float(abs(r1-r0).max())
        require(residual_change <= 1e-10, 'Derivative-only option changed residual values')
        ids = np.asarray(fore.vsol.Is[0]); lam = ids[fore.vsol.turb[ids] == 0]
        p = mf.build_param(fore, 0)
        active = [int(j) for j in lam if mf.get_damp(fore.glob.U[:, j], p)[0] > 1e-4]
        selected = sorted(set([max(active, key=lambda j: mf.get_Hk(fore.glob.U[:, j], p)[0]), int(lam[-1])]))
        for node in selected:
            require(node not in fore.isol.Istag, 'Probe must avoid moving stagnation coordinate')
            for variable in range(4):
                direction = np.zeros(4*solver.nsys)
                direction[4*node+variable] = max(abs(saved[0][variable, node]), 1e-6)
                expected = j1@direction; old = j0@direction
                for step in (1e-4, 1e-5):
                    try:
                        fore.glob.U[variable, node] = saved[0][variable, node]+step*direction[4*node+variable]
                        plus = solver._system()[0]
                        fore.glob.U[variable, node] = saved[0][variable, node]-step*direction[4*node+variable]
                        minus = solver._system()[0]
                    finally:
                        fore.glob.U[:] = saved[0]
                    numerical = (plus-minus)/(2*step)
                    # Entire coupled directional derivative, including the transition cell
                    # and velocity coupling; no row is omitted after inspecting an error.
                    scale = max(float(abs(expected).max()), 1e-10)
                    rows.append(dict(case=name, node=node, x=float(fore.foil.x[0, node]),
                        variable=('theta', 'delta_star', 'N', 'ue')[variable], relative_step=step,
                        corrected_relative_error=float(abs(numerical-expected).max()/scale),
                        legacy_relative_error=float(abs(numerical-old).max()/scale),
                        direction_scale=float(direction[4*node+variable]), derivative_scale=scale))
        configure(solver.models, LEGACY); back, jb = solver._system()
        require(np.array_equal(back, r0) and np.array_equal(jb, j0), 'Legacy system not restored exactly')
        require(all(np.array_equal(m.glob.U, u) and np.array_equal(m.vsol.turb, f)
                    for m, u, f in zip(solver.models, saved, flags)), 'Saved states changed')
        cases.append(dict(case=name, residual_value_maximum_change=residual_change,
            legacy_system_restored_exactly=True, saved_arrays_restored_exactly=True, selected_nodes=selected))
        print(name, 'max corrected full-system relative error', max(r['corrected_relative_error'] for r in rows if r['case'] == name), flush=True)
    csv_write(out/'assembled_derivatives.csv', rows)
    for name in ('amplification_jacobian.py', 'validation/check_s414_amplification_jacobian.py',
                 'validation/report_s414_wake_transfer.py', 'test_amplification_jacobian.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(files=inputs.files))
    passed = all(r['corrected_relative_error'] < 2e-5 for r in rows)
    summary = dict(complete=True, passed=passed, criterion=2e-5, cases=cases, comparisons=len(rows),
        maximum_corrected_relative_error=max(r['corrected_relative_error'] for r in rows),
        maximum_legacy_relative_error=max(r['legacy_relative_error'] for r in rows),
        new_nonlinear_flow_solves=0, runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')}))
    write(out/'summary.json', summary)
    require(passed, 'Assembled corrected Jacobian check failed; preserve receipt and do not launch')


if __name__ == '__main__':
    main()
