#!/usr/bin/env python3
"""Evaluate archived S414 states with a FROZEN pre-integration source snapshot.

Run only by verify_s414_default_integration.py, in a separate process, so the
frozen modules never share an interpreter with the current ones. Geometry and
saved arrays arrive as data; nothing is imported from the current repository.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys

snapshot, exchange = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
# Drop this script's own directory: only the snapshot may supply solver modules.
sys.path[:] = [str(snapshot)]+[p for p in sys.path[1:] if Path(p or '.').resolve() != Path(__file__).resolve().parent]
import numpy as np
import viscous, multifoil, paneling
from viscous import CoupledSolver, mf, te_info

loaded = {name: str(Path(module.__file__).resolve()) for name, module in
          (('viscous', viscous), ('multifoil', multifoil), ('paneling', paneling), ('vendor.mfoil', mf))}
for name, path in loaded.items():
    if not Path(path).is_relative_to(snapshot):
        raise ValueError('Frozen child imported a current module: '+name)
if 'amplification_jacobian' in sys.modules:
    raise ValueError('Optional correction imported before the plain reference')


def restore(solver, alpha, arrays):
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = arrays[model.geom.name]
        if np.max(abs(model.foil.x-a['xy'])) > 2e-12 or np.max(abs(model.wake.x-a['wake_xy'])) > 2e-11:
            raise ValueError('Frozen rebuild differs from saved geometry')
        model.glob.U = a['U'].copy(); model.glob.Nsys = a['U'].shape[1]
        model.vsol.turb = a['turbulent'].copy()
        model.isol.Istag = [int(a['lower'][0]), int(a['upper'][0])]
        model.isol.sgnue = np.ones(model.foil.N); model.isol.sgnue[a['lower']] = -1
        model.isol.sstag = float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model); mf.stagpoint_move(model)
        if np.max(abs(model.isol.xi-a['xi'])) > 2e-12:
            raise ValueError('Frozen restore differs from saved stagnation distances')
        if abs(te_info(model.foil.x)[1]) < 1e-12: model.vsol.wgap = np.zeros(model.wake.N)
        else: mf.set_wake_gap(model)
    solver._build_source_coupling()
    return solver


cases = json.loads((exchange/'cases.json').read_text())
settings = cases['settings']; record = dict(loaded=loaded, cases={})


def build(kind, case, raw, arrays):
    kwargs = dict(reynolds=settings['reynolds'], mach=settings['mach'], ncrit=settings['ncrit'],
                  iterations=settings['iterations'], tolerance=settings['residual_tolerance'],
                  wake_length=case['wake_length'])
    solver = kind(raw, **kwargs)
    return restore(solver, case['alpha_deg'], arrays)


inputs = {}
for case in cases['cases']:
    with np.load(exchange/(case['name']+'_input.npz')) as data:
        raw = [(name, data['raw_'+name].copy()) for name in ('fore', 'aft')]
        arrays = {name: {key[len(name)+1:]: data[key].copy() for key in data.files if key.startswith(name+'_')}
                  for name in ('fore', 'aft')}
    inputs[case['name']] = (case, raw, arrays)

# 1. Historical solver with no correction module loaded: the original mathematics.
for name, (case, raw, arrays) in inputs.items():
    if mf.get_damp.__module__ != mf.__name__:
        raise ValueError('Plain reference no longer uses vendored get_damp')
    residual, jacobian = build(CoupledSolver, case, raw, arrays)._system()
    np.save(exchange/(name+'_plain_R.npy'), residual); np.save(exchange/(name+'_plain_J.npy'), jacobian)

# 2. The frozen optional implementation, as verified in the earlier derivative study.
import amplification_jacobian as option
if not Path(option.__file__).resolve().is_relative_to(snapshot):
    raise ValueError('Frozen child imported the current correction module')
record['optional_default'] = option.AmplificationJacobianSolver.__init__.__kwdefaults__['amplification_jacobian']
for name, (case, raw, arrays) in inputs.items():
    for mode in (option.CORRECTED, option.LEGACY):
        solver = build(lambda *a, **k: option.AmplificationJacobianSolver(*a, amplification_jacobian=mode, **k), case, raw, arrays)
        if any(m.param.amplification_jacobian != mode for m in solver.models):
            raise ValueError('Frozen optional mode not applied')
        residual, jacobian = solver._system()
        np.save(exchange/(name+'_'+mode+'_R.npy'), residual); np.save(exchange/(name+'_'+mode+'_J.npy'), jacobian)
    record['cases'][name] = dict(nsys=int(solver.nsys))
(exchange/'frozen_record.json').write_text(json.dumps(record, indent=2)+'\n')
