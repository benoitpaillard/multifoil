#!/usr/bin/env python3
"""Run/resume the TN808 slotted-flap reconstruction at fixed lift.

Every attempted incidence, including failure, is archived. Geometry/transition
interpretations are assumptions to assess, not parameters fitted to experiment.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from multifoil import repanel,validate_geometry,save_results
from viscous import CoupledSolver
from optimization import solve_fixed_lift,LiftTargetError
from tn808_geometry import geometry,geometry_audit
import numpy as np

SETTINGS={
    'coarse':dict(panels=[180,120]),
    'base':dict(panels=[260,180]),
    'fine':dict(panels=[340,240]),
    'ncrit6':dict(panels=[260,180],ncrit=6.),
    'ncrit3':dict(panels=[260,180],ncrit=3.),
    'table_axis':dict(panels=[260,180],datum='table_axis'),
    'entry_plus':dict(panels=[260,180],entry_adjustment=.001),
    'gap_minus':dict(panels=[260,180],gap=.019),
    'gap_plus':dict(panels=[260,180],gap=.021),
    'long_wake':dict(panels=[260,180],wake_length=2.),
}
BOUNDS={10:[0.,8.],20:[-4.,4.],30:[-8.,0.]}


def specification(settings):
    files=['multifoil.py','viscous.py','vendor/mfoil.py','optimization.py',
           'validation/tn808_geometry.py','validation/tn808_benchmark.py',
           'validation/data/naca23012_tn808/source.json']
    params={s:dict(ncrit=9.,gap=.02,datum='nose_to_te',entry_adjustment=0.,
                   wake_length=1.,**{}) for s in settings}
    for s in params:params[s].update(SETTINGS[s])
    return dict(schema_version=1,reference='NACA TN808 Figure8: flap1-b, gap0.02c',
        target_cl=1.5,reynolds=2190000.,mach=.105,reference_chord=1.,
        moment_reference=[.25,0.],cl_tolerance=2e-5,residual_tolerance=1e-7,
        iterations=100,alpha_bounds=BOUNDS,settings=params,
        jobs=[dict(deflection=d,setting=s) for s in settings for d in BOUNDS],
        geometry_audit=geometry_audit(),
        source_qualifications=['Re is physical, not historical effective Re.',
            'Mach0.105 and natural-transition Ncrit are assumptions.',
            'Table-based geometry reconstruction and flap datum require qualification.',
            'Reference drag includes unquantified fitting contribution.'],
        sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files})


def run_job(out_string,job,spec):
    d,s=job['deflection'],job['setting']
    out=Path(out_string)/f'd{d:02d}__{s}'
    out.mkdir(parents=True,exist_ok=True)
    path=out/'result.json'
    if path.exists():return json.loads(path.read_text())
    cfg=spec['settings'][s]
    start=time.perf_counter()
    records=[]
    result=dict(deflection=d,setting=s,feasible=False,CL=None,CD=None,CM=None,
                alpha_deg=None,failure_reason='',flow_evaluations=0)
    try:
        raw,audit=geometry(d,cfg['gap'],cfg['datum'],cfg['entry_adjustment'])
        (out/'geometry_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        for name,xy in raw:
            np.savetxt(out/f'input_{name}.dat',xy,fmt='%.16g',header=name,comments='')
        elements=[(name,repanel(xy,n)) for (name,xy),n in zip(raw,cfg['panels'])]
        validate_geometry(elements)
        bounds=spec['alpha_bounds'][str(d)]
        solver=CoupledSolver(elements,reynolds=spec['reynolds'],mach=spec['mach'],
            ncrit=cfg['ncrit'],iterations=spec['iterations'],
            tolerance=spec['residual_tolerance'],wake_length=cfg['wake_length'],
            start_alpha=.5*sum(bounds))
        def checkpoint(solution):
            records.append(solution)
            save_results(out/'flow',records,elements,dict(study=spec,job=job),'viscous')
        answer,_=solve_fixed_lift(solver,spec['target_cl'],bounds,
                                  spec['cl_tolerance'],on_evaluation=checkpoint)
        result.update(feasible=True,CL=answer.cl,CD=answer.cd,CM=answer.cm,
                      alpha_deg=answer.alpha,lift_error=answer.cl-spec['target_cl'],
                      residual=answer.residual,
                      local_supersonic=any(answer.diagnostics.get('local_supersonic_detected',{}).values()),
                      max_separated_arclength_fraction=max(
                        side['separated_arclength_fraction']
                        for bl in answer.diagnostics['boundary_layers'].values() for side in bl.values()))
    except (ValueError,RuntimeError,FloatingPointError,LiftTargetError) as exc:
        result['failure_reason']=f'{type(exc).__name__}: {exc}'
    result.update(flow_evaluations=len(records),wall_seconds=time.perf_counter()-start)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',default='base')
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/tn808_base')
    args=parser.parse_args()
    settings=args.settings.split(',')
    if not 1<=args.workers<=4 or len(set(settings))!=len(settings) or any(s not in SETTINGS for s in settings):
        parser.error('Use 1..4 workers and unique known settings: '+','.join(SETTINGS))
    spec=json.loads(json.dumps(specification(settings)))
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=args.out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text())!=spec:
            parser.error('Code, geometry, inputs or requested jobs changed; choose a new output directory.')
    else:manifest.write_text(json.dumps(spec,indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_job,str(args.out),j,spec) for j in spec['jobs']]
        for future in as_completed(futures):print(json.dumps(future.result()),flush=True)
    records=[json.loads((args.out/f"d{j['deflection']:02d}__{j['setting']}"/'result.json').read_text())
             for j in spec['jobs']]
    summary=dict(requested=len(records),converged=sum(r['feasible'] for r in records),
                 new_independent_validation_claim=False,records=records)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return 0 if all(r['feasible'] for r in records) else 2


if __name__=='__main__':
    raise SystemExit(main())
