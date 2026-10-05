"""Retry one failed TN808 job with explicit, archived incidence initialization."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from tn808_benchmark import ROOT,geometry,repanel,validate_geometry,CoupledSolver,save_results
from optimization import solve_fixed_lift,LiftTargetError
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',type=Path,required=True)
    parser.add_argument('--deflection',type=int,required=True)
    parser.add_argument('--setting',required=True)
    parser.add_argument('--bounds',type=float,nargs=2,required=True)
    parser.add_argument('--start-alpha',type=float,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    args.parent=args.parent.resolve();args.out=args.out.resolve()
    spec=json.loads((args.parent/'study.json').read_text())
    for name,digest in spec['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            parser.error('Parent numerical sources changed')
    parent_result=args.parent/f'd{args.deflection:02d}__{args.setting}'/'result.json'
    parent=json.loads(parent_result.read_text())
    if parent['feasible']:parser.error('This helper retains a retry of a failed job only')
    spec['jobs']=[dict(deflection=args.deflection,setting=args.setting)]
    spec['alpha_bounds']={str(args.deflection):args.bounds}
    spec['retry_of']=dict(path=str(parent_result.relative_to(ROOT)),sha256=hashlib.sha256(parent_result.read_bytes()).hexdigest(),reason=parent['failure_reason'])
    spec['start_alpha']=args.start_alpha;spec['max_evaluations']=24
    spec['sha256']['validation/tn808_retry.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=args.out/'study.json'
    if manifest.exists() and json.loads(manifest.read_text())!=spec:parser.error('Study changed')
    manifest.write_text(json.dumps(spec,indent=2)+'\n')
    out=args.out/f'd{args.deflection:02d}__{args.setting}';out.mkdir(exist_ok=True)
    result_path=out/'result.json'
    if result_path.exists():print(result_path.read_text());return
    cfg=spec['settings'][args.setting]
    raw,audit=geometry(args.deflection,cfg['gap'],cfg['datum'],cfg['entry_adjustment'])
    (out/'geometry_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    for name,xy in raw:np.savetxt(out/f'input_{name}.dat',xy,fmt='%.16g',header=name,comments='')
    elements=[(name,repanel(xy,n)) for (name,xy),n in zip(raw,cfg['panels'])]
    records=[];start=time.perf_counter()
    result=dict(deflection=args.deflection,setting=args.setting,feasible=False,CL=None,CD=None,CM=None,
                alpha_deg=None,failure_reason='',flow_evaluations=0)
    def checkpoint(solution):
        records.append(solution);save_results(out/'flow',records,elements,dict(study=spec),'viscous')
    try:
        validate_geometry(elements)
        solver=CoupledSolver(elements,reynolds=spec['reynolds'],mach=spec['mach'],ncrit=cfg['ncrit'],
            iterations=spec['iterations'],tolerance=spec['residual_tolerance'],
            wake_length=cfg['wake_length'],start_alpha=args.start_alpha)
        answer,_=solve_fixed_lift(solver,spec['target_cl'],args.bounds,spec['cl_tolerance'],
                                 max_evaluations=24,on_evaluation=checkpoint)
        result.update(feasible=True,CL=answer.cl,CD=answer.cd,CM=answer.cm,alpha_deg=answer.alpha,
            lift_error=answer.cl-spec['target_cl'],residual=answer.residual,
            local_supersonic=any(answer.diagnostics['local_supersonic_detected'].values()),
            max_separated_arclength_fraction=max(side['separated_arclength_fraction']
                for bl in answer.diagnostics['boundary_layers'].values() for side in bl.values()))
    except (ValueError,RuntimeError,FloatingPointError,LiftTargetError) as exc:
        result['failure_reason']=f'{type(exc).__name__}: {exc}'
    result.update(flow_evaluations=len(records),wall_seconds=time.perf_counter()-start)
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    (args.out/'summary.json').write_text(json.dumps(dict(requested=1,converged=int(result['feasible']),records=[result]),indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return 0 if result['feasible'] else 2


if __name__=='__main__':raise SystemExit(main())
