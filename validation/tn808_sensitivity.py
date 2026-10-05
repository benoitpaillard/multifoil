"""TN808 follow-up with narrower incidence brackets, retaining the first run.

The brackets use the first converged preliminary lift points, not drag fitting.
All settings in this follow-up use the same incidence procedure.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
from tn808_benchmark import ROOT,SETTINGS,specification,run_job


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',default='base,coarse,ncrit6,table_axis')
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/tn808_sensitivity')
    args=parser.parse_args()
    settings=args.settings.split(',')
    if not 1<=args.workers<=4 or len(set(settings))!=len(settings) or any(s not in SETTINGS for s in settings):
        parser.error('Use 1..4 workers and unique known settings')
    spec=json.loads(json.dumps(specification(settings)))
    spec['alpha_bounds']={'10':[0.,6.],'20':[-4.,1.],'30':[-8.,-3.]}
    spec['procedure']='Narrowed brackets after preliminary guarded run: CL at alpha0/-4/-8 '
    spec['procedure']+='was 1.090/1.275/1.242 for deflections10/20/30. No drag calibration.'
    spec['sha256']['validation/tn808_sensitivity.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=args.out/'study.json'
    if manifest.exists() and json.loads(manifest.read_text())!=spec:
        parser.error('Study changed; use a new output directory')
    manifest.write_text(json.dumps(spec,indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_job,str(args.out),j,spec) for j in spec['jobs']]
        for future in as_completed(futures):print(json.dumps(future.result()),flush=True)
    records=[json.loads((args.out/f"d{j['deflection']:02d}__{j['setting']}"/'result.json').read_text())
             for j in spec['jobs']]
    summary=dict(requested=len(records),converged=sum(r['feasible'] for r in records),records=records)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return 0 if all(r['feasible'] for r in records) else 2


if __name__=='__main__':raise SystemExit(main())
