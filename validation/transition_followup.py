#!/usr/bin/env python3
"""Fine-grid follow-up after the pilot's Ncrit-dependent ranking reversal.

This adaptive check was selected after seeing the pilot; it is not a holdout.
It preserves the original 26-job study and reuses its unchanged flow driver.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from validation.optimization_pilot import specification, one_job
from validate import verify_data, write_csv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/transition_followup')
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    if not 1<=args.workers<=8:
        parser.error('workers must be between 1 and 8')
    verify_data()
    spec=specification()
    spec['settings']={'ncrit6_fine':dict(panels=[340,240],wake_length=1.,ncrit=6,start_alpha=4.)}
    spec['jobs']=[(g,'ncrit6_fine') for g in ('baseline','camber_minus','camber_plus')]
    spec['reason']='Adaptive numerical check of the Ncrit-dependent camber ranking reversal; not a holdout.'
    spec['sha256']['validation/transition_followup.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    spec=json.loads(json.dumps(spec))
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    manifest=out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text())!=spec:
            raise ValueError('Changed study/code: use a new output directory')
    elif any(out.iterdir()):
        raise ValueError('Nonempty folder lacks a study manifest')
    else:
        manifest.write_text(json.dumps(spec,indent=2)+'\n')
    records=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(one_job,str(out),g,s,spec) for g,s in spec['jobs']]
        for future in as_completed(jobs):
            r=future.result();records.append(r)
            print(f"{r['geometry']}: feasible={r['feasible']} CL={r['CL']} CD={r['CD']}",flush=True)
    baseline=next(r for r in records if r['geometry']=='baseline')
    for r in records:
        r['delta_cd_counts']=(1e4*(r['CD']-baseline['CD']) if r['feasible'] and baseline['feasible'] else None)
    write_csv(out/'results.csv',records,list(dict.fromkeys(k for r in records for k in r)))
    summary=dict(requested=3,completed=len(records),feasible=sum(r['feasible'] for r in records),
                 reason=spec['reason'],records=records)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return 0 if summary['feasible']==3 else 2


if __name__=='__main__':
    raise SystemExit(main())
