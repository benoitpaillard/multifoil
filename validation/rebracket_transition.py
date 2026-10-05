#!/usr/bin/env python3
"""Retry the failed fine-grid baseline inside the original incidence domain.

The earlier attempts failed at 4 degrees, before reaching target lift. This
adaptive retry uses [5, 6] within the original [4, 7] domain; no physics changes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from validation.optimization_pilot import one_job


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=ROOT/'validation/results/transition_followup/study.json')
    p.add_argument('--out',type=Path,default=ROOT/'validation/results/transition_rebracket')
    a=p.parse_args()
    source=a.study.read_bytes();spec=json.loads(source)
    for f,h in spec['sha256'].items():
        if hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h:
            raise ValueError(f'Changed source code/input: {f}')
    if not spec['alpha_bounds'][0]<=5<6<=spec['alpha_bounds'][1]:
        raise ValueError('Retry must stay inside the original domain')
    spec['jobs']=[['baseline','ncrit6_fine']]
    spec['settings']={'ncrit6_fine':dict(spec['settings']['ncrit6_fine'],start_alpha=5.)}
    spec['alpha_bounds']=[5.,6.]
    spec['source_study_sha256']=hashlib.sha256(source).hexdigest()
    spec['reason']='Adaptive retry with [5,6] bracket after failure at the original lower endpoint 4; physical conditions unchanged.'
    spec['sha256']['validation/rebracket_transition.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    manifest=out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text())!=spec:
            raise ValueError('Changed study: use a new directory')
    elif any(out.iterdir()):
        raise ValueError('Nonempty directory lacks a study manifest')
    else:
        manifest.write_text(json.dumps(spec,indent=2)+'\n')
    r=one_job(str(out),'baseline','ncrit6_fine',spec)
    (out/'summary.json').write_text(json.dumps(dict(requested=1,completed=1,feasible=int(r['feasible']),record=r),indent=2)+'\n')
    print(json.dumps(r,indent=2))
    return 0 if r['feasible'] else 2


if __name__=='__main__':
    raise SystemExit(main())
