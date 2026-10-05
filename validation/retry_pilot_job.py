#!/usr/bin/env python3
"""Retry one declared pilot job with a different initialization in a new folder.

Physical settings are copied unchanged. The original result remains intact.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from validation.optimization_pilot import one_job
from validate import verify_data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True)
    p.add_argument('--geometry',required=True)
    p.add_argument('--setting',required=True)
    p.add_argument('--start-alpha',type=float,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    verify_data()
    source=a.study.read_bytes();spec=json.loads(source)
    if [a.geometry,a.setting] not in spec['jobs']:
        raise ValueError('Requested job is not part of the source study')
    for f,h in spec['sha256'].items():
        if hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h:
            raise ValueError(f'Changed source code/input: {f}')
    spec['jobs']=[[a.geometry,a.setting]]
    spec['settings']={a.setting:dict(spec['settings'][a.setting],start_alpha=a.start_alpha)}
    spec['reason']='Retry with different initialization only; original attempt retained separately.'
    spec['source_study_sha256']=hashlib.sha256(source).hexdigest()
    spec['sha256']['validation/retry_pilot_job.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    manifest=out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text())!=spec:
            raise ValueError('Changed retry study: use a new directory')
    elif any(out.iterdir()):
        raise ValueError('Nonempty directory lacks a study manifest')
    else:
        manifest.write_text(json.dumps(spec,indent=2)+'\n')
    r=one_job(str(out),a.geometry,a.setting,spec)
    (out/'summary.json').write_text(json.dumps(dict(requested=1,completed=1,feasible=int(r['feasible']),record=r),indent=2)+'\n')
    print(json.dumps(r,indent=2))
    return 0 if r['feasible'] else 2


if __name__=='__main__':
    raise SystemExit(main())
