#!/usr/bin/env python3
"""Check the N=9 fine-mesh target without the wide [-2,1] incidence excursion.

The [-0.8,-0.2] bracket surrounds the previously computed base-mesh target
alpha=-0.491 degrees. It is not chosen from experimental drag. This calculation
is additional to the retained wide-bracket attempt; it never overwrites it.
"""
import argparse
import hashlib
import json
from pathlib import Path

from s702_transition_study import ROOT, run_job


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_transition_fine_local')
    args = ap.parse_args()
    job = dict(name='d+10_cl0.8_local', flap=10, target_cl=.8, scenario='both',
               ncrit_by_element=dict(fore=9., aft=9.), bounds=[-.8, -.2])
    paths = ['multifoil.py', 'viscous.py', 'vendor/mfoil.py', 'paneling.py',
        'optimization.py', 'validation/s702_transition_study.py',
        'validation/s702_local_bracket_check.py', 'validation/s702_geometry.py',
        'validation/s702_benchmark.py', 'validation/data/s702/fore.dat',
        'validation/data/s702/aft_d10.dat',
        'validation/results/s702_transition_fine_n9/study.json',
        'validation/results/s702_mach_check/flap10_mach0.1/result.json']
    spec = dict(schema=1, reynolds=1e6, mach=.1, panels=[340, 240], wake_length=1.,
        target_cl=.8, cl_tolerance=2e-5, residual_tolerance=1e-7,
        iterations=100, start_alpha=-1., jobs=[job], validation_claim=False,
        purpose='Additional continuation-path check around the pre-existing numerical target',
        sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out/'study.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != spec:
            raise ValueError('Study/source changed; select a new output directory')
    elif any(args.out.iterdir()):
        raise ValueError('Nonempty output directory lacks a manifest')
    else:
        manifest.write_text(json.dumps(spec, indent=2)+'\n')
    result = run_job(str(args.out), job, spec)
    (args.out/'summary.json').write_text(json.dumps([result], indent=2)+'\n')
    print({k: result[k] for k in ('name', 'feasible', 'CD', 'error')}, flush=True)
    return 0 if result['feasible'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
