#!/usr/bin/env python3
"""Independent Cp comparison at S702 flap 0, CL=1.04, Re=1e6, M=0.1."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path

from s702_transition_study import ROOT, run_job


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_pressure_check')
    ap.add_argument('--workers', type=int, default=2)
    args = ap.parse_args()
    if not 1 <= args.workers <= 2:
        ap.error('Use one or two workers')
    jobs = [dict(name=f'ncrit{n:g}', flap=0, target_cl=1.04, scenario='both',
                 ncrit_by_element=dict(fore=n, aft=n), bounds=[4.5, 6.]) for n in (9., 6.)]
    paths = ['multifoil.py', 'viscous.py', 'vendor/mfoil.py', 'paneling.py',
        'optimization.py', 'validation/s702_transition_study.py',
        'validation/s702_pressure_check.py', 'validation/s702_geometry.py',
        'validation/s702_benchmark.py', 'validation/extract_s702_pressure.py',
        'validation/data/s702/fore.dat', 'validation/data/s702/aft_d0.dat']
    paths += [str(p.relative_to(ROOT)) for p in
              sorted((ROOT/'validation/data/s702/pressure_cl1.04').glob('*')) if p.is_file()]
    spec = dict(schema=1, reynolds=1e6, mach=.1, panels=[260, 180], wake_length=1.,
        target_cl=1.04, cl_tolerance=2e-5, residual_tolerance=1e-7,
        iterations=100, start_alpha=-1., jobs=jobs, validation_claim=False,
        purpose='Additional physical observable; Cp is not used to choose Ncrit',
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
    with ProcessPoolExecutor(args.workers) as pool:
        for fut in as_completed([pool.submit(run_job, str(args.out), j, spec) for j in jobs]):
            r = fut.result()
            print({k: r[k] for k in ('name', 'feasible', 'CD', 'error')}, flush=True)
    rows = [json.loads((args.out/j['name']/'result.json').read_text()) for j in jobs]
    (args.out/'summary.json').write_text(json.dumps(rows, indent=2)+'\n')
    return 0 if all(r['feasible'] for r in rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
