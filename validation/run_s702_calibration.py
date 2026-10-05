#!/usr/bin/env python3
"""Reproduce the complete S702 sweep and frozen-parameter validation study.

Run from the distributed project, which includes the original N=6/9 reference
runs and extracted experimental data. Choose a fresh output directory inside
the project. Numerical failures are retained; setup/report errors stop the run.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/reproduced_s702_calibration')
    parser.add_argument('--workers',type=int,default=8)
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(ROOT) or not 1 <= args.workers <= 8:
        parser.error('Use an output directory inside this project and 1-8 workers')
    if out.exists() and any(out.iterdir()):
        parser.error('Choose a fresh output directory; see resume_s702_sweeps.py for interrupted studies')
    out.mkdir(parents=True,exist_ok=True)
    commands = []

    def invoke(script, flags, numerical=False):
        command = [sys.executable,str(ROOT/'validation'/script)]+[str(v) for v in flags]
        commands.append(command)
        (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        print('Running '+script+' '+str(flags),flush=True)
        result = subprocess.run(command,cwd=ROOT)
        if result.returncode not in ((0,2) if numerical else (0,)):
            raise RuntimeError(f'{script} failed with exit status {result.returncode}')

    def sweep(name,ncrit,extra=()):
        folder = out/name
        invoke('s702_ncrit_sweep.py',['--out',folder,'--workers',args.workers,
            '--ncrit',ncrit]+list(extra),numerical=True)
        return folder

    coarse = sweep('coarse','3,4,5,6,7,8,9,10,11,12')
    refined = [sweep('refined','5.5,5.75,6.25,6.5')]
    analysis = out/'analysis'

    def report():
        invoke('report_s702_ncrit_sweep.py',['--coarse',coarse,'--refined',*refined,'--out',analysis])
        return read(analysis/'summary.json')

    summary = report()
    if summary['best']['ncrit'] == 6.5:
        refined.append(sweep('refined_upper','6.75'))
        summary = report()
    selected = summary['best']['ncrit']
    frozen = dict(selected_ncrit=selected,
        frozen_utc=datetime.now(timezone.utc).isoformat(),
        criterion='Equal-point MAPE on 12 Re=1e6, M=.1 cases; common Ncrit on both elements',
        coarse_grid=summary['coarse_ncrit_values'],refined_grid=summary['refined_ncrit_values'],
        calibration_points_sha256=digest(analysis/'points.csv'),
        calibration_summary_sha256=digest(analysis/'summary.json'),
        holdout_reference_sha256=digest(ROOT/'validation/data/s702/re700k_holdout/source.json'),
        holdout_used_in_selection=False,
        qualification='Effective S702 calibration; no claim of measured disturbance level or validated shape derivatives')
    (analysis/'frozen_selection.json').write_text(json.dumps(frozen,indent=2)+'\n')
    holdout_flags = ['--flaps','0,5,10','--reynolds','700000','--mach','.07','--no-reuse']
    h9 = sweep('holdout_n9','9',holdout_flags)
    hs = sweep('holdout_selected',f'{selected:g}',holdout_flags)
    from s702_ncrit_sweep import local_seed
    retries = []
    for source in sorted(hs.glob('n*/cl*/result.json')):
        result = read(source)
        if result['feasible']:
            continue
        try:
            local_seed(result)
        except ValueError as exc:
            print(f'Local retry unavailable for {source}: {exc}',flush=True)
            continue
        flags = ['--flaps',str(result['flap']),'--targets',str(result['target_cl']),
                 '--reynolds','700000','--mach','.07','--no-reuse','--seed-from',source]
        cold = sweep(f'holdout_local_d{result["flap"]}_cl{result["target_cl"]:g}',f'{selected:g}',flags)
        retries.append(cold)
        cold_target = read(next(cold.glob('n*/cl*/result.json')))
        if not cold_target['feasible'] and result['records'][0]['admissible']:
            retries.append(sweep(f'holdout_warm_d{result["flap"]}_cl{result["target_cl"]:g}',
                                 f'{selected:g}',flags+['--warm-seed']))
    control = sweep('control_n6','6',['--no-reuse'])
    fine = [sweep('fine_n6','6',['--panels','340,240','--no-reuse'])]
    if selected != 6.:
        fine.append(sweep('fine_selected',f'{selected:g}',['--panels','340,240','--no-reuse']))
    mach = sweep('mach_zero',f'{selected:g}',['--flaps','5,10','--targets','.8','--mach','0','--no-reuse'])
    pressure = sweep('pressure',f'{selected:g}',['--flaps','0','--targets','1.04','--no-reuse'])
    invoke('check_s702_ncrit_transfer.py',['--analysis',analysis,'--holdout',h9,hs,
        '--fine',*fine,'--mach0',mach,'--pressure',pressure,'--control',control,'--retry',*retries])
    invoke('write_s702_ncrit_report.py',['--analysis',analysis])
    print('Completed report: '+str(analysis/'NCRIT_SWEEP_RESULTS.md'),flush=True)


if __name__ == '__main__':
    main()
