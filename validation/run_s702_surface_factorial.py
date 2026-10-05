#!/usr/bin/env python3
"""Restartable fixed-source/fixed-wake fore/aft sampling study; five new targets."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import fcntl
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from run_s702_discretization import ROOT, VERSIONS, prepare, read, sha, stamp, write, worker, numerical_setup
import numpy as np

STUDY = 'validation/results/s702_surface_factorial_20261001'
PARENT = 'validation/results/s702_discretization_20261001'
AUDIT = 'validation/results/s702_discretization_analysis'
PARENT_HASH = '032ca093ea6cd2dff5007cd4177fb49c392fafc7b63896604529aa350119e9fd'
AUDIT_HASH = '2ec00402f775909582a22050fe07e22d507bb4b72270a8c37519581c90934c66'
INPUTS_HASH = '3bf58a151cc91f3fbf32ea7d71547c656c04fc562aa7cc9bb0dabcc7d4655d0f'
CELLS = {'B':[340,240], 'F':[440,240], 'A':[340,320], 'C':[440,320]}
COMPARISON = 'validation/results/s702_drag_mechanism/s702_flap_5_to_10.csv'


def specification():
    for name, digest in ((PARENT+'/study.json', PARENT_HASH), (AUDIT+'/summary.json', AUDIT_HASH),
                         (AUDIT+'/input_manifest.json', INPUTS_HASH)):
        if sha(ROOT/name) != digest:
            raise ValueError('Frozen parent/audit changed: '+name)
    prior = read(ROOT/PARENT/'study.json'); audit = read(ROOT/AUDIT/'summary.json')
    if not (audit['full_state_audit_pass'] and audit['factor_controls_pass']
            and audit['control_saved_arrays_exactly_equal'] and audit['endpoints_shared_arrays_exactly_equal']):
        raise ValueError('Parent full-state/control audit failed')
    for name, digest in prior['sha256'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('Frozen source/input changed: '+name)
    pinned = read(ROOT/AUDIT/'input_manifest.json')['files']
    references, jobs = [], []
    files = list(prior['sha256']) + [PARENT+'/study.json', PARENT+'/runtime.json',
        AUDIT+'/summary.json', AUDIT+'/input_manifest.json', AUDIT+'/FINDINGS.md',
        'validation/run_s702_surface_factorial.py', 'validation/S702_SURFACE_FACTORIAL_PROTOCOL.md',
        'test_s702_surface_factorial.py', COMPARISON]
    for flap in (5,10):
        endpoints = {}
        for cell, stage in (('B','contour'), ('C','surface')):
            original = next(j for j in prior['jobs'] if j['name'] == f'd{flap}_{stage}')
            folder = Path(PARENT)/original['name']; path = str(folder/'result.json')
            endpoints[cell] = path
            references.append(dict(flap=flap, cell=cell, result=path, job=original))
            for name, digest in pinned.items():
                if Path(name).is_relative_to(folder):
                    if sha(ROOT/name) != digest:
                        raise ValueError('Audited reference changed: '+name)
                    files.append(name)
        base = next(j for j in prior['jobs'] if j['name'] == f'd{flap}_surface')
        for cell, stage in (('F','fore_only'), ('A','aft_only')):
            jobs.append(dict(base, name=f'd{flap}_{stage}', stage=stage, cell=cell,
                panels=CELLS[cell], control=False, repeat_control=False,
                reference_result=endpoints['B'], surface_reference_result=endpoints['C']))
    base = next(j for j in prior['jobs'] if j['name'] == 'd10_surface')
    jobs.append(dict(base, name='d10_surface_control', cell='C', control=False, repeat_control=True,
        reference_result=PARENT+'/d10_surface/result.json', surface_reference_result=PARENT+'/d10_surface/result.json'))
    if sha(ROOT/COMPARISON) != pinned[COMPARISON]:
        raise ValueError('Experimental/MSES comparison context changed')
    inherited = {k:v for k,v in prior.items() if k not in ('jobs','sha256','stages','parent_study')}
    return dict(inherited, experiment='S702 fixed-source/wake fore/aft surface factorial',
        jobs=jobs, references=references, parent_study=PARENT, cells=CELLS,
        reference_comparison=COMPARISON, workers_limit=4,
        qualification='Conditional surface sampling effects and interaction; full audit pending. No independent validation or error bound.',
        sha256={name:sha(ROOT/name) for name in sorted(set(files))})


def effects(values):
    """Full factorial contrasts; refuse an incomplete/nonfinite four-cell design."""
    if set(values) != set(CELLS) or any(v is None or not np.isfinite(v) for v in values.values()):
        raise ValueError('Complete finite B/F/A/C values required')
    b, f, a, c = (values[k] for k in ('B','F','A','C'))
    answer = dict(fore_at_low_aft=f-b, fore_at_high_aft=c-a,
                  aft_at_low_fore=a-b, aft_at_high_fore=c-f,
                  interaction=c-f-a+b, total=c-b)
    if abs(answer['fore_at_low_aft']+answer['aft_at_low_fore']+answer['interaction']-answer['total']) > 1e-10:
        raise ValueError('Factor contrasts do not telescope')
    return answer


def summarize(out, spec):
    results = {j['name']:read(out/j['name']/'result.json') for j in spec['jobs']}
    states = {(j['flap'],j['cell']):results[j['name']] for j in spec['jobs'] if not j['repeat_control']}
    states.update({(ref['flap'],ref['cell']):read(ROOT/ref['result']) for ref in spec['references']})
    rows = []
    for job in spec['jobs']:
        r = results[job['name']]
        if not r['feasible'] and any(r[key] is not None for key in ('CD','CL','CM','alpha_deg')):
            raise ValueError('Failed target has finite force/angle')
        admissible = bool(r['feasible'] and all(r[k] is not None and np.isfinite(r[k]) for k in ('CD','CL','CM','alpha_deg')))
        if admissible:
            selected = r['records'][r['selected_flow_index']]
            admissible = bool(selected['admissible'] and selected['residual'] <= spec['residual_tolerance']
                              and abs(r['CL']-job['target_cl']) <= spec['cl_tolerance'])
        rows.append(dict(job=job['name'], cell=job['cell'], flap=job['flap'], repeat_control=job['repeat_control'],
            feasible=r['feasible'], force_tolerances_pass=admissible, CD=r['CD'], CL=r['CL'],
            selected_path=r.get('selected_path'), evaluations=len(r['records']),
            failed_evaluations=sum(not v['admissible'] for v in r['records']), wall_seconds=r['wall_seconds']))
    pairs = []
    for cell in CELLS:
        a, b = states[5,cell], states[10,cell]
        complete = bool(a['feasible'] and b['feasible'] and all(r['CD'] is not None and np.isfinite(r['CD']) for r in (a,b)))
        pairs.append(dict(cell=cell, panels=CELLS[cell], complete=complete,
            flap_delta_CD_counts=(b['CD']-a['CD'])*1e4 if complete else None))
    repeat = results['d10_surface_control']; original = states[10,'C']
    delta = (repeat['CD']-original['CD'])*1e4 if repeat['feasible'] and original['feasible'] else None
    exact = delta is not None
    if exact:
        for e in ('fore','aft'):
            with np.load(out/'d10_surface_control'/f'{e}_state.npz', allow_pickle=False) as a, np.load(ROOT/PARENT/'d10_surface'/f'{e}_state.npz', allow_pickle=False) as b:
                exact = exact and set(a.files)==set(b.files) and all(np.array_equal(a[k],b[k]) for k in a.files)
    coverage = bool(all(v['force_tolerances_pass'] for v in rows) and all(v['complete'] for v in pairs))
    control_pass = bool(delta is not None and abs(delta) <= spec['acceptance_control_CD_counts'] and exact)
    eligible = coverage and control_pass
    individual, paired = {}, None
    if eligible:
        for flap in (5,10):
            individual[str(flap)] = dict(total_drag=effects({cell:states[flap,cell]['CD']*1e4 for cell in CELLS}))
            for e in ('fore','aft'):
                individual[str(flap)][e+'_drag'] = effects({cell:states[flap,cell]['element_forces'][e]['CD_wake']*1e4 for cell in CELLS})
        paired = dict(total_drag=effects({v['cell']:v['flap_delta_CD_counts'] for v in pairs}))
        for e in ('fore','aft'):
            paired[e+'_drag'] = effects({cell:(states[10,cell]['element_forces'][e]['CD_wake']-states[5,cell]['element_forces'][e]['CD_wake'])*1e4 for cell in CELLS})
        for term in paired['total_drag']:
            if abs(paired['total_drag'][term]-(individual['10']['total_drag'][term]-individual['5']['total_drag'][term])) > 1e-10:
                raise ValueError('Paired/individual factorial arithmetic differs')
    with (ROOT/spec['reference_comparison']).open() as stream:
        context = next(v for v in csv.DictReader(stream) if float(v['reynolds'])==spec['reynolds']
                       and float(v['mach'])==spec['mach'] and float(v['target_cl'])==.8)
    summary = dict(completed_utc=stamp(), required_jobs=5, feasible_jobs=sum(v['feasible'] for v in rows),
        required_factorial_cells=8, coverage_complete=coverage, control_pass=control_pass,
        control_CD_difference_counts=delta, control_saved_arrays_exactly_equal=bool(exact),
        execution_checks_pass=eligible, jobs=rows, paired_cells=pairs,
        individual_effects_counts=individual, paired_effects_counts=paired,
        reference_context=dict(experiment_delta_counts=float(context['experiment_delta_10_minus_5_counts']),
            published_mses_ncrit9_delta_counts=float(context['mses_n9_delta_10_minus_5_counts']),
            multifoil_ncrit=spec['ncrit'], mses_is_published_curve=True,
            transition_settings_matched=False, source=spec['reference_comparison']),
        qualification='Provisional contrasts only; full selected-state/source/geometry audit required. Experiment and MSES are context, not acceptance targets.')
    write(out/'summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/STUDY)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 4:
        parser.error('Use a results directory and 1–4 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out,spec)
        if args.prepare_only:
            print(json.dumps(dict(state='prepared',jobs=len(spec['jobs']),study_sha256=sha(out/'study.json')))); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k:importlib.metadata.version(k) for k in VERSIONS},
            thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        if runtime['versions'] != VERSIONS or sys.version_info[:2] != (3,12):
            raise ValueError('Use pinned Python 3.12 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json',runtime)
        done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state,**extra):
            write(out/'status.json',dict(state=state,updated_utc=stamp(),pid=os.getpid(),completed_jobs=done,total_jobs=5,**extra))
        status('running')
        try:
            pending = []
            for job in spec['jobs']:
                folder = out/job['name']
                if (folder/'result.json').exists(): continue
                if folder.exists():
                    archived = out/'interrupted_attempts'/stamp().replace(':','')/job['name']
                    archived.parent.mkdir(parents=True,exist_ok=True); shutil.move(folder,archived)
                pending.append(job)
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(worker,out,j,spec):j for j in pending}
                for future in as_completed(futures):
                    r=future.result(); done+=1
                    status('running',last_job=futures[future]['name'],last_job_feasible=r['feasible'])
                    print(f"{done}/5 {futures[future]['name']} feasible={r['feasible']}",flush=True)
            summary=summarize(out,spec); status('finished',execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed',error=f'{type(exc).__name__}: {exc}'); raise


if __name__ == '__main__':
    raise SystemExit(main())
