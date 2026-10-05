#!/usr/bin/env python3
"""V2 tranche 2: refinement and cold checks through the S414 transition-boundary interval.

Eighteen independent evaluations. Continuation starts from audited saved states
restored into a new solver; two restore controls must first reproduce tranche-1
live results within the declared screens. Protocol:
validation/S414_TRANSITION_BOUNDARY_T2_PROTOCOL.md.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import fcntl
import importlib.metadata
import os
from pathlib import Path
import platform
import shutil
import sys
import time
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_s414_transition_boundary as deck
from run_s414_transition_boundary import read, sha, stamp, write, finite, prepare, VERSIONS, PYTHON, THREADS
from multifoil import save_results
from optimization import admissible
from viscous import CoupledSolver, mf, te_info

STUDY = 'validation/results/s414_transition_boundary_t2_20261004'
T1 = deck.STUDY
T1_RECEIPT = 'validation/results/s414_transition_boundary_launch_20261004/retrieved_files_sha256.json'
ANALYSIS = 'validation/results/s414_transition_boundary_analysis_20261004'
ANCHOR = 'validation/results/s414_incidence_20261003'
ANCHOR_JOB = 'cl0056_measured_alpha'
ANCHOR_RECEIPT = 'validation/results/s414_incidence_launch_20261003/retrieved_files_sha256.json'
ANCHOR_AUDIT = 'validation/results/s414_incidence_analysis_20261003/residual_recheck.json'
PROTOCOL = 'validation/S414_TRANSITION_BOUNDARY_T2_PROTOCOL.md'
TRANSITION_CHANGE = .05
MODE = 'corrected_af'


def transition_x(record):
    tr = record['boundary_layers']['fore']['lower']['transition']
    return None if tr is None else tr['xy_over_cref'][0]


def tranche1_base(scenario, direction, alpha, files):
    """Nearest audited, converged tranche-1 state on that direction's own path, on the upstream side."""
    job = f'{scenario}_corrected_{direction}'; result = read(ROOT/T1/job/'result.json')
    side = [r for r in result['records'] if r['role'] == 'path' and r['converged']
            and (r['requested_alpha_deg'] < alpha if direction == 'increasing' else r['requested_alpha_deg'] > alpha)]
    if not side:
        raise ValueError('No converged tranche-1 state on the required side')
    rec = (max if direction == 'increasing' else min)(side, key=lambda r: r['requested_alpha_deg'])
    names = {e: f'{T1}/{job}/states/point_{rec["index"]:02d}_{e}.npz' for e in deck.ELEMENTS}
    files.update(names.values()); files.add(f'{T1}/{job}/result.json')
    return dict(study=T1, job=job, index=rec['index'], alpha_deg=rec['requested_alpha_deg'], state_files=names,
                fore_lower_transition_x=transition_x(rec))


def specification():
    verification = read(ROOT/ANALYSIS/'verification.json'); audit = read(ROOT/ANALYSIS/'residual_recheck.json')
    if not (verification['complete'] and audit['complete'] and len(audit['states']) == audit['saved_states'] == 25):
        raise ValueError('Tranche 1 is not completely audited; tranche 2 is not admitted')
    # The two analysis outputs this deck depends on must be the verified ones.
    for name in ('residual_recheck.json', 'tranche2_selection.json'):
        if sha(ROOT/ANALYSIS/name) != verification['output_sha256'][name]:
            raise ValueError('Audited tranche-1 analysis changed: '+name)
    t1 = read(ROOT/T1/'study.json'); angles = t1['angles_deg']; selection = read(ROOT/ANALYSIS/'tranche2_selection.json')
    files = {PROTOCOL, T1+'/study.json', T1_RECEIPT, ANCHOR_RECEIPT, ANCHOR_AUDIT, ANCHOR+'/study.json',
             ANALYSIS+'/verification.json', ANALYSIS+'/residual_recheck.json', ANALYSIS+'/tranche2_selection.json',
             'vendor/mfoil.py', 'vendor/__init__.py', deck.QUAL+'/report_test_fore.dat', deck.QUAL+'/report_test_aft.dat',
             deck.QUAL+'/verification.json'} | set(deck.TEST_DATA)
    jobs = []
    common = dict(amplification_jacobian=MODE, wake_length=1.)
    # Restore controls: replay one audited tranche-1 continuation step from a restored state.
    for scenario in deck.ARCHIVED_JOBS:
        base = tranche1_base(scenario, 'increasing', angles[1], files)
        if base['alpha_deg'] != angles[0]:
            raise ValueError('Restore control must start at the lower bracket end')
        reference = dict(study=T1, job=f'{scenario}_corrected_increasing', index=1, alpha_deg=angles[1])
        files.update(f'{T1}/{reference["job"]}/flow/surface_001_{e}.csv' for e in deck.ELEMENTS)
        jobs.append(dict(common, name=scenario+'_restore_control', scenario=scenario, kind='restore_control',
            direction='increasing', alpha_deg=angles[1], base=base, reference=reference))
    # Cubic geometry: the four incidences selected by the declared rule from audited tranche 1.
    chosen = selection['fore_source_cubic']
    if len(chosen['new_locations']) != 4 or chosen['fallback'] or selection['nominal']['new_locations'] or not selection['nominal']['fallback']:
        raise ValueError('Tranche-1 selection differs from the frozen protocol')
    cubic = sorted(l['alpha_deg'] for l in chosen['new_locations'])
    for k, alpha in enumerate(cubic):
        for direction in ('increasing', 'decreasing'):
            jobs.append(dict(common, name=f'fore_source_cubic_{direction}_{k}', scenario='fore_source_cubic', kind='continuation',
                direction=direction, alpha_deg=alpha, base=tranche1_base('fore_source_cubic', direction, alpha, files)))
    # Nominal geometry: declared fallback to the bracket from the archived history.
    anchor_spec = read(ROOT/ANCHOR/'study.json'); anchor_job = next(j for j in anchor_spec['jobs'] if j['name'] == ANCHOR_JOB)
    anchor = read(ROOT/ANCHOR/ANCHOR_JOB/'result.json'); receipt = read(ROOT/ANCHOR_RECEIPT)
    recheck = next(s for s in read(ROOT/ANCHOR_AUDIT)['states'] if s['job'] == ANCHOR_JOB)
    for key in ('reynolds', 'mach', 'ncrit', 'panels', 'iterations', 'residual_tolerance'):
        if anchor_spec[key] != t1[key]:
            raise ValueError('Archived anchor uses different settings: '+key)
    if not (anchor['feasible'] and anchor_job['scenario'] == 'nominal' and anchor_job['amplification_jacobian'] == MODE
            and anchor_job['wake_length'] == 1. and recheck['max_residual'] <= t1['residual_tolerance']):
        raise ValueError('Archived anchor state is not an audited corrected nominal state')
    names = {e: f'{ANCHOR}/{ANCHOR_JOB}/{e}_state.npz' for e in deck.ELEMENTS}
    for name in [*names.values(), f'{ANCHOR}/{ANCHOR_JOB}/result.json']:
        if sha(ROOT/name) != receipt[name[len(ANCHOR)+1:]]:
            raise ValueError('Archived anchor changed: '+name)
        files.add(name)
    lower = dict(study=ANCHOR, job=ANCHOR_JOB, index=anchor['selected_flow_index'], alpha_deg=anchor['alpha_deg'],
                 state_files=names, fore_lower_transition_x=transition_x(anchor))
    upper = tranche1_base('nominal', 'decreasing', angles[0]-1e-9, files)
    if upper['alpha_deg'] != angles[0] or not lower['alpha_deg'] < upper['alpha_deg']:
        raise ValueError('Nominal fallback bracket is not ordered')
    change = abs(upper['fore_lower_transition_x']-lower['fore_lower_transition_x'])
    if change <= TRANSITION_CHANGE:
        raise ValueError('Archived nominal bracket does not contain the transition change')
    width = upper['alpha_deg']-lower['alpha_deg']
    # One eligible interval: its midpoint, then its quarter points (lower first), as the rule states.
    nominal = sorted(lower['alpha_deg']+f*width for f in (.5, .25, .75))
    for k, alpha in enumerate(nominal):
        for direction, base in (('increasing', lower), ('decreasing', upper)):
            jobs.append(dict(common, name=f'nominal_{direction}_{k}', scenario='nominal', kind='continuation',
                direction=direction, alpha_deg=alpha, base=base))
    # Cold checks. Cubic: approved deviation from the declared angle; nominal: midpoint of its bracket.
    jobs.append(dict(common, name='fore_source_cubic_cold_check', scenario='fore_source_cubic', kind='cold_check',
        direction='cold', alpha_deg=angles[1], reference=dict(study=T1, job='fore_source_cubic_corrected_increasing', index=1, alpha_deg=angles[1]),
        deviation=dict(declared_alpha_deg=chosen['cold_check']['alpha_deg'], approved_by='user, 4 October 2026',
                       reason='Three cold starts already failed at the declared incidence in tranche 1')))
    jobs.append(dict(common, name='nominal_cold_check', scenario='nominal', kind='cold_check', direction='cold',
        alpha_deg=lower['alpha_deg']+.5*width))
    files.update(f'{T1}/fore_source_cubic_corrected_increasing/flow/surface_001_{e}.csv' for e in deck.ELEMENTS)
    receipt1 = read(ROOT/T1_RECEIPT)
    for name in sorted(files):
        if name.startswith(T1+'/') and sha(ROOT/name) != receipt1[name[len(T1)+1:]]:
            raise ValueError('Audited tranche-1 file changed: '+name)
    for folder in (ROOT, ROOT/'validation'):
        files.update(str(p.relative_to(ROOT)) for p in folder.glob('*.py'))
    settings = {k: t1[k] for k in ('reynolds', 'mach', 'ncrit', 'panels', 'reference_chord', 'model_chord_mm', 'aft_shift_c',
                                   'iterations', 'residual_tolerance', 'continuation_step', 'min_step', 'minimum_te_gap',
                                   'equivalence_screens', 'expected_versions', 'expected_python')}
    return dict(schema=1, experiment='S414 fixed-incidence transition-boundary paths, V2 tranche 2', **settings,
        tranche1_angles_deg=angles, cubic_new_incidences_deg=cubic, nominal_new_incidences_deg=nominal,
        nominal_fallback_bracket=dict(lower_alpha_deg=lower['alpha_deg'], upper_alpha_deg=upper['alpha_deg'],
            fore_lower_transition_change=change), jobs=jobs, workers_limit=2,
        requested_points=len(jobs), tranche1_requested_points=24, allowance=42,
        gate='Both restore controls reproduce the tranche-1 live result within the declared screens before any new evaluation starts',
        qualification='Numerical path diagnostic at fixed incidence; not matched-lift drag validation, not physical '
                      'hysteresis, not solution uniqueness; audit required. V2 ends with this tranche',
        sha256={name: sha(ROOT/name) for name in sorted(files)})


def restore(solver, alpha, arrays, tolerance):
    """Put an audited saved state into a new solver and prove it is the same converged state."""
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))]); worst = dict(geometry=0., wake=0., stations=0.)
    def same(a, b, limit, key, label):
        difference = float(np.max(abs(a-b))) if a.shape == b.shape else np.inf
        worst[key] = max(worst[key], difference)
        if difference > limit:
            raise ValueError('Restored base state differs: '+label)
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = arrays[model.geom.name]
        same(model.foil.x, a['xy'], 2e-12, 'geometry', 'foil nodes')
        for field, saved in ((model.wake.x, a['wake_xy']), (model.wake.t, a['wake_t']), (model.wake.s, a['wake_s'])):
            same(field, saved, 2e-11, 'wake', 'inviscid wake')
        model.glob.U = a['U'].copy(); model.glob.Nsys = a['U'].shape[1]
        model.vsol.turb = a['turbulent'].copy()
        model.isol.Istag = [int(a['lower'][0]), int(a['upper'][0])]
        model.isol.sgnue = np.ones(model.foil.N); model.isol.sgnue[a['lower']] = -1
        model.isol.sstag = float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model); mf.stagpoint_move(model)
        same(model.isol.xi, a['xi'], 2e-12, 'stations', 'station coordinates')
        if not all(np.array_equal(model.vsol.Is[k], a[label]) for k, label in enumerate(('lower', 'upper'))):
            raise ValueError('Restored base state differs: surface indices')
        if abs(te_info(model.foil.x)[1]) < 1e-12: model.vsol.wgap = np.zeros(model.wake.N)
        else: mf.set_wake_gap(model)
        model.initialization = 'restored audited saved state'
    solver._build_source_coupling()
    residual, _ = solver._system()
    blnorm = float(np.linalg.norm(residual[:3*solver.nsys])); coupling = float(np.max(abs(residual[3*solver.nsys:])))
    if not max(blnorm, coupling) <= tolerance:
        raise ValueError('Restored base state fails the residual tolerance')
    for model in solver.models:
        if not np.array_equal(model.glob.U, arrays[model.geom.name]['U']):
            raise ValueError('Restore changed the saved solution')
    solver.initialized = True; solver.last_alpha = float(alpha); solver.last_solution = None
    return dict(BL_residual_l2=blnorm, coupling_residual_max=coupling, maximum_geometry_difference=worst['geometry'],
                maximum_wake_difference=worst['wake'], maximum_station_difference=worst['stations'])


def worker(out, job, spec):
    folder = out/job['name']; folder.mkdir(exist_ok=True)
    with (folder/'run.log').open('a', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        begun = time.perf_counter(); alpha = job['alpha_deg']; cold = job['kind'] == 'cold_check'
        elements = deck.geometry_elements(job['scenario'], spec['panels'])
        solver = CoupledSolver(elements, reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
            iterations=spec['iterations'], tolerance=spec['residual_tolerance'], wake_length=job['wake_length'],
            minimum_te_gap=spec['minimum_te_gap'], continuation_step=spec['continuation_step'], min_step=spec['min_step'],
            start_alpha=alpha if cold else None, amplification_jacobian=job['amplification_jacobian'])
        if solver.amplification_jacobian != MODE:
            raise ValueError('Solver mode differs from declared job')
        result = dict(job=job['name'], scenario=job['scenario'], kind=job['kind'], direction=job['direction'],
            amplification_jacobian=MODE, ncrit=spec['ncrit'], complete=False, restore=None, records=[], error='',
            flow_directory=str((folder/'flow').relative_to(ROOT)))
        (folder/'states').mkdir(exist_ok=True)
        try:
            if not cold:
                arrays = {}
                for element, name in job['base']['state_files'].items():
                    with np.load(ROOT/name, allow_pickle=False) as saved:
                        arrays[element] = {k: saved[k].copy() for k in saved.files}
                result['restore'] = restore(solver, job['base']['alpha_deg'], arrays, spec['residual_tolerance'])
        except ValueError as exc:
            # No evaluation is attempted from a base state that cannot be proven identical.
            result.update(error=f'Base state not admitted: {exc}', complete=True, wall_seconds=time.perf_counter()-begun)
            write(folder/'result.json', result); return result
        sol = solver.solve(alpha)
        sol.diagnostics['ncrit_by_element'] = dict(fore=spec['ncrit'], aft=spec['ncrit']); diag = sol.diagnostics
        trail = [dict(alpha_deg=t['alpha_deg'], converged=bool(t['converged']), CL=finite(t['CL'])) for t in diag['continuation']]
        entry = dict(index=0, role=job['kind'], requested_alpha_deg=alpha, alpha_deg=sol.alpha, CL=finite(sol.cl), CD=finite(sol.cd),
            CM=finite(sol.cm), converged=bool(sol.converged), admissible=bool(admissible(sol)), residual=finite(sol.residual),
            failure_reason=diag.get('failure_reason', ''), iterations_final_step=diag['iterations'],
            reached_by=(deck.reached_by(False, None, alpha, trail) if cold else
                        f"continuation from restored audited state {job['base']['job']}[{job['base']['index']}] "
                        f"at incidence {job['base']['alpha_deg']!r}"),
            approach='cold' if cold else job['direction'], continuation=trail, wall_seconds=diag['wall_seconds'], state_files=[])
        if sol.converged:
            if abs(solver.last_alpha-alpha) > 1e-12 or solver.last_solution.cl != sol.cl:
                raise ValueError('Solver state does not match the returned solution')
            entry.update(boundary_layers=diag['boundary_layers'], element_forces=diag['element_forces'],
                transition_audit=deck.transition_audit(solver), wake_discretization={})
            stream = np.array([np.cos(np.deg2rad(sol.alpha)), np.sin(np.deg2rad(sol.alpha))])
            for model in solver.models:
                name = f'states/point_00_{model.geom.name}.npz'
                np.savez_compressed(folder/name, U=model.glob.U, xy=model.foil.x, xi=model.isol.xi, turbulent=model.vsol.turb,
                    wake_xy=model.wake.x, lower=np.asarray(model.vsol.Is[0]), upper=np.asarray(model.vsol.Is[1]),
                    wake_s=model.wake.s, wake_t=model.wake.t)
                entry['state_files'].append(name)
                spacing = np.linalg.norm(np.diff(model.wake.x, axis=1), axis=0)
                entry['wake_discretization'][model.geom.name] = dict(nodes=int(model.wake.N), arclength=float(np.sum(spacing)),
                    first_spacing=float(spacing[0]), last_spacing=float(spacing[-1]),
                    streamwise_start=float(stream@model.wake.x[:, 0]), streamwise_end=float(stream@model.wake.x[:, -1]))
        elif any(entry[k] is not None for k in ('CL', 'CD', 'CM')):
            raise ValueError('Failed evaluation carries a finite force')
        result['records'].append(entry)
        save_results(folder/'flow', [sol], elements, dict(study=spec, job=job), 'viscous')
        print(f"alpha={alpha!r} converged={sol.converged} CL={entry['CL']} CD={entry['CD']} {entry['reached_by']}", flush=True)
        result.update(complete=True, wall_seconds=time.perf_counter()-begun)
        write(folder/'result.json', result)
        return result


def surface(path):
    return np.genfromtxt(path, delimiter=',', names=True, dtype=None, encoding='utf8')


def compare_reference(folder, job, screens):
    """Declared screens between this evaluation and the tranche-1 live state at the same incidence."""
    result = read(folder/'result.json'); ref = job['reference']
    old = next(r for r in read(ROOT/ref['study']/ref['job']/'result.json')['records'] if r['index'] == ref['index'])
    if not result['records'] or not result['records'][0]['converged'] or not old['converged']:
        return dict(passed=False, reason=result['error'] or 'not converged', reference=ref)
    new = result['records'][0]; cp = 0.; shift = 0.; topology = True
    for element in deck.ELEMENTS:
        a = surface(folder/f'flow/surface_000_{element}.csv'); b = surface(ROOT/ref['study']/ref['job']/f'flow/surface_{ref["index"]:03d}_{element}.csv')
        cp = max(cp, float(np.max(abs(a['Cp']-b['Cp'])))); topology = topology and bool(np.array_equal(a['turbulent'], b['turbulent']))
        for side in ('upper', 'lower'):
            la, lb = new['boundary_layers'][element][side], old['boundary_layers'][element][side]
            ta, tb = la['transition'], lb['transition']
            topology = topology and (ta is None) == (tb is None) and len(la['separation_regions']) == len(lb['separation_regions'])
            if ta is not None and tb is not None: shift = max(shift, abs(ta['xy_over_cref'][0]-tb['xy_over_cref'][0]))
    out = dict(reference=ref, delta_CL=new['CL']-old['CL'], delta_CD_counts=(new['CD']-old['CD'])*1e4, max_abs_delta_Cp=cp,
               max_transition_shift_x=shift, same_topology=bool(topology))
    out['passed'] = bool(abs(out['delta_CL']) <= screens['CL'] and abs(out['delta_CD_counts']) <= screens['CD_counts']
                         and cp <= screens['Cp'] and shift <= screens['transition_x_over_c'] and topology)
    return out


def controls(out, spec):
    result = {j['name']: compare_reference(out/j['name'], j, spec['equivalence_screens']) for j in spec['jobs'] if j['kind'] == 'restore_control'}
    return dict(passed=all(r['passed'] for r in result.values()), jobs=result)


def summarize(out, spec):
    rows = []
    for job in spec['jobs']:
        r = read(out/job['name']/'result.json'); rec = r['records'][0] if r['records'] else None; ok = True
        if rec:
            meta = read(out/job['name']/'flow/metadata.json')
            if meta['case']['study'] != spec or meta['case']['job'] != job:
                raise ValueError('Exported study/job differs')
            diag = meta['diagnostics'][0]
            if diag['amplification_jacobian'] != MODE or diag['ncrit_by_element'] != {'fore': 9., 'aft': 9.}:
                raise ValueError('Exported mode or Ncrit differs')
            if rec['converged']:
                ok = bool(rec['admissible'] and rec['residual'] <= spec['residual_tolerance']
                    and abs(rec['alpha_deg']-job['alpha_deg']) <= 1e-12
                    and abs(sum(v['CD_wake'] for v in rec['element_forces'].values())-rec['CD']) < 1e-12
                    and len(rec['state_files']) == 2 and all((out/job['name']/n).exists() for n in rec['state_files'])
                    and all(t['inside_cell'] and abs(t['amplification_residual']) <= 1e-7 and t['station_residual_max'] <= 1e-7
                            for t in rec['transition_audit']))
            elif any(rec[k] is not None for k in ('CL', 'CD', 'CM')):
                raise ValueError('Failed evaluation contains a finite force')
        rows.append(dict(job=job['name'], scenario=job['scenario'], kind=job['kind'], direction=job['direction'],
            alpha_deg=job['alpha_deg'], evaluated=rec is not None, converged=bool(rec and rec['converged']),
            execution_check_pass=ok, CL=rec['CL'] if rec else None, CD=rec['CD'] if rec else None,
            reached_by=rec['reached_by'] if rec else 'not evaluated', iterations_final_step=rec['iterations_final_step'] if rec else None,
            failure_reason=(rec['failure_reason'] if rec else r['error'])))
    control = controls(out, spec); new = [r for r in rows if r['kind'] != 'restore_control']
    result = dict(required_jobs=len(spec['jobs']), completed_jobs=len(rows),
        requested_path_points=len(new), converged_path_points=sum(r['converged'] for r in new),
        failed_path_points=sum(not r['converged'] for r in new), control_points=len(rows)-len(new),
        execution_checks_pass=bool(control['passed'] and all(r['execution_check_pass'] for r in rows)),
        controls=control, points=rows,
        qualification='Execution only; nonconverged points are findings, not forces. Full audit and interpretation pending')
    write(out/'summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/STUDY)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(ROOT/'validation/results') or not 1 <= args.workers <= 2:
        parser.error('Use a results directory and 1–2 workers')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        spec = specification(); prepare(out, spec)
        if args.prepare_only:
            print(f'Prepared {len(spec["jobs"])} tranche-2 evaluations; study sha256='+sha(out/'study.json')); return 0
        runtime = dict(python=platform.python_version(), hostname=platform.node(), executable=sys.executable,
            versions={k: importlib.metadata.version(k) for k in VERSIONS}, thread_limits={k: os.environ[k] for k in THREADS})
        if runtime['versions'] != VERSIONS or runtime['python'] != PYTHON:
            raise ValueError('Use pinned Python 3.12.13 numerical environment')
        if (out/'runtime.json').exists() and read(out/'runtime.json') != runtime:
            raise ValueError('Runtime changed; use a new study directory')
        write(out/'runtime.json', runtime)
        total = len(spec['jobs']); done = sum((out/j['name']/'result.json').exists() for j in spec['jobs'])
        def status(state, **extra):
            write(out/'status.json', dict(state=state, updated_utc=stamp(), pid=os.getpid(), completed_jobs=done, total_jobs=total, **extra))
        try:
            for phase in ('restore_controls', 'new_evaluations'):
                if phase == 'new_evaluations':
                    control = controls(out, spec); write(out/'controls.json', control)
                    if not control['passed']:
                        status('execution_failed', error='Restore controls differ from tranche-1 live results; new evaluations not launched'); return 2
                pending = []
                for job in spec['jobs']:
                    if (job['kind'] == 'restore_control') != (phase == 'restore_controls'): continue
                    folder = out/job['name']
                    if (folder/'result.json').exists(): continue
                    if folder.exists():
                        dest = out/'interrupted_attempts'/stamp().replace(':', '')/job['name']
                        dest.parent.mkdir(parents=True, exist_ok=True); shutil.move(folder, dest)
                    pending.append(job)
                status('running', phase=phase)
                with ProcessPoolExecutor(args.workers) as pool:
                    futures = {pool.submit(worker, out, j, spec): j for j in pending}
                    for future in as_completed(futures):
                        result = future.result(); done += 1
                        converged = bool(result['records'] and result['records'][0]['converged'])
                        status('running', phase=phase, last_job=futures[future]['name'], last_job_converged=converged)
                        print(f'{done}/{total}: {futures[future]["name"]} converged={converged}', flush=True)
            summary = summarize(out, spec); status('finished', execution_checks_pass=summary['execution_checks_pass'])
            return 0 if summary['execution_checks_pass'] else 2
        except Exception as exc:
            status('execution_failed', error=f'{type(exc).__name__}: {exc}'); raise


if __name__ == '__main__':
    raise SystemExit(main())
