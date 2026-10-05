#!/usr/bin/env python3
"""Independent audit of every saved V2 tranche-2 state; no nonlinear flow solve.

Each converged evaluation gets the export/array consistency check and a fresh
residual, force, Cp/Cf, station, transition-root and wake-drag evaluation in a
newly built solver (the tranche-1 routine, unchanged). Base states, restore
records and the restore-control gate are re-derived from the retrieved files.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, Inputs, sha, write
from report_s702_drag_mechanism import require, close, audit_surface, side_data
from report_s702_wake_extent import thermo, array_close
from viscous import mf
from audit_s414_transition_boundary import recompute
import run_s414_transition_boundary as deck
import run_s414_transition_boundary_t2 as t2

STUDY = ROOT/t2.STUDY
LAUNCH = ROOT/'validation/results/s414_transition_boundary_t2_launch_20261004'
OUT = ROOT/'validation/results/s414_transition_boundary_t2_analysis_20261004'
STUDY_HASH = 'e0d69a92e4bc755d0dfa319bbfd41eec56657ef3e350a2ae5b81f3cc8b35ac18'
AUDIT_CODE = ['validation/audit_s414_transition_boundary_t2.py', 'validation/audit_s414_transition_boundary.py',
              'validation/report_s702_wake_extent.py', 'validation/report_s702_drag_mechanism.py',
              'validation/report_design_initialization.py']


def load(inputs):
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Unexpected study')
    launch = inputs.json(LAUNCH/'launch_receipt.json'); preparation = inputs.json(LAUNCH/'preparation.json')
    require(launch['study_sha256'] == preparation['study_sha256'] == STUDY_HASH, 'Launched study differs')
    spec = inputs.json(STUDY/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    receipt = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in receipt.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved output changed: '+name)
    retrieval = inputs.json(LAUNCH/'retrieval.json')
    require(retrieval['all_files_match'] and retrieval['files_verified'] == len(receipt), 'Retrieval coverage changed')
    runtime = inputs.json(STUDY/'runtime.json'); status = inputs.json(STUDY/'status.json')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == len(spec['jobs']) == 18,
            'Incomplete execution')
    require(runtime['python'] == spec['expected_python'] and runtime['versions'] == spec['expected_versions'],
            'Solver environment changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread limits changed')
    automatic = inputs.json(STUDY/'summary.json')
    # Pure recomputation of the restore-control gate from the retrieved exports.
    controls = t2.controls(STUDY, spec)
    require(controls == automatic['controls'] == inputs.json(STUDY/'controls.json') and controls['passed'],
            'Restore-control gate differs')
    for name in AUDIT_CODE: inputs.add(ROOT/name)
    return spec, runtime, controls, len(receipt)


def audit_exports(folder, job, spec, runtime, inputs):
    """One evaluation: record, polar, diagnostics, base state and saved arrays must agree."""
    result = inputs.json(folder/'result.json'); param = thermo(spec)
    require(result['complete'] and result['job'] == job['name'] and result['kind'] == job['kind']
            and result['amplification_jacobian'] == job['amplification_jacobian'] == 'corrected_af', 'Result does not belong to the declared job')
    cold = job['kind'] == 'cold_check'
    if cold:
        require(result['restore'] is None, 'Cold check carries a restored base')
    else:
        require(not result['error'].startswith('Base state not admitted'), 'Base state was not admitted: '+result['error'])
        restore = result['restore']
        require(max(restore['BL_residual_l2'], restore['coupling_residual_max']) <= spec['residual_tolerance']
                and restore['maximum_geometry_difference'] <= 2e-12 and restore['maximum_wake_difference'] <= 2e-11
                and restore['maximum_station_difference'] <= 2e-12, 'Restored base state was not the audited state')
    require(len(result['records']) == 1, 'Exactly one evaluation is expected')
    rec = result['records'][0]; meta = inputs.json(folder/'flow/metadata.json'); polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ from manifest')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(polar) == len(meta['diagnostics']) == 1, 'Evaluation counts differ')
    p, diag = polar[0], meta['diagnostics'][0]
    require(rec['role'] == job['kind'] and rec['requested_alpha_deg'] == job['alpha_deg'], 'Evaluation differs from the declared job')
    require(rec['converged'] == bool(int(p['converged'])) == (diag['status'] == 'converged'), 'Convergence differs')
    for field in ('CL', 'CD', 'CM', 'residual'):
        if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null evaluation has finite export')
        else: close(rec[field], float(p[field]), context='record/export '+field)
    close(float(p['alpha_deg']), job['alpha_deg'], context='exported incidence')
    for key in ('mach', 'reynolds', 'ncrit'): close(diag[key], spec[key], context=key)
    require(diag['ncrit_by_element'] == dict(fore=spec['ncrit'], aft=spec['ncrit']), 'Per-element Ncrit differs')
    require(diag['amplification_jacobian'] == 'corrected_af', 'Recorded derivative mode differs')
    require(rec['iterations_final_step'] == diag['iterations'] == len(diag['history']) <= spec['iterations'], 'Iteration count differs')
    require([(t['alpha_deg'], t['converged'], t['CL']) for t in rec['continuation']]
            == [(t['alpha_deg'], t['converged'], t['CL']) for t in diag['continuation']], 'Continuation trail differs')
    # The recorded start must be the declared one: a restored warm start is never called cold.
    if cold:
        require(rec['approach'] == 'cold' and rec['reached_by'] == 'cold solve at requested incidence'
                and len(rec['continuation']) == 1 and rec['continuation'][0]['alpha_deg'] == job['alpha_deg'],
                'Cold check was not a single cold solve at the requested incidence')
        require(set(diag['initialization'].values()) <= {'inviscid marching', 'gentle-profile seed; exact global equations restored', 'not reached'},
                'Cold check did not start from a cold initialization')
    else:
        base = job['base']
        require(rec['approach'] == job['direction'] and base['job'] in rec['reached_by'] and repr(base['alpha_deg']) in rec['reached_by'],
                'Recorded base state differs from the declared one')
        require((job['alpha_deg'] > base['alpha_deg']) == (job['direction'] == 'increasing'), 'Direction differs from the base side')
        require(set(diag['initialization'].values()) == {'restored audited saved state'}, 'Evaluation did not start from the restored state')
        steps = [t['alpha_deg'] for t in rec['continuation']]
        require(all((s > base['alpha_deg']) == (job['direction'] == 'increasing') for s in steps), 'Internal step left the declared side')
    if not rec['converged']:
        require(all(rec[k] is None for k in ('CL', 'CD', 'CM')) and not rec['state_files'] and diag['failure_reason'],
                'Failed evaluation has a force or a saved state')
        require(not any((folder/f'flow/{kind}_000_{e}.csv').exists() for kind in ('surface', 'wake') for e in deck.ELEMENTS),
                'Failed evaluation has surface exports')
        return result, None
    require(rec['admissible'] and not diag['failure_reason'], 'Converged record is inadmissible')
    residual = float(p['residual']); require(residual <= spec['residual_tolerance'], 'Tolerance failure')
    close(residual, max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
    for key in ('element_forces', 'boundary_layers'): require(rec[key] == diag[key], 'Recorded diagnostics differ')
    for field, term in (('CL', 'CL'), ('CD', 'CD_wake'), ('CM', 'CM')):
        close(rec[field], sum(f[term] for f in rec['element_forces'].values()), context='element sum '+field)
    for row in rec['transition_audit']:
        require(row['inside_cell'] and abs(row['amplification_residual']) <= 1e-7 and row['station_residual_max'] <= 1e-7,
                'Subcell transition audit failed')
    arrays, surfaces = {}, {}
    for element, panels in zip(deck.ELEMENTS, spec['panels']):
        surface = inputs.array(folder/f'flow/surface_000_{element}.csv'); wake = inputs.array(folder/f'flow/wake_000_{element}.csv')
        name = f'states/point_00_{element}.npz'; require(name in rec['state_files'], 'State file not recorded')
        with np.load(inputs.add(folder/name), allow_pickle=False) as saved:
            a = {k: saved[k].copy() for k in saved.files}
        n = len(surface)
        require(n == panels+1 and a['U'].shape == (4, n+len(wake)), 'State dimensions differ')
        require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite Newton state')
        for k, field in ((0, 'theta'), (1, 'delta_star'), (2, 'amplification_or_shear')):
            array_close(a['U'][k, :n][::-1], surface[field], 'Newton/surface '+field)
        for k, field in ((0, 'theta'), (1, 'delta_star')):
            array_close(a['U'][k, n:], wake[field], 'Newton/wake '+field)
        uk = mf.get_uk(a['U'][3], param)[0]
        array_close(uk[:n][::-1], surface['ue_over_vinf'], 'Newton/surface corrected velocity')
        array_close(uk[n:], wake['ue_over_vinf'], 'Newton/wake corrected velocity')
        hk = np.array([mf.get_Hk(a['U'][:, j], param)[0] for j in range(n+len(wake))])
        array_close(hk[:n][::-1], surface['Hk'], 'Newton/surface Hk'); array_close(hk[n:], wake['Hk'], 'Newton/wake Hk')
        require(np.array_equal(a['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence mismatch')
        labels = np.full(n, 'upper', dtype='<U5'); labels[a['lower']] = 'lower'
        require(np.array_equal(labels[::-1], surface['side']), 'Surface indices differ')
        require(set(a['lower']) | set(a['upper']) == set(range(n)) and not set(a['lower']) & set(a['upper']), 'Surface partition differs')
        for k, field in enumerate(('x_over_cref', 'y_over_cref')):
            array_close(a['xy'][k][::-1], surface[field], 'Surface geometry differs')
            array_close(a['wake_xy'][k], wake[field], 'Wake geometry differs')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
        for side in ('upper', 'lower'): audit_surface(side_data(surface, side), rec['boundary_layers'][element][side])
        ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0); arclength = np.r_[0., np.cumsum(ds)]
        array_close(a['wake_s']-a['wake_s'][0], arclength, 'Wake arclength differs')
        array_close(np.linalg.norm(a['wake_t'], axis=0), np.ones(len(wake)), 'Wake tangent normalization')
        info = rec['wake_discretization'][element]; require(info['nodes'] == len(wake), 'Wake node count differs')
        for k, value in (('arclength', arclength[-1]), ('first_spacing', ds[0]), ('last_spacing', ds[-1])):
            close(info[k], value, context=k)
        h = a['U'][1, n:]/a['U'][0, n:]; sy = 2*a['U'][0, n:]*(uk[n:]/param.Vinf)**((5+h)/2)
        close(sy[-1], rec['element_forces'][element]['CD_wake'], context='terminal Squire–Young drag')
        arrays[element], surfaces[element] = a, surface
    return result, dict(record=rec, residual=residual, arrays=arrays, surfaces=surfaces)


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs()
    spec, runtime, controls, retrieved = load(inputs)
    states, failed, restores = {}, [], []
    for job in spec['jobs']:
        result, state = audit_exports(STUDY/job['name'], job, spec, runtime, inputs)
        if state: states[job['name']] = state
        else: failed.append(dict(job=job['name'], kind=job['kind'], direction=job['direction'], alpha_deg=job['alpha_deg'],
            failure_reason=result['records'][0]['failure_reason'], continuation=result['records'][0]['continuation']))
        if result['restore']: restores.append(dict(job=job['name'], base_job=job['base']['job'], base_alpha_deg=job['base']['alpha_deg'], **result['restore']))
    environment = dict(python=platform.python_version(),
        versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})
    receipt = dict(input_sha256=inputs.files, audit_runtime=environment, states=[], complete=False, new_flow_solves=0,
        retrieved_files_verified=retrieved, controls=controls, failed_evaluations=failed, restored_base_states=restores,
        saved_states=len(states),
        method='Full export/array audit of every evaluation and fresh coupled residual, force, Cp/Cf, station, '
               'transition-root and wake-drag evaluation of every saved state; no nonlinear solve')
    path = OUT/'residual_recheck.json'
    if path.exists():
        previous = json.loads(path.read_text())
        require(previous['input_sha256'] == receipt['input_sha256'] and previous['audit_runtime'] == environment,
                'Audit changed; preserve previous receipt')
        receipt = previous
    for job in spec['jobs']:
        if job['name'] not in states or any(r['job'] == job['name'] for r in receipt['states']): continue
        row = recompute(job, spec, states[job['name']]); receipt['states'].append(row); write(path, receipt)
        print(f"{len(receipt['states'])}/{len(states)} {job['name']}: residual {row['max_residual']:.3g}", flush=True)
    receipt['complete'] = True; write(path, receipt)
    print('audit complete:', len(states), 'saved states;', len(failed), 'failed evaluation(s) without force')


if __name__ == '__main__':
    main()
