#!/usr/bin/env python3
"""Audit measured-incidence cases without imposing an experimental lift target."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import platform
import importlib.metadata
import numpy as np
from audit_s414_jacobian import recompute
from report_s702_wake_extent import ROOT, Inputs, sha, audit, require, close, array_close, selected_index, thermo, ELEMENTS, audit_surface, side_data, mf
from run_s414_incidence import controls
from s702_ncrit_sweep import write

STUDY = ROOT/'validation/results/s414_incidence_20261003'
OUT = ROOT/'validation/results/s414_incidence_analysis_20261003'
LAUNCH = ROOT/'validation/results/s414_incidence_launch_20261003'
STUDY_HASH = '996466c4fcdc4a49590cc5b57df111d3cc14b278b4bed547d30e1b92ca99b725'


def check_fixed_angle_identity(result, job, spec):
    require(job['comparison_basis'] == result['comparison_basis'] == 'fixed_alpha', 'Wrong comparison basis')
    require('target_cl' not in job and 'target_cl' not in result, 'Fixed-angle case silently given a lift target')
    require(result['target_alpha'] == job['target_alpha'] and
        result['reference_experimental_cl'] == job['reference_experimental_cl'], 'Measured condition changed')
    close(result['alpha_deg'], job['target_alpha'], tolerance=1e-12, context='Requested incidence')
    require(len(result['records']) == 1 and result['selected_flow_index'] == 0, 'Fixed-angle coverage differs')
    require(result['selected_path'] == 'cold_measured_incidence', 'Unqualified warm path')
    require(result['records'][0]['initialization'] == dict(index=0, path='cold_measured_incidence', phase='fixed_alpha'), 'Wrong fixed-angle phase')
    require(result['attempts'] == [dict(path='cold_measured_incidence', start_alpha=job['target_alpha'],
        requested_alpha=job['target_alpha'], feasible=True)], 'Fixed-angle attempt differs')


def load_states(inputs):
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Study identity changed')
    spec = inputs.json(STUDY/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    retrieved = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in retrieved.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved file changed: '+name)
    receipt = inputs.json(LAUNCH/'retrieval.json')
    require(receipt['all_files_match'] and receipt['files_verified'] == len(retrieved) == 1977, 'Retrieval incomplete')
    runtime = inputs.json(STUDY/'runtime.json'); status = inputs.json(STUDY/'status.json')
    automatic = inputs.json(STUDY/'summary.json')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == 4, 'Run incomplete')
    require(automatic['execution_checks_pass'] and automatic['required_jobs'] == automatic['converged_jobs'] == 4, 'Coverage incomplete')
    require(runtime['python'] == '3.12.13' and runtime['versions'] == spec['expected_versions'], 'Runtime changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread settings changed')
    control = controls(STUDY, spec)
    require(control['passed'] and control == automatic['controls'] == inputs.json(STUDY/'controls.json'), 'Exact control differs')
    states = {}
    for job in spec['jobs']:
        function = audit if job['comparison_basis'] == 'fixed_lift' else audit_fixed_angle
        state = function(STUDY/job['name'], job, spec, runtime, inputs)
        require(state['result']['feasible'], 'Incomplete scientific comparison')
        require(state['result']['amplification_jacobian'] == job['amplification_jacobian'] == 'corrected_af', 'Wrong derivative mode')
        require(all(d['amplification_jacobian'] == 'corrected_af' for d in state['metadata']['diagnostics']), 'Evaluation mode changed')
        states[job['name']] = state
    for name in ('validation/audit_s414_incidence.py', 'validation/audit_s414_jacobian.py', 'validation/report_s702_wake_extent.py'):
        inputs.add(ROOT/name)
    return spec, states, control


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs(); spec, states, control = load_states(inputs)
    receipt = dict(complete=False, input_sha256=inputs.files.copy(), states=[], controls=control,
        new_nonlinear_flow_solves=0, runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')}))
    path = OUT/'residual_recheck.json'
    if path.exists():
        import json
        old = json.loads(path.read_text())
        require(old['input_sha256'] == receipt['input_sha256'] and old['runtime'] == receipt['runtime'], 'Preserve changed audit')
        receipt = old
    for job in spec['jobs']:
        if any(r['job'] == job['name'] for r in receipt['states']): continue
        row = recompute(job, spec, states[job['name']])
        row.update(amplification_jacobian=job['amplification_jacobian'], comparison_basis=job['comparison_basis'])
        receipt['states'].append(row); write(path, receipt)
        print(job['name'], 'fresh residual', row['max_residual'], flush=True)
    receipt['complete'] = True; write(path, receipt)


# The full export audit below is adapted from the historical auditor. Only
# fixed-lift admissibility and initialization checks are replaced by the
# explicit fixed-angle checks above. Saved files and metadata are not changed.

def audit_fixed_angle(folder, job, spec, runtime, inputs):
    """Check all exports plus saved arrays, including nonzero-Mach velocity conversion."""
    r = inputs.json(folder/'result.json')
    meta = inputs.json(folder/'flow/metadata.json'); polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ from manifest')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(r['records']) == len(polar) == len(meta['diagnostics']), 'Evaluation counts differ')
    for rec, p, diag in zip(r['records'], polar, meta['diagnostics']):
        require(rec['converged'] == bool(int(p['converged'])) == (diag['status']=='converged'), 'Convergence differs')
        for field in ('alpha_deg', 'CL', 'CD', 'CM', 'residual'):
            if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null evaluation has finite export')
            else: close(rec[field], float(p[field]), context='record/export '+field)
        for key in ('mach', 'reynolds', 'ncrit'): close(diag[key], spec[key], context=key)
        require(diag['ncrit_by_element'] == dict(fore=spec['ncrit'], aft=spec['ncrit']), 'Per-element Ncrit differs')
    if not r['feasible']:
        require(all(r[k] is None for k in ('CL','CD','CM','alpha_deg')) and 'selected_flow_index' not in r,
                'Failed target has selected force/state')
        return dict(result=r, metadata=meta, arrays={}, surfaces={}, wakes={})
    i = selected_index(r, polar); diag = meta['diagnostics'][i]
    require(r['records'][i]['admissible'] and not diag['failure_reason'], 'Selected record is inadmissible')
    residual = float(polar[i]['residual'])
    require(residual <= spec['residual_tolerance'], 'Residual tolerance failure')
    close(residual, max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
    for key in ('element_forces', 'boundary_layers'): require(r[key] == diag[key], 'Selected diagnostics differ')
    for field, term in (('CL','CL'), ('CD','CD_wake'), ('CM','CM')):
        close(r[field], sum(f[term] for f in r['element_forces'].values()), context='element sum '+field)
    require(r['selected_path'] == r['records'][i]['initialization']['path'], 'Selected path differs')
    check_fixed_angle_identity(r, job, spec)
    for row in r['transition_audit']:
        require(row['inside_cell'] and abs(row['amplification_residual']) <= 1e-7 and row['station_residual_max'] <= 1e-7,
                'Subcell transition audit failed')
    param = thermo(spec); arrays, surfaces, wakes, stations = {}, {}, {}, {}
    for element, panels in zip(ELEMENTS, spec['panels']):
        inputs.add(folder/f'flow/geometry_{element}.dat')
        surface = inputs.array(folder/f'flow/surface_{i:03d}_{element}.csv')
        wake = inputs.array(folder/f'flow/wake_{i:03d}_{element}.csv')
        with np.load(inputs.add(folder/f'{element}_state.npz'), allow_pickle=False) as saved:
            a = {k:saved[k].copy() for k in saved.files}
        n = len(surface)
        require(n == panels+1 and a['U'].shape == (4, n+len(wake)), 'State dimensions differ')
        require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite Newton state')
        for k, field in ((0,'theta'), (1,'delta_star'), (2,'amplification_or_shear')):
            array_close(a['U'][k,:n][::-1], surface[field], 'Newton/surface '+field)
        for k, field in ((0,'theta'), (1,'delta_star')):
            array_close(a['U'][k,n:], wake[field], 'Newton/wake '+field)
        # U[3] is uncorrected speed; exports contain Karman–Tsien corrected speed.
        uk = mf.get_uk(a['U'][3], param)[0]
        array_close(uk[:n][::-1], surface['ue_over_vinf'], 'Newton/surface corrected velocity')
        array_close(uk[n:], wake['ue_over_vinf'], 'Newton/wake corrected velocity')
        hk = np.array([mf.get_Hk(a['U'][:,j], param)[0] for j in range(n+len(wake))])
        array_close(hk[:n][::-1], surface['Hk'], 'Newton/surface Hk')
        array_close(hk[n:], wake['Hk'], 'Newton/wake Hk')
        require(np.array_equal(a['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence mismatch')
        labels = np.full(n, 'upper', dtype='<U5'); labels[a['lower']] = 'lower'
        require(np.array_equal(labels[::-1], surface['side']), 'Surface indices differ')
        require(set(a['lower']) | set(a['upper']) == set(range(n)) and not set(a['lower']) & set(a['upper']), 'Surface index partition differs')
        for k, field in enumerate(('x_over_cref','y_over_cref')):
            array_close(a['xy'][k][::-1], surface[field], 'Surface geometry differs')
            array_close(a['wake_xy'][k], wake[field], 'Wake geometry differs')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
        for side in ('upper','lower'): audit_surface(side_data(surface,side), r['boundary_layers'][element][side])
        ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0)
        require(np.all(ds > 0), 'Nonpositive wake spacing')
        arclength = np.r_[0., np.cumsum(ds)]
        if 'wake_s' in a: array_close(a['wake_s']-a['wake_s'][0], arclength, 'Wake arclength differs')
        if 'wake_t' in a: array_close(np.linalg.norm(a['wake_t'], axis=0), np.ones(len(wake)), 'Wake tangent normalization')
        if 'wake_discretization' in r:
            info = r['wake_discretization'][element]
            require(info['nodes'] == len(wake), 'Wake node count differs')
            for k, value in (('arclength',arclength[-1]),('first_spacing',ds[0]),('last_spacing',ds[-1])):
                close(info[k],value,context=k)
        # Same Squire–Young estimator evaluated downstream; bookkeeping, not independent drag.
        h = a['U'][1,n:]/a['U'][0,n:]
        sy = 2* a['U'][0,n:] * (uk[n:]/param.Vinf)**((5+h)/2)
        close(sy[-1], r['element_forces'][element]['CD_wake'], context='terminal Squire–Young drag')
        stations[element] = dict(s=arclength, sy=sy)
        arrays[element], surfaces[element], wakes[element] = a, surface, wake
    return dict(result=r, metadata=meta, residual=residual, arrays=arrays, surfaces=surfaces, wakes=wakes, stations=stations)

if __name__ == '__main__':
    main()
