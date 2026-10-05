#!/usr/bin/env python3
"""Audit completed Mach-zero integration outputs, including the failed target."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_drag_mechanism import audit_surface, side_data, selected_index, close, require

ROOT = Path(__file__).resolve().parents[1]
STUDY_HASH = 'f95ec243e739b69d883062ae9ecb37133b0e71910d62e352d6cdafd103a201f6'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


class Inputs:
    def __init__(self): self.files = {}

    def add(self, path):
        path = Path(path).resolve()
        self.files[str(path.relative_to(ROOT))] = sha(path)
        return path

    def json(self, path): return json.loads(self.add(path).read_text())

    def csv(self, path):
        with self.add(path).open(newline='') as stream: return list(csv.DictReader(stream))

    def array(self, path):
        values = np.atleast_1d(np.genfromtxt(self.add(path), delimiter=',', names=True, dtype=None, encoding='utf-8'))
        for field in values.dtype.names:
            if values[field].dtype.kind in 'fiu':
                require(np.isfinite(values[field]).all(), 'Nonfinite selected '+str(path))
        return values


def csv_write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def audit(folder, job, spec, runtime, inputs):
    result = inputs.json(folder/'result.json')
    config = inputs.json(folder/'configuration.json')
    meta = inputs.json(folder/'flow/metadata.json')
    polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['config'] == config, 'Export/config mismatch')
    require(meta['case']['point'] == config['points'][0], 'Point mismatch')
    require(meta['case']['panels'] == spec['config']['panels'], 'Panel settings mismatch')
    require(config['constraints'] == spec['config']['constraints'], 'Design restrictions changed')
    require(config['mach'] == 0 and config['reynolds'] == 700000., 'Physical conditions mismatch')
    require(config['scenarios'] == [dict(name='ncrit6p25', ncrit=6.25)], 'Transition scenario changed')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime mismatch')
    require(len(polar) == len(result['records']) == len(meta['diagnostics']) == result['flow_evaluations'], 'Work count mismatch')
    require(meta['case']['initialization_history'] == [r['initialization'] for r in result['records']], 'Initialization history mismatch')
    require(sum(r['admissible'] for r in result['records']) == result['admissible_evaluations'], 'Admissibility count mismatch')
    for record, row, diag in zip(result['records'], polar, meta['diagnostics']):
        require(record['admissible'] == (diag['status'] == 'converged') == bool(int(row['converged'])), 'Convergence mismatch')
        for field in ('alpha_deg', 'CL', 'CD', 'CM', 'residual'):
            if record[field] is None:
                require(not np.isfinite(float(row[field])), 'Null record has finite export')
            else: close(record[field], float(row[field]), context='record/export '+field)
        close(diag['mach'], 0.); close(diag['reynolds'], 700000.); close(diag['ncrit'], 6.25)
    if not result['feasible']:
        require(all(result[k] is None for k in ('CL', 'CD', 'CM', 'alpha_deg')), 'Failed flow has final forces')
        require('selected_flow_index' not in result, 'Failed flow has selected state')
        return dict(result=result, metadata=meta, surfaces={}, wakes={})
    index = selected_index(result, polar)
    diag = meta['diagnostics'][index]
    require(result['records'][index]['admissible'] and not diag['failure_reason'], 'Selected flow is invalid')
    require(abs(result['CL']-spec['target_cl']) <= 2e-5 and result['residual'] <= 1e-7, 'Selected tolerances failed')
    close(result['residual'], max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
    for key in ('boundary_layers', 'element_forces'):
        require(result[key] == diag[key], 'Selected diagnostics mismatch: '+key)
    for key, term in (('CD', 'CD_wake'), ('CL', 'CL'), ('CM', 'CM')):
        close(result[key], sum(v[term] for v in diag['element_forces'].values()), context='element sum '+key)
    require(result['selected_path'] == result['records'][index]['initialization']['path'], 'Selected path mismatch')
    surfaces, wakes = {}, {}
    for element, panels in zip(('fore', 'aft'), spec['config']['panels']):
        inputs.add(folder/f'input_{element}.dat'); inputs.add(folder/f'flow/geometry_{element}.dat')
        surface = inputs.array(folder/f'flow/surface_{index:03d}_{element}.csv')
        wake = inputs.array(folder/f'flow/wake_{index:03d}_{element}.csv')
        require(len(surface) == panels+1, 'Selected panel count mismatch')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake state')
        for side in ('upper', 'lower'):
            audit_surface(side_data(surface, side), diag['boundary_layers'][element][side])
        if job['mode'] == 'policy':
            with np.load(inputs.add(folder/f'{element}_state.npz'), allow_pickle=False) as saved:
                require(all(np.isfinite(saved[k]).all() for k in saved.files), 'Nonfinite Newton state')
                n = len(surface)
                for k, field in ((0, 'theta'), (1, 'delta_star'), (2, 'amplification_or_shear'), (3, 'ue_over_vinf')):
                    require(np.max(abs(saved['U'][k, :n][::-1]-surface[field])) < 1e-12, 'Newton/surface mismatch')
                for k, field in ((0, 'theta'), (1, 'delta_star'), (3, 'ue_over_vinf')):
                    require(np.max(abs(saved['U'][k, n:]-wake[field])) < 1e-12, 'Newton/wake mismatch')
                require(np.array_equal(saved['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence mismatch')
                labels = np.full(n, 'upper', dtype='<U5'); labels[saved['lower']] = 'lower'
                require(np.array_equal(labels[::-1], surface['side']), 'Surface index mismatch')
                for k, field in enumerate(('x_over_cref', 'y_over_cref')):
                    require(np.max(abs(saved['xy'][k][::-1]-surface[field])) < 1e-12, 'Surface geometry mismatch')
                    require(np.max(abs(saved['wake_xy'][k]-wake[field])) < 1e-12, 'Wake geometry mismatch')
        surfaces[element], wakes[element] = surface, wake
    separation = max(v['separated_arclength_fraction'] for layers in diag['boundary_layers'].values() for v in layers.values())
    close(result['max_separated_fraction'], separation)
    require(separation <= config['constraints']['max_separated_fraction'], 'Design separation restriction failed')
    return dict(result=result, metadata=meta, surfaces=surfaces, wakes=wakes)


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected integration study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Source snapshot mismatch: '+name)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Runtime mismatch')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Study is not complete')
    original = inputs.json(study/'summary.json')
    states = {j['name']: audit(study/j['name'], j, spec, runtime, inputs) for j in spec['jobs']}
    rows, comparisons, transitions, failures = [], [], [], []
    for job in spec['jobs']:
        state = states[job['name']]; r = state['result']
        rows.append(dict(job=job['name'], flap=job['flap'], amplitude=job['amplitude'], mode=job['mode'],
            feasible=r['feasible'], CL=r['CL'], CD=r['CD'], CM=r['CM'], alpha_deg=r['alpha_deg'],
            lift_error=r['CL']-spec['target_cl'] if r['feasible'] else None, residual=r.get('residual'),
            path=r.get('selected_path'), requested_evaluations=r['flow_evaluations'],
            failed_requested_evaluations=sum(not a['admissible'] for a in r['records']), wall_seconds=r['wall_seconds']))
        for index, record in enumerate(r['records']):
            if not record['admissible']:
                failures.append(dict(job=job['name'], record=record, diagnostics=state['metadata']['diagnostics'][index]))
        if r['feasible']:
            for element, layers in r['boundary_layers'].items():
                for side, layer in layers.items():
                    transitions.append(dict(job=job['name'], element=element, side=side,
                        transition_x=layer['transition']['x_over_local_chord'] if layer['transition'] else None,
                        separated_fraction=layer['separated_arclength_fraction'], max_Hk=layer['max_Hk']))
    pairs = []
    for name in sorted({j['pair'] for j in spec['jobs']}):
        a, b = states[name+'_legacy'], states[name+'_policy']
        complete = a['result']['feasible'] and b['result']['feasible']
        delta = (b['result']['CD']-a['result']['CD'])*1e4 if complete else None
        pairs.append(dict(pair=name, complete=complete, difference_CD_counts=delta,
                          passed=bool(complete and abs(delta) <= spec['acceptance_CD_counts'])))
        if not complete: continue
        for element in ('fore', 'aft'):
            require(sha(study/(name+'_legacy')/f'input_{element}.dat') == sha(study/(name+'_policy')/f'input_{element}.dat'), 'Pair has different input geometry')
            for kind in ('surfaces', 'wakes'):
                x, y = a[kind][element], b[kind][element]
                require(x.shape == y.shape, 'Pair export dimensions differ')
                for field in x.dtype.names:
                    exact = bool(np.array_equal(x[field], y[field]))
                    difference = float(np.max(abs(x[field]-y[field]))) if x[field].dtype.kind in 'fiu' else None
                    if field in ('side', 'turbulent'): require(exact, 'Pair discrete topology changed')
                    comparisons.append(dict(pair=name, element=element, kind=kind, field=field,
                                            exactly_equal=exact, max_abs_difference=difference))
    require(sum(r['feasible'] for r in rows) == original['converged_jobs'], 'Coverage mismatch')
    require(sum(r['requested_evaluations'] for r in rows) == original['requested_evaluations'], 'Work summary mismatch')
    for pair, old in zip(pairs, original['pairs']):
        require(all(pair[k] == old[k] for k in pair), 'Paired summary mismatch')
    report = dict(completed_jobs=len(rows), feasible_jobs=sum(r['feasible'] for r in rows),
        selected_states_audited=sum(r['feasible'] for r in rows), saved_policy_states_audited=sum(r['feasible'] and r['mode']=='policy' for r in rows),
        failed_jobs=[r['job'] for r in rows if not r['feasible']], passed_pairs=sum(p['passed'] for p in pairs), required_pairs=len(pairs),
        requested_evaluations=sum(r['requested_evaluations'] for r in rows), failed_requested_evaluations=len(failures),
        max_absolute_complete_pair_CD_difference_counts=max(abs(p['difference_CD_counts']) for p in pairs if p['complete']),
        successful_state_consistency_pass=True, complete_pair_discrete_topology_identical=True,
        regression_pass=all(p['passed'] for p in pairs),
        qualification='Incomplete baseline prevents acceptance or policy-based ranking at flap 10; numerical agreement is not independent aerodynamic validation.')
    out.mkdir(parents=True, exist_ok=True)
    for name, values in (('targets', rows), ('pairs', pairs), ('state_differences', comparisons), ('transitions', transitions)):
        csv_write(out/(name+'.csv'), values)
    write(out/'summary.json', report); write(out/'failed_evaluations.json', failures)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    for ax, flap in zip(axes, (5, 10)):
        for mode, style in (('legacy', '-o'), ('policy', '--s')):
            values = sorted([r for r in rows if r['flap']==flap and r['mode']==mode], key=lambda r:r['amplitude'])
            ax.plot([r['amplitude']*1000 for r in values],
                    [r['CD']*1e4 if r['feasible'] else np.nan for r in values], style, label=mode)
        ax.set(title=f'Flap {flap}°', xlabel='Fore-camber amplitude [0.001 chord]', ylabel='CD [counts]')
        ax.grid(alpha=.2); ax.legend()
    axes[1].text(.5, .52, 'Policy baseline failed\n(no force assigned)', transform=axes[1].transAxes, ha='center', color='#a32d30')
    fig.suptitle('Mach 0, CL=0.8: 11/12 feasible; 5/6 complete comparisons')
    fig.tight_layout()
    for ext in ('png', 'svg'): fig.savefig(out/f'integration_audit.{ext}', dpi=170)
    plt.close(fig)
    (out/'FINDINGS.md').write_text(
        '# Mach-zero design integration: incomplete\n\n'
        'All 12 jobs finished. Eleven are feasible; five of six paired comparisons pass. '
        'The flap-10° baseline fails with both declared policy paths, while the legacy control converges. '
        'The original acceptance criterion therefore fails and the policy cannot supply a complete flap-10° design ranking.\n\n'
        f"The largest absolute drag difference among complete pairs is {report['max_absolute_complete_pair_CD_difference_counts']:.9g} counts "
        '(limit 0.01). All 11 selected exports passed force, residual, lift, surface/wake and transition checks; '
        'the five selected policy Newton states match their exports. Successful pairs have identical turbulence flags '
        'and surface labels. Small pointwise differences are tabulated without claiming identical continuous states.\n\n'
        'There are 55 requested evaluations and three inadmissible requests. Both baseline policy attempts fail '
        'at −0.216226115° with Newton residual 0.0147453. Their internal starts (−0.166226115° and −1°) also fail. '
        'Thus the nominal warm-anchor path never obtains a converged starting state at Mach zero. It then repeats '
        'the same unsuccessful requested incidence. The legacy path also fails internally at −1°, but its first '
        'requested incidence, −0.566226115°, converges (CL=0.7592255), allowing continuation to the target. '
        'This explains why declaring a warm starting angle alone is insufficient; it does not yet prove which '
        'additional path will reproduce the selected target.\n\n'
        'Next diagnostic: one explicit seed at the legacy lower endpoint (center−0.4°), then the unchanged local '
        'fixed-lift policy. Check the baseline and both neighboring flap-10° contours with the same declared path. '
        'Use a fresh study and retain this failed matrix. No closure, tolerance, Ncrit or geometry fitting is proposed.\n\n'
        'Reproduce: `python3 validation/report_design_initialization.py`. '
        'Sources, configurations, all requested failures and internal continuation histories remain traceable in '
        'the study and input_manifest.json. This is implementation evidence, not independent aerodynamic validation.\n')
    for name in ('validation/report_design_initialization.py', 'validation/report_s702_drag_mechanism.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return report


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, default=ROOT/'validation/results/design_initialization_20260921')
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/design_initialization_analysis')
    args = ap.parse_args()
    print(json.dumps(report(args.study.resolve(), args.out.resolve()), indent=2))
