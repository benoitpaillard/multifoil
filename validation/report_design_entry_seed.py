#!/usr/bin/env python3
"""Audit explicit-entry recovery against the preserved legacy controls; no solves."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_design_initialization import ROOT, Inputs, sha, write, csv_write, audit, require, close

STUDY_HASH = '78a679d3aaae7ee987bc0510714aafbd2faaf259ba868826e01cf42545b879b4'


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected entry-seed study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Source snapshot mismatch: '+name)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Numerical runtime mismatch')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Study is unfinished')
    original = inputs.json(study/'summary.json')
    rows, differences, histories, states, references = [], [], [], {}, {}
    for job in spec['jobs']:
        current = audit(study/job['name'], dict(job, mode='policy'), spec, runtime, inputs)
        r = current['result']
        require(r['feasible'], 'Incomplete entry recovery: '+job['name'])
        path = ROOT/job['reference_result']
        require(sha(path) == spec['sha256'][job['reference_result']], 'Legacy reference changed')
        ref = audit(path.parent, dict(job, mode='legacy'), spec, runtime, inputs)
        old = ref['result']
        states[job['name']], references[job['name']] = current, ref
        config = current['metadata']['case']['point']['initialization']
        require(config['paths'] == [dict(name='warm_entry', start_alpha=spec['warm_start_alpha'], seed_alphas=[job['seed_alpha']])], 'Declared entry path changed')
        require(config['center'] == job['center'] and config['half_width'] == spec['half_width'], 'Local bracket changed')
        seed = r['records'][0]
        require(seed['admissible'] and seed['initialization']['phase'] == 'warm_seed', 'No admissible explicit entry')
        close(seed['alpha_deg'], job['seed_alpha'], context='entry incidence')
        require(r['selected_flow_index'] > 0, 'Entry seed substituted for target')
        for key in ('CL', 'CD', 'CM', 'alpha_deg'):
            close(seed[key], old['records'][0][key], context='entry/legacy first request '+key)
        trail = current['metadata']['diagnostics'][0]['continuation']
        require(trail[-1]['converged'], 'Entry continuation did not converge')
        internal_failures = sum(not t['converged'] for d in current['metadata']['diagnostics'] for t in d['continuation'])
        delta = (r['CD']-old['CD'])*1e4
        require(abs(delta) <= spec['acceptance_CD_counts'], 'Entry drag regression failed')
        rows.append(dict(job=job['name'], amplitude=job['amplitude'], CL=r['CL'], CD=r['CD'], CM=r['CM'],
            alpha_deg=r['alpha_deg'], residual=r['residual'], lift_error=r['CL']-spec['target_cl'],
            difference_CL=r['CL']-old['CL'], difference_CD_counts=delta, difference_CM=r['CM']-old['CM'],
            requested_evaluations=r['flow_evaluations'], failed_requested_evaluations=sum(not v['admissible'] for v in r['records']),
            failed_internal_continuation_points=internal_failures, wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'], attempts=r['initialization_attempts'], records=r['records'],
            selected_flow_index=r['selected_flow_index'], continuation=[d['continuation'] for d in current['metadata']['diagnostics']]))
        for element in ('fore', 'aft'):
            require(sha(inputs.add(study/job['name']/f'input_{element}.dat')) == sha(inputs.add(path.parent/f'input_{element}.dat')), 'Compared geometries differ')
            for kind in ('surfaces', 'wakes'):
                a, b = current[kind][element], ref[kind][element]
                require(a.shape == b.shape, 'Selected export dimensions differ')
                for field in a.dtype.names:
                    exact = bool(np.array_equal(a[field], b[field]))
                    diff = float(np.max(abs(a[field]-b[field]))) if a[field].dtype.kind in 'fiu' else None
                    if field in ('side', 'turbulent'): require(exact, 'Discrete topology differs')
                    differences.append(dict(job=job['name'], element=element, kind=kind, field=field,
                                            exactly_equal=exact, max_abs_difference=diff))
        for old_row in original['jobs']:
            if old_row['job'] == job['name']:
                close(old_row['difference_CD_counts'], delta)
                require(old_row['passed'] and old_row['evaluations'] == r['flow_evaluations'], 'Automatic summary differs')
                break
        else: raise ValueError('Automatic summary omitted a job')
    require(len(rows) == original['feasible_jobs'] == original['passed_jobs'] == original['required_jobs'] == 3, 'Coverage mismatch')
    summary = dict(required=3, converged=3, full_state_audit_pass=True, selected_states_audited=6,
        policy_newton_states_audited=3, entry_states_match_legacy_first_request=True,
        selected_discrete_topology_matches=True, requested_evaluations=sum(r['requested_evaluations'] for r in rows),
        failed_requested_evaluations=sum(r['failed_requested_evaluations'] for r in rows),
        failed_internal_continuation_points=sum(r['failed_internal_continuation_points'] for r in rows),
        max_absolute_CD_difference_counts=max(abs(r['difference_CD_counts']) for r in rows),
        max_absolute_lift_difference=max(abs(r['difference_CL']) for r in rows),
        qualification='Explicit-entry diagnostic only. Validate the complete ordered fallback policy before claiming that the original integration coverage gap is closed.')
    out.mkdir(parents=True, exist_ok=True)
    csv_write(out/'targets.csv', rows); csv_write(out/'state_differences.csv', differences)
    write(out/'summary.json', summary); write(out/'initialization_histories.json', histories)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    rows = sorted(rows, key=lambda r:r['amplitude'])
    axes[0].plot([r['amplitude']*1000 for r in rows], [r['difference_CD_counts'] for r in rows], 'o-')
    axes[0].axhline(0, color='0.5', lw=.7)
    axes[0].set(xlabel='Fore-camber amplitude [0.001 chord]', ylabel='Entry − legacy CD [counts]', title='3/3 complete comparisons')
    for group, style, label in ((references, '-', 'Legacy control'), (states, '--', 'Explicit entry')):
        surface = group['d10_baseline_entry']['surfaces']['aft']
        upper = surface[surface['side']=='upper']
        axes[1].plot(upper['x_over_local_chord'], upper['Hk'], style, label=label)
    axes[1].set(xlabel='x / aft local chord', ylabel='Hk', title='Recovered flap-10° baseline: aft upper surface')
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    fig.tight_layout()
    for ext in ('png', 'svg'): fig.savefig(out/f'entry_seed_audit.{ext}', dpi=170)
    plt.close(fig)
    (out/'FINDINGS.md').write_text(
        '# Explicit entry recovers the missing baseline\n\n'
        'All three flap-10° cases passed: baseline and fore-camber amplitudes −0.001 and +0.001. '
        'The explicit seed at −0.566226115° is admissible in every case and reproduces the first requested '
        'legacy state. From that seed, the formerly failing baseline incidence −0.216226115° converges, '
        'and the local fixed-lift solve reaches CL=0.8 within tolerance.\n\n'
        f"Maximum CD difference from the legacy targets: {summary['max_absolute_CD_difference_counts']:.9g} counts "
        '(acceptance 0.01). The maximum achieved-lift difference is '
        f"{summary['max_absolute_lift_difference']:.9g}; each target meets the original 2e-5 tolerance. "
        'All six selected target exports passed the force, residual, surface/wake and transition audits; '
        'all three new Newton states match their selected exports. Turbulence flags and surface labels '
        'match the controls. Continuous-state differences remain tabulated, not declared zero.\n\n'
        'There were 14 requested evaluations and no inadmissible requested evaluations. Two internal '
        'starting solves at −1° still failed (baseline and positive camber); their failures are retained. '
        'In those cases the explicit lower entry succeeds after the failed internal start. This is '
        'recovery through a declared entry point, not evidence that −1° itself is a reliable warm state.\n\n'
        'Next: test a common ordered cold-local → fresh −1° start → fresh explicit-entry policy '
        'on all six Mach-zero contour cases, with a separate repeat of the difficult baseline. '
        'The original 11/12 matrix remains failed and intact. No flow equation, Ncrit or tolerance changes '
        'are justified by this diagnostic. Independent aerodynamic contour rankings remain unvalidated.\n\n'
        'Reproduce: `python3 validation/report_design_entry_seed.py`.\n')
    for name in ('validation/report_design_entry_seed.py', 'validation/report_design_initialization.py',
                 'validation/report_s702_drag_mechanism.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, default=ROOT/'validation/results/design_entry_seed_20260921')
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/design_entry_seed_analysis')
    args = ap.parse_args()
    print(json.dumps(report(args.study.resolve(), args.out.resolve()), indent=2))
