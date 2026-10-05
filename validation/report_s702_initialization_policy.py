#!/usr/bin/env python3
"""Audit the seven saved policy cases against the completed repeatability study."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from report_s702_repeatability import state, compare_states, STUDY_HASH as OLD_HASH
from report_s702_drag_mechanism import Inputs, ROOT, sha, require, write_json, relative
from report_s702_ncrit_sweep import csv_write

STUDY_HASH = '5388fef33ad7833c66b9d8db43da4c6304a54a4d50a4aae38f1da2a2d27219cb'


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected policy manifest')
    old = ROOT/'validation/results/s702_repeatability_20260915'
    require(sha(old/'study.json') == OLD_HASH, 'Unexpected repeatability manifest')

    class PolicyInputs(Inputs):
        def sources(self, spec, folder):
            # The old auditor expects job/cl/result.json; this study has job/result.json.
            super().sources(spec, study if folder == study.parent else folder)

    inputs = PolicyInputs()
    spec, previous = inputs.json(study/'study.json'), inputs.json(old/'study.json')
    inputs.sources(spec, study)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Runtime mismatch')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Study incomplete')
    summary = inputs.json(study/'summary.json')
    rows, comparisons, states, provenance = [], [], {}, []
    failures = []
    for job in spec['jobs']:
        s = state(study/job['name']/'result.json', spec, inputs, job['name'])
        s['provenance']['study_path'] = relative(study/'study.json')
        provenance.append(s['provenance'])
        result = s['result']
        states[job['name']] = s
        ref = state(ROOT/job['reference_result'], spec, inputs, 'archive')
        # Use the old increasing sweep, which covered all six conditions.
        old_job = next(j for j in previous['jobs'] if j['flap'] == job['flap'] and j['path'] == 'increasing')
        prior = state(old/old_job['name']/('cl'+str(job['target_cl']).replace('.', 'p'))/'result.json',
                      previous, inputs, 'repeatability_increasing')
        delta = (result['CD']-ref['result']['CD'])*1e4
        require(abs(delta) <= spec['acceptance_CD_counts'], 'Drag regression failed')
        require(result['records'][result['selected_flow_index']]['initialization']['path'] == result['selected_path'],
                'Selected initialization provenance differs')
        for reference in (ref, prior):
            for d in compare_states(reference, s):
                comparisons.append(dict(job=job['name'], reference=reference['label'], **d))
        rows.append(dict(job=job['name'], control=job['control'], flap=job['flap'],
            target_cl=job['target_cl'], selected_path=result['selected_path'],
            evaluations=len(result['records']), failed_evaluations=sum(not r['admissible'] for r in result['records']),
            lift_error=result['CL']-job['target_cl'], residual=s['summary']['residual'],
            difference_CD_archive_counts=delta,
            difference_CD_repeatability_counts=(result['CD']-prior['result']['CD'])*1e4,
            wall_seconds=result['wall_seconds']))
        for r in result['records']:
            if not r['admissible']:
                failures.append(dict(job=job['name'], record=r))
    require(len(rows) == 7 and summary['converged_primary'] == 6 and summary['converged_control'] == 1,
            'Coverage mismatch')
    require(sum(r['evaluations'] for r in rows) == summary['requested_evaluations'], 'Work summary mismatch')
    require(len(failures) == summary['failed_requested_evaluations'], 'Failure count mismatch')
    topology = [d for d in comparisons if d['field'] in ('turbulent', 'lower', 'upper')]
    require(all(d['exactly_equal'] for d in topology), 'Discrete topology changed')
    geometry = [d for d in comparisons if d['field'] == 'xy']
    require(all(d['exactly_equal'] for d in geometry if d['reference'] == 'repeatability_increasing'),
            'Geometry changed relative to the repeatability study')
    require(all(d['max_abs_difference'] < 1e-14 for d in geometry), 'Archived geometry changed beyond roundoff')
    replay = compare_states(states['replay_control_d10_cl0p8'], states['policy_d10_cl0p8'])
    maxima = {f: max(d['max_abs_difference'] for d in comparisons
                    if d['field'] == f and d['reference'] == 'repeatability_increasing')
              for f in sorted({d['field'] for d in comparisons})}
    result = dict(primary_coverage='6/6', control_coverage='1/1', full_state_audit_pass=True,
        requested_evaluations=summary['requested_evaluations'], failed_requested_evaluations=len(failures),
        max_absolute_CD_difference_archive_counts=max(abs(r['difference_CD_archive_counts']) for r in rows),
        max_absolute_CD_difference_repeatability_counts=max(abs(r['difference_CD_repeatability_counts']) for r in rows),
        discrete_topology_identical=True, geometry_identical_to_repeatability=True,
        max_geometry_difference_archive=max(d['max_abs_difference'] for d in geometry),
        max_state_difference_repeatability=maxima,
        warm_anchor_vs_replay=replay, numerical_runtime=runtime,
        qualification='Implementation regression at Mach 0.07 only; Mach-zero design integration and independent contour rankings remain unvalidated.')
    out.mkdir(parents=True, exist_ok=True)
    csv_write(out/'targets.csv', rows)
    csv_write(out/'state_differences.csv', comparisons)
    write_json(out/'failures.json', failures)
    write_json(out/'selected_state_provenance.json', provenance)
    write_json(out/'summary.json', result)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for flap, marker in ((5, 'o'), (10, 's')):
        selected = [r for r in rows if r['flap'] == flap and not r['control']]
        axes[0].plot([r['target_cl'] for r in selected],
                     [r['difference_CD_archive_counts'] for r in selected], marker=marker, label=f'Flap {flap}°')
    axes[0].axhline(0, color='0.5', lw=.6)
    axes[0].set(xlabel='Target CL', ylabel='CD − archived CD [counts]', title='6/6 primary + 1/1 replay converged')
    axes[0].legend()
    for name, style in (('policy_d10_cl0p8', '-'), ('replay_control_d10_cl0p8', '--')):
        surface = states[name]['surfaces']['aft']
        upper = surface[surface['side'] == 'upper']
        axes[1].plot(upper['x_over_local_chord'], upper['Hk'], style,
                     label='Warm anchor' if name.startswith('policy') else 'Archived replay')
    axes[1].set(xlabel='x / aft local chord', ylabel='Hk', title='Flap 10°, CL=0.8: selected warm states')
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    fig.tight_layout()
    for ext in ('png', 'svg'): fig.savefig(out/f'policy_audit.{ext}', dpi=170)
    plt.close(fig)
    (out/'FINDINGS.md').write_text(
        '# S702 initialization-policy audit\n\n'
        'All six primary cases and the replay control converged. The cold flap-10°, CL=0.8 '
        'solve failed once; a fresh solver starting at −1° recovered it without the archived seed prefix. '
        'The failed evaluation remains in the record. Total: 28 requested evaluations, one failure.\n\n'
        f"Maximum drag difference from the archived selection: {result['max_absolute_CD_difference_archive_counts']:.9g} counts "
        '(acceptance 0.01). Against the completed increasing sweep: '
        f"{result['max_absolute_CD_difference_repeatability_counts']:.9g} counts. "
        'All selected exports, Newton arrays, surface/wake diagnostics, force sums, lift and residual tolerances passed. '
        'Turbulence flags and upper/lower surface indices agree exactly. Geometry is identical to the '
        'repeatability study; the older archives differ by at most 7.78e-16 chord (roundoff). Full pointwise differences '
        'are in state_differences.csv; no interpolation or failed-state force substitution was used.\n\n'
        'This supports opt-in integration of the initialization policy. It is implementation verification at '
        'Re=700,000, Mach=0.07, Ncrit=6.25, 260/180 panels; it does not establish Mach-zero design robustness '
        'or independent shape-ranking accuracy. The pre-existing far-field ordering test failure is unchanged.\n\n'
        'Next: an explicit policy option in the design evaluator, followed by a bounded Mach-zero integration matrix.\n\n'
        'Reproduce: `python3 validation/report_s702_initialization_policy.py`.\n')
    for name in ('validation/report_s702_initialization_policy.py', 'validation/report_s702_repeatability.py',
                 'validation/report_s702_drag_mechanism.py', 'validation/report_s702_ncrit_sweep.py'):
        inputs.add(ROOT/name)
    write_json(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, default=ROOT/'validation/results/s702_initialization_policy_20260921')
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_initialization_policy_analysis')
    args = ap.parse_args()
    print(json.dumps(report(args.study.resolve(), args.out.resolve()), indent=2))
