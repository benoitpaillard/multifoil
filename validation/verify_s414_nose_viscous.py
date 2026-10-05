#!/usr/bin/env python3
"""Independent check of the V4 convex-nose viscous analysis from the raw polar exports; no flow solve."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, sha, write
from report_s702_drag_mechanism import require

STUDY = ROOT/'validation/results/s414_nose_viscous_20261004'
OUT = ROOT/'validation/results/s414_nose_viscous_analysis_20261004'
ARCHIVE = ROOT/'validation/results/s414_jacobian_20261003'
LOWER = ROOT/'validation/results/s414_lower_bucket_20261003'
REFERENCES = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/targets.csv'
TESTS = ['test_s414_nose_viscous', 'test_s414_nose_geometry']


def polar(folder):
    with (folder/'flow/polar.csv').open(newline='') as stream:
        return [dict(alpha=float(r['alpha_deg']), converged=r['converged'] == '1', CL=float(r['CL']), CD=float(r['CD'])) for r in csv.DictReader(stream)]


def main():
    require(not (OUT/'verification.json').exists(), 'Preserve completed verification')
    audit = json.loads((OUT/'residual_recheck.json').read_text()); summary = json.loads((OUT/'summary.json').read_text())
    spec = json.loads((STUDY/'study.json').read_text()); screens = spec['equivalence_screens']; jobs = {j['name']: j for j in spec['jobs']}
    require(audit['complete'] and len(audit['states']) == audit['saved_states'] == 19, 'Audit coverage differs')
    require(max(s['max_residual'] for s in audit['states']) <= spec['residual_tolerance'], 'Fresh residual exceeds tolerance')
    changed = [n for n, d in audit['input_sha256'].items() if sha(ROOT/n) != d]; require(not changed, 'Audited inputs changed: '+str(changed[:3]))
    # Selected fixed-lift states from the raw polar rows.
    selected = {}
    for name in ('nominal_cl0056_control', 'convex_nose_cl0056', 'convex_nose_cl0175', 'convex_nose_cl0284'):
        rows = polar(STUDY/name); hit = [r for r in rows if r['converged'] and abs(r['CL']-jobs[name]['target_cl']) <= spec['cl_tolerance']]
        require(len(hit) == 1, 'Selected state is not unique in the raw polar: '+name); selected[name] = hit[0]
    archived = [r for r in polar(ARCHIVE/'cl0056_corrected_af') if abs(r['CL']-.056) <= spec['cl_tolerance']][0]
    require(selected['nominal_cl0056_control'] == archived, 'Nominal control differs from the archived polar row')
    nominal_0175 = [r for r in polar(LOWER/'cl0175_nominal') if r['converged'] and abs(r['CL']-.175) <= spec['cl_tolerance']][0]
    with REFERENCES.open(newline='') as stream: reference = {r['case']: r for r in csv.DictReader(stream)}
    experiment = {c: float(reference[c]['experiment_CD'])*1e4 for c in ('cl0056', 'cl0175', 'cl028')}
    errors = dict(nominal_cl0056=selected['nominal_cl0056_control']['CD']*1e4-experiment['cl0056'],
        convex_cl0056=selected['convex_nose_cl0056']['CD']*1e4-experiment['cl0056'],
        nominal_cl0175=nominal_0175['CD']*1e4-experiment['cl0175'], convex_cl0175=selected['convex_nose_cl0175']['CD']*1e4-experiment['cl0175'],
        convex_cl0284=selected['convex_nose_cl0284']['CD']*1e4-experiment['cl028'])
    reported = {(t['geometry'], t['target_cl']): t for t in summary['targets']}
    for key, geometry, cl in (('nominal_cl0056', 'nominal', .056), ('convex_cl0056', 'convex_nose', .056), ('nominal_cl0175', 'nominal', .175),
                              ('convex_cl0175', 'convex_nose', .175), ('convex_cl0284', 'convex_nose', .284)):
        require(abs(errors[key]-reported[(geometry, cl)]['error_CD_counts']) < 1e-9, 'Drag error differs from the raw polar: '+key)
    paired = dict(nominal=(nominal_0175['CD']-selected['nominal_cl0056_control']['CD'])*1e4,
                  convex_nose=(selected['convex_nose_cl0175']['CD']-selected['convex_nose_cl0056']['CD'])*1e4,
                  experiment=experiment['cl0175']-experiment['cl0056'])
    for row in summary['paired_drag_change']:
        if row['series'] in paired: require(abs(row['CD_change_counts']-paired[row['series']]) < 1e-9, 'Paired change differs: '+row['series'])
    # Paths from the raw polar rows.
    up, down = polar(STUDY/'convex_nose_increasing'), polar(STUDY/'convex_nose_decreasing')[::-1]
    require([r['alpha'] for r in up] == [r['alpha'] for r in down] == spec['path_angles_deg'], 'Path incidences differ')
    pairs = [dict(alpha_deg=a['alpha'], delta_CL=a['CL']-b['CL'], delta_CD_counts=(a['CD']-b['CD'])*1e4) for a, b in zip(up, down) if a['converged'] and b['converged']]
    one_sided = [a['alpha'] for a, b in zip(up, down) if a['converged'] != b['converged']]
    require(len(pairs) == 7 == summary['direction_pairs_within_screens'] and one_sided == [spec['path_angles_deg'][5]], 'Path coverage differs')
    require(all(abs(p['delta_CL']) <= screens['CL'] and abs(p['delta_CD_counts']) <= screens['CD_counts'] for p in pairs), 'Raw polars violate a screen')
    for a, b in zip(up, down):
        if not a['converged']: require(math.isnan(a['CL']) and math.isnan(a['CD']), 'Failed point carries a force')
    lift = [(a if a['converged'] else b)['CL'] for a, b in zip(up, down)]; drag = [(a if a['converged'] else b)['CD']*1e4 for a, b in zip(up, down)]
    require(all(y > x for x, y in zip(lift, lift[1:])) and sum((x-.175)*(y-.175) < 0 for x, y in zip(lift, lift[1:])) == 1
            and summary['convex_CL0175_admitted'], 'CL 0.175 admission differs from the raw polars')
    require(min(drag) > 105, 'Convex-nose path drag differs from the finding')
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TESTS], cwd=ROOT, text=True, capture_output=True)
    (OUT/'tests.log').write_text(tests.stdout+tests.stderr); require(tests.returncode == 0, 'Tests failed')
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, stage='V4 viscous check of the convex-nose reconstruction',
        retrieved_files_verified=audit['retrieved_files_verified'], saved_states_audited=audit['saved_states'],
        maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        nominal_control_identical_to_archive=audit['controls']['passed'],
        control_arrays_compared=sum(len(j['array_comparisons']) for j in audit['controls']['jobs'].values()),
        drag_errors_from_raw_polars_counts=errors, paired_drag_change_counts=paired,
        path_lift=lift, path_drag_counts=drag, direction_pairs_rechecked=pairs, one_sided_failure_alpha_deg=one_sided,
        failed_path_points=audit['failed_path_points'], failed_targets=audit['failed_targets'],
        tests=dict(modules=TESTS, passed=True), figures_reviewed=['drag_against_lift.png', 'incidence_response.png'],
        figure_review='Rendered and inspected: experiment and published MSES shown as reference markers; unreached states absent; '
                      'no value drawn where no route converged',
        new_flow_solves=0, output_sha256={p.name: sha(p) for p in outputs},
        qualification='Numerical records verified. Reconstruction sensitivity only; neither contour is the tested one.')
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('saved_states_audited', 'maximum_fresh_residual', 'nominal_control_identical_to_archive',
        'control_arrays_compared', 'drag_errors_from_raw_polars_counts', 'paired_drag_change_counts', 'one_sided_failure_alpha_deg')}, indent=1))


if __name__ == '__main__':
    main()
