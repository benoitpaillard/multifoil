#!/usr/bin/env python3
"""Verify the diagnostic receipt without altering prior studies or running flow."""
import ast
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_fore_nose_20261003'
PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
PRIOR_HASH = '3325127f24a971f345f9857529aa250ed6e7b05d10e588e93b2472720e03489f'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def rows(name):
    with (OUT/name).open() as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = read(OUT/'input_manifest.json')
    require(manifest['prior_verification_sha256'] == PRIOR_HASH == sha(PRIOR/'verification.json'),
            'Prior audit receipt changed')
    for name, digest in manifest['files'].items():
        require(sha(ROOT/name) == digest, 'Frozen diagnostic input changed: '+name)
    old_outputs = 0
    for study in ('s414_lower_bucket_analysis_20261003', 's414_reconstruction_analysis_20261003',
                  's414_wake_transfer_20261003', 's414_trip_qualification_20261003'):
        folder = ROOT/'validation/results'/study
        for name, digest in read(folder/'verification.json')['output_sha256'].items():
            require(sha(folder/name) == digest, 'Previous audited output changed: '+str(folder/name))
            old_outputs += 1

    summary = read(OUT/'summary.json')
    require(summary['complete'] and [c['case'] for c in summary['cases']] == ['cl0056', 'cl0175', 'cl028'],
            'Incomplete declared cases')
    require(summary['ncrit'] == 9 and summary['prior_exact_repeat'], 'Baseline changed')
    require(summary['new_nonlinear_flow_solves'] == 0 and summary['unchanged_inviscid_operator_rebuilds'] == 1
            and summary['archived_incidence_inviscid_evaluations'] == 3, 'Scope changed')
    require(not summary['exact_tested_geometry_established'] and not summary['numerical_convergence_established'],
            'Unsupported validation claim')
    require(summary['amplification_balances_and_transition_roots_pass'], 'Balance check failed')
    for case in summary['cases']:
        require(case['state_and_flags_restored_exactly'] and case['frozen_march_same_transition_cell'],
                'State restoration or frozen march failed')
        require(case['maximum_amplification_balance'] <= 1e-7
                and case['frozen_march_maximum_N_difference'] <= 1e-6, 'Amplification check failed')
        if case['transition'] is not None:
            t = case['transition']
            require(abs(t['xi']-t['archived_xi']) <= 2e-11 and abs(t['root_residual']) <= 1e-10
                    and 0 < t['fraction_of_cell'] < 1, 'Transition root check failed')
        iv = [r for r in rows('amplification_intervals.csv') if r['case'] == case['case']]
        for column, key in [('regularization_trapezoid', 'total_regularization_N'),
                            ('base_trapezoid', 'total_base_N'),
                            ('gauss_minus_trapezoid', 'frozen_state_gauss_minus_trapezoid')]:
            require(abs(sum(float(r[column]) for r in iv)-case[key]) < 1e-12, 'Growth accounting differs')
        require(max(abs(float(r['balance'])) for r in iv) == case['maximum_amplification_balance'],
                'Balance summary differs')
    require(summary['cases'][0]['N_at_x0034'] is None and summary['cases'][2]['transition'] is None,
            'Missing values silently filled')
    for r in rows('fore_lower_stations.csv'):
        laminar = r['laminar'] == 'True'
        require(bool(r['N']) == laminar and bool(r['shear']) != laminar
                and bool(r['dN_dxi']) == laminar, 'Shear interpreted as amplification')

    derivative = rows('growth_derivative_checks.csv')
    require(len(derivative) == 72 and {r['relative_step'] for r in derivative} == {'1e-05', '1e-06'},
            'Derivative coverage changed')
    original_errors = [float(r['scaled_error']) for r in derivative]
    corrected_errors = [float(r['corrected_scaled_error']) for r in derivative]
    for r in derivative:
        for field, error in [('analytic', 'scaled_error'), ('diagnostic_corrected', 'corrected_scaled_error')]:
            a, fd = float(r[field]), float(r['finite_difference'])
            require(abs(abs(fd-a)/max(1, abs(a))-float(r[error])) < 1e-14, 'Derivative metric changed')
    failed = sum(e >= 2e-5 for e in original_errors)
    require(failed == 18 and not summary['original_growth_derivative_check_pass'], 'Baseline failure lost')
    require(max(original_errors) == summary['maximum_original_growth_derivative_error'], 'Failure magnitude differs')
    require(max(corrected_errors) == summary['maximum_corrected_growth_derivative_error'] < 2e-5
            and summary['diagnostic_corrected_derivative_check_pass'], 'Diagnostic correction failed')
    failure = read(OUT/'initial_derivative_check_failed/failure.json')
    require(not failure['passed'] and failure['criterion'] == 2e-5 and not failure['production_code_changed'],
            'Original failure record changed')
    require((OUT/'initial_derivative_check_failed/report_s414_fore_nose.py').is_file(), 'Failed script missing')

    # Compare the recorded author expressions to the local AST, not to an unavailable author file.
    source_check = read(OUT/'sources/author_source_check.json')
    routine = next(n for n in ast.parse((ROOT/'vendor/mfoil.py').read_text()).body
                   if isinstance(n, ast.FunctionDef) and n.name == 'get_damp')
    expressions = {t.id: n.value for n in ast.walk(routine) if isinstance(n, ast.Assign)
                   for t in n.targets if isinstance(t, ast.Name)}
    for name, expression in source_check['checked_expressions'].items():
        require(ast.dump(expressions[name]) == ast.dump(ast.parse(expression, mode='eval').body),
                'Vendored expression differs from author expression record')
    require(not source_check['whole_file_hash_available'], 'Unsupported whole-file comparison')
    pressure = rows('pressure_at_taps.csv')
    require(len(pressure) == 210 and len(rows('archived_cubic_control.csv')) == 70, 'Pressure coverage changed')
    require({r['coordinate'] for r in pressure} == {'table', 'plotted'}, 'Coordinate sensitivity lost')

    test_log = (OUT/'tests.log').read_text()
    require('Ran 23 tests in ' in test_log and test_log.rstrip().endswith('OK'), 'Targeted tests not passed')
    visual = read(OUT/'visual_review.json')
    require(visual['passed'] and len(visual['images']) == 3, 'Visual review incomplete')
    for name, digest in visual['images'].items():
        require(sha(OUT/name) == digest, 'Reviewed image changed: '+name)
    links = 0
    for target in re.findall(r'\]\(([^)]+)\)', (OUT/'FINDINGS.md').read_text()):
        if target.startswith(('https:', 'http:', '#')) or target == 'verification.json':
            continue
        require((OUT/target).resolve().exists(), 'Broken report link: '+target)
        links += 1
    output_hashes = {str(p.relative_to(OUT)): sha(p) for p in sorted(OUT.rglob('*'))
                     if p.is_file() and p != OUT/'verification.json'}
    code = [Path(__file__), ROOT/'test_s414_fore_nose.py']
    receipt = dict(verified_utc=datetime.now(timezone.utc).isoformat(), complete=True,
        input_files_verified=len(manifest['files']), prior_output_files_unchanged=old_outputs,
        original_growth_derivative_check_pass=False, original_derivative_failed_components=failed,
        derivative_component_step_checks=len(derivative), diagnostic_corrected_derivatives_pass=True,
        amplification_balances_and_transition_roots_pass=True, targeted_tests_passed=23,
        local_report_links_checked=links, visual_review_pass=True, new_nonlinear_flow_solves=0,
        code_sha256={str(p.relative_to(ROOT)): sha(p) for p in code}, output_sha256=output_hashes)
    (OUT/'verification.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print({k: v for k, v in receipt.items() if k not in ('code_sha256', 'output_sha256')})


if __name__ == '__main__':
    main()
