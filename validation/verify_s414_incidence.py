#!/usr/bin/env python3
"""Check the measured-incidence analysis, coverage and comparison bases."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_incidence_analysis_20261003'
STUDY = ROOT/'validation/results/s414_incidence_20261003'


def read(path):
    return json.loads(path.read_text())


def rows(name):
    with (OUT/name).open() as stream: return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = read(OUT/'input_manifest.json')['files']
    for name, digest in manifest.items():
        if sha(ROOT/name) != digest: raise ValueError('Frozen input changed: '+name)
    summary = read(OUT/'summary.json'); residual = read(OUT/'residual_recheck.json')
    assert summary['complete'] and summary['full_state_audit_pass'] and residual['complete']
    assert summary['cases'] == len(residual['states']) == 4 and summary['requests'] == 12 and summary['failed_requests'] == 0
    assert max(r['max_residual'] for r in residual['states']) == summary['maximum_fresh_residual'] <= 1e-7
    assert summary['controls']['passed'] and sum(len(r['array_comparisons']) for r in summary['controls']['jobs'].values()) == 36
    for job in read(STUDY/'study.json')['jobs']:
        actual = read(STUDY/job['name']/'result.json')
        target = next(t for t in summary['targets'] if t['job'] == job['name'])
        assert target['CL'] == actual['CL'] and target['CD_context'] == actual['CD']
        assert target['CL_minus_caption'] == actual['CL']-target['reference_caption_CL']
        if job['comparison_basis'] == 'fixed_alpha':
            assert abs(actual['alpha_deg']-job['target_alpha']) <= 1e-12
            assert target['matched_lift_CD_error_counts'] is None and 'target_cl' not in job
        else:
            assert abs(actual['CL']-job['target_cl']) <= 2e-5
            assert target['matched_lift_CD_error_counts'] == (actual['CD']-target['matched_lift_experiment_CD'])*1e4
    pressure = rows('pressure_at_taps.csv')
    included = [r for r in pressure if r['included'] == 'True']
    assert len(pressure) == summary['pressure_comparison_rows'] == 752
    assert len(included) == summary['scored_pressure_rows'] == 696
    assert len({(r['case'], r['surface'], r['table_index']) for r in included}) == summary['unique_supported_pressure_observations'] == 174
    for row in rows('pressure_metrics.csv'):
        selected = [r for r in included if r['job']==row['job'] and r['coordinate']==row['coordinate']
            and (row['surface']=='all_surfaces' or r['surface']==row['surface'])
            and (row['region']=='all' or r['surface'].startswith('fore_') and float(r['x_over_c'])<=.06)]
        assert len(selected) == int(row['scored_taps'])
        errors = [float(r['model_Cp'])-float(r['experiment_Cp']) for r in selected]
        assert abs((sum(e*e for e in errors)/len(errors))**.5-float(row['rmse'])) < 1e-13
    assert summary['mses_pressure_cases'] == ['cl028']
    assert not summary['exact_geometry_established'] and not summary['numerical_convergence_established'] and not summary['solution_uniqueness_established']
    tests = (OUT/'tests.log').read_text(); assert 'Ran 18 tests in ' in tests and tests.rstrip().endswith('OK')
    visual = read(OUT/'visual_review.json'); assert visual['passed'] and len(visual['images']) == 3
    for name, digest in visual['images'].items(): assert sha(OUT/name) == digest
    links = 0
    for name in re.findall(r'\]\(([^)]+)\)', (OUT/'FINDINGS.md').read_text()):
        if name.startswith(('http:', 'https:')) or name == 'verification.json': continue
        assert (OUT/name).resolve().exists(), name
        links += 1
    previous = ROOT/'validation/results/s414_jacobian_analysis_20261003'
    prior = read(previous/'verification.json')
    for name, digest in prior['output_sha256'].items(): assert sha(previous/name) == digest
    code = [Path(__file__), ROOT/'test_s414_incidence_analysis.py']
    receipt = dict(complete=True, verified_utc=datetime.now(timezone.utc).isoformat(),
        input_files_verified=len(manifest), frozen_study_inputs_verified=1879, retrieved_files_verified=1977,
        selected_states_audited=4, exact_control_arrays=36, tests_passed=18,
        pressure_metrics_independently_recomputed=True, previous_audited_outputs_unchanged=len(prior['output_sha256']),
        visual_review_pass=True, local_report_links_checked=links, new_nonlinear_flow_solves=0,
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in code},
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'})
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print({k:v for k,v in receipt.items() if k not in ('code_sha256','output_sha256')})


if __name__ == '__main__': main()
