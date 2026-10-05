#!/usr/bin/env python3
"""Verify archived inputs, independent drag arithmetic and completed report artifacts."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
SOURCE = ROOT/'validation/results/s414_trip_qualification_20261003'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def linear(points, target):
    require(all(b[0] > a[0] for a, b in zip(points, points[1:])), 'Folded reference')
    if not points[0][0] <= target <= points[-1][0]:
        return None
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= target <= x1:
            return y0+(target-x0)/(x1-x0)*(y1-y0)
    raise ValueError('Uncovered reference')


def main():
    manifest = read(OUT/'input_manifest.json')
    for name, digest in manifest['files'].items():
        require(sha(ROOT/name) == digest, 'Analysis input changed: '+name)
    summary = read(OUT/'summary.json')
    residual = read(OUT/'residual_recheck.json')
    require(summary['converged_targets'] == summary['required_targets'] == 4 and residual['complete'], 'Incomplete audit')
    require(len(residual['states']) == 4 and all(r['max_residual'] <= 1e-7 for r in residual['states']), 'Residual failure')
    require(summary['new_flow_solves'] == residual['new_flow_solves'] == 0, 'Scope changed')
    require(summary['supported_physical_pressure_observations'] == 348 and summary['raw_pressure_markers'] == 376, 'Pressure coverage changed')
    require(summary['mses_pressure_cases'] == ['cl028'] and len(summary['mses_drag_cases']) == 3, 'MSES coverage changed')
    require(not summary['exact_tested_geometry_established'] and not summary['numerical_convergence_established'], 'Unsupported qualification')
    source_rows = [r for r in rows(SOURCE/'experimental_tables.csv') if r['reynolds'] == '1000000'
                   and r['configuration'] == 'free' and r['prestall_branch'] == 'True' and r['drag_available'] == 'True']
    # Read original printed tokens independently of the float fields used by the report.
    experiment = [(float(t[1]), float(t[2])) for t in (r['raw_tokens'].split() for r in source_rows)]
    theory = [(float(r['CL']), float(r['CD'])) for r in rows(SOURCE/'published_mses_drag.csv')
              if r['reynolds'] == '1000000' and r['configuration'] == 'free']
    maximum_error = 0.
    for target in summary['targets']:
        for method, curve in [('experiment', experiment), ('mses', theory)]:
            expected = linear(curve, target['target_cl'])
            actual = target[method+'_CD']
            require((expected is None) == (actual is None), 'Missing reference silently filled')
            if expected is not None:
                maximum_error = max(maximum_error, abs(expected-actual))
        require(abs((target['CD']-target['experiment_CD'])*1e4-target['error_CD_counts']) < 1e-10, 'Drag subtraction differs')
    require(maximum_error < 1e-14, 'Independent drag interpolation differs')
    old_outputs = 0
    for study in ('s414_reconstruction_analysis_20261003', 's414_wake_transfer_20261003', 's414_trip_qualification_20261003'):
        folder = ROOT/'validation/results'/study
        receipt = read(folder/'verification.json')
        for name, digest in receipt['output_sha256'].items():
            require(sha(folder/name) == digest, 'Previous audited output changed: '+str(folder/name))
            old_outputs += 1
    require(not read(SOURCE/'summary.json')['strict_repeated_pressure_copy_check_pass'], 'Source precision failure lost')
    test_log = (OUT/'tests.log').read_text()
    require('Ran 38 tests in ' in test_log and test_log.rstrip().endswith('OK'), 'Targeted tests not passed')
    visual = read(OUT/'visual_review.json')
    require(visual['passed'] and len(visual['images']) == 6, 'Visual review incomplete')
    for name, digest in visual['images'].items():
        require(sha(OUT/name) == digest, 'Reviewed image changed')
    links = 0
    for target in re.findall(r'\]\(([^)]+)\)', (OUT/'FINDINGS.md').read_text()):
        if target.startswith(('https:', 'http:', '#')):
            continue
        if target == 'verification.json':
            continue  # This receipt is produced below.
        require((OUT/target).resolve().exists(), 'Broken report link: '+target)
        links += 1
    artifact_hashes = {p.name: sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'verification.json'}
    code = [Path(__file__), ROOT/'test_s414_lower_bucket_analysis.py']
    receipt = dict(verified_utc=datetime.now(timezone.utc).isoformat(), complete=True,
        input_files_verified=len(manifest['files']), prior_output_files_unchanged=old_outputs,
        independent_drag_maximum_difference=maximum_error, targeted_tests_passed=38,
        local_report_links_checked=links, visual_review_pass=True, new_flow_solves=0,
        code_sha256={str(p.relative_to(ROOT)): sha(p) for p in code}, output_sha256=artifact_hashes)
    (OUT/'verification.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print({k: v for k, v in receipt.items() if k not in ('code_sha256', 'output_sha256')})


if __name__ == '__main__':
    main()
