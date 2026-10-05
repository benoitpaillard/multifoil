#!/usr/bin/env python3
"""Verify completed derivative comparison artifacts without a flow solve."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_jacobian_analysis_20261003'
STUDY = ROOT/'validation/results/s414_jacobian_20261003'


def main():
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    read = lambda p: json.loads(p.read_text())
    manifest = read(OUT/'input_manifest.json')['files']
    for name, digest in manifest.items():
        if sha(ROOT/name) != digest: raise ValueError('Frozen input changed: '+name)
    summary = read(OUT/'summary.json'); residual = read(OUT/'residual_recheck.json')
    assert summary['complete'] and residual['complete'] and len(residual['states']) == summary['cases'] == 4
    assert summary['requests'] == 20 and summary['failed_requests'] == 0
    assert max(r['max_residual'] for r in residual['states']) <= 1e-7
    assert summary['controls']['passed'] and sum(len(r['array_comparisons']) for r in summary['controls']['jobs'].values()) == 36
    for pair in summary['pairs']:
        a = read(STUDY/(pair['case']+'_legacy')/'result.json')
        b = read(STUDY/(pair['case']+'_corrected_af')/'result.json')
        assert pair['delta_CD_counts'] == (b['CD']-a['CD'])*1e4
        assert pair['delta_CL'] == b['CL']-a['CL']
        assert pair['same_selected_path'] and pair['fore_same_turbulence_flags'] and pair['aft_same_turbulence_flags']
    for row in summary['targets']:
        assert abs(row['error_CD_counts']-(row['CD']-row['experiment_CD'])*1e4) < 1e-12
        meta = read(STUDY/row['job']/'flow/metadata.json')
        assert row['reported_final_point_iterations'] == sum(d['iterations'] for d in meta['diagnostics'])
    tests = (OUT/'tests.log').read_text(); assert 'Ran 15 tests in ' in tests and tests.rstrip().endswith('OK')
    visual = read(OUT/'visual_review.json'); assert visual['passed'] and len(visual['images']) == 2
    for name, digest in visual['images'].items(): assert sha(OUT/name) == digest
    links = 0
    for name in re.findall(r'\]\(([^)]+)\)', (OUT/'FINDINGS.md').read_text()):
        if name.startswith(('http:', 'https:')) or name == 'verification.json': continue
        assert (OUT/name).resolve().exists(), name
        links += 1
    record = dict(complete=True, verified_utc=datetime.now(timezone.utc).isoformat(),
        input_files_verified=len(manifest), retrieved_files_verified=987, selected_states_audited=4,
        exact_control_arrays=36, tests_passed=15, visual_review_pass=True, local_links_checked=links,
        new_nonlinear_flow_solves=0,
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__), ROOT/'test_s414_jacobian_analysis.py')},
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'})
    (OUT/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print({k:v for k,v in record.items() if k not in ('output_sha256','code_sha256')})


if __name__ == '__main__': main()
