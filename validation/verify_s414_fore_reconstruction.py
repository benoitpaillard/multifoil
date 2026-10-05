#!/usr/bin/env python3
"""Independent coverage, reference arithmetic, wake drag and artifact checks."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import re
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_fore_reconstruction_analysis_20261004'
STUDY = ROOT/'validation/results/s414_fore_reconstruction_20261004'
PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
SOURCE = ROOT/'validation/results/s414_trip_qualification_20261003'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def rows(p):
    with p.open() as stream: return list(csv.DictReader(stream))


def linear(points, target):
    assert all(b[0] > a[0] for a, b in zip(points, points[1:]))
    assert points[0][0] <= target <= points[-1][0]
    for (x0,y0),(x1,y1) in zip(points,points[1:]):
        if x0 <= target <= x1: return y0+(target-x0)*(y1-y0)/(x1-x0)
    raise ValueError('Missing reference support')


def main():
    manifest = read(OUT/'input_manifest.json')['files']
    for name, expected in manifest.items(): assert sha(ROOT/name) == expected, name
    summary = read(OUT/'summary.json'); audit = read(OUT/'residual_recheck.json')
    assert summary['analysis_complete'] and audit['complete']
    assert summary['converged_targets'] == summary['selected_states_freshly_audited'] == len(audit['states']) == 3
    assert summary['required_targets'] == 4 and not summary['target_coverage_complete']
    assert summary['requests'] == 28 and summary['failed_requests'] == 2
    assert summary['paired_drag_changes_counts']['fore_source_cubic'] is None
    assert summary['exact_control_arrays'] == 36 and summary['controls']['passed']
    assert max(r['max_residual'] for r in audit['states']) < 1e-7
    assert summary['maximum_growth_balance'] < 1e-7
    assert not any(summary[k] for k in ('exact_geometry_established','numerical_convergence_established','solution_uniqueness_established'))
    prior = read(PRIOR/'verification.json')
    assert sha(PRIOR/'input_manifest.json') == prior['output_sha256']['input_manifest.json']
    prior_inputs = read(PRIOR/'input_manifest.json')['files']
    independent_inputs = {}
    for name in ('experimental_tables.csv','published_mses_drag.csv'):
        path = SOURCE/name; relative = str(path.relative_to(ROOT))
        assert sha(path) == prior_inputs[relative]
        independent_inputs[relative] = sha(path)
    experiment = [(float(t[1]),float(t[2])) for t in (r['raw_tokens'].split() for r in rows(SOURCE/'experimental_tables.csv')
                  if r['reynolds']=='1000000' and r['configuration']=='free' and r['prestall_branch']=='True' and r['drag_available']=='True')]
    mses = [(float(r['CL']),float(r['CD'])) for r in rows(SOURCE/'published_mses_drag.csv')
            if r['reynolds']=='1000000' and r['configuration']=='free']
    force_checks = 0
    for row in summary['targets']:
        result = read(STUDY/row['job']/'result.json')
        assert row['feasible'] == result['feasible']
        for key, curve in (('experiment_CD',experiment),('mses_CD',mses)):
            assert abs(row[key]-linear(curve,row['target_cl'])) < 1e-14
        if not row['feasible']:
            assert all(row[k] is None for k in ('CL','CD','alpha_deg','error_CD_counts','fore_transition_x','fore_CD','aft_CD'))
            assert 'selected_flow_index' not in result
            continue
        assert row['CD'] == result['CD']
        assert abs(row['error_CD_counts']-(result['CD']-row['experiment_CD'])*1e4) < 1e-10
        force_checks += wake_checks(row['job'], result['selected_flow_index'], result['element_forces'])
    failed = read(STUDY/'cl0175_fore_source_cubic/result.json')
    metadata = read(STUDY/'cl0175_fore_source_cubic/flow/metadata.json')['diagnostics']
    requests = rows(OUT/'failed_target_requests.csv')
    assert len(requests)==11 and sum(r['admissible']=='True' for r in requests)==9
    for row in requests:
        i=int(row['evaluation']); original=failed['records'][i]
        if not original['admissible']:
            assert row['CL']==row['CD']==row['fore_transition_x']==''
            assert all(original[k] is None for k in ('CL','CD','CM'))
            continue
        assert abs(float(row['CD'])-original['CD']) < 1e-15
        force_checks += wake_checks('cl0175_fore_source_cubic', i, metadata[i]['element_forces'])
    nominal = read(STUDY/'cl0175_nominal/result.json')
    nommeta = read(STUDY/'cl0175_nominal/flow/metadata.json')['diagnostics']
    for row in summary['same_angle_endpoints']:
        i=row['evaluation']; a=nominal['records'][i]; b=failed['records'][i]
        assert a['alpha_deg']==b['alpha_deg']==row['alpha_deg']
        assert row['comparison_basis']=='same_angle_not_matched_lift'
        assert row['nominal_CL']==a['CL'] and row['cubic_CL']==b['CL']
        assert abs(row['delta_CD_counts']-(b['CD']-a['CD'])*1e4)<1e-10
        assert abs(row['delta_CD_counts']-row['fore_delta_CD_counts']-row['aft_delta_CD_counts'])<1e-10
        force_checks += wake_checks('cl0175_nominal',i,nommeta[i]['element_forces'])
    old=read(STUDY/'cl0056_nominal/result.json'); new=read(STUDY/'cl0056_fore_source_cubic/result.json')
    assert abs(summary['low_lift_comparison']['delta_CD_counts']-(new['CD']-old['CD'])*1e4)<1e-10
    assert abs(summary['paired_drag_changes_counts']['nominal']-(nominal['CD']-old['CD'])*1e4)<1e-10
    for label,curve in (('experiment',experiment),('published_mses',mses)):
        assert abs(summary['paired_drag_changes_counts'][label]-(linear(curve,.175)-linear(curve,.056))*1e4)<1e-10
    old_outputs=0
    for name,digest in prior['output_sha256'].items():
        assert sha(PRIOR/name)==digest, 'Changed old result: '+name
        old_outputs+=1
    tests=(OUT/'tests.log').read_text(); assert 'Ran 13 tests in ' in tests and tests.rstrip().endswith('OK')
    visual=read(OUT/'visual_review.json'); assert visual['passed'] and len(visual['images'])==3
    for name,digest in visual['images'].items(): assert sha(OUT/name)==digest
    links=0
    for name in re.findall(r'\]\(([^)]+)\)',(OUT/'FINDINGS.md').read_text()):
        if name.startswith(('https:','http:')) or name=='verification.json': continue
        assert (OUT/name).resolve().exists(),name
        links+=1
    code=[Path(__file__),ROOT/'test_s414_fore_reconstruction_analysis.py']
    receipt=dict(verified_utc=datetime.now(timezone.utc).isoformat(),analysis_complete=True,
        target_coverage_complete=False,inputs_verified=len(manifest),retrieved_files_verified=430,
        selected_states_freshly_audited=3,additional_request_exports_audited=11,
        independent_terminal_wake_force_checks=force_checks,prior_output_files_unchanged=old_outputs,
        exact_control_arrays=36,tests_passed=13,visual_review_pass=True,local_links_checked=links,
        new_nonlinear_flow_solves=0,independent_reference_inputs=independent_inputs,
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in code},
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'})
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print({k:v for k,v in receipt.items() if k not in ('output_sha256','code_sha256','independent_reference_inputs')})


def wake_checks(job,index,forces):
    for element in ('fore','aft'):
        last=rows(STUDY/job/f'flow/wake_{index:03d}_{element}.csv')[-1]
        theta=float(last['theta']); h=float(last['delta_star'])/theta
        cd=2*theta*float(last['ue_over_vinf'])**((5+h)/2)
        assert math.isfinite(cd) and abs(cd-forces[element]['CD_wake'])<1e-12
    return 2


if __name__=='__main__': main()
