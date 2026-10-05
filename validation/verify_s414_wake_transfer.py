#!/usr/bin/env python3
"""Verify saved diagnostic evidence and produce a receipt; never solve a flow."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'validation/results/s414_wake_transfer_20261003'
STUDY=ROOT/'validation/results/s414_reconstruction_20261002'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(test, message):
    if not test:raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def verify():
    manifest=read(OUT/'input_manifest.json');operators=read(OUT/'operator_checks.json')
    summary=read(OUT/'summary.json')
    for name,digest in manifest['files'].items():
        require(sha(ROOT/name)==digest,'Input changed: '+name)
    require(operators['input_sha256']==manifest['files'],'Operator/analysis inputs differ')
    require(operators['runtime']==manifest['runtime'],'Operator/analysis runtime differs')
    require(operators['complete'] and len(operators['cases'])==2,'Incomplete operator audit')
    require(summary['complete'] and summary['operator_states_checked']==2 and
            summary['geometry_states_checked']==8 and summary['central_difference_checks']==16,
            'Incomplete diagnostic coverage')
    require(summary['new_flow_solves']==operators['new_flow_solves']==0,'Unexpected flow solve')
    require(summary['ncrit']==9 and not summary['missing_channel_proven_cause'],
            'Scientific scope changed')
    initial=read(OUT/'initial_guard_receipt/operator_checks.json')
    require(initial['cases']==operators['cases'],'Dependency-guard rerun changed numerical evidence')
    require(initial['complete'],'Initial numerical evidence was incomplete')
    spec=read(STUDY/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(ROOT/name)==sha(STUDY/'source_snapshot'/name)==digest,'Frozen study changed: '+name)
    retrieved=read(ROOT/'validation/results/s414_reconstruction_launch_20261002/retrieved_files_sha256.json')
    for name,digest in retrieved.items():
        require(sha(STUDY/name)==digest,'Archived result changed: '+name)
    # The helper was already used and pinned in the completed pressure audit.
    prior=read(ROOT/'validation/results/s414_reconstruction_analysis_20261003/input_manifest.json')
    helper='validation/report_s702_matched_pressure.py'
    require(manifest['files'][helper]==prior['files'][helper],'Pressure/branch helper changed')
    proximity=list(csv.DictReader((OUT/'wake_proximity.csv').open()))
    contexts=list(csv.DictReader((OUT/'fore_layer_and_wake_context.csv').open()))
    require(len(contexts)==8 and len(set(r['job'] for r in contexts))==8,'Geometry coverage differs')
    require(all(r['all_fore_wake_nodes_flagged_turbulent']=='True' for r in contexts),'Wake classification differs')
    diagnostics=[];max_sum_error=0.;max_identity_error=0.
    for case in operators['cases']:
        require(case['states_restored_exactly'] and case['baseline_residual']<=1e-7,'Bad restored state')
        require(case['new_flow_solves']==0 and case['ncrit']==9.,'Operator scope changed')
        require(len(case['finite_differences'])==8 and len(case['blocks'])==4,'Missing perturbations')
        for block in case['blocks']:
            require(block['aft_BL_block_max_abs']==0,'Unexpected cross-element BL block')
            active=block['variable'] in ('delta_star','edge_velocity')
            require((block['aft_velocity_block_max_abs']>0)==active,'Mean-flow Jacobian channel differs')
        require({(r['variable'],r['relative_step']) for r in case['finite_differences']}==
                {(v,h) for v in ('theta','delta_star','root_shear','edge_velocity') for h in (.001,.0001)},
                'Variable/step coverage differs')
        for row in case['finite_differences']:
            require(row['max_aft_BL_residual_change']==0 and row['max_fore_BL_residual_change']>1e-10,
                    'Zero channel or positive control fails')
            active=row['variable'] in ('delta_star','edge_velocity')
            require((row['max_aft_velocity_residual_change']>0)==active,'Finite-difference channel differs')
            require(max(row['aft_BL_Jacobian_relative_discrepancy'],
                        row['aft_velocity_Jacobian_relative_discrepancy'])<1e-7,'Directional derivative fails')
        columns=('fore_surface','fore_wake','aft_surface','aft_wake')
        rows=sorted(case['source_columns'],key=lambda r:r['x_over_c'])
        for row in rows:
            correction=sum(row[k] for k in columns)
            max_sum_error=max(max_sum_error,abs(correction-row['summed_displacement_velocity']))
            max_identity_error=max(max_identity_error,abs(row['edge_velocity_state']-
                row['inviscid_edge_velocity']-correction-row['coupling_residual']))
        x=[r['x_over_c'] for r in rows];target=.89313
        require(x[0]<=target<=x[-1],'Source-value extrapolation')
        parts={k:float(np.interp(target,x,[r[k] for r in rows])) for k in
               (*columns,'summed_displacement_velocity','edge_velocity_state','inviscid_edge_velocity')}
        scales=[]
        for job in sorted({r['job'] for r in proximity if r['case']==case['case']}):
            rows=sorted([r for r in proximity if r['job']==job],key=lambda r:float(r['x_over_c']))
            x=[float(r['x_over_c']) for r in rows]
            require(x[0]<=target<=x[-1],'Proximity extrapolation')
            scales.append(dict(job=job,**{k:float(np.interp(target,x,[float(r[k]) for r in rows]))
                for k in ('centerline_distance_c','fore_wake_delta_star_c','fore_wake_theta_c',
                          'aft_delta_star_c','distance_over_displacement_sum','fore_wake_root_shear')}))
        require(len(scales)==4,'Scenario coverage differs')
        diagnostics.append(dict(case=case['case'],x_over_c=target,nominal_source_contributions=parts,
                                geometry_scenarios=scales,method='Piecewise linear interpolation of saved tables'))
    require(max_sum_error<=2e-12 and max_identity_error<=2e-12,'Source-column identities fail')
    tests=(OUT/'tests.log').read_text()
    require(re.search(r'Ran 25 tests in ',tests) and tests.rstrip().endswith('OK'),'Targeted tests not passed')
    links=0
    for doc in (OUT/'FINDINGS.md',ROOT/'validation/S414_WAKE_TRANSFER_PROTOCOL.md'):
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
            if target.startswith(('http:','https:','#')):continue
            if (doc.parent/target).resolve()==OUT/'verification.json':continue  # Produced after these checks.
            require((doc.parent/target.split('#')[0]).exists(),'Broken document link: '+target)
            links+=1
    (OUT/'diagnostic_values.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    receipt=dict(verified_utc=datetime.now(timezone.utc).isoformat(),
        checks=dict(input_files=len(manifest['files']),all_input_hashes_match=True,
            frozen_study_files=len(spec['sha256']),frozen_study_and_snapshots_match=True,
            retrieved_files=len(retrieved),all_retrieved_hashes_match=True,
            restored_nominal_states=2,geometry_states=8,finite_difference_checks=16,
            maximum_receiving_Jacobian_relative_discrepancy=summary['largest_receiving_Jacobian_relative_discrepancy'],
            source_column_sum_error=max_sum_error,source_residual_identity_error=max_identity_error,
            dependency_guard_rerun_exact=True,states_restored_exactly=True,targeted_tests=25,local_document_links=links),
        scientific_limits=dict(new_flow_solves=0,ncrit=9.,causal_error_attribution=False,
            turbulence_intensity_recovered=False,measured_wake_overlap=False,
            numerical_convergence_established=False,exact_tested_geometry=False,new_MSES_runs=0),
        source_access='Reused checksum-verified public source PDFs; no new acquisition',
        primary_sources=[dict(url='https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf',
            local='validation/results/s414_geometry_qualification_20261002/sources/s414_report_2010.pdf',
            scope='2010 report printed pp.9–11; Figure 7; free/fore-only/all-fixed configurations'),
            dict(url='https://www.ae.metu.edu.tr/tuncer/ae546/multi/s414-1.pdf',
            local='validation/results/reference_review_20260924/sources/s414_article.pdf',scope='2014 article section VI')],
        diagnostic_runtime=manifest['runtime'])
    # Final artifact hashes make content/visual review traceable without altering older receipts.
    receipt['output_sha256']={p.name:sha(p) for p in sorted(OUT.iterdir())
        if p.is_file() and p.name!='verification.json'}
    receipt['documentation_and_test_sha256']={str(p.relative_to(ROOT)):sha(p) for p in
        (Path(__file__).resolve(),ROOT/'test_s414_wake_transfer.py',ROOT/'AGENTS.md',ROOT/'DEVELOPMENT_PLAN.md')}
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt['checks'],indent=2))


if __name__=='__main__':verify()
