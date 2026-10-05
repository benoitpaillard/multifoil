#!/usr/bin/env python3
"""Check the completed S414 source qualification and preserve an audit receipt."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'validation/results/s414_trip_qualification_20261003'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def require(condition,message):
    if not condition:raise ValueError(message)


def main():
    manifest=read(OUT/'source_manifest.json');summary=read(OUT/'summary.json')
    for name,digest in manifest['input_sha256'].items():require(sha(ROOT/name)==digest,'Changed input: '+name)
    for name,digest in manifest['output_sha256'].items():require(sha(OUT/name)==digest,'Changed extracted output: '+name)
    require(summary['complete'] and summary['new_flow_solves']==0,'Invalid completion/scope')
    require((summary['experiment_table_rows'],summary['independent_plot_table_checks'],
             summary['free_pressure_markers'],summary['repeated_pressure_checks'])==(169,108,376,94),'Incomplete coverage')
    require(not summary['strict_repeated_pressure_copy_check_pass'],'Original precision failure lost')
    require(not summary['causal_wake_transfer_benchmark_admitted'],'Causal admission incorrectly promoted')
    evidence=read(OUT/'evidence_qualification.json')
    require(not evidence['pressure_coverage']['fore_only_fixed_Cp_found_in_report'],'Pressure availability changed')
    require(evidence['mses']['transition_method']=='2010 envelope e^N','MSES methods conflated')
    with (OUT/'matched_lift_targets.csv').open() as stream:targets=list(csv.DictReader(stream))
    with (OUT/'experimental_tables.csv').open() as stream:rows=list(csv.DictReader(stream))
    require(len(targets)==8,'Incomplete target coverage')
    max_error=0.
    # Independent interpolation from raw printed tokens, not stored float fields.
    for target in targets:
        re_value=int(target['reynolds']);cl=float(target['CL']);values={}
        for state in ('free','fore_fixed','both_fixed'):
            branch=[r for r in rows if int(r['reynolds'])==re_value and r['configuration']==state and r['prestall_branch']=='True']
            tokens=[r['raw_tokens'].split() for r in branch]
            x=np.array([float(r[1]) for r in tokens]);y=np.array([float(r[2]) for r in tokens])
            require(np.all(np.diff(x)>0) and x[0]<=cl<=x[-1],'Invalid independent interpolation support')
            values[state]=float(np.interp(cl,x,y))
            max_error=max(max_error,abs(values[state]-float(target['experiment_'+state+'_CD'])))
        require(abs(1e4*(values['fore_fixed']-values['free'])-float(target['experiment_fore_minus_free_counts']))<1e-10,
                'Trip delta does not reproduce from original table tokens')
    require(max_error<1e-14,'Printed-token interpolation differs')
    old_outputs=0
    for name in ('s414_reconstruction_analysis_20261003','s414_wake_transfer_20261003'):
        folder=ROOT/'validation/results'/name;receipt=read(folder/'verification.json')
        for name,digest in receipt['output_sha256'].items():
            require(sha(folder/name)==digest,'Prior audited output changed: '+str(folder/name));old_outputs+=1
    study=ROOT/'validation/results/s414_reconstruction_20261002';spec=read(study/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(ROOT/name)==sha(study/'source_snapshot'/name)==digest,'Frozen solver/input changed: '+name)
    retrieved=read(ROOT/'validation/results/s414_reconstruction_launch_20261002/retrieved_files_sha256.json')
    for name,digest in retrieved.items():require(sha(study/name)==digest,'Old simulation output changed: '+name)
    tests=(OUT/'tests.log').read_text()
    require('Ran 19 tests in ' in tests and tests.rstrip().endswith('OK'),'Targeted tests not passed')
    links=0
    for target in re.findall(r'\]\(([^)]+)\)',(OUT/'FINDINGS.md').read_text()):
        if target.startswith(('https:','http:','#')):continue
        if target=='verification.json':continue
        require((OUT/target.split('#')[0]).exists(),'Broken findings link: '+target);links+=1
    # These figures were inspected by the project agent before this receipt.
    visual=read(OUT/'visual_review.json')
    for name,digest in visual['reviewed_png_sha256'].items():require(sha(OUT/name)==digest,'Reviewed figure changed')
    receipt=dict(verified_utc=datetime.now(timezone.utc).isoformat(),
        checks=dict(input_hashes=len(manifest['input_sha256']),generated_artifact_hashes=len(manifest['output_sha256']),
            independent_table_rows=169,table_marker_checks=108,pressure_markers=376,duplicate_pressure_checks=94,
            matched_lift_targets=8,independent_table_token_interpolation_error=max_error,targeted_tests=19,
            prior_audited_output_hashes=old_outputs,frozen_solver_and_snapshot_files=len(spec['sha256']),
            retrieved_simulation_hashes=len(retrieved),local_findings_links=links),
        retained_failed_checks=dict(strict_pressure_copy=False,maximum_Cp_offset=summary['maximum_repeated_pressure_Cp_difference'],
            maximum_x_offset=summary['maximum_repeated_pressure_x_difference']),
        admission=evidence['admission'],new_flow_solves=0,reviewed_png_sha256=visual['reviewed_png_sha256'],
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'},
        current_documentation_and_test_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__).resolve(),ROOT/'test_s414_trips.py',ROOT/'AGENTS.md',ROOT/'DEVELOPMENT_PLAN.md',ROOT/'validation/REFERENCES.md')})
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt['checks'],indent=2))


if __name__=='__main__':main()
