#!/usr/bin/env python3
"""Audit the frozen four-case S414 sequence without a nonlinear flow solve."""
import json
import platform
import importlib.metadata
from audit_s414_state import recompute, AUDIT_CODE as BASE_CODE
from report_s702_wake_extent import ROOT, Inputs, sha, audit, require
from run_s414_lower_bucket import compare_control
from s702_ncrit_sweep import write

STUDY = ROOT/'validation/results/s414_lower_bucket_20261003'
OUT = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
LAUNCH = ROOT/'validation/results/s414_lower_bucket_launch_20261003'
STUDY_HASH = '2c1d57eba8101dcbc8283f23032447bd18e656dc3f18812a8d0b4b17c82f8cec'
AUDIT_CODE = BASE_CODE + ['validation/audit_s414_lower_bucket.py']


def load_states(inputs):
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Unexpected study')
    spec = inputs.json(STUDY/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    receipt = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in receipt.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved output changed: '+name)
    retrieval = inputs.json(LAUNCH/'retrieval.json')
    require(retrieval['all_files_match'] and retrieval['files_verified'] == len(receipt) == 238,
            'Retrieval coverage changed')
    runtime = inputs.json(STUDY/'runtime.json')
    status = inputs.json(STUDY/'status.json')
    automatic = inputs.json(STUDY/'summary.json')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == 4,
            'Incomplete execution')
    require(automatic['converged_jobs'] == automatic['required_jobs'] == 4
            and automatic['execution_checks_pass'], 'Incomplete coverage')
    require(runtime['python'] == '3.12.13' and runtime['versions'] == spec['expected_versions'],
            'Solver environment changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread limits changed')
    # Pure comparison: never rewrite controls.json in the archived study.
    controls = compare_control(STUDY/'cl028_control', ROOT/spec['reference_control'],
                               spec['acceptance_control_CD_counts'])
    require(controls['passed'] and controls == automatic['controls'] == inputs.json(STUDY/'controls.json'),
            'Exact control gate differs')
    states = {j['name']: audit(STUDY/j['name'], j, spec, runtime, inputs) for j in spec['jobs']}
    require(all(s['result']['feasible'] for s in states.values()), 'Incomplete comparison')
    for name in AUDIT_CODE:
        inputs.add(ROOT/name)
    return spec, states, controls


def main():
    OUT.mkdir(exist_ok=True)
    inputs = Inputs()
    spec, states, controls = load_states(inputs)
    environment = dict(python=platform.python_version(),
        versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')})
    receipt = dict(input_sha256=inputs.files, audit_runtime=environment, states=[], complete=False,
        new_flow_solves=0, controls=controls,
        method='Four full export/state audits and fresh coupled residual/force/Cp/Cf/BL evaluations; no nonlinear solve')
    path = OUT/'residual_recheck.json'
    if path.exists():
        previous = json.loads(path.read_text())
        require(previous['input_sha256'] == receipt['input_sha256']
                and previous['audit_runtime'] == environment, 'Audit changed; preserve previous receipt')
        receipt = previous
    for job in spec['jobs']:
        if any(row['job'] == job['name'] for row in receipt['states']):
            continue
        row = recompute(job, spec, states[job['name']])
        receipt['states'].append(row)
        write(path, receipt)
        print(f"{len(receipt['states'])}/4 {job['name']}: residual {row['max_residual']:.3g}", flush=True)
    receipt['complete'] = True
    write(path, receipt)


if __name__ == '__main__':
    main()
