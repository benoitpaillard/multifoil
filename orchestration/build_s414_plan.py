#!/usr/bin/env python3
"""Write the reviewable S414 pilot plan; does not freeze, test or launch it."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_nose_resolution as runner


def build():
    required = ['validation/audit_s414_nose_resolution.py', 'validation/report_s414_nose_resolution.py',
                'test_s414_nose_resolution_analysis.py']
    for name in required:
        if not (ROOT/name).is_file():
            raise ValueError('Finish the analysis dependency before freezing: '+name)
    numerical = runner.specification()  # File/hash checks only; no solve.
    files = set(numerical['sha256']) | {
        'LICENSE', 'THIRD_PARTY.md', 'AGENTS.md',
        'CODEX_HANDOVER.md', 'DEVELOPMENT_PLAN.md', 'validation/S414_VERIFICATION_STATUS.md',
        'validation/S414_V5_NOTES.md',
        'orchestration/pilot.py', 'orchestration/test_pilot.py',
        'orchestration/build_s414_plan.py', 'orchestration/README.md',
        'validation/results/s414_post_verification_review_20261005/FINDINGS.md'}
    study, launch = runner.STUDY, runner.LAUNCH
    analysis = 'validation/results/s414_nose_resolution_analysis_20261005'
    python = sys.executable
    script = 'validation/s414_nose_resolution_remote.py'
    commands = {
        'prepare': {'argv': [python, 'validation/prepare_s414_nose_resolution_payload.py', '--result', '{result}'],
                    'timeout_seconds': 300},
        'submit': {'argv': [python, script, 'submit', '--job-id', '{job_id}', '--workers', '{workers}',
                           '--payload', '{snapshot}/'+launch+'/s414_nose_resolution_payload.tar.gz',
                           '--result', '{result}'], 'timeout_seconds': 180},
        'poll': {'argv': [python, script, 'poll', '--job-id', '{job_id}', '--result', '{result}'],
                 'timeout_seconds': 60},
        'collect': {'argv': [python, script, 'collect', '--job-id', '{job_id}', '--result', '{result}'],
                    'timeout_seconds': 600},
        'audit': {'argv': [python, 'validation/report_s414_nose_resolution.py', '--audit',
                          '--study', '{snapshot}/'+study, '--launch', '{snapshot}/'+launch,
                          '--out', '{snapshot}/'+analysis, '--result', '{result}'],
                  'timeout_seconds': 1800}}
    descriptions = {
        'execution_complete': 'The bounded six-case supervisor terminates without execution/budget failure; scientific solve failures remain explicit.',
        'evidence_integrity': 'Input/transfer hashes, all converged saved states, residuals <=1e-7, forces and exported fields pass fresh checks.',
        'archived_controls': 'Both factor-1 cases exactly reproduce the archived corrected-AF forces and 18 state arrays before refinement.',
        'selected_coverage': 'All six selected states are available and converged; no missing case is omitted from the denominator.',
        'numerical_screen': 'Both fixed-incidence cases meet every protocol tolerance, contraction, topology and transition-cell screen.',
        'isolated_nose_refinement': 'Frozen polygon, metric and non-nose controls are preserved; wake-coordinate movement <=1e-5c supports nose-only attribution.'}
    return dict(schema=1, id='s414-nose-resolution-20261005-r1', project_root=str(ROOT),
        objective='Determine whether early fore-lower transition and excess drag survive local coupled refinement on the same polygon and station metric.',
        authorization='User approved the bounded six-case pilot and independent prelaunch/postanalysis review. No physical changes, extra cases, retries or outbound data request.',
        resources=dict(remote_workers=2, max_final_cases=6, max_requested_evaluations=30,
            max_remote_wall_seconds=5400, max_case_seconds=1800, attempts_per_case=1,
            thread_limits=1, monitor_max_hours=3),
        reviewer=dict(effort='high', timeout_seconds=900),
        files=sorted(files), writable_snapshot_paths=[study, launch, analysis], commands=commands,
        criteria=[dict(id=key, description=value) for key, value in descriptions.items()],
        protocol=runner.PROTOCOL,
        scientific_decision='A failed numerical criterion is a valid result, not a command failure or retry trigger. Improved fixed-angle drag is not matched-lift accuracy. No correction follows without independent review.',
        stages=['freeze', 'prepare', 'prelaunch_review', 'submit', 'poll', 'collect', 'audit', 'postanalysis_review'],
        stop_conditions=['Changed source/input or transfer hashes', 'Rejected independent prelaunch review',
            'Unknown remote identity or lost submission reply: reconcile only', 'Budget/worker limit exceeded',
            'Fresh state audit failure', 'Protocol, closure, Ncrit or geometry change required'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'orchestration/s414_nose_resolution.yaml')
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Preserve the existing plan; choose a new revision filename')
    plan = build()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(plan, indent=2)+'\n')
    print(args.out)


if __name__ == '__main__':
    main()
