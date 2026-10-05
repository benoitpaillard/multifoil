#!/usr/bin/env python3
"""Resume the frozen S702 study after verifying its exact known baseline failure.

This does not waive arbitrary test failures. Both completed test logs, the
original baseline receipt, all payload hashes and pinned runtime must match.
The seven unaffected viscous checks are reused from that verified preflight;
all 38 fast checks run again. The 100-chord order test remains failed/unmodified.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT/'handover/policy_runtime'
STUDY = 'validation/results/s702_initialization_policy_20260921'
FAILURE = 'test_multiple_elements_order_and_far_limit (test_viscous.ViscousPhysics.test_multiple_elements_order_and_far_limit)'
ASSERTION = 'AssertionError: 1.1415574480699853e-06 not less than 5e-07'
FAST = ['test_fixed_lift_policy', 'test_s702_initialization_policy', 'test_s702_repeatability',
        'test_optimization', 'test_s702_sweep', 'test_s702', 'test_s702_pressure']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, indent=2)+'\n')
    temp.replace(path)


def check_failure(log, count):
    if (re.findall(r'^FAIL: (.+)$', log, re.M) != [FAILURE] or
            re.search(r'^ERROR:', log, re.M) or ASSERTION not in log or
            not re.search(rf'^Ran {count} tests? in ', log, re.M) or
            'FAILED (failures=1)' not in log):
        raise ValueError('Unexpected baseline/preflight outcome; execution is stopped')


def verify_evidence():
    continuation = json.loads((ROOT/'handover/policy_continuation_payload_20260921.json').read_text())
    if sha(ROOT/'handover/policy_payload_20260921.json') != continuation['original_payload_sha256']:
        raise ValueError('Original payload manifest changed')
    for name, expected in continuation['sha256'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Continuation tooling changed: '+name)
    payload = json.loads((ROOT/'handover/policy_payload_20260921.json').read_text())
    for name, expected in payload['sha256'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Frozen payload changed: '+name)
    receipt_path = RUNTIME/'original_far_field.json'
    receipt = json.loads(receipt_path.read_text())
    baseline = Path(receipt['cwd'])
    if receipt['exit_code'] != 1:
        raise ValueError('Baseline did not reproduce the recorded failure')
    for name, expected in receipt['sha256'].items():
        if sha(ROOT/name) != expected or sha(baseline/name) != expected:
            raise ValueError('Current and original baseline differ: '+name)
    # Freeze the original stopped preflight before the append-only launch log
    # receives this continuation's output. On resume, verify the frozen evidence.
    original_log = RUNTIME/'preflight_attempt_001.log'
    raw = original_log.read_text() if original_log.exists() else (RUNTIME/'launch.log').read_text()
    check_failure(raw, 46)
    baseline_log = RUNTIME/'original_far_field.log'
    check_failure(baseline_log.read_text(), 1)
    if not original_log.exists():
        original_log.write_text(raw)
    if sys.version_info[:2] != (3, 12):
        raise ValueError('Use the pinned Python 3.12 environment')
    expected_versions = {'numpy': '2.3.5', 'scipy': '1.17.0', 'matplotlib': '3.10.8'}
    actual_versions = {n: importlib.metadata.version(n) for n in expected_versions}
    if actual_versions != expected_versions:
        raise ValueError('Pinned dependencies changed')
    if (ROOT/'.venv').resolve() != (baseline/'.venv').resolve():
        raise ValueError('The preflight and baseline must use the same pinned environment')
    decision = dict(classification='Reproduced pre-existing far-field element-order failure',
        decision='Proceed only with the declared S702 initialization experiment; retain failed test',
        failed_test=FAILURE, assertion=ASSERTION, preflight_passed=45, preflight_total=46,
        unchanged_core_and_test_sha256=receipt['sha256'],
        evidence_sha256={p.name: sha(p) for p in (original_log, baseline_log, receipt_path)},
        python=platform.python_version(), versions=actual_versions,
        scope='No claim that the full physics suite passes; full-state S702 review remains required')
    frozen = RUNTIME/'baseline_classification.json'
    if frozen.exists() and json.loads(frozen.read_text()) != decision:
        raise ValueError('Previously classified baseline evidence changed')
    write(frozen, decision)
    return decision


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=7, choices=range(1, 8))
    args = ap.parse_args()
    RUNTIME.mkdir(parents=True, exist_ok=True)
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ.update(MPLBACKEND='Agg', PYTHONUNBUFFERED='1')
    with (RUNTIME/'launch.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        started = stamp()
        def phase(name, code=None, **extra):
            data = dict(phase=name, exit_code=code, launcher_pid=os.getpid(),
                        hostname=platform.node(), started_utc=started, updated_utc=stamp(), **extra)
            write(RUNTIME/'launch_status.json', data)
            print(json.dumps(data), flush=True)
        try:
            phase('verifying_baseline_evidence')
            verify_evidence()
            phase('fast_preflight', known_baseline_failure=FAILURE)
            subprocess.run([sys.executable, '-m', 'unittest', '-v', *FAST], cwd=ROOT, check=True)
            subprocess.run([sys.executable, 'validation/run_s702_initialization_policy.py',
                            '--out', STUDY, '--prepare-only'], cwd=ROOT, check=True)
            phase('running_study', known_baseline_failure=FAILURE)
            result = subprocess.run([sys.executable, 'validation/run_s702_initialization_policy.py',
                                     '--out', STUDY, '--workers', str(args.workers)], cwd=ROOT)
            name = {0: 'finished', 2: 'finished_with_regression_failures'}.get(result.returncode, 'execution_failed')
            phase(name, result.returncode, known_baseline_failure=FAILURE,
                  followup='Review selected full states before adopting the policy in optimization')
            return result.returncode
        except Exception as exc:
            phase('execution_failed', 1, error=f'{type(exc).__name__}: {exc}')
            raise


if __name__ == '__main__':
    raise SystemExit(main())
