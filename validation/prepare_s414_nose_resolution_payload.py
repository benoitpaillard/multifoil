#!/usr/bin/env python3
"""Freeze a self-contained numerical payload; short contract checks only, no flow solves."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import run_s414_nose_resolution as run

TESTS = ('test_s414_nose_resolution', 'test_s414_nose_resolution_runner', 'test_s414_nose_resolution_analysis')
SHELL = ('#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\n'
    'exec 9>launch.lock\nflock -n 9\n'
    'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1\n'
    'export PYTHONDONTWRITEBYTECODE=1\n'
    'exec .venv/bin/python validation/run_s414_nose_resolution_launch.py\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path)
    args = parser.parse_args(); launch = run.ROOT/run.LAUNCH; launch.mkdir(parents=True, exist_ok=True)
    if (launch/'launch_receipt.json').exists():
        raise ValueError('Preserve launched payload')
    spec = run.specification(); run.prepare(run.ROOT/run.STUDY, spec)
    archive = launch/'s414_nose_resolution_payload.tar.gz'
    if (launch/'preparation.json').exists():
        saved = run.read(launch/'preparation.json')
        if saved['study_sha256'] != run.sha(run.ROOT/run.STUDY/'study.json') or saved['archive_sha256'] != run.sha(archive):
            raise ValueError('Prepared payload changed; use a new immutable revision')
        if args.result: run.write(args.result, saved)
        print(archive); return 0
    # The coordinator remains local. Its fake-command recovery tests run in the
    # frozen ROOT snapshot, without adding orchestration code to the remote payload.
    check = subprocess.run([sys.executable, '-m', 'unittest', '-v', 'orchestration/test_pilot.py'],
        cwd=run.ROOT, text=True, capture_output=True, timeout=120)
    orchestration_log = launch/'orchestration_tests.log'
    orchestration_log.write_text(check.stdout+check.stderr)
    if check.returncode: raise RuntimeError('Frozen orchestration recovery tests failed')
    orchestration_receipt = launch/'orchestration_tests.json'
    run.write(orchestration_receipt, dict(passed=True, exit_code=check.returncode,
        command=[sys.executable, '-m', 'unittest', '-v', 'orchestration/test_pilot.py'],
        cwd=str(run.ROOT), log_sha256=run.sha(orchestration_log), nonlinear_solves=0,
        source_sha256={n: run.sha(run.ROOT/n) for n in ('orchestration/pilot.py', 'orchestration/test_pilot.py')}))
    names = set(spec['sha256']) | {'LICENSE', 'THIRD_PARTY.md', run.STUDY+'/study.json'}
    files = {name: run.ROOT/name for name in names}
    shell = launch/'run_s414_nose_resolution.sh'; shell.write_text(SHELL); files[shell.name] = shell
    manifest = {name: run.sha(path) for name, path in sorted(files.items())}
    run.write(launch/'payload_manifest.json', manifest)
    with tempfile.TemporaryDirectory(prefix='s414_nose_resolution_payload_') as tmp:
        root = Path(tmp)
        for name, path in files.items():
            destination = root/name; destination.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, destination)
        check = subprocess.run([sys.executable, 'validation/run_s414_nose_resolution.py', '--prepare-only'],
            cwd=root, text=True, capture_output=True, timeout=60)
        (launch/'isolated_payload_check.log').write_text(check.stdout+check.stderr)
        if check.returncode or run.sha(root/run.STUDY/'study.json') != run.sha(run.ROOT/run.STUDY/'study.json'):
            raise RuntimeError('Isolated payload preparation failed')
        check = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TESTS], cwd=root,
            text=True, capture_output=True, timeout=120)
        (launch/'isolated_payload_tests.log').write_text(check.stdout+check.stderr)
        if check.returncode:
            raise RuntimeError('Isolated payload contract/geometry checks failed')
        check = subprocess.run([sys.executable, 'validation/check_s414_nose_resolution.py', '--output', 'linear_preflight.json'],
            cwd=root, text=True, capture_output=True, timeout=240)
        (launch/'isolated_linear_preflight.log').write_text(check.stdout+check.stderr)
        if check.returncode:
            raise RuntimeError('Isolated linear preflight failed')
        shutil.copyfile(root/'linear_preflight.json', launch/'linear_preflight.json')
        shutil.copyfile(root/'linear_preflight.npz', launch/'linear_preflight.npz')
    with tarfile.open(archive, 'w:gz') as tar:
        for name, path in sorted(files.items()): tar.add(path, arcname=name)
        tar.add(launch/'payload_manifest.json', arcname='payload_manifest.json')
    receipt = dict(status='prepared', files=len(files), archive=str(archive), archive_sha256=run.sha(archive),
        archive_bytes=archive.stat().st_size, study_sha256=run.sha(run.ROOT/run.STUDY/'study.json'),
        isolated_preparation_pass=True, isolated_short_tests_pass=True, isolated_linear_preflight_pass=True,
        orchestration_recovery_tests_pass=True, launch_budget_scope='complete_remote_launch_including_preflight',
        nonlinear_solves_run=0, test_modules=list(TESTS),
        workers_limit=2, evidence=[dict(path=str(p), sha256=run.sha(p)) for p in
            (launch/'payload_manifest.json', archive, run.ROOT/run.STUDY/'study.json',
             launch/'linear_preflight.json', launch/'linear_preflight.npz', orchestration_log,
             orchestration_receipt, shell, launch/'isolated_payload_tests.log', launch/'isolated_linear_preflight.log')])
    run.write(launch/'preparation.json', receipt)
    if args.result: run.write(args.result, receipt)
    print(archive)
    return 0


if __name__ == '__main__': raise SystemExit(main())
