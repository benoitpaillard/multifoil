#!/usr/bin/env python3
"""Package and test the guarded S414 transition-boundary run (V4 convex-nose viscous check) in isolation."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import run_s414_nose_viscous as run

LAUNCH = run.ROOT/'validation/results/s414_nose_viscous_launch_20261004'
TESTS = ['test_s414_nose_viscous', 'test_s414_nose_geometry', 'test_fixed_lift_policy', 'test_default_integration.DefaultAndDispatch',
    'test_default_integration.GrowthDerivativeBranches',
    'test_default_integration.PropagationAndProvenance.test_cli_config_and_failed_output_record_mode',
    'test_default_integration.PropagationAndProvenance.test_cli_rejects_bad_or_misplaced_mode_without_output',
    'test_default_integration.PropagationAndProvenance.test_optimizer_flow_job_passes_configured_mode',
    'test_default_integration.PropagationAndProvenance.test_optimizer_manifest_hashes_correction_and_guards_resume',
    'test_amplification_jacobian', 'test_paneling',
    'test_viscous.ViscousPhysics.test_boundary_layer_residual_and_jacobian_match_upstream',
    'test_viscous.ViscousPhysics.test_failure_is_not_a_force_or_a_stall_prediction']
SHELL = ('#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\nexec 9>launch.lock\nflock -n 9\n'
    'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1\n'
    '.venv/bin/python -m unittest -v '+' '.join(TESTS)+' > preflight.log 2>&1\n'
    'exec .venv/bin/python validation/run_s414_nose_viscous.py --workers 2\n')
SMOKE = r'''
import json, sys
sys.path.insert(0, 'validation')
import run_s414_nose_viscous as run
spec = run.specification(); out = run.ROOT/run.STUDY; job = spec['jobs'][0]
result = run.worker(out, job, spec); gate = run.compare_control(out/job['name'], job, spec['acceptance_control_CD_counts'])
print(json.dumps(dict(job=job['name'], feasible=result['feasible'], alpha_deg=result['alpha_deg'], CL=result['CL'], CD=result['CD'],
    selected_path=result.get('selected_path'), evaluations=len(result['records']), control_passed=gate['passed'],
    CD_difference_counts=gate['CD_difference_counts'], arrays_equal=sum(a['exactly_equal'] for a in gate['array_comparisons']),
    arrays_compared=len(gate['array_comparisons']), wall_seconds=result['wall_seconds'])))
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true', help='Also run the real nominal control in the isolated copy (local runtime)')
    args = parser.parse_args()
    LAUNCH.mkdir(exist_ok=True)
    if (LAUNCH/'launch_receipt.json').exists(): raise ValueError('Preserve already launched payload')
    spec = run.specification(); run.prepare(run.ROOT/run.STUDY, spec)
    names = set(spec['sha256']) | {'LICENSE', 'THIRD_PARTY.md', run.STUDY+'/study.json'}
    files = {name: run.ROOT/name for name in names}
    shell = LAUNCH/'run_s414_nose_viscous.sh'; shell.write_text(SHELL); files[shell.name] = shell
    manifest = {name: run.sha(p) for name, p in sorted(files.items())}; run.write(LAUNCH/'payload_manifest.json', manifest)
    smoke = None
    with tempfile.TemporaryDirectory(prefix='s414_nose_viscous_payload_') as tmp:
        root = Path(tmp)
        for name, path in files.items():
            dest = root/name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, dest)
        result = subprocess.run([sys.executable, 'validation/run_s414_nose_viscous.py', '--prepare-only'], cwd=root, text=True, capture_output=True)
        (LAUNCH/'isolated_payload_check.log').write_text(result.stdout+result.stderr)
        if result.returncode: raise RuntimeError('Isolated preparation failed')
        expected = run.sha(run.ROOT/run.STUDY/'study.json')
        if run.sha(root/run.STUDY/'study.json') != expected: raise ValueError('Isolated study manifest differs')
        tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TESTS], cwd=root, text=True, capture_output=True)
        (LAUNCH/'isolated_payload_tests.log').write_text(tests.stdout+tests.stderr)
        if tests.returncode: raise RuntimeError('Isolated payload tests failed')
        if args.smoke:
            done = subprocess.run([sys.executable, '-c', SMOKE], cwd=root, text=True, capture_output=True)
            (LAUNCH/'isolated_payload_smoke.log').write_text(done.stdout+done.stderr)
            if done.returncode: raise RuntimeError('Isolated smoke job failed')
            smoke = json.loads(done.stdout.strip().splitlines()[-1])
    archive = LAUNCH/'s414_nose_viscous_payload.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for name, path in sorted(files.items()): tar.add(path, arcname=name)
        tar.add(LAUNCH/'payload_manifest.json', arcname='payload_manifest.json')
    receipt = dict(files=len(files), archive_bytes=archive.stat().st_size, archive_sha256=run.sha(archive),
        study_sha256=expected, isolated_payload_preparation_pass=True, isolated_payload_tests_pass=True,
        test_modules=TESTS, workers_limit=2, isolated_smoke_job=smoke,
        smoke_note='Local-runtime exercise of the real fixed-lift worker on the nominal control; not the pinned-runtime gate')
    run.write(LAUNCH/'preparation.json', receipt); print(json.dumps(receipt, indent=2))


if __name__ == '__main__': main()
