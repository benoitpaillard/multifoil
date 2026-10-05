#!/usr/bin/env python3
"""Package and test the guarded four-case S414 fore-reconstruction run in isolation."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import run_s414_fore_reconstruction as run

LAUNCH=run.ROOT/'validation/results/s414_fore_reconstruction_launch_20261004'
TESTS=['test_s414_fore_reconstruction', 'test_amplification_jacobian', 'test_fixed_lift_policy', 'test_optimization', 'test_paneling', 'test_viscous.ViscousPhysics.test_boundary_layer_residual_and_jacobian_match_upstream', 'test_viscous.ViscousPhysics.test_failure_is_not_a_force_or_a_stall_prediction']
SHELL='#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\nexec 9>launch.lock\nflock -n 9\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1\n.venv/bin/python -m unittest -v test_s414_fore_reconstruction test_amplification_jacobian test_fixed_lift_policy test_optimization test_paneling test_viscous.ViscousPhysics.test_boundary_layer_residual_and_jacobian_match_upstream test_viscous.ViscousPhysics.test_failure_is_not_a_force_or_a_stall_prediction > preflight.log 2>&1\nexec .venv/bin/python validation/run_s414_fore_reconstruction.py --workers 2\n'


def main():
    LAUNCH.mkdir(exist_ok=True)
    if (LAUNCH/'launch_receipt.json').exists():raise ValueError('Preserve already launched payload')
    spec=run.specification();run.prepare(run.ROOT/run.STUDY,spec)
    names=set(spec['sha256'])|{'vendor/__init__.py','LICENSE','THIRD_PARTY.md',
        'validation/tn808_geometry.py','validation/data/nlr7301/NLR7301MainECARP.dat',
        'validation/data/nlr7301/NLR7301FlapECARP.dat','validation/data/nlr7301/LICENSE',
        'validation/monitor_s414_fore_reconstruction.py','validation/retrieve_s414_fore_reconstruction.py',
        'validation/launch_s414_fore_reconstruction.py','validation/prepare_s414_fore_reconstruction_payload.py',run.STUDY+'/study.json'}
    files={name:run.ROOT/name for name in names}
    shell=LAUNCH/'run_s414_fore_reconstruction.sh';shell.write_text(SHELL);files[shell.name]=shell
    manifest={name:run.sha(p) for name,p in sorted(files.items())};run.write(LAUNCH/'payload_manifest.json',manifest)
    with tempfile.TemporaryDirectory(prefix='s414_fore_reconstruction_payload_') as tmp:
        root=Path(tmp)
        for name,path in files.items():
            dest=root/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        result=subprocess.run([sys.executable,'validation/run_s414_fore_reconstruction.py','--prepare-only'],cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_check.log').write_text(result.stdout+result.stderr)
        if result.returncode:raise RuntimeError('Isolated preparation failed')
        expected=run.sha(run.ROOT/run.STUDY/'study.json')
        if run.sha(root/run.STUDY/'study.json')!=expected:raise ValueError('Isolated study manifest differs')
        tests=subprocess.run([sys.executable,'-m','unittest','-v',*TESTS],cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_tests.log').write_text(tests.stdout+tests.stderr)
        if tests.returncode:raise RuntimeError('Isolated payload tests failed')
    archive=LAUNCH/'s414_fore_reconstruction_payload.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for name,path in sorted(files.items()):tar.add(path,arcname=name)
        tar.add(LAUNCH/'payload_manifest.json',arcname='payload_manifest.json')
    receipt=dict(files=len(files),archive_bytes=archive.stat().st_size,archive_sha256=run.sha(archive),
        study_sha256=expected,isolated_payload_preparation_pass=True,isolated_payload_tests_pass=True,
        test_modules=TESTS,workers_limit=2)
    run.write(LAUNCH/'preparation.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
