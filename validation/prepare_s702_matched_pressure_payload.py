#!/usr/bin/env python3
"""Package all inputs and verify preparation in an isolated checkout."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import run_s702_matched_pressure as run

LAUNCH=run.ROOT/'validation/results/s702_matched_pressure_launch_20261002'
SHELL='''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec 9>launch.lock
flock -n 9
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -m unittest -v test_s702_matched_pressure test_fixed_lift_policy test_optimization test_paneling > preflight.log 2>&1
exec .venv/bin/python validation/run_s702_matched_pressure.py --workers 4
'''


def main():
    LAUNCH.mkdir(exist_ok=True)
    if (LAUNCH/'launch_receipt.json').exists():
        raise ValueError('Preserve launched payload; receipt already exists')
    spec=run.specification();run.prepare(run.ROOT/run.STUDY,spec)
    names=list(spec['sha256'])+['vendor/__init__.py','LICENSE','THIRD_PARTY.md',
        'validation/tn808_geometry.py','validation/data/nlr7301/NLR7301MainECARP.dat',
        'validation/data/nlr7301/NLR7301FlapECARP.dat','validation/data/nlr7301/LICENSE',
        'validation/monitor_s702_matched_pressure.py',
        'validation/prepare_s702_matched_pressure_payload.py',
        run.QUAL+'/FINDINGS.md',run.STUDY+'/study.json']
    files={name:run.ROOT/name for name in names}
    shell=LAUNCH/'run_matched_pressure.sh';shell.write_text(SHELL)
    files['run_matched_pressure.sh']=shell
    manifest={name:run.sha(path) for name,path in sorted(files.items())}
    run.write(LAUNCH/'payload_manifest.json',manifest)
    with tempfile.TemporaryDirectory(prefix='s702_matched_pressure_payload_') as tmp:
        root=Path(tmp)
        for name,path in files.items():
            dest=root/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        check=subprocess.run([sys.executable,'validation/run_s702_matched_pressure.py','--prepare-only'],
            cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_check.log').write_text(check.stdout+check.stderr)
        if check.returncode:raise RuntimeError('Isolated payload preparation failed')
        expected=run.sha(run.ROOT/run.STUDY/'study.json')
        if run.sha(root/run.STUDY/'study.json')!=expected:raise ValueError('Isolated manifest differs')
        tests=subprocess.run([sys.executable,'-m','unittest','-v','test_s702_matched_pressure',
            'test_fixed_lift_policy','test_optimization','test_paneling'],
            cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_tests.log').write_text(tests.stdout+tests.stderr)
        if tests.returncode:raise RuntimeError('Isolated payload tests failed')
    archive=LAUNCH/'s702_matched_pressure_payload.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for name,path in sorted(files.items()):tar.add(path,arcname=name)
        tar.add(LAUNCH/'payload_manifest.json',arcname='payload_manifest.json')
    receipt=dict(files=len(files),archive_bytes=archive.stat().st_size,archive_sha256=run.sha(archive),
        study_sha256=expected,isolated_payload_preparation_pass=True,isolated_payload_tests_pass=True)
    run.write(LAUNCH/'preparation.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
