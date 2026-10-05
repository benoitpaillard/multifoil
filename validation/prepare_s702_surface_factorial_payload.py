#!/usr/bin/env python3
"""Build and prepare-test a complete isolated payload for the five-target study."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import run_s702_surface_factorial as run

LAUNCH = run.ROOT/'validation/results/s702_surface_factorial_launch_20261001'
SHELL = '''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec 9>launch.lock
flock -n 9
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -m unittest -v test_s702_surface_factorial test_diagnostic_discretization test_s702_discretization test_fixed_lift_policy test_optimization test_paneling > preflight.log 2>&1
exec .venv/bin/python validation/run_s702_surface_factorial.py --workers 4
'''


def main():
    LAUNCH.mkdir(exist_ok=True)
    spec=run.specification(); run.prepare(run.ROOT/run.STUDY,spec)
    if (LAUNCH/'launch_receipt.json').exists():
        raise ValueError('Launch receipt exists; preserve delivered archive instead of rebuilding it')
    old=run.read(run.ROOT/'validation/results/s702_discretization_launch_20261001/payload_manifest.json')
    files={}
    for name,digest in old.items():
        if name=='run_discretization.sh':continue
        if run.sha(run.ROOT/name)!=digest:raise ValueError('Prior payload input changed: '+name)
        files[name]=run.ROOT/name
    for name in list(spec['sha256'])+['validation/monitor_s702_surface_factorial.py',
                                     'validation/prepare_s702_surface_factorial_payload.py',run.STUDY+'/study.json']:
        files[name]=run.ROOT/name
    shell=LAUNCH/'run_surface_factorial.sh';shell.write_text(SHELL)
    files['run_surface_factorial.sh']=shell
    manifest={name:run.sha(path) for name,path in sorted(files.items())}
    run.write(LAUNCH/'payload_manifest.json',manifest)
    with tempfile.TemporaryDirectory(prefix='s702_surface_payload_') as tmp:
        root=Path(tmp)
        for name,path in files.items():
            target=root/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        check=subprocess.run([sys.executable,'validation/run_s702_surface_factorial.py','--prepare-only'],
                             cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_check.log').write_text(check.stdout+check.stderr)
        if check.returncode:raise RuntimeError('Isolated preparation failed; inspect isolated_payload_check.log')
        expected=run.sha(run.ROOT/run.STUDY/'study.json')
        if run.sha(root/run.STUDY/'study.json')!=expected:raise ValueError('Isolated manifest differs')
    archive=LAUNCH/'s702_surface_factorial_payload.tar.gz'
    with tarfile.open(archive,'w:gz') as stream:
        for name,path in sorted(files.items()):stream.add(path,arcname=name)
        stream.add(LAUNCH/'payload_manifest.json',arcname='payload_manifest.json')
    info=dict(files=len(files),archive_bytes=archive.stat().st_size,archive_sha256=run.sha(archive),
              study_sha256=expected,isolated_payload_preparation_pass=True)
    run.write(LAUNCH/'preparation.json',info);print(json.dumps(info,indent=2))


if __name__=='__main__':main()
