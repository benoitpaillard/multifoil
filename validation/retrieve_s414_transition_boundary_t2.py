#!/usr/bin/env python3
"""Retrieve a completed S414 study and verify every file against remote SHA-256."""
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
STUDY='validation/results/s414_transition_boundary_t2_20261004'
LAUNCH=ROOT/'validation/results/s414_transition_boundary_t2_launch_20261004'
REMOTE='/home/jony/multifoil-runs/s414-transition-boundary-t2-20261004'
OPTIONS=['-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15']
PROBE=r'''
import hashlib,json
from pathlib import Path
root=Path('/home/jony/multifoil-runs/s414-transition-boundary-t2-20261004')
out=root/'validation/results/s414_transition_boundary_t2_20261004'
status=json.loads((out/'status.json').read_text())
if status['state']!='finished':raise ValueError('Study has not finished')
files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted(out.rglob('*')) if p.is_file()}
print(json.dumps(dict(status=status,files=files,launch=json.loads((root/'launch_receipt.json').read_text()))))
'''


def main():
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    remote=subprocess.run(['ssh',*OPTIONS,'jony4-anywhere','python3 -'],input=PROBE,
        text=True,capture_output=True,check=True,timeout=30)
    evidence=json.loads(remote.stdout);files=evidence['files']
    launch=json.loads((LAUNCH/'launch_receipt.json').read_text())
    if launch!=evidence['launch'] or files['study.json']!=launch['study_sha256']:
        raise ValueError('Remote launch/study identity changed')
    for name,digest in files.items():
        path=ROOT/STUDY/name
        if path.exists() and sha(path)!=digest:
            raise ValueError('Preserve differing local file before retrieval: '+name)
    manifest=LAUNCH/'retrieved_files_sha256.json'
    if manifest.exists() and json.loads(manifest.read_text())!=files:
        raise ValueError('Previously retrieved remote result changed')
    manifest.write_text(json.dumps(files,indent=2)+'\n')
    subprocess.run(['rsync','-a','--partial','-e','ssh '+' '.join(OPTIONS),
        'jony4-anywhere:'+REMOTE+'/'+STUDY+'/',str(ROOT/STUDY)+'/'],check=True)
    for name,digest in files.items():
        if sha(ROOT/STUDY/name)!=digest:raise ValueError('Retrieved file differs: '+name)
    receipt=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),files_verified=len(files),
        study_sha256=files['study.json'],remote_status=evidence['status'],all_files_match=True)
    (LAUNCH/'retrieval.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
