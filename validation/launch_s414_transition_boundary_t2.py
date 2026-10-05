#!/usr/bin/env python3
"""Deploy the verified S414 payload once and detach its guarded runner."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

ROOT=Path(__file__).resolve().parents[1]
LAUNCH=ROOT/'validation/results/s414_transition_boundary_t2_launch_20261004'
REMOTE='/home/jony/multifoil-runs/s414-transition-boundary-t2-20261004'
HOST='jony4-anywhere'
OPTIONS=['-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15']
DEPLOY=r'''
import datetime,hashlib,json,os,subprocess,sys,tarfile
from pathlib import Path
root,archive=Path(sys.argv[1]),Path(sys.argv[2]);expected=sys.argv[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if sha(archive)!=expected:raise ValueError('Transferred archive hash mismatch')
root.mkdir(exist_ok=False)
with tarfile.open(archive,'r:gz') as tar:
    for member in tar.getmembers():
        if not (root/member.name).resolve().is_relative_to(root) or not (member.isfile() or member.isdir()):
            raise ValueError('Unsafe archive member')
    tar.extractall(root)
manifest=json.loads((root/'payload_manifest.json').read_text())
for name,digest in manifest.items():
    if sha(root/name)!=digest:raise ValueError('Payload file changed: '+name)
env=Path('/home/jony/multifoil-runs/s702-repeatability-20260916/.venv')
if not (env/'bin/python').exists():raise ValueError('Pinned environment missing')
(root/'.venv').symlink_to(env,target_is_directory=True)
with (root/'launch.log').open('ab',buffering=0) as log:
    process=subprocess.Popen(['bash','run_s414_transition_boundary_t2.sh'],cwd=root,
        stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
        start_new_session=True,close_fds=True)
receipt=dict(launched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    host=os.uname().nodename,root=str(root),pid=process.pid,archive_sha256=expected,
    study_sha256=sha(root/'validation/results/s414_transition_boundary_t2_20261004/study.json'),
    payload_files_verified=len(manifest),environment=str(env),workers_limit=2,
    command='bash run_s414_transition_boundary_t2.sh',detached=True)
(root/'launch_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
'''


def main():
    receipt=LAUNCH/'launch_receipt.json'
    if receipt.exists():raise ValueError('Launch already recorded; use monitor/resume instructions')
    preparation=json.loads((LAUNCH/'preparation.json').read_text())
    archive=LAUNCH/'s414_transition_boundary_t2_payload.tar.gz'
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=preparation['archive_sha256']:raise ValueError('Local payload changed')
    remote_archive=REMOTE+'.tar.gz'
    subprocess.run(['ssh',*OPTIONS,HOST,'test ! -e '+shlex.quote(REMOTE)+' && test ! -e '+shlex.quote(remote_archive)],check=True)
    subprocess.run(['scp',*OPTIONS,str(archive),HOST+':'+remote_archive],check=True)
    command='python3 - '+' '.join(map(shlex.quote,[REMOTE,remote_archive,digest]))
    response=subprocess.run(['ssh',*OPTIONS,HOST,command],input=DEPLOY,text=True,capture_output=True,check=True)
    data=json.loads(response.stdout)
    if data['study_sha256']!=preparation['study_sha256']:raise ValueError('Launched study differs')
    receipt.write_text(json.dumps(data,indent=2)+'\n');print(receipt.read_text())


if __name__=='__main__':main()
