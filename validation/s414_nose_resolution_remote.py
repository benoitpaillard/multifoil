#!/usr/bin/env python3
"""Idempotent remote submission, read-only polling and verified retrieval for the S414 pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STUDY = 'validation/results/s414_nose_resolution_20261005'
LAUNCH = ROOT/'validation/results/s414_nose_resolution_launch_20261005'
HOST = 'jony4-anywhere'
REMOTE_BASE = '/home/jony/multifoil-runs'
OPTIONS = ['-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes', '-o', 'ConnectTimeout=15']
ENVIRONMENT = '/home/jony/multifoil-runs/s702-repeatability-20260916/.venv'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+'.tmp'); temp.write_text(json.dumps(data, indent=2)+'\n'); temp.replace(path)


COMMON = r'''
import datetime,fcntl,hashlib,json,os,subprocess,sys,tarfile,time
from pathlib import Path
root=Path(sys.argv[1]);job_id=sys.argv[2]
study='validation/results/s414_nose_resolution_20261005'
out=root/study
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text()) if p.exists() else {}
def state():
    receipt=read(root/'launch_receipt.json'); detail=read(out/'status.json');alive=False
    for pid in set((detail.get('pid'),receipt.get('pid'))):
        try:
            proc=Path('/proc')/str(pid)
            cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
            alive=alive or ((proc/'cwd').resolve()==root and 'run_s414_nose_resolution' in cmd and (proc/'stat').read_text().split()[2]!='Z')
        except OSError:pass
    supervisor_alive=alive;jobs=[]
    for j in read(out/'study.json').get('jobs',[]):
        r=read(out/j['name']/'result.json');p=r or read(out/j['name']/'progress.json')
        attempt=read(out/j['name']/'attempt.json');worker_alive=False
        try:
            proc=Path('/proc')/str(attempt.get('pid'))
            cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
            worker_alive=(proc/'cwd').resolve()==root and '--worker '+j['name'] in cmd and (proc/'stat').read_text().split()[2]!='Z'
        except OSError:pass
        alive=alive or worker_alive
        jobs.append(dict(job=j['name'],complete=bool(r),alive=worker_alive,feasible=r.get('feasible'),evaluations=len(p.get('records',[])),error=p.get('error','')))
    status='running' if alive else 'complete' if detail.get('state')=='finished' else 'failed' if receipt else 'unknown'
    preflight=root/'preflight.log'
    return dict(job_id=job_id,status=status,alive=alive,supervisor_alive=supervisor_alive,remote_root=str(root),detail=detail,jobs=jobs,
        controls=read(out/'controls.json'),summary=read(out/'summary.json'),launch=receipt,
        preflight_tail=preflight.read_text().splitlines()[-8:] if preflight.exists() else [],evidence=[])
'''

SUBMIT = COMMON+r'''
archive=Path(sys.argv[3]);expected=sys.argv[4];expected_study=sys.argv[5];env=Path(sys.argv[6])
if sha(archive)!=expected:raise ValueError('Transferred archive hash mismatch')
root.parent.mkdir(parents=True,exist_ok=True)
with (root.parent/(root.name+'.submit.lock')).open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    if (root/'launch_receipt.json').exists():
        old=read(root/'launch_receipt.json')
        if old['job_id']!=job_id or old['archive_sha256']!=expected or old['study_sha256']!=expected_study:
            raise ValueError('Job identity already used for a different immutable payload')
        print(json.dumps(state()));sys.exit(0)
    if root.exists():raise ValueError('Interrupted deployment retained; no automatic resubmission')
    root.mkdir()
    with tarfile.open(archive,'r:gz') as tar:
        for m in tar.getmembers():
            if not (root/m.name).resolve().is_relative_to(root) or not (m.isfile() or m.isdir()):
                raise ValueError('Unsafe archive member')
        tar.extractall(root)
    manifest=read(root/'payload_manifest.json')
    for name,digest in manifest.items():
        if sha(root/name)!=digest:raise ValueError('Payload file differs: '+name)
    if sha(out/'study.json')!=expected_study:raise ValueError('Study identity differs')
    if not (env/'bin/python').is_file():raise ValueError('Pinned environment missing')
    (root/'.venv').symlink_to(env,target_is_directory=True)
    started=time.time();started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
    spec=read(out/'study.json')
    budget=dict(started_utc=started_utc,started_unix=started,deadline_unix=started+spec['study_timeout_seconds'],
        limit_seconds=spec['study_timeout_seconds'],study_sha256=expected_study,
        scope='complete_remote_launch_including_preflight',cleanup_reserve_seconds=10)
    temp=out/'budget.json.tmp';temp.write_text(json.dumps(budget,indent=2)+'\n');temp.replace(out/'budget.json')
    with (root/'launch.log').open('ab',buffering=0) as log:
        process=subprocess.Popen(['bash','run_s414_nose_resolution.sh'],cwd=root,stdin=subprocess.DEVNULL,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    receipt=dict(job_id=job_id,launched_utc=started_utc,
        host=os.uname().nodename,root=str(root),pid=process.pid,archive_sha256=expected,
        study_sha256=expected_study,payload_files_verified=len(manifest),environment=str(env),workers_limit=2,
        command='bash run_s414_nose_resolution.sh',detached=True,launch_deadline_unix=budget['deadline_unix'],
        budget_scope=budget['scope'])
    tmp=root/'launch_receipt.json.tmp';tmp.write_text(json.dumps(receipt,indent=2)+'\n');tmp.replace(root/'launch_receipt.json')
    print(json.dumps(state()))
'''

POLL = COMMON+'\nprint(json.dumps(state()))\n'
COLLECT = COMMON+r'''
data=state()
if data['status'] not in ('complete','failed') or data['alive']:raise ValueError('Remote job is not quiescent')
if not data['launch']:raise ValueError('No launch identity')
data['files']={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
data['launch_files']={n:sha(root/n) for n in ('launch_receipt.json','launch.log','supervisor.log','preflight.log','linear_preflight.json','linear_preflight.npz') if (root/n).is_file()}
print(json.dumps(data))
'''


def remote(script, root, job_id, *values):
    command='python3 - '+' '.join(shlex.quote(str(v)) for v in (root, job_id, *values))
    response=subprocess.run(['ssh', *OPTIONS, HOST, command],input=script,text=True,capture_output=True,check=True,timeout=45)
    return json.loads(response.stdout)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('submit','poll','collect'))
    parser.add_argument('--job-id',default='s414-nose-resolution-20261005')
    parser.add_argument('--result',type=Path)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--payload',type=Path,default=LAUNCH/'s414_nose_resolution_payload.tar.gz')
    args=parser.parse_args()
    if not re.fullmatch(r'[a-z0-9][a-z0-9-]{0,80}',args.job_id) or args.workers!=2:
        parser.error('Require a safe deterministic job id and the frozen two-worker limit')
    root=REMOTE_BASE+'/'+args.job_id
    if args.operation=='submit':
        preparation=read(LAUNCH/'preparation.json');digest=sha(args.payload)
        if digest!=preparation['archive_sha256']:raise ValueError('Prepared payload changed')
        remote_archive=root+'.tar.gz'
        subprocess.run(['scp',*OPTIONS,str(args.payload),HOST+':'+remote_archive],check=True,timeout=90)
        data=remote(SUBMIT,root,args.job_id,remote_archive,digest,preparation['study_sha256'],ENVIRONMENT)
        write(LAUNCH/'launch_receipt.json',data['launch'])
    elif args.operation=='poll':
        data=remote(POLL,root,args.job_id)
        if data.get('launch'):
            preparation=read(LAUNCH/'preparation.json');receipt=data['launch']
            if (receipt['job_id']!=args.job_id or receipt['archive_sha256']!=preparation['archive_sha256']
                    or receipt['study_sha256']!=preparation['study_sha256']):
                raise ValueError('Remote launch does not match the prepared identity')
            local=LAUNCH/'launch_receipt.json'
            if local.exists() and read(local)!=receipt:
                raise ValueError('Previously recorded launch identity differs')
            write(local,receipt)
    else:
        data=remote(COLLECT,root,args.job_id);launch=read(LAUNCH/'launch_receipt.json')
        if data['launch']!=launch or data['files'].get('study.json')!=launch['study_sha256']:
            raise ValueError('Remote launch/study identity differs')
        for name,digest in data['files'].items():
            path=ROOT/STUDY/name
            if path.exists() and sha(path)!=digest:
                raise ValueError('Preserve different local result before retrieval: '+name)
        manifest=LAUNCH/'retrieved_files_sha256.json'
        if manifest.exists() and read(manifest)!=data['files']:
            raise ValueError('Previously retrieved results changed')
        subprocess.run(['rsync','-a','--partial','-e','ssh '+' '.join(OPTIONS),HOST+':'+root+'/'+STUDY+'/',str(ROOT/STUDY)+'/'],check=True)
        for name,digest in data['files'].items():
            if sha(ROOT/STUDY/name)!=digest:raise ValueError('Retrieved file differs: '+name)
        write(manifest,data['files'])
        for name,digest in data['launch_files'].items():
            destination=LAUNCH/('remote_'+name)
            subprocess.run(['scp',*OPTIONS,HOST+':'+root+'/'+name,str(destination)],check=True,timeout=60)
            if sha(destination)!=digest:raise ValueError('Retrieved launch evidence differs: '+name)
        data=dict(job_id=args.job_id,status='collected',remote_status=data['status'],files_verified=len(data['files']),
            all_files_match=True,study_sha256=launch['study_sha256'],
            launch_file_sha256=data['launch_files'], evidence=[dict(path=str(p),sha256=sha(p)) for p in
                [manifest, *(LAUNCH/('remote_'+name) for name in data['launch_files'])]])
        write(LAUNCH/'retrieval.json',data)
    if args.result:write(args.result,data)
    print(json.dumps(data,indent=2))
    return 0


if __name__=='__main__':raise SystemExit(main())
