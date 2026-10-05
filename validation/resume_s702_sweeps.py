#!/usr/bin/env python3
"""Resume interrupted studies using their exact archived execution versions.

Linux/fork is required to inherit the verified in-memory worker modules.
Completed targets are immutable; partial target folders are moved into an
interruption archive before restarting. All input hashes must match either the
current input or an exact source snapshot, which is actually executed.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import importlib.util
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import sys
import types
from datetime import datetime, timezone

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'validation'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('studies',type=Path,nargs='+')
    ap.add_argument('--workers',type=int,default=8)
    args = ap.parse_args()
    if not 1 <= args.workers <= 8:
        ap.error('Use one to eight processes')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    stages = []
    for index, raw in enumerate(args.studies):
        folder = raw.resolve()
        spec = json.loads((folder/'study.json').read_text())
        runner = 'validation/s702_ncrit_sweep.py'
        for relative, expected in spec['sha256'].items():
            source = ROOT/relative
            if digest(source) != expected:
                snapshot = folder/'source_snapshot'/relative
                # Only the orchestration runner may differ. Imported numerical
                # and geometry modules must still match the actual runtime.
                if relative != runner or not snapshot.exists() or digest(snapshot) != expected:
                    raise ValueError('Changed execution dependency: '+relative)
        source = folder/'source_snapshot'/runner
        if not source.exists():
            source = ROOT/runner
        if digest(source) != spec['sha256'][runner]:
            raise ValueError('Cannot recover the executed worker version')
        name = f's702_archived_worker_{index}'
        module = types.ModuleType(name)
        module.__file__ = str(ROOT/runner)
        sys.modules[name] = module
        exec(compile(source.read_bytes(),str(ROOT/runner),'exec'),module.__dict__)
        interruptions = []
        for job in spec['jobs']:
            group = folder/job['name']
            if (group/'result.json').exists():
                continue
            for cl in spec['targets']:
                target = group/f'cl{module.tag(cl)}'
                if target.exists() and not (target/'result.json').exists():
                    archive = folder/'interrupted_attempts'/stamp/job['name']/target.name
                    archive.parent.mkdir(parents=True,exist_ok=True)
                    shutil.move(str(target),str(archive))
                    interruptions.append(dict(job=job['name'],target_cl=cl,
                        archive=str(archive.relative_to(ROOT)),
                        reason='Execution service disconnected; no target result was finalized'))
        dump(folder/f'resume_{stamp}.json',dict(resumed_utc=stamp,
             worker_source_sha256=digest(source), workers=args.workers,
             interrupted_attempts=interruptions))
        jobs = [j for j in spec['jobs'] if not (folder/j['name']/'result.json').exists()]
        stages.append((folder,spec,module,jobs))
    # Round-robin stages keeps one slow transition case from blocking the other
    # studies, without oversubscribing the machine with multiple large pools.
    work = []
    for jobs in itertools.zip_longest(*(s[3] for s in stages)):
        for stage,job in zip(stages,jobs):
            if job is not None:
                work.append((stage,job))
    with ProcessPoolExecutor(args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
        futures = {pool.submit(stage[2].run_group,stage[0],job,stage[1]):(stage[0],job) for stage,job in work}
        for future in as_completed(futures):
            folder,job = futures[future]
            result = future.result()
            print(dict(study=folder.name,job=job['name'],success=result['complete_success']),flush=True)
    successful = True
    for folder,spec,_,_ in stages:
        groups = [json.loads((folder/j['name']/'result.json').read_text()) for j in spec['jobs']]
        dump(folder/'summary.json',groups)
        successful = successful and all(g['complete_success'] for g in groups)
    return 0 if successful else 2


if __name__ == '__main__':
    raise SystemExit(main())
