#!/usr/bin/env python3
"""Bound the complete detached launch, including preflight and worker cleanup."""
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import run_s414_nose_resolution as run

TESTS = ('test_s414_nose_resolution', 'test_s414_nose_resolution_runner', 'test_s414_nose_resolution_analysis')
CLEANUP_RESERVE_SECONDS = 10


def worker_groups(root):
    """Return only recorded workers still belonging to this isolated launch."""
    groups = []
    for path in (root/run.STUDY).glob('*/attempt.json'):
        attempt = run.read(path); pid = attempt.get('pid')
        if not isinstance(pid, int): continue
        try:
            proc = Path('/proc')/str(pid)
            argv = (proc/'cmdline').read_bytes().split(b'\0')
            if ((proc/'cwd').resolve() == root.resolve() and b'--worker' in argv
                    and attempt['job'].encode() in argv
                    and any(a.endswith(b'run_s414_nose_resolution.py') for a in argv)
                    and os.getpgid(pid) == pid):
                groups.append(pid)
        except (OSError, ProcessLookupError): pass
    return groups


def kill_group(pid, sig):
    try: os.killpg(pid, sig)
    except ProcessLookupError: pass


def process_identity(pid):
    try:
        fields = (Path('/proc')/str(pid)/'stat').read_text().rsplit(')', 1)[1].split()
        return int(fields[1]), fields[19]  # Parent PID and kernel start tick.
    except (OSError, IndexError, ValueError): return None


def descendants(pid):
    """Capture detached preflight children too; start ticks prevent PID-reuse mistakes."""
    processes = {int(p.name): process_identity(int(p.name)) for p in Path('/proc').iterdir() if p.name.isdigit()}
    found = {}; parents = {pid}
    while parents:
        children = {p: identity for p, identity in processes.items() if identity and identity[0] in parents}
        found.update(children); parents = set(children)
    return found


def stop(process, root, grace_seconds=5):
    """Ask the supervisor to clean up, then kill any detached workers left behind."""
    children = descendants(process.pid)
    if process.poll() is None: kill_group(process.pid, signal.SIGTERM)
    try: process.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        kill_group(process.pid, signal.SIGKILL); process.wait(timeout=1)
    for pid in worker_groups(root): kill_group(pid, signal.SIGKILL)
    for pid, identity in children.items():
        current = process_identity(pid)
        if current and current[1] == identity[1]:
            try: os.kill(pid, signal.SIGKILL)
            except ProcessLookupError: pass


def execute(command, root, log, timeout_seconds, grace_seconds=5):
    if timeout_seconds <= 0: raise TimeoutError('Complete launch budget exhausted before next phase')
    process = subprocess.Popen(command, cwd=root, stdin=subprocess.DEVNULL, stdout=log,
        stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)
    try:
        code = process.wait(timeout=timeout_seconds)
        if code: raise RuntimeError(f'Launch phase failed with exit {code}: {command[1:]}')
    except BaseException:
        stop(process, root, grace_seconds); raise


def signal_stop(signum, frame):
    raise TimeoutError(f'Launch interrupted by signal {signum}; cleaning up all workers')


def main():
    root = run.ROOT; out = root/run.STUDY; spec = run.read(out/'study.json')
    budget_path = out/'budget.json'
    if budget_path.exists():
        budget = run.read(budget_path)
        if (budget['study_sha256'] != run.sha(out/'study.json') or
                budget['limit_seconds'] != spec['study_timeout_seconds'] or
                budget.get('scope') != 'complete_remote_launch_including_preflight' or
                budget.get('cleanup_reserve_seconds') != CLEANUP_RESERVE_SECONDS):
            raise ValueError('Recorded launch budget differs')
    else:
        started = time.time()
        budget = dict(started_utc=run.stamp(), started_unix=started,
            deadline_unix=started+spec['study_timeout_seconds'], limit_seconds=spec['study_timeout_seconds'],
            study_sha256=run.sha(out/'study.json'), scope='complete_remote_launch_including_preflight',
            cleanup_reserve_seconds=CLEANUP_RESERVE_SECONDS)
        run.write(budget_path, budget)
    signal.signal(signal.SIGTERM, signal_stop); signal.signal(signal.SIGINT, signal_stop)
    deadline = time.monotonic()+max(0., budget['deadline_unix']-time.time())-CLEANUP_RESERVE_SECONDS
    try:
        run.check_inputs(spec)
        commands = [
            [sys.executable, '-m', 'unittest', '-v', *TESTS],
            [sys.executable, 'validation/check_s414_nose_resolution.py', '--output', 'linear_preflight.json'],
            [sys.executable, 'validation/run_s414_nose_resolution.py', '--workers', '2'],
        ]
        for index, command in enumerate(commands):
            log_path = root/('preflight.log' if index < 2 else 'supervisor.log')
            with log_path.open('ab', buffering=0) as log:
                execute(command, root, log, deadline-time.monotonic())
        return 0
    except BaseException as exc:
        # During preflight there may be no numerical supervisor/status yet.
        run.write(out/'status.json', dict(state='execution_failed', updated_utc=run.stamp(), pid=os.getpid(),
            total_jobs=6, completed_jobs=sum((out/j['name']/'result.json').exists() for j in spec['jobs']),
            error=f'{type(exc).__name__}: {exc}', budget_scope=budget['scope']))
        raise
    finally:
        for pid in worker_groups(root): kill_group(pid, signal.SIGKILL)


if __name__ == '__main__': raise SystemExit(main())
