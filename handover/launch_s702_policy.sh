#!/usr/bin/env bash
# Restartable, detached follow-up; uses the existing pinned Python environment.
set -Eeuo pipefail
multifoil_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd -- "$multifoil_root"
multifoil_workers=${1:-7}
case "$multifoil_workers" in [1-7]) ;; *) exit 64 ;; esac
multifoil_runtime="$multifoil_root/handover/policy_runtime"
mkdir -p -- "$multifoil_runtime"
exec 9>"$multifoil_runtime/launch.lock"
flock -n 9 || { printf 'Another policy launcher is active\n' >&2; exit 73; }
exec >>"$multifoil_runtime/launch.log" 2>&1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1

phase() {
    python3 - "$multifoil_runtime/launch_status.json" "$1" "${2:-0}" "$$" <<'PY'
from datetime import datetime, timezone
import json, pathlib, socket, sys
p=pathlib.Path(sys.argv[1])
data=dict(phase=sys.argv[2], exit_code=int(sys.argv[3]), launcher_pid=int(sys.argv[4]),
          hostname=socket.gethostname(), updated_utc=datetime.now(timezone.utc).isoformat())
t=p.with_suffix('.tmp'); t.write_text(json.dumps(data,indent=2)+'\n'); t.replace(p)
print(json.dumps(data),flush=True)
PY
}
trap 'multifoil_exit=$?; phase execution_failed "$multifoil_exit"; exit "$multifoil_exit"' ERR
phase verifying_payload
python3 - <<'PY'
import hashlib,json,pathlib
r=pathlib.Path.cwd()
m=json.loads((r/'handover/policy_payload_20260921.json').read_text())
for name, expected in m['sha256'].items():
    if hashlib.sha256((r/name).read_bytes()).hexdigest()!=expected:
        raise ValueError('Remote payload changed: '+name)
print('Verified payload files:',len(m['sha256']),flush=True)
PY
phase preflight_tests
.venv/bin/python -m unittest -v test_fixed_lift_policy test_s702_initialization_policy \
    test_s702_repeatability test_optimization test_s702_sweep test_s702 test_s702_pressure test_viscous
.venv/bin/python validation/run_s702_initialization_policy.py \
    --out validation/results/s702_initialization_policy_20260921 --prepare-only
phase running_study
if .venv/bin/python validation/run_s702_initialization_policy.py \
    --out validation/results/s702_initialization_policy_20260921 --workers "$multifoil_workers"; then
    phase finished
else
    multifoil_exit=$?
    if [[ "$multifoil_exit" == 2 ]]; then
        phase finished_with_regression_failures 2
    else
        phase execution_failed "$multifoil_exit"
    fi
    exit "$multifoil_exit"
fi
