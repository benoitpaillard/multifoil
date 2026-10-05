#!/usr/bin/env bash
# Detached, restartable environment setup and execution for the prepared matrix.
set -Eeuo pipefail
multifoil_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd -- "$multifoil_root"
multifoil_workers=${1:-8}
case "$multifoil_workers" in [1-8]) ;; *) exit 64 ;; esac
multifoil_runtime="$multifoil_root/handover/remote_runtime"
mkdir -p -- "$multifoil_runtime"
exec 9>"$multifoil_runtime/launch.lock"
flock -n 9 || { printf 'Another launcher is active\n' >&2; exit 73; }
exec >>"$multifoil_runtime/launch.log" 2>&1

phase() {
    python3 - "$multifoil_runtime/launch_status.json" "$1" "${2:-0}" "$$" <<'PY'
from datetime import datetime, timezone
import json, pathlib, socket, sys
path=pathlib.Path(sys.argv[1])
value=dict(phase=sys.argv[2],exit_code=int(sys.argv[3]),launcher_pid=int(sys.argv[4]),
           hostname=socket.gethostname(),updated_utc=datetime.now(timezone.utc).isoformat())
temp=path.with_suffix('.tmp'); temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)
print(json.dumps(value),flush=True)
PY
}
trap 'multifoil_exit=$?; phase failed "$multifoil_exit"; exit "$multifoil_exit"' ERR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
export UV_PYTHON_INSTALL_DIR="$multifoil_root/.python"
export UV_CACHE_DIR="$multifoil_root/.cache/uv"
multifoil_uv="$multifoil_root/.tools/uv"
phase verifying_payload
python3 - <<'PY'
import hashlib,json,pathlib
root=pathlib.Path.cwd()
manifest=json.loads((root/'handover/remote_payload_20260916.json').read_text())
for name, expected in manifest['sha256'].items():
    actual=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if actual!=expected:
        raise ValueError('Remote payload changed: '+name)
print('Verified payload files:',len(manifest['sha256']),flush=True)
PY
phase installing_python
# 3.12.14 is absent from this uv catalogue; 3.12.13 is also used in archived
# S702 metadata. Keep the exact numerical dependency pins and record runtime.
"$multifoil_uv" python install 3.12.13 --no-bin
if [[ ! -x "$multifoil_root/.venv/bin/python" ]]; then
    "$multifoil_uv" venv --python 3.12.13 "$multifoil_root/.venv"
fi
phase installing_dependencies
"$multifoil_uv" pip install --python "$multifoil_root/.venv/bin/python" -r requirements-handover.txt
phase preflight_tests
.venv/bin/python -m unittest -v test_s702_repeatability test_optimization test_s702_sweep test_s702 test_s702_pressure
.venv/bin/python validation/run_s702_repeatability.py \
    --out validation/results/s702_repeatability_20260915 --prepare-only
phase running_study
if .venv/bin/python validation/run_s702_repeatability.py \
    --out validation/results/s702_repeatability_20260915 --workers "$multifoil_workers"; then
    phase finished
else
    multifoil_exit=$?
    if [[ "$multifoil_exit" == 2 ]]; then
        phase finished_with_numerical_failures 2
    else
        phase failed "$multifoil_exit"
    fi
    exit "$multifoil_exit"
fi
