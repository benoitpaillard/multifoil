#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python - <<'PY'
import hashlib,json
from pathlib import Path
p=json.loads(Path('validation/design_entry_seed_payload.json').read_text())
for name,digest in p['sha256'].items():
    if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
        raise ValueError('Deployment changed: '+name)
print('Verified deployment:',len(p['sha256']),'files',flush=True)
PY
.venv/bin/python -m unittest -v test_shape_initialization test_shape_optimization \
  test_fixed_lift_policy test_optimization test_paneling > preflight.log 2>&1
exec .venv/bin/python -u validation/run_design_entry_seed.py \
  --out validation/results/design_entry_seed_20260921 --workers 3
