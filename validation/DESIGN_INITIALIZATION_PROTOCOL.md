# Mach-zero design initialization integration

Prerequisite: the seven-case Mach-0.07 policy study and saved-state audit passed.
The flow equations, closure, Ncrit and policy helper remain unchanged. This
matrix exercises the explicit policy through `shape_optimization.flow_job`.

There are twelve independent jobs: flap 5° and 10°, fore-camber amplitudes
−0.001, 0 and +0.001, each with a legacy warm-start fixed-lift control and the
new cold-local then fresh-warm policy. All use CL=0.8, Mach=0, Re=700,000,
Ncrit=6.25, 260/180 panels, unit reference chord/wake length, 100 Newton
iterations, residual tolerance 1e-7 and lift tolerance 2e-5. Shape constraints
are the existing design example's restrictions. This is an integration matrix,
not an optimization or independent contour-ranking validation.

Centers are frozen numerical incidences from the earlier baseline study,
identical across amplitudes. The legacy control starts at −1° and searches
center ±0.4°. Policy paths start at the center then −1°, use center ±0.05°,
and permit one lift-only recenter per path inside [−6°, 8°]. Each bracket
allows 20 requested evaluations. No experimental force selects an incidence,
initialization path or geometry. Internal continuation evaluations are extra.

Acceptance requires 12/12 feasible design points, 6/6 complete pairs,
the declared lift/residual tolerances, and absolute paired CD difference
≤0.01 count. Failed points cannot be omitted. On completion, review the saved
selected flow exports, transition diagnostics and restored policy Newton arrays;
the automatic summary checks only coverage, forces and tolerances. Retain any
failure without changing this protocol to make it pass.

Each job saves its configuration, exact deformed coordinates, every requested
flow and initialization context, selected index, timing and failure reason.
Sources and inputs are hashed and snapshotted. Runtime versions are pinned to
Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Six workers use
one BLAS/OpenMP thread each. Completed jobs resume unchanged; partial attempts
are archived before a fresh solver reruns them. Changes require a new directory.

Run or resume from the project root:

```bash
.venv/bin/python validation/run_design_initialization.py \
  --out validation/results/design_initialization_20260921 --workers 6
```

The remote copy is `/home/jony/multifoil-runs/design-initialization-20260921`
on jony4. `launch.log` and `preflight.log` are in that directory. From the
local workspace, `python3 validation/study_status.py --watch` follows completion.
The SSH monitor can disconnect without stopping the job. An interrupted server
run resumes with the command above in its remote directory.

The unchanged 100-chord far-field ordering test has a separately reproduced
baseline failure. This integration does not fix or waive that limitation.
