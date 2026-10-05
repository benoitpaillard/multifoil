# Complete initialization sequence: numerical acceptance matrix

Prerequisite: the explicit-entry diagnostic passed 3/3 cases and saved-state
audit. This fresh study evaluates the complete ordered policy through the
unchanged `shape_optimization.flow_job`, rather than combining results from
different paths after the fact. The original 11/12 study remains failed.

Run six primary jobs: flap 5° and 10°, each with fore-camber amplitudes
−0.001, 0 and +0.001. A seventh job independently repeats the flap-10° baseline.
Each job owns fresh solvers and applies the same ordered rule:

1. Cold local start at the frozen numerical center.
2. Fresh solver starting at −1°, without an explicit seed.
3. Fresh solver starting at −1°, requesting center−0.4° before the local solve.

The first successful path wins; drag never selects a path. The extra entry is
a configuration of the already tested policy helper, not a flow-core change.
All requested failures and internal continuation failures remain in exports.

Keep CL=0.8, Mach=0, Re=700,000, Ncrit=6.25, 260/180 panels, unit reference
chord/wake length, the original design restrictions, 100 Newton iterations,
residual tolerance 1e-7 and lift tolerance 2e-5. Centers are frozen from the
earlier numerical study and identical across amplitudes for each flap setting.
All starts/seeds stay in [−6°, 8°]. Local half-width is 0.05°, with at most two
brackets (one lift-based recenter), 20 requested evaluations per bracket/path.
Worst-case requested budget: 121/job, 847 total; internal continuation is extra.

Acceptance requires 7/7 feasible jobs, absolute CD difference ≤0.01 count
against each contour's archived legacy control, and the original lift/residual
tolerances. The independent repeated baseline must also agree within 0.01 count;
report its actual force and full-state differences rather than inferring a
universal noise bound. After completion, audit all selected surfaces, transition,
Newton arrays and initialization history. Report paired contour changes against
each flap's own baseline. Independent aerodynamic rankings remain unvalidated.

The source/geometry/reference hashes, runtime pin, incremental outputs, lock and
partial-attempt archival follow the entry diagnostic runner. Up to six workers
use one BLAS/OpenMP thread each. All seven jobs are new solves. Completed original
study directories and their source snapshots remain intact; changed runner hashes
deliberately prevent reusing old output directories with the new runner.

Remote root on jony4: `/home/jony/multifoil-runs/design-combined-policy-20260921`.
Launch/resume there:

```bash
.venv/bin/python -u validation/run_design_entry_seed.py --combined \
  --out validation/results/design_combined_policy_20260921 --workers 6
```

Monitor locally:

```bash
python3 validation/run_design_entry_seed.py --combined --remote-status --watch
```

The remote root contains `launch.log` and `preflight.log`. The study directory
contains per-job progress, results, configurations and full selected states.
