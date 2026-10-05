# S702 transition calibration protocol

The objective is to determine whether a single effective natural-transition
amplification threshold improves drag predictions across operating points and
flap settings. Geometry, viscous closures, reference chord and the solver's
default Ncrit are unchanged. Ncrit is common to both elements.

The predeclared coarse grid is 3 through 12 at integer spacing. Calibration uses
flaps -2, 0, 5 and 10 degrees at CL=0.5, 0.8 and 1.2, Re=1,000,000 and M=0.1.
The primary score is the equally weighted mean absolute relative drag error
(MAPE) across those 12 points. RMSE in drag counts, maximum relative error and
all pairwise flap drag differences are reported separately. One drag count is
0.0001. No experimental standard deviations are assumed or invented.

A candidate needs complete converged coverage to qualify for the global score.
Failed targets have no finite drag and cannot be silently omitted. Every
requested flow evaluation, internal continuation metadata, attempted bracket,
surface distribution and converged target state is retained. Failure of the
Newton solve is not classified as physical stall.

Base panel counts are 260 fore / 180 aft, wake length is one assembly chord,
the nonlinear tolerance is 1e-7 and the fixed-lift tolerance is 2e-5. The previous
Ncrit=6 and 9 points are reused only at matching physical conditions, mesh and
convergence criteria. Their earlier incidence paths remain explicit. New runs
use one solver per flap/Ncrit pair, advancing through increasing lift targets;
initial brackets use only prior numerical incidences, never experimental drag.
All bracket retries are preserved. The selected flow state is restored before
transition residual checks, including when a bracket endpoint is the solution.

Leave-one-flap-out checks select Ncrit using the other three flaps and then
evaluate the excluded flap. Only the predeclared integer grid enters those
checks: a globally chosen refinement grid would otherwise leak the held-out
flap errors into candidate selection. The existing flap data were already
inspected in earlier studies, so these are grouped calibration checks rather
than a blind experiment.

Refine the complete coarse-grid minimum within its two neighboring integer
values at half-unit spacing, and use quarter-unit spacing if the minimum or a
transition jump warrants it. Sampled values within one percentage point of the
best MAPE define an engineering sensitivity range, not a statistical confidence
interval. No interpolation through failed values establishes an admissible
range. Finer panels check the promising values and the sensitive drag ranking.

The reserved Reynolds-number check is Re=700,000, M=0.07, flaps 0, 5 and 10,
at the same three lift targets. NASA Figures 23(b), 25(a) and 26(b) provide
previously unused measurements and published MSES Ncrit=9 curves. There is no
published -2-degree case at this Reynolds number. Freeze the globally selected
Ncrit before running/scoring the reserved cases; compare with a fresh Ncrit=9
baseline and do not retune using the reserved errors. Independent pressure
data at flap 0, CL=1.04 supply an additional physical observable. Mach-zero
checks assess the intended incompressible application; they do not constitute
matched-Mach experimental validation.

The existing code represents natural e^N transition. Sweeping Ncrit does not
add forced trips, bypass transition, wake-induced transition, wake merging,
as-built geometry or the experimental hinge-cover film. A fitted threshold may
compensate for those omissions. Published MSES reference curves remain Ncrit=9;
no MSES executable is run during this study.

For shape optimization, freeze the calibrated physical assumptions while
varying geometry and evaluate improvements over the credible transition range.
Agreement on the baseline polars does not establish the accuracy of local
contour-design gradients.

## Reproduction

From the project root, run with Python, NumPy, SciPy and Matplotlib installed:

The complete driver reproduces calibration, freezes the selection, then runs
the reserved Reynolds-number, mesh, continuation, Mach-zero and pressure
checks before writing the report:

```bash
python validation/run_s702_calibration.py --out validation/results/reproduced_s702_calibration --workers 8
```

Choose a fresh output directory inside the project. The full study can take
hours. Existing Ncrit=6/9 calibration results are reused only after checking
their conditions and source hashes. Numerical failures remain in the results.
To run or inspect the individual stages instead:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python validation/s702_ncrit_sweep.py --workers 8 --out validation/results/my_ncrit_coarse
python validation/report_s702_ncrit_sweep.py --coarse validation/results/my_ncrit_coarse --out validation/results/my_ncrit_analysis
python -m unittest test_s702_sweep test_s702 test_s702_pressure
```

Use `--ncrit`, `--flaps`, `--targets`, `--panels`, `--mach`, `--reynolds`,
`--no-reuse` and a new `--out` directory to run the documented follow-up checks.
The completed results report records their exact selections. A manifest guards
each output directory against mixing changed source or inputs; completed
targets resume only under the original manifest. The first coarse runner's
exact source is archived under its `source_snapshot` directory. Its legacy
Ncrit=9 `flow_directory` label omitted the old per-CL folder; the report resolves
that path from the authoritative source result and lift target. Later runners
write the corrected metadata path. Numerical outputs are unaffected.

Reference: Somers and Maughmer (2022),
[NASA/CR-20220014964](https://ntrs.nasa.gov/citations/20220014964).
Transition interpretation and resolution guidance:
[XFOIL documentation](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt).

For an interrupted execution, preserve partial attempts and use the original
archived worker version with the Linux resume helper:

```bash
python validation/resume_s702_sweeps.py validation/results/s702_ncrit_coarse validation/results/s702_ncrit_refined --workers 8
```

This helper verifies every numerical/geometry dependency, executes the exact
recorded runner source, archives unfinished target folders, and reuses completed
targets. It does not replace the aerodynamic equations or relax convergence
criteria. Each resume event is recorded separately from numerical failures.

## Separate local initialization follow-up

The reserved check exposed an incidence-bracketing limitation: the original
seeds come from Re=1,000,000, and a symmetric bracket expansion at Re=700,000
can traverse an inadmissible lower-lift region. After a target has failed the
original protocol, `--seed-from PATH/TO/FAILED/result.json` enables a separate
study using the nearest two admissible computed lift values to estimate the
target incidence. It cold-starts there with an initial half-width of 0.05°.
It requires the same Ncrit, flap, lift, Reynolds number, Mach and panel counts;
experimental drag never enters this estimate. Original target files are
immutable, and the retry manifest hashes both the failed result and its study.

This numerical follow-up was added after the reserved run revealed the
bracketing problem. It is not part of the predeclared selection protocol.
The full driver now performs eligible reserved retries automatically. The
report keeps the original coverage and scores, then reports any recovery
separately; hatched bars identify successful local retries. No new Ncrit is
selected from these results.

If a cold local initialization fails, `--warm-seed` replays up to the first two
consecutive admissible source evaluations, starting from the usual -1° initial
incidence, before entering the same local bracket. Warm evaluations are labeled
`warm_seed` in the records. The complete driver tries cold and then warm; the
first successful initialization supplies the follow-up score. Failed attempts
and any later initialization checks remain visible. Drag error does not select
the initialization path. In the executed study, the warm CL=0.8 diagnostic was
started after the cold initializer failed, while its remaining retries finished.
