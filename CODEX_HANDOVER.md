# MultiFoil: Codex handover

Updated 21 September 2026. The archived-state diagnosis and initialization
repeatability audit are complete. A reusable initialization policy is implemented
as an explicit addition. The seven-job regression is launched on jony4 after
classifying an exactly reproduced pre-existing far-field test failure. The
original numerical core and shape-optimization behavior are unchanged; the
14 September archive manifest and verification receipt are retained.

## Latest completed task and next step

Development priority 1 is complete. Read
[the drag mechanism findings](validation/results/s702_drag_mechanism/FINDINGS.md)
and [paired plots/tables](validation/results/s702_drag_mechanism/README.md).
`validation/report_s702_drag_mechanism.py` audits all 42 selected states (21
conditions, two Ncrit scenarios), including all 16 primary states and the two
separately recovered reserved cases. It reproduces all six flap-delta rows and
verifies selected polar/metadata/surface/wake consistency. All 24 targeted tests
pass; the core is unchanged. Reporting used the available Python 3.10 environment;
original solver versions remain in the provenance. See
`handover/drag_mechanism_verification.json` for the new receipt.

At CL=0.8, the calibrated 5-to-10-degree fore wake drag changes are -31.54 and
-21.72 counts at Re=1,000,000 and 700,000; aft penalties of +10.37 and +20.38
counts offset those benefits. Fore upper transition is much later at 10 degrees,
while the aft upper layer thickens and aft lower separation appears. At CL=1.2,
the smaller fore benefit cannot outweigh the aft penalty. These are modeled
associations, not an experimental allocation of the error to either element.
At Re=1,000,000/CL=0.8, calibration changes the ranking mostly by increasing
the 5-degree drag (+27.52 counts), rather than a large change at 10 degrees
(-2.97 counts). This does not establish correct geometry sensitivities.

**Repeatability audit completed 21 September:** the full original study and
runtime logs have been retrieved. Read
[the findings](validation/results/s702_repeatability_analysis/FINDINGS.md).
All 18 successful targets converged on their first attempt; the independent cold
flap-10/CL=0.8 target failed all three attempts. Increasing, decreasing and
archived replay paths solve that condition with matching selected states.
The largest CD span over five fully covered conditions is 2.77e-7 count;
the largest difference from the older archived selected results is 0.0027721
count. The audit checks 24 selected states, including six references, their
surface/wake exports and saved Newton arrays. Turbulence flags and surface
index lists match across all compared successful states. The missing cold path
prevents a complete score at flap 10/CL=0.8; this is not universal path independence.

**Common policy implemented, numerical regression launched:**
`fixed_lift_policy.py` adds ordered fresh initialization paths, one bounded
computed-lift recenter when justified, complete failed-evaluation provenance,
and restoration of the complete selected solver in memory. This is not a
portable checkpoint. The original `optimization.solve_fixed_lift`, viscous
equations and sweep worker remain unchanged. All 38 fast contract/data tests
pass locally. The original snapshot audit matches 11,053 files and reports only
the four previously edited documentation/status files.

`validation/S702_INITIALIZATION_POLICY_PROTOCOL.md` declares six primary
conditions with cold-local then warm-at-minus-one-degree initialization, and a
seventh exact historical warm-replay control. A simple warm anchor is a hypothesis
until this run is reviewed. Its isolated remote root is
`/home/jony/multifoil-runs/s702-policy-20260921`; the existing pinned Python 3.12.13
environment is reused without package changes. Seven single-threaded workers;
the original preflight passed 45/46 tests. The unchanged far-field element-order
test fails identically in the untouched baseline: CL difference
1.1415574480699853e-6 exceeds 5e-7 at 100-chord separation. This existing
limitation remains recorded; no threshold or physics was changed. The new
continuation verifies both logs, all frozen sources and the shared pinned
environment, then reruns 38 fast checks before the seven S702 cases. The seven
unaffected viscous checks are reused from the verified preflight, not claimed
as a fresh full-suite pass. The local study remains preparation only until
remote outputs are retrieved.

**Follow completion:** `python3 handover/policy_status.py --watch` from this
project, or omit `--watch` for one check. It exits at completion/failure; Ctrl-C
only stops watching. Latest launch receipt:
`handover/policy_resume_launch_20260921.json`. See
[remote run instructions](handover/S702_REMOTE_RUN.md) for logs and retrieval.
**Next session:** retrieve the completed follow-up, then audit coverage and selected full states before promoting
the policy into shape optimization. Do not relaunch the completed original matrix or retune
Ncrit. Independent aerodynamic rankings remain unvalidated. SSH works unattended
through `jony4-anywhere`; setup notes are in `/home/jony/remoteConnection/README.md`.

## Objective and current decision

Develop an open, reproducible multi-element analogue of XFOIL, ultimately useful
for incompressible shape optimization. Viscosity, natural transition, limited
separation, fixed-lift evaluation and a bounded contour optimizer are implemented.
The immediate development objective is to improve the reliability of predicted
drag changes between geometries.

Use the current model for constrained design screening. A fitted transition
threshold improves baseline drag but does not establish correct shape gradients.
The most recent incompressible shape trial retained the baseline. No experimentally
validated contour optimum has been accepted.

Read this file, then [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). Start from the
archived data rather than repeating the completed parameter sweep. A ready-to-use
initial instruction is in [CODEX_START_PROMPT.md](CODEX_START_PROMPT.md).

## What was learned

### Transition explains part of the S702 drag discrepancy

The S702 study samples 14 common Ncrit values: 3, 4, 5, 5.5, 5.75, 6, 6.25,
6.5, 7, 8, 9, 10, 11 and 12. Calibration covers four flap settings (-2, 0, 5,
10 degrees), three lift targets (0.5, 0.8, 1.2), Re = 1,000,000 and M = 0.1.

| Comparison | MultiFoil Ncrit 9 | MultiFoil Ncrit 6.25 |
| --- | ---: | ---: |
| Calibration mean absolute relative drag error, 12 points | 8.84% | 4.85% |
| Reserved Re 700,000 / M 0.07 mean error, 9 points | 10.15% | 5.70% |
| Worst reserved point error | 25.88% | 24.04% |

Ncrit = 6.25 minimizes the predeclared relative-error objective. Ncrit = 6
minimizes absolute drag RMSE and flap-difference RMSE; its MAPE is 4.95%.
The difference in MAPE is only 0.10 percentage points. Finer panels retain the
same broad minimum. Treat approximately 6-7 as a useful sensitivity region for
this case, not a measured tunnel Ncrit or a universal physical correction.
The wider sampled set within one percentage point of the best score is
5.75, 6, 6.25, 6.5 and 7; that is an engineering band, not a confidence interval.

There are 168 calibration entries: 165 converged and 24 reused from matching
prior studies. Three candidates have missing targets and are not assigned a
complete score. The recorded 649 new requested flow evaluations include 11
failed evaluations; these are different counts from failed final targets and
do not include internal continuation steps.

Leave-one-flap-out checks use only the predeclared integer grid. Combined
excluded-flap MAPE is 5.63%, but the excluded -2-degree flap gets worse
(5.93% at Ncrit 9 versus 7.60% at the selected Ncrit 6). These data had already
been inspected, so this is an internal grouped check, not blind validation.

### The flap-ranking discrepancy remains important

The table reports CD(10 degrees) - CD(5 degrees) at equal lift. Negative means
that the 10-degree setting has lower drag. One count is 0.0001 CD.

| Re | CL | Experiment | Published MSES Ncrit 9 | MultiFoil Ncrit 9 | MultiFoil Ncrit 6.25 |
| --- | --- | ---: | ---: | ---: | ---: |
| 1,000,000 | 0.5 | +8.10 | +14.64 | +15.20 | +7.99 |
| 1,000,000 | 0.8 | -29.86 | +4.05 | +9.32 | -21.17 |
| 1,000,000 | 1.2 | -6.45 | +2.49 | +2.96 | +0.74 |
| 700,000 | 0.5 | +12.72 | +25.45 | +27.44 | +16.66 |
| 700,000 | 0.8 | -20.72 | +12.65 | +15.17 | -1.34 |
| 700,000 | 1.2 | -8.76 | +9.38 | +11.79 | +5.79 |

Calibration corrects the raw 5-to-10-degree ordering at CL = 0.8 but severely
understates the benefit at Re = 700,000. At CL = 1.2 all numerical variants
still predict the opposite direction from experiment. Improving mean error is
therefore insufficient for optimization. The latest
[total-drag and baseline-relative plots](validation/results/s702_flap_effect/README.md)
make these differences explicit.

At Re = 700,000, flap 10 degrees and CL = 0.8, calibrated MultiFoil predicts
CD = 0.0100257243 versus experiment 0.0080824909: +24.04%. Published MSES
Ncrit 9 gives 0.0099748848. MSES was not run in this project; the reference
values are extracted from a published comparison.

### Numerical initialization matters independently of the physics

The reserved Ncrit 6.25 run initially completed 7 of 9 targets. A separate
follow-up recovered flap 10 / CL 0.5 with a local cold start, and flap 10 /
CL 0.8 with a warm replay after the cold retry failed. Ncrit, equations, panels
and tolerances were unchanged. The two replayed starting CL/CD states matched
the original states exactly. The 5.70% score uses all nine points after those
recoveries; it must not be described as nine successes under the original
protocol. The threshold was frozen before the reserved calculations and was
not retuned using their errors. Reserved measurements had been extracted and
inspected before freezing; they were excluded from the selection loss, not blind.

The first successful initialization in the documented cold-then-warm order is
used. Numerical lift selects the incidence seed; experimental drag does not.
Original failures and unsuccessful retries remain archived. The 15 execution-
service interruptions across calibration and reserved runs are also retained
separately from numerical nonconvergence. A resumed run is not a restored full
Newton checkpoint.

### Mach zero and finer panels did not remove the experimental discrepancy

The original Ncrit 9 Mach-zero check changes the CL = 0.8 ranking margin by
only +0.036 count. The later Ncrit 6.25 check changes individual CD by at most
0.169 count at the two examined settings. These are sensitivity calculations;
experimental/MSES comparisons remain at M = 0.1 or 0.07.

At Ncrit 6.25, refining 260/180 to 340/240 panels shifts CD by at most about
0.55 count over the checked calibration points. At CL = 1.2 the 10-minus-5
difference changes from +0.74 to +0.50 count; experiment gives -6.45 counts.
These observations weaken simple compressibility and panel-density explanations
for this discrepancy. They are not formal numerical uncertainty bounds.

The separate Cp comparison is at flap 0, CL = 1.04, Re = 1,000,000, M = 0.1.
Calibrating Ncrit improves some pressure differences but leaves substantial
ones. Those taps do not directly validate the difficult 5/10-degree states.

### Shape optimization exists, but physical trust is still limited

`shape_optimization.py` implements a bounded derivative-free pattern search
over smooth camber modes, with paired fixed-lift scenarios, geometry/moment/
off-design restrictions and finer-grid screening of a qualifying candidate.
The executed S702 incompressible trial achieved only 0.19 and 0.47 counts of
mean improvement under Ncrit 9 and 6. Both did not exceed its predeclared
0.5-count requirement, so the baseline was retained and no candidate fine-grid
flow check was triggered. That threshold belongs to this screening protocol;
it is not an experimental precision claim.

The earlier NLR7301 pilot showed transition-dependent camber preferences and
combined panel/initialization sensitivity. A single-element native-mfoil check
matched the sign of a local sensitivity but differed by 8.34% in magnitude.
It shares boundary-layer closures with MultiFoil and is implementation
verification, not independent aerodynamic validation.

## Code and evidence map

| Location | Role |
| --- | --- |
| `multifoil.py` | Geometry handling, inviscid backend, CLI and exports |
| `viscous.py` | Coupled multi-body viscous solve, wakes, Newton system and diagnostics |
| `vendor/mfoil.py` | Bundled unchanged MIT-licensed boundary-layer closures |
| `paneling.py` | Conditional monotone panel-spacing recovery; keeps valid upstream spacing |
| `optimization.py` | Fixed-lift solve, admissibility, camber perturbation and clearance |
| `fixed_lift_policy.py` | Explicit initialization and selected-state restoration; opt-in pending regression review |
| `validation/report_s702_repeatability.py` | Reproducible audit and plots of the completed original matrix |
| `validation/run_s702_initialization_policy.py` | Restartable seven-job follow-up for the common policy |
| `shape_optimization.py` | Multipoint/scenario design screening and pattern search |
| `validation/s702_ncrit_sweep.py` | Audited sweep, reuse, local/warm retry drivers |
| `validation/report_s702_ncrit_sweep.py` | Source verification, reference interpolation, state-path resolution and scores |
| `validation/check_s702_ncrit_transfer.py` | Reserved, fine-panel, Mach-zero, Cp and retry comparisons |
| `validation/plot_s702_flap_effect.py` | Latest flap-angle plots; no solver execution |
| `handover/status.json` | Machine-readable findings, conditions, hashes and priorities |

`validation/results/s702_ncrit_analysis/` is the principal results entry point.
Start with `NCRIT_SWEEP_RESULTS.md`, `summary.json`, `transfer_summary.json`,
`frozen_selection.json`, `points.csv` and `holdout_after_retries.csv`.
`holdout_points.csv` preserves original coverage. `local_retry_points.csv`
distinguishes the successful and unsuccessful follow-ups.

All original studies are included, not just the S702 summary. Earlier benchmarks
cover NACA 0012, S809, NLR7301, TN808 and an unscored TR664 candidate. The
NLR7301 result near 13 degrees remains about 51% high in drag. S809 transition
and stall comparisons and NACA 0012's abrupt lift loss expose limitations.
TN808 exposed a panel-spacing fold that was fixed, but its sharp-cove drag
disagreement and about one-count continuation sensitivity remain unresolved.
See `validation/README.md`, `VALIDATION.md` and the case reports before
generalizing S702 improvements.

The S702 geometry is reconstructed from vector Figure 1 of
[Somers and Maughmer, NASA/CR-20220014964](https://ntrs.nasa.gov/citations/20220014964).
Exact author coordinates are unavailable. The small plain flap is on the aft
element; its deflection is not a rotation of that entire element. Preserve the
assembly reference chord, sign conventions and actual exported geometry.

## Model facts to preserve during diagnosis

Natural transition uses a simplified e^N envelope. Forced trips, roughness,
wake-induced bypass transition, wake merging and boundary-layer/wake confluence
are not implemented. Wakes follow the global inviscid field and are not relaxed
with viscous displacement. Large separation, dynamic stall and shocks are beyond
the model's demonstrated scope. Missing physics are hypotheses for investigation,
not established causes of the S702 discrepancy.

Total profile CD is the sum of per-element Squire-Young wake estimates.
`diagnostics.element_forces` contains `CD_wake` and `CD_friction`; friction is a
diagnostic component, not an extra term to add to wake CD. Inviscid pressure drag
is a discretization diagnostic. There is no three-dimensional induced drag here.
Surface `H` is Hk; `amplification_or_shear` changes meaning after transition.
The report loader's `flow_directory()` resolves the old reused Ncrit 9 path
label correctly. Use it instead of assuming all legacy result paths are uniform.

## Start and verify

Use a Python 3.12 interpreter for the pinned environment (`python3` below must
refer to Python 3.12). From the extracted `multifoil/` directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-handover.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python handover/verify_snapshot.py
python -m unittest -v test_s702_sweep test_s702 test_s702_pressure
python validation/plot_s702_flap_effect.py --out validation/results/handover_plot_check
```

The pinned environment records Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8 and optional extraction dependency PyMuPDF 1.26.6. Core
requirements remain in `requirements.txt`. Installation needs access to the
packages; no virtual environment or dependency wheels are bundled.

The snapshot checker verifies every delivered file listed in its manifest.
The 13 fast S702 audit tests check extraction, provenance and calibration rules;
they do not rerun the flow solver or certify aerodynamic accuracy. The exact
handover checks are recorded in `handover/verification.json` and its test log.
After numerical-core changes, also run the relevant existing physics tests;
the complete suite command is `python -m unittest discover -p 'test_*.py' -v`.

The optional full numerical reproduction is:

```bash
python validation/run_s702_calibration.py --out validation/results/reproduced_s702_calibration --workers 4
```

Use a fresh directory inside the project; choose 1-8 workers appropriate to RAM.
This can take hours and reuses eligible archived Ncrit 6/9 calibration points.
It is not a fully cold rerun of every historical solve. Use the individual sweep
driver's `--no-reuse` for an explicitly fresh numerical experiment, under a new
study manifest. Read `validation/NCRIT_SWEEP_PROTOCOL.md` before replay/resume.

The archive includes code, tests, examples, raw/completed/failed study outputs,
source snapshots, available source documents, extraction audits, figures and
licenses. It excludes only runtime caches. Historical absolute paths inside
logs identify the old execution location; current entry points derive the
project root from their files. Additional literature downloads outside the
project are not dependencies of the bundled calculations.
