# TN808: experimental flap-angle challenge at fixed lift

This addition runs a **qualified reconstruction comparison**, not a certificate
of optimization accuracy. It also fixes a real panel-distribution failure found
while building the case. Numerical results and plots are generated in
`results/tn808_comparison/TN808_RESULTS.md` by `report_tn808.py`.

## Independent source and selected conditions

John G. Lowry (1941), *Wind-Tunnel Investigation of an NACA 23012 Airfoil with
Several Arrangements of Slotted Flaps with Extended Lips*, NACA-TN-808.
[Original NASA report](https://ntrs.nasa.gov/api/citations/19930085051/downloads/19930085051.pdf).
The original PDF, inspected figures, coordinate tables, digitization clicks and
calibration are distributed in `data/naca23012_tn808/`.

Figure8 (PDF18) gives individual polars for slotted flap1-b at a fixed gap of
0.02 original wing chords. We select flap angles10,20,30 degrees and CL1.5.
Each plotted marker belongs to a fixed configuration; the report's envelope
comparisons are not substituted for these data. This tests flap-angle ranking,
not contour-shape gradients or an optimum over all possible geometries.

| Quantity | Selected value or treatment |
|---|---|
| Reference chord | Original NACA23012 wing chord, 3ft; coordinates normalized by that chord |
| Flap extent | Table I's 0.30c modified Clark Y flap |
| Main lip | x/c0.90; upper/lower ordinates0.0168/0.0159 |
| Nominal flap reference station | Nose x/c0.8975, from Figure1 and text |
| Gap interpretation | Vertical distance from lip underside to upper flap at x/c0.90, 0.02c |
| Physical Reynolds number | 2.19 million |
| Historical effective Reynolds number | Approximately3.5 million; **not** substituted for physical Re |
| Mach | Assumed0.105, consistent with the reported approximately80mph at standard conditions |
| Transition | Ncrit9 nominal; Ncrit6 is a sensitivity assumption, not measured tunnel turbulence or a forced trip |
| Objective | Profile CD at CL1.5, fixed-lift tolerance2e-5 |
| Main/flap panels | Base260/180; coarse180/120; transition/datum checks260/180 |

The report does not quantify a correction for hinge fittings. None is subtracted.
Its drag coefficients use the original wing chord, including the retracted flap,
not the deployed total length. MultiFoil retains that reference throughout.

## Geometry reconstruction

`tn808_geometry.py` embeds visually checked Table I ordinates. A standard NACA
five-digit construction resolves the small nose overhang, then smooth ordinate
corrections match the main-airfoil table. The uncorrected analytical formula
differs from the rounded published table by at most0.000300c; the dense corrected
curve interpolates those ordinates within1e-7c. This is a reconstruction check,
not metrology of the original model or a bound on aerodynamic error.

The lower main surface joins the tabulated slot entrance through a prescribed
Hermite segment over x/c0.65–0.6996. The entrance has two ordinates at the same
x station, producing a vertical step. The remaining slot uses monotone coordinate
interpolation in arc length. The modified Clark Y flap uses interpolation in
sqrt(x), imposing its published leading-edge radius0.0045c. Finite trailing-edge
thicknesses are retained. All generated and actual solver geometries are saved.

The tabulated flap nose is (0,0.0105)c and its trailing-edge midpoint is
(0.30,0.00125)c. Their connecting chord is inclined−1.76606 degrees relative to
the table's horizontal axis. The nominal interpretation sets the physical
nose-to-TE chord to the stated deflection. A separate `table_axis` calculation
rotates the tabulated axis by the stated angle instead. Both keep the prescribed
vertical gap. No convention is selected by fitting experimental drag.

Gap/reference-point interpretation, interpolation, the small main-surface join,
unknown transition and unmeasured fitting drag qualify the comparison. The
minimum Euclidean element clearance is not the specified vertical gap: it is
about0.0169,0.0181,0.0191c at the three nominal flap angles. The current study
does not exhaustively bound all geometric ambiguities.

## Panel redistribution defect and fix

The first three attempts were rejected before solving flow: the valid input
main element became self-intersecting inside the solver. At the slot-entry
curvature change, mfoil's cubic interpolation of inverse cumulative spacing
returned nonmonotone arc-length coordinates, folding the panel sequence.

`paneling.py` first evaluates the original algorithm. It returns that result
unchanged when the coordinates are finite and strictly ordered. Otherwise it
reconstructs the same curvature measure and inverts it with a monotone PCHIP
interpolator. It does not replace boundary-layer equations, smooth away the
source step, or suppress the geometry validator. Diagnostics identify each
element that needed recovery.

The regression reproduces the original fold and checks that the recovered
contour has ordered arc lengths, preserved endpoints and no intersections.
Smooth NACA0012 and NLR7301 paneling arrays are bit-for-bit unchanged in the
tested inputs. The full46-test suite passes, including existing aerodynamic
checks. That verifies the fix's tested behavior; it does not validate cove-flow
physics.

## Reference extraction and comparison allowance

`extract_tn808_reference.py` stores9 manually identified marker centres and four
grid-calibration intersections in the archived high-resolution crop. An affine
pixel-to-(CL,CD) mapping accounts for scan skew. Linear interpolation between
recorded bounding markers gives the three values at CL1.5; no extrapolation or
fitting to MultiFoil is used. The source-image overlay makes the clicks reviewable.

| Flap angle | Digitized/interpolated experimental CD at CL1.5 |
|---|---:|
| 10 degrees | 0.03528 |
| 20 degrees | 0.03275 |
| 30 degrees | 0.04246 |

The report states CD repeatability limits at minimum drag, CL1 and CL2.5, but
does not give one at CL1.5. For this challenge we assume the larger stated
limit0.002 applies at CL1.5 and add0.0012 for clicking, calibration, lift mapping
and interpolation. The resulting **±0.0032 per point is a working comparison
allowance**, not a measured total uncertainty or confidence interval. In
particular it does not bound unknown configuration-dependent fitting bias.

On that allowance, the two-point difference allowance is±64 drag counts.
The measured10–20 difference is about25 counts and remains unresolved; the
30–20 difference is about97 counts and retains its sign. This distinction
prevents a small nominal ordering from becoming a false validation success.

Model sensitivity ranges in the plotted comparison are also observed ranges,
not bounds on unknown physical-model error. The strict independent-ranking
comparator is therefore not used to declare a fully audited validation pass.

## Reproduction

From the project root, with the supplied NumPy/SciPy/Matplotlib dependencies:

```bash
python validation/extract_tn808_reference.py
OPENBLAS_NUM_THREADS=1 python validation/tn808_benchmark.py \
  --settings base --workers 3 --out validation/results/tn808_guarded
OPENBLAS_NUM_THREADS=1 python validation/tn808_sensitivity.py \
  --settings coarse,table_axis,ncrit6 --workers 3 \
  --out validation/results/tn808_sensitivity
python -m unittest -v test_multifoil.py test_viscous.py test_paneling.py \
  test_optimization.py test_validation.py test_design_rankings.py
```

Study manifests prevent resuming runs after numerical code or input changes.
Choose a new output folder for changed settings. Failed geometry/flow attempts
retain null forces, their reason and all completed requested flow evaluations.
An exit code2 indicates at least one rejected attempt. It is not a stall label.
Continuation points are retained in each flow's diagnostics.

The original broad-bracket run and the narrower sensitivity procedure are
separate studies. The narrower brackets were chosen from the first converged
model lift points, without fitting drag. Any additional recovery run must be
retained separately and identified in the comparison report.

## Earlier results and exact historical reproduction

The earlier NLR7301 optimization pilot was produced before this fix and retains
its original source hashes. Those archived results are not relabeled as reruns
of the new solver. `audit/before_paneling_guard/` preserves the previous
`viscous.py` and affected driver files with SHA256 hashes. To reproduce an older
study exactly, restore the matching files from that overlay in an isolated copy
of the project. A current-code rerun belongs in a new output folder. New pilot
manifests include `paneling.py` among their numerical dependencies.

## Retained recovery runs

The archived retries used the original helper preserved in the historical
overlay. The current helper additionally preserves exit code2 when returning
a cached failure. For an exact historical rerun restore that helper too; for
current-code retries select new output directories.

The initial coarse10 target search and base30 flow attempt failed, followed by
an unsuccessful fresh initialization for base30. The following commands preserve
separate recovery attempts. They do not overwrite the original failures.

```bash
OPENBLAS_NUM_THREADS=1 python validation/tn808_retry.py \
  --parent validation/results/tn808_guarded --deflection 30 --setting base \
  --bounds -6.05 -5.85 --start-alpha -6.05 --out validation/results/tn808_retry_base30
OPENBLAS_NUM_THREADS=1 python validation/tn808_retry.py \
  --parent validation/results/tn808_sensitivity --deflection 10 --setting coarse \
  --bounds 3.04 3.09 --start-alpha 3.04 --out validation/results/tn808_retry_coarse10
OPENBLAS_NUM_THREADS=1 python validation/tn808_retry.py \
  --parent validation/results/tn808_retry_base30 --deflection 30 --setting base \
  --bounds -6.05 -5.85 --start-alpha -4 --out validation/results/tn808_retry_base30_warm
OPENBLAS_NUM_THREADS=1 python validation/tn808_retry.py \
  --parent validation/results/tn808_sensitivity --deflection 20 --setting ncrit6 \
  --bounds -2.8 -1.8 --start-alpha 0 --out validation/results/tn808_retry_ncrit6_20
OPENBLAS_NUM_THREADS=1 python validation/tn808_retry.py \
  --parent validation/results/tn808_sensitivity --deflection 30 --setting table_axis \
  --bounds -6.55 -6.05 --start-alpha -8 --out validation/results/tn808_retry_datum30
python validation/report_tn808.py
```

A recovered target is not proof of branch uniqueness. In the initial coarse10
trace, two near-identical incidence values differ in CD by about1.01 counts
and in CL by0.000653 following different continuation histories. This diagnostic
is recorded in `results/tn808_sensitivity/continuation_sensitivity.json`.
It concerns this case and is not a measured noise bound for the earlier NLR pilot.

Here `feasible` in a run record means an admissible numerical solution at target
lift. It does not mean that the thin-layer assumptions or the optimization
objective have been experimentally validated.
