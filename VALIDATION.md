# Validation and numerical verification record

For the consolidated current interpretation and remaining work, read
[CODEX_HANDOVER.md](CODEX_HANDOVER.md). This file retains the historical
verification stages; the completed Ncrit calibration and the latest flap-angle
plots are documented in `validation/results/s702_ncrit_analysis/` and
`validation/results/s702_flap_effect/`. The 14 September handover adds no new
aerodynamic calculations; its checks are recorded in `handover/verification.json`.

Executed 2026-09-08, Linux x86_64, Python 3.12.13, NumPy 2.3.5,
SciPy 1.17.0, Matplotlib 3.10.8.

The original numerical checks below are now supplemented by published
wind-tunnel comparisons in `validation/`. The demo geometries discussed below
remain synthetic. No independent multi-element RANS validation has been
performed, and the prototype is not a validated engineering design tool.

See `validation/README.md` for data provenance and qualifications and
`validation/results/REPORT.md` for executed error metrics, coverage and plots.
The comparison suite does not tune aerodynamic closures to the reference data.

## Incompressible contour-design workflow added

`shape_optimization.py` adds fixed-lift, Mach-zero contour optimization with
smooth camber modes, transition-scenario comparisons, design restrictions and
a mandatory finer-grid check before exporting an accepted changed shape.
See `validation/INCOMPRESSIBLE_DESIGN.md` and the executed
`validation/results/s702_incompressible_design/INCOMPRESSIBLE_RESULTS.md`.

The S702 Mach check repeats the 5/10-degree flap comparison at CL=0.8.
Changing Mach 0.1 to zero changes the ranking margin only from 9.31836 to
9.35407 counts and does not reverse the preference for 5 degrees. Published
MSES and experiment remain at Mach 0.1; the incompressible calculations are
a separate sensitivity study.

The new workflow passed 20 checks covering shape invariants, constrained
objectives, budgets, fixed lift, paneling and the actual final-geometry export
decision under a finer-grid reversal. A three-shape NACA0012 comparison against
native mfoil at Mach zero and fixed CL=0.5 gave drag differences no greater than
0.102 counts and the same central sensitivity sign, with an 8.34% difference in
magnitude. Native Python mfoil required documented array-layout adapters in
the comparison runner. This is numerical consistency with shared closures,
not independent physical validation. The core MultiFoil equations, vendored
source and prior experimental/MSES results are unchanged.

## Published wind-tunnel comparisons added

The three experimental benchmark families requested 44 operating points; 42 converged.
Six comparison-integrity tests passed at that stage. The current suite and its
results are recorded in `validation/test_results.txt`. The
previous 20 physics checks remain recorded below; the aerodynamic solver was
not modified for this validation addition.

* **NACA 0012:** all 18 points converged. CL RMSE is 0.0262 for the seven
  observations with |alpha|<=10 degrees. At 18.65 degrees the model gives
  CL=1.751 versus measured 1.005: it misses the abrupt lift loss. The overall
  drag RMSE is 0.001167 (11.7 drag counts) over 16 reported CD values.
* **S809:** 20/22 points converged; -0.01 and 6.16 degrees failed. The model
  gives CL=1.107 at 9.22 degrees versus digitized 1.014 near the report's usable
  maximum, and continues rising. Its sampled peak occurs at 16.24 degrees.
  Ten of fourteen qualified transition/attached-turbulence observations could
  be compared; two lie outside the measured interval plus reading allowance.
  The inferred angle mapping and differing measurement observable limit that
  conclusion. Later lift observations are marked as suspect separated-flow data.
* **NLR7301:** all four requested states around 6 and 13.1 degrees converged.
  CL errors are about +0.139 and -0.184. CD is 0.023220 versus 0.022967 near
  6 degrees, and 0.067405 versus 0.044592 near 13.1 degrees (about 51% high).
  Nominal geometry, uncertain Ncrit, unresolved cove/wake physics and remaining
  discretization error must be investigated. These two incidences do not test
  multi-element stall. The optional full-polar case is included but not executed.

These are observed discrepancies, not a successful validation certificate.
No grid-convergence or Ncrit sensitivity study has yet separated numerical,
model and condition errors for these experimental cases. Source data, hashes,
angle conventions, missing-value handling and extraction audits are included.

## Published MSES comparison added

`nlr7301_mses` adds a code-to-code benchmark against the green MSES CL curve
in Neto et al., ICAS2006 Figure6 (NLR7301 nominal 2.6% gap, Re 2.51e6, Mach 0.185).
All ten points from 0 to 15 degrees converged. CL RMSE against the digitized
MSES curve is 0.07283, MAE 0.06629 and maximum absolute difference 0.12137.
The sampled maxima occur at 12 degrees (MultiFoil, CL 3.0374) and 13 degrees
(published MSES curve, CL approximately 3.093). These are coarse-grid sampled
peaks, not resolved physical stall angles.

See `validation/MSES.md` for source audits and qualifications. The original
MSES deck, exact grid, Ncrit/trips and version are unavailable. The comparison
therefore matches nominal conditions, with explicit assumptions. No MSES
software was run, and this NLR case does not score CD/CM/Cp/transition against MSES.
The reference classifier and deterministic curve extraction add two passing
integrity tests, bringing that suite to eight. The full default benchmark now
contained 54 requested operating points at that stage, of which 52 converged.

## Published MSES drag comparison added

`naca0012_mses_drag` adds eight CD(alpha) values from COBEM2005 Figure3,
NACA0012 at Re6e6 and Mach0.3. All eight points converge. At 0–10 degrees,
CD RMSE is 0.000560 (5.60 counts), within the graph's ±12-count reading allowance.
It does not establish accuracy to a few drag counts. Local supersonic flow is
flagged at12,14,16,18 degrees; the solver cannot resolve shocks there. At16 degrees,
CD is0.034685 versus the graphical MSES0.058507, an error of−238.22 counts.

See `validation/DRAG.md` for the extraction audit, drag definition, limitations,
and full comparison. The aerodynamic solver remains unchanged. The comparison
suite now has10 passing integrity tests; the default five cases request62 flow
points, of which60 converge. The two S809 failures and four unavailable transition
observations remain reported. Neither MSES installation nor a new refinement or
wake-length sensitivity study was used for this addition.

## Automated checks

All **20 tests passed**; see `test_results.txt`. Reproduce with:

```bash
python -m unittest -v test_multifoil.py test_viscous.py
```

The 12 original checks exercise inviscid symmetry, exact Joukowski lift,
refinement, element ordering and interaction, the far-separated limit,
reflection, chord normalization, off-body velocity, collision rejection,
coordinate order, blunt TEs and supplied geometries.

Eight viscous checks exercise symmetric lift/positive drag/residuals; upstream
BL residual and Jacobian agreement; comparison with native mfoil; independent
source velocity/streamfunction consistency; Reynolds/Ncrit sensitivity;
multi-element ordering and the far limit; local-chord/Reynolds scaling; and
failed-force/stall handling.

At a separation of 100 reference chords, reversing the element order changed
CL by 6.89e-8 and CD by 2.43e-10. The far-separated result is within the test's
0.2% bound of twice the single-element result. The large coordinate offsets
reduce the conditioning of the streamfunction system; equality is checked to
a numerical tolerance, not bit for bit.

## Single-element comparison with the upstream solver

Same effective 140-panel NACA 0012 geometry, alpha 2 degrees, Re 1e6,
Mach 0, Ncrit 9. The reference runs the original single-body mfoil Newton
solver, with its sparse allocation compatibility issue bypassed.

| Solver | CL | CD |
|---|---:|---:|
| New coupling | 0.197352588 | 0.005690917 |
| Native mfoil | 0.197340021 | 0.005688067 |

Differences are approximately 0.0064% in CL and 0.050% in CD. Wake tracing and
half-panel velocity weights differ as documented in `THIRD_PARTY.md`.
This comparison checks the extension against the model it uses; it is not
independent experimental validation.

The dense BL residual and Jacobian reproduce the original sparse equations
within the tolerances in `test_viscous.py`. Integrating analytic normal velocity
along a separate closed contour reproduces cross-body source streamfunction
differences to about 1.1e-13 maximum absolute error in the diagnostic run.
That check exposed and removed cross-body branch-cut jumps and missing wake
half-panel velocity terms.

Reducing Ncrit from 9 to 5 moves upper-surface transition forward and increases
CD in the tested NACA 0012 case. Increasing Re from 1e6 to 2e6 reduces its CD.
These are qualitative sensitivity tests, not universal monotonicity claims.

## Single-element separation and sampled lift peak

`examples/stall_single.json`: 200-panel NACA 0012, Re 1e6, Mach 0, Ncrit 9.
The full saved sweep covers 0 to 20 degrees in 2-degree requested steps with
1-degree internal continuation. All requested points converge. Selected rows:

| Alpha [deg] | CL | CD | Converged |
|---:|---:|---:|---|
| 0 | 0.000000 | 0.005321 | yes |
| 4 | 0.416218 | 0.007219 | yes |
| 8 | 0.932384 | 0.012118 | yes |
| 12 | 1.221284 | 0.018838 | yes |
| 14 | 1.316693 | 0.025669 | yes |
| 16 | 1.363331 | 0.040237 | yes |
| 18 | 1.263634 | 0.083291 | yes |
| 20 | 1.039481 | 0.156060 | yes |

The sampled maximum is at 16 degrees, bracketed by converged points at 14 and
18 degrees. This is a model result, not a validated stall angle. At 20 degrees,
the model reports a large separated upper-surface region; the thin-layer
assumption is no longer credible for quantitative post-stall forces.

## Interacting main element and flap

`examples/stall_two_element.json`: main NACA 2412 with 200 panels; NACA 0012
flap with 140 panels, chord 0.30, deflection -20 degrees and translation
[0.93, -0.055]. Re 1e6, Mach 0, Ncrit 9, all coefficients based on reference
chord 1.0. Extended, independent wakes are used.

| Alpha [deg] | CL | CD | Converged |
|---:|---:|---:|---|
| 0 | 1.320992 | 0.025619 | yes |
| 2 | 1.645643 | 0.028102 | yes |
| 4 | 1.899177 | 0.033065 | yes |
| 6 | 2.134336 | 0.038558 | yes |
| 8 | 2.337679 | 0.045729 | yes |
| 10 | 2.532673 | 0.055598 | yes |
| 12 | 2.640042 | 0.072637 | yes |
| 14 | nan | nan | no |
| 16 | 2.054125 | 0.244736 | yes |
| 18 | nan | nan | no |

Automatic stall diagnostic: **not_bracketed**. Failed angles: [14.0, 18.0].
The failure rows remain NaN; later recovery does not turn a failed point into
a stall event. A missing neighbor prevents a peak bracket even if subsequent
converged points show lower lift. Continuation history is retained in metadata.

## Panel refinement of the two-element case

The same coordinates and operating conditions were rerun with **300 main and
210 flap panels** (50% more on each element). This is an initial sensitivity
check, not a grid-convergence proof.

| Alpha | CL, 200/140 | CL, 300/210 | CD, 200/140 | CD, 300/210 |
|---:|---:|---:|---:|---:|
| 0 | 1.320992 | 1.355267 | 0.025619 | 0.024036 |
| 4 | 1.899177 | 1.927720 | 0.033065 | 0.032229 |
| 10 | 2.532673 | 2.574116 | 0.055598 | 0.055446 |
| 12 | 2.640042 | 2.716697 | 0.072637 | 0.073282 |
| 14 | nan | 2.595819 | nan | 0.118035 |

The refined case brackets a sampled maximum at 12 degrees (CL approximately 2.7167) between converged 10- and 14-degree points. The coarser run has no bracket because its 14-degree point failed.

Lift, drag and transition retain material panel sensitivity. Near stall,
convergence can also change with resolution and continuation history. Refine
panels, alpha spacing, wake length, and TE gap before drawing design conclusions.

## Three-element checks

The three-foil stack in `examples/three_foil_stack.json` contains three NACA
0012 contours at y = -0.7, 0 and +0.7, with 100 panels each, Re 1e6.

| Alpha [deg] | CL | CD | Converged |
|---:|---:|---:|---|
| 0 | 0.000000 | 0.018255 | yes |
| 2 | 0.418772 | 0.019028 | yes |
| 4 | 0.850474 | 0.021984 | yes |

The exploratory slat/main/flap case in `examples/three_element_viscous.json`
uses an unrotated 0.18-chord slat with translation [-0.25, 0.16], plus the
same main and flap geometry as above. It starts at -8 degrees and allows
100 Newton iterations per attempted point. Re is 1e6.

| Alpha [deg] | CL | CD | Converged |
|---:|---:|---:|---|
| -8 | 0.154358 | 0.026552 | yes |
| -6 | 0.433899 | 0.027609 | yes |
| -4 | 0.782129 | 0.025025 | yes |
| -2 | 1.068302 | 0.027014 | yes |
| 0 | 1.382233 | 0.031191 | yes |
| 2 | 1.698255 | 0.061466 | yes |
| 4 | 1.988029 | 0.092486 | yes |
| 6 | 2.251606 | 0.127052 | yes |

This demonstrates coupled three-body execution. It does not validate slat
performance, bypass transition, wake confluence or the very large drag values
that may occur in separated branches. The original tightly spaced
`three_element.json` geometry had poor viscous convergence during development;
it is retained as an inviscid example, not presented as a successful viscous run.

## Compressibility check

A NACA 0012 case at alpha 2 degrees, Re 1e6, Ncrit 9 and Mach 0.3 was executed:

| Alpha [deg] | CL | CD | Converged |
|---:|---:|---:|---|
| 2 | 0.209315 | 0.005949 | yes |

The case converged with no local-supersonic flag. No shock model is implemented.

## Retained inviscid verification

For the Joukowski mapping of a circle of radius 1.08 centered at -0.08, the
360-panel inviscid solver at 5 degrees gives CL = 0.588161625, compared with
exact CL = 0.588179807, a 0.003091% error. Pressure integration of CL is less
accurate (about 0.25% error on that grid). These checks apply to the panel
backend, not the viscous stall model.

The prior Viiflow Basic 1.9.4 single-element integration outputs are retained in
`demo_results/viiflow_single`. Its multi-element Pro adapter was not exercised.
The archived inviscid examples are explicitly named `inviscid_two_element` and
`inviscid_three_element`.

## Remaining limits

The NLR7301 experimental comparisons cover two nominal incidences and use secondary
digitizations. Formal near-stall grid convergence and validation of summed
wake drag under strong interaction remain incomplete. Wakes follow the inviscid velocity and do not merge. No dynamic
stall, hysteresis, forced trips, bypass transition, shocks, cavitation or 3-D
physics are modeled. Numerical convergence alone cannot establish validity.
The reported peak, large separated regions, and post-peak force coefficients
must be interpreted within these limits.

## TN808 fixed-lift challenge (8 September 2026)

The [executed TN808 comparison](validation/results/tn808_comparison/TN808_RESULTS.md)
adds 12 setting/angle solutions from 20 retained attempts and 87 requested flow
evaluations. The initial geometry failures exposed a nonmonotone panel-spacing
inverse, now corrected by a conditional monotone fallback. All 46 software tests
passed, including the existing aerodynamic tests and the new regression.
A further CLI check confirms that cached retry failures retain exit code2.

Under the documented reconstruction assumptions, nominal drag is approximately
59–60% below the digitized experimental data. The larger-angle ordering agrees
within the working reference allowance, while the smaller10–20degree comparison
is unresolved experimentally. Initialization failures and about1-count numerical
path sensitivity are retained. These results do not establish drag-optimization
accuracy or a physical stall boundary.

Historical source overlays preserve the exact code hashes of earlier studies.
See `validation/TN808_BENCHMARK.md` before resuming them with changed code.

## S702 MSES and experimental trial (8 September 2026)

The [S702 benchmark](validation/S702_BENCHMARK.md) adds four published flap
configurations at Re=1e6, Mach0.10 and free transition with computational
Ncrit9. Geometry and both reference series are recovered from vector paths in
Somers and Maughmer's NASA report. The
[executed comparison](validation/results/s702_comparison/S702_RESULTS.md)
contains polars, fixed-lift drag, panel sensitivity and retained failures.

The numerical solver, paneling and fixed-lift helper are byte-identical to the
preceding complete archive. New tests check the flap/hinge geometry, printed
lift values and consistency between independently drawn experimental markers.
The source audit regenerates all 32 extraction files identically. Qualification
as a figure reconstruction is retained: original author coordinates, complete
experimental uncertainty, forced-transition comparisons and near-stall testing
remain outside this trial.

All 24 fixed-lift targets converged. Across the 12 base-resolution targets,
mean absolute CD differences are 1.46% from published MSES and 8.84% from
experiment; the largest experimental discrepancy is 30.4% at CL=0.8 with
flap deflection 5 degrees. Both numerical resolutions favor 5 over 10 degrees
at CL=0.8, while the nominal experimental ordering favors 10 degrees. The
largest panel-resolution CD change at common lift is 1.409 counts. The 192
requested comparison evaluations contain one retained coarse-grid failure at
alpha=-4 degrees; the separate initialization probe also converged. All 13
selected geometry/reference integrity tests passed.
