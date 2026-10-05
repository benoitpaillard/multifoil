# Multifoil — viscous multi-element airfoils in Python

For the complete **14 September 2026 Codex handover**, start with
[CODEX_HANDOVER.md](CODEX_HANDOVER.md), then
[DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). They consolidate the executed
studies, remaining drag-ranking discrepancies, reproducibility commands and
the next development priorities. The latest
[flap-angle plots](validation/results/s702_flap_effect/README.md) compare
MultiFoil, published MSES and experiment at both tested Reynolds numbers.

The default backend now couples **all elements' boundary layers and wakes to
one global potential-flow solution**. It includes Reynolds-dependent viscous
drag, natural transition, limited separation, and diagnostics for a sampled
lift peak. The open solver is included; a Viiflow license is not required.

This is a research prototype using the integral boundary-layer models in
Krzysztof Fidkowski's MIT-licensed [mfoil](https://websites.umich.edu/~kfid/codes.html).
It is not a validated replacement for XFOIL, MSES, or RANS. In particular,
converged results in massive separation are not reliable post-stall predictions.

## Run

Use Python 3.10+; tested on Python 3.12 with NumPy, SciPy, and Matplotlib.
From the unzipped `multifoil` directory:

```bash
python -m pip install -r requirements.txt
python multifoil.py examples/two_element.json --out results
```

For a bounded **incompressible contour-shape search** at fixed lift:

```bash
python shape_optimization.py examples/s702_incompressible_design.json \
  --out design_run --workers 4
```

This changes the fore-element camber using smooth modes, evaluates two lift
targets and two natural-transition assumptions, and checks a selected change
on a finer grid. It enforces geometry, moment and off-design drag restrictions.
Read [INCOMPRESSIBLE_DESIGN.md](validation/INCOMPRESSIBLE_DESIGN.md) for the
objective, acceptance rules and limits. `summary.json` distinguishes a changed
geometry that passes numerical screening from retention of the baseline.
The aerodynamic closures remain unchanged. A separate S702 study calibrates
the effective Ncrit assumption against baseline drag data; a numerical
screening result is not independent validation of shape improvement.

A sweep intended to investigate separation and a lift peak:

```bash
python multifoil.py examples/refined_two_element.json --out results_stall
python multifoil.py examples/stall_single.json --out results_single
```

Set operating conditions directly:

```bash
python multifoil.py examples/two_element.json --re 3000000 --ncrit 5 \
  --alpha 0 2 4 6 8 --out results_re3m
```

The command returns **0** if every requested point converged, **2** if some
failed (results are still saved), and **1** for an input/setup error. Output
directories must be new or empty. Failed force coefficients are NaN in CSV and
Python, and null in JSON. A failed solve never receives an extrapolated CL or CD.

## What the upgrade adds

| Quantity or behavior | `viscous` backend |
|---|---|
| Element interaction | Global potential-flow and displacement-source coupling |
| Viscosity | Integral momentum/energy equations, laminar and turbulent closures |
| Transition | Simplified e^N amplification envelope; adjustable Ncrit |
| Turbulent development | mfoil shear-lag closure |
| Separation | Negative skin friction and separated intervals on each surface |
| Lift and moment | Integration of viscously modified pressure |
| Profile drag | Sum of Squire–Young estimates from individual wakes |
| Stall diagnostic | Sampled lift-peak brackets using converged neighboring points |
| Compressibility | Kármán–Tsien correction; subsonic, shock-free flow only |
| Continuation | Reuses converged states, subdivides failed angle steps |

There is no imposed CL limit or empirical post-stall curve. A lift roll-off is
an outcome of the coupled equations. `stall.json` distinguishes a sampled
maximum from a peak bracketed on both sides. A peak remains a **candidate**:
repeat with finer panels and angle spacing before interpreting its location.
A failed point, or the last successful point before failure, is not a stall angle.

Natural transition is modeled on both sides of every element. Lower Ncrit
usually moves transition forward; it is an amplification threshold, not a
specified transition position. Forced trips, surface roughness, crossflow
transition, and wake-induced bypass transition are not implemented.

## Output

| File | Contents |
|---|---|
| `polar.csv` | Alpha, convergence, CL, CD, CM, residual |
| `surface_000_main.csv`, etc. | Actual nodes, Cp, Cf, theta, delta-star, Hk, transition flag, edge speed |
| `wake_000_main.csv`, etc. | Wake coordinates, momentum/displacement thickness and edge speed |
| `transition_separation.csv` | Transition x/c, separated arclength fraction, minimum Cf, maximum Hk |
| `stall.json` | Peak candidates, brackets, sampled maximum, failed angles |
| `metadata.json` | Effective settings, convergence histories, per-element forces, separation intervals |
| `geometry_main.dat`, etc. | Actual normalized, repaneled and TE-regularized geometry |
| `pressure.png`, `polar.png`, `boundary_layers.png` | Scientific plots |

Pressure is plotted at the middle converged requested point. The boundary-layer
plot uses the largest converged requested angle. Surface and wake CSVs contain
every converged requested point. Internal continuation points are recorded in
metadata; they do not replace requested points in the polar.

Cf means wall shear stress divided by freestream dynamic pressure. Negative Cf
indicates reversed skin friction within the integral model. Thicknesses are in
reference-chord units; Hk is the kinematic shape factor. The Python surface key
`H` contains Hk for compatibility with the original adapter. The third state,
`amplification_or_shear`, is N in laminar flow and the square root of the modeled
shear-stress coefficient in turbulent flow; interpret it using `turbulent`.
A blank transition coordinate means no transition before the trailing edge.

## Geometry and reference quantities

Use one Selig-format `.dat` contour per element: optional title, then x/y columns,
from upper trailing edge around the leading edge to lower trailing edge. Reversed
order is accepted. A duplicated sharp-TE endpoint is accepted. Do not add a
separate closing edge across a blunt TE. Lednicer point-count headers need
conversion first. Alternatively specify a four-digit NACA code.

```json
{
  "reference_chord": 1.0,
  "moment_reference": [0.25, 0.0],
  "reynolds": 1000000,
  "mach": 0.0,
  "ncrit": 9,
  "alpha_deg": [0, 2, 4, 6],
  "elements": [
    {"name": "main", "file": "main.dat", "panels": 200},
    {"name": "flap", "file": "flap.dat", "panels": 160,
     "scale": 0.30, "angle_deg": -20,
     "pivot": [0, 0], "translation": [0.93, -0.055]}
  ]
}
```

File paths are relative to the case JSON. The example assumes unit-chord input
profiles. For coordinates already assembled, leave scale at 1, rotation at 0,
and translation at [0, 0]. Transform order is scale, rotate about the scaled
pivot, translate. Positive rotation is counterclockwise; a downward flap uses
a negative angle. Coordinates are x downstream and y upward.

All elements share one reference chord. Reynolds number is U-infinity times
reference chord divided by kinematic viscosity. CL and CD use that chord; CM
uses its square and is nose-up positive about the stated reference point.
The local Reynolds number follows automatically from each element's size.
No element is silently rescaled to unit chord. Freestream direction is
(cos(alpha), sin(alpha)); positive alpha gives positive lift on a symmetric foil.

The viscous solver uses curvature-based repaneling in each element's local
frame. A sharp TE is smoothly opened to a default gap of **1e-5 local chords**
to avoid the singular viscous source limit. The added gap is recorded in metadata
and the effective geometry is exported. Test sensitivity using `minimum_te_gap`.
Coarse source coordinates and poorly resolved slots cannot be repaired by the solver.

## Solver controls

These optional JSON settings are also accepted as constructor arguments:

| Setting | Default | Meaning |
|---|---:|---|
| `iterations` | 80 | Maximum Newton iterations per attempted point |
| `tolerance` | 1e-7 | Both BL residual L2 norm and maximum coupling residual must pass |
| `continuation_step` | 1.0 | Maximum internal alpha step, degrees |
| `min_step` | 0.25 | Smallest retry step, degrees |
| `start_alpha` | null | Explicit initial angle; otherwise start near zero |
| `wake_length` | 1.0 | Approximate downstream extension past aftmost TE, in reference chords |
| `minimum_te_gap` | 1e-5 | Minimum TE gap divided by local chord |
| `amplification_jacobian` | `corrected_af` | Corrected amplification derivative; `legacy` replays the historical derivative |

The standard viscous solver includes the missing factor of two in the
amplification Jacobian. The amplification growth value is unchanged. Select
`amplification_jacobian="legacy"` in Python, `"amplification_jacobian": "legacy"`
in case/design JSON, or `--amplification-jacobian legacy` on the main CLI for
comparisons using the historical derivative. Each solution records the selected
mode in its diagnostics, including failed solves. The vendored source is unchanged.
Existing studies retain their source hashes; replay archived studies with their
frozen source snapshots instead of bypassing manifest guards after this change.

For difficult configurations, begin at an easier angle and sweep upward.
`--start-alpha`, `--continuation-step`, `--iterations`, and `--wake-length`
override the corresponding settings. `--verbose` prints Newton residuals.
The last converged state is restored after a failed continuation attempt.
A gentle velocity profile is used as an initialization fallback if the one-way
boundary-layer march breaks down; the final residuals always use the exact
specified operating conditions and full element interaction.

Wakes follow the **global inviscid** velocity field and are rebuilt at each
attempted alpha. They are not relaxed with viscous displacement. Every wake is
extended approximately beyond the aftmost element; its number of nodes stays
fixed during continuation. Increase wake length to assess drag sensitivity.
Wake/body intersections are rejected. Wake merging and boundary-layer/wake
confluence are not modeled, even if the geometric curves do not intersect.

Dense work scales approximately as the cube of the total station count. The
limits are 1200 surface nodes and 1800 surface-plus-wake nodes. Start with about
200 panels on the main element and 140–200 on smaller elements. Narrow gaps,
transition bubbles, and stall investigations require refinement. Python speed
is much lower than compiled XFOIL. The included viscous sweeps took about one
to ten minutes in the test environment; difficult points can take much longer.

## Python API

```python
from multifoil import load_case, save_results
from viscous import CoupledSolver, analyze_stall

cfg, elements, reference = load_case("examples/two_element.json")
solver = CoupledSolver(elements, reference, reynolds=1e6, ncrit=9)
results = [solver.solve(a) for a in range(0, 17)]

for r in results:
    print(r.alpha, r.converged, r.cl, r.cd)
print(analyze_stall(results))
# Export solver.elements: these include the actual viscous repaneling/TE gap.
save_results("api_results", results, solver.elements, cfg, "viscous")
```

Rebuild the solver to change geometry, Reynolds number, Mach number, or Ncrit.
`Solution` contains forces, surface arrays and detailed diagnostics. Check
`converged` before using coefficients; also check model applicability.

## Other backends and examples

The original inviscid backend remains available:

```bash
python multifoil.py examples/two_element.json --backend panel --flow --out inviscid
```

It uses linear vortex panels for sharp TEs and Hess–Smith panels when a blunt TE
is present. It has no Reynolds effects, viscous drag, transition or stall. Its
pressure drag is a discretization diagnostic, not a usable profile CD.

The original `--backend viiflow` adapter is retained. It requires a separate
compatible installation; the tested Basic 1.9.4 wheel permits one element,
whereas multi-element use needs a suitable Pro license. No proprietary binary
is bundled and the Pro integration has not been tested here.

`three_foil_stack.json` exercises three interacting viscous foils.
`three_element_viscous.json` is an exploratory slat/main/flap case with increased
slat clearance and continuation from -8 degrees. `three_element.json` retains
the original tightly spaced inviscid geometry; its viscous convergence is poor.
The examples above are synthetic. The separate `validation/` directory now
contains published NACA 0012, S809 and NLR7301 benchmarks, with traced reference
data and executed comparisons. Included outputs retain failed rows.

## Published validation cases

```bash
python validate.py --out validation_run
python validate.py --compare-only --out validation/results
```

See `validation/README.md` for sources, conditions, digitization qualifications,
refinement controls and individual cases. The new `nlr7301_mses` case compares
CL directly against published MSES results; see `validation/MSES.md`. Run it
with `python validate.py --case nlr7301_mses --out mses_benchmark`.
The `naca0012_mses_drag` case adds a published MSES CD(alpha) comparison at
Re=6 million and Mach=0.3. Run it with
`python validate.py --case naca0012_mses_drag --out mses_drag_benchmark`.
See `validation/DRAG.md` for drag definitions, results and the remaining
multi-element drag discrepancy.
Run `python validation/plot_three_way.py` to reproduce the combined MultiFoil,
published MSES and experimental overlays. See `validation/THREE_WAY.md` for
the figure and data sources, including the newly added single official NLR MSES
drag point (exact input geometry unconfirmed).
For published evidence and a proposed path to trustworthy shape optimization,
see `validation/OPTIMIZATION_READINESS.md`. It also records a newly located
official MSES NLR7301 example with a printed drag value at one incidence.

The first fixed-lift design-screening pilot is described in
[`validation/OPTIMIZATION_PILOT.md`](validation/OPTIMIZATION_PILOT.md), with
executed results in [`validation/results/optimization_pilot/REPORT.md`](validation/results/optimization_pilot/REPORT.md).
It tests smooth camber changes and flap placement against panel, wake,
transition-parameter and initialization sensitivity. Independent geometry
rankings remain unvalidated; it is not a continuous shape optimizer.
The archived comparisons show
material transition, drag and stall errors; convergence does not establish
experimental accuracy. Read `validation/results/REPORT.md` before using results.

## Verification and limitations

```bash
python -m unittest -v test_multifoil.py test_viscous.py
```

See `VALIDATION.md` for executed results and their scope. Checks include an
exact inviscid solution, upstream single-element equations, independent source
velocity/streamfunction consistency, Reynolds/Ncrit effects, scaling,
multiple-element interaction, ordering, and failure handling.

The equations are steady, two-dimensional, thin-layer approximations. They do
not resolve massive separation, unsteady or dynamic stall, hysteresis, merging
wakes, shocks, 3-D induced drag, or cavitation. Local supersonic flow is flagged
when using nonzero Mach. A small Newton residual verifies the equations were
solved; it does not validate the flow model or establish experimental accuracy.

See `THIRD_PARTY.md` for source provenance and licensing. The mfoil boundary-layer
closures are reused unchanged. Cross-body source branches, global interaction,
continuation, diagnostics, and corrected wake half-panel velocity weights are
implemented in `viscous.py`; the upstream self-surface near-TE approximation is
retained. This extension has not been reviewed or endorsed by mfoil's author.

## Independent design-ranking checks

See [`validation/INDEPENDENT_RANKING.md`](validation/INDEPENDENT_RANKING.md) for the
NACA TR664 source audit and uncertainty-aware fixed-lift comparison command.
The candidate is explicitly blocked pending geometry, setting and uncertainty
audits; no new independent ranking or aerodynamic solve is claimed.

## TN808 experimental challenge and paneling fix

The [TN808 benchmark](validation/TN808_BENCHMARK.md) reconstructs the published
slotted-flap geometry and compares drag at common lift across flap angles.
It also reproduces and fixes a panel-spacing fold at the slot entrance.
The fix preserves valid upstream paneling and records when monotone recovery
was used. The challenge exposes large drag discrepancies and numerical path
sensitivity; it does not establish drag-optimization accuracy.

## S702: MultiFoil, published MSES and experiment

The [transition sweep](validation/results/s702_ncrit_analysis/NCRIT_SWEEP_RESULTS.md)
compares 14 Ncrit values from 3 to 12, including fractional refinement, across
four flaps and three lift targets. The lowest complete relative drag error is
4.85% at Ncrit=6.25, versus 8.84% at 9. Absolute drag and flap-difference RMSE
both prefer Ncrit=6. Finer panels retain the broad minimum, with errors of
4.75% at 6.25 and 4.81% at 6. This supports an effective calibration region,
not a measured disturbance threshold or validated contour-design gradients.
The report includes grouped flap checks, a reserved Re=700,000 comparison,
pressure checks and Mach-zero sensitivity. All three failed calibration
targets remain visible; incomplete candidates cannot win by losing cases.

At the reserved Re=700,000 conditions, seven of nine targets converged under
the original incidence protocol. Separate local/warm initialization retries
recovered the other two with the frozen Ncrit=6.25, reducing nine-point mean
relative drag error from 10.15% at Ncrit=9 to 5.70%. The worst point still has
24.04% error. The driver records the original failures and uses computed lift
to guide the retries; it does not use experimental drag to choose incidence.

Reproduce the full study, including frozen-parameter validation, from a fresh
directory inside the project:

```bash
python validation/run_s702_calibration.py --out validation/results/reproduced_s702_calibration --workers 8
```

The full study can take hours; completed results are included in this archive.
See [NCRIT_SWEEP_PROTOCOL.md](validation/NCRIT_SWEEP_PROTOCOL.md) for individual
commands, scoring rules and continuation provenance. The solver default stays
at Ncrit=9; use `--ncrit 6.25` to apply the selected S702 assumption explicitly.

The [discrepancy investigation](validation/results/s702_discrepancy/DISCREPANCY_RESULTS.md)
adds controlled per-element transition tests, a 12-point flap/lift transfer
check, finer panels, Mach-zero checks, and independently extracted experimental
pressure data. Earlier fore-element transition removes the 5°/10° ranking
reversal at CL=0.8, but a blanket Ncrit=6 setting makes some other points worse.
This is a partial diagnostic remedy, not a validated replacement of the
transition physics. The default remains Ncrit=9. Reproduce the completed report
with `python validation/report_s702_discrepancy.py`; full solve commands and
qualifications are in the report.

The [S702 benchmark](validation/S702_BENCHMARK.md) extracts geometry and
published comparison curves directly from the vector NASA report. It includes
four small-flap settings, matched free-transition conditions, fixed-lift drag
comparisons and a second panel resolution. The
[executed results](validation/results/s702_comparison/S702_RESULTS.md) distinguish
agreement with MSES from experimental configuration rankings. The original
author geometry deck remains unavailable; the distributed coordinates are a
documented figure reconstruction. The numerical flow solver is unchanged.
