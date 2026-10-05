# Incompressible contour design with MultiFoil

This upgrade adds a working bounded shape search at **Mach zero**. It corrects
the design formulation: compare drag at prescribed lift, keep the assembly
reference chord fixed, use smooth shape modes, and check paired changes under
different transition assumptions. It does not alter the boundary-layer closures
or fit drag to S702 experimental data. Existing S702 Mach-0.1 results remain
unchanged and separate from this incompressible study.

## Run the supplied example

From the project directory, with the normal requirements installed:

```bash
python shape_optimization.py examples/s702_incompressible_design.json \
  --out validation/results/s702_incompressible_design --workers 4
```

The included completed directory can be read or resumed with identical code and
configuration. Use a new output directory after changing the configuration,
geometry or source code. Hashes prevent accidental reuse of old evaluations.
The run may take many minutes. Workers parallelize independent operating points.

The example changes the **fore-element contour**, keeping the S702 small flap
at its published 5-degree setting. Its source coordinates come from vector
Figure 1 of NASA/CR–20220014964, not an original author coordinate deck. See
[S702_BENCHMARK.md](S702_BENCHMARK.md) for provenance and qualifications.

Every objective evaluation uses Mach=0, Re=1,000,000 and fixed lift at CL=0.5
and 0.8, with equal weights. Ncrit=9 is the nominal natural-transition model;
Ncrit=6 is an explicitly chosen sensitivity scenario. It is neither a forced
trip nor a calibrated representation of this tunnel. Each shape, lift target
and transition scenario receives a fresh solver with the same declared
incidence bounds and starting angle; state is reused only within that
fixed-lift solve. This makes evaluation order between geometries irrelevant.

## Deformation and geometry restrictions

Each element uses a chord frame derived once from its **baseline** coordinates.
The midpoint of its two trailing-edge ends and the source vertex farthest from
that midpoint define the local chord. This is a deformation convention for
conventional airfoil contours; it does not change their actual incidence or
position. The assembly reference remains one for all elements and shapes.

The camber mode adds the following displacement along the frozen chord normal:

\[
\Delta y = 64\,a\,c_e\,t^3(1-t)^3,\qquad
t=\mathrm{clip}\left(\frac{x/c_e-x_0}{x_1-x_0},0,1\right).
\]

Here a is the maximum displacement as a fraction of the local element chord.
Value, first derivative and second derivative vanish at the support ends.
Both surfaces receive the same displacement at a given chordwise coordinate,
preserving the continuous thickness distribution and section area. In the
example, support is 5–95% of the fore-element chord, amplitude bounds are
±0.002, and the initial search step is 0.001. The leading edge, trailing edge
and placement stay fixed. Additional broad camber modes may be specified on
any element; no adjoint or analytic aerodynamic shape gradient is claimed.

Input polygons, outer repaneling and actual solver surfaces are checked for
intersection and element clearance. Raw and actual surface area checks also
limit discretization drift: evaluating a smooth deformation on a finite
polygon does not preserve its polygonal area exactly. Coordinates use fractions
of the common assembly chord; no element is independently normalized.

The example's engineering restrictions are a minimum clearance of 0.015
assembly chords, relative polygon-area tolerance 0.002, moment change at each
operating point no greater than 0.015, and separated arclength fraction on each
surface no greater than 0.25. These are declared **design restrictions**, not
validated boundaries of the aerodynamic model. In particular, convergence and
passing a separated-length restriction do not establish accurate bubble drag.

## Objective, acceptance and verification

For each transition scenario s, compute a weighted mean of the drag change
relative to that scenario's own baseline at identical lift:

\[
D_s(a)=\frac{\sum_j w_j[C_D(a,C_{L,j},s)-C_D(0,C_{L,j},s)]}{\sum_jw_j},
\qquad J(a)=\max_s D_s(a).
\]

This prevents an improvement under one transition assumption from masking a
worsening under another. In addition, no individual example operating point
may worsen by more than one drag count. One drag count is 0.0001 in CD.
Failed flow solutions, unmatched lift, missing points and moment violations
cannot become an objective. All requested flow evaluations are retained.

Bounded coordinate polling tries positive and negative mode changes, accepts
an objective improvement larger than 0.5 counts, and halves the step after an
unsuccessful complete poll. The example allows five unique designs and a
minimum step of half the initial step. A budget stop is explicitly reported;
it is not convergence to an optimum. Nonmonotone lift or exhausted fixed-lift
iterations reject the point. No finite penalty is substituted for a failed
flow solution.

A changed candidate is recomputed on 340/240 panels after the 260/180-panel
search. Every fine candidate is compared with its **fine-grid baseline**.
Acceptance for numerical screening requires:

- All fine-grid points satisfy the same flow and design restrictions.
- Each transition scenario's mean gain exceeds the 0.5-count threshold.
- The gain on each grid exceeds the observed difference between the paired
  base- and fine-grid gains, plus 0.5 counts.

This two-grid check is an observed numerical sensitivity check, not a physical
error bound. `selected_geometry/` receives the changed shape only if this check
passes; otherwise it receives the baseline. Both rejected and accepted trial
coordinates remain in their individual folders. A numerical-screen pass does
not validate the predicted improvement experimentally.

## Outputs

An optional `initialization` object on each lift point enables the audited
local-bracket policy. Existing configurations retain the legacy fixed bracket.
For example, a point with `alpha_bounds: [-6, 8]` can declare:

```json
"initialization": {
  "center": 1.3689680656734002,
  "half_width": 0.05,
  "max_brackets": 2,
  "max_evaluations": 20,
  "paths": [
    {"name": "cold_local", "start_alpha": 1.3689680656734002},
    {"name": "warm_anchor", "start_alpha": -1}
  ]
}
```

Choose and freeze numerical seeds for the declared operating conditions. Every
start and optional `seed_alphas` entry must lie inside `alpha_bounds`. With the
option enabled, each path's `start_alpha` replaces the legacy point start; the
first successful path wins. The same ordered policy applies independently to
every geometry, including its baseline. One bounded recenter may use computed
lift; drag does not choose the seed. No state is shared between geometries.

`records` and flow metadata retain failed trials and initialization context.
`selected_flow_index` identifies the returned flow, and policy runs save the
restored selected Newton arrays as `<element>_state.npz`. These arrays are audit
artifacts, not portable solver restart files. Changing code or policy invalidates
the study cache and requires a new output directory.

The Mach-0.07 policy audit passed all seven cases; see
[its findings](results/s702_initialization_policy_analysis/FINDINGS.md).
The bounded [Mach-zero integration protocol](DESIGN_INITIALIZATION_PROTOCOL.md)
must still be assessed before extending the policy to a broader design campaign.

| File or folder | Meaning |
|---|---|
| `study.json` | Conditions, constraints, source and input hashes |
| `design_NNN/fore.dat`, `aft.dat` | Exact input contours of every trial |
| `design_NNN/result.json` | Paired objective, feasibility and geometry measures |
| `design_NNN/scenario__point/flow/` | Every requested solve, surfaces, wakes and transition diagnostics |
| `search.json` | Search history, selected trial and stopping reason |
| `verification_baseline/`, `verification_candidate/` | Independent finer-grid evaluations |
| `verification.json` | Explicit acceptance decision and differential grid sensitivity |
| `selected_geometry/` | Numerically screened geometry, or retained baseline |
| `summary.json` | Machine-readable outcome, with physical validation explicitly false |

The CLI returns 0 for a completed search, including a completed search that
retains the baseline; 2 for an infeasible baseline; and 1 for an uncaught setup
error. Check `summary.json` for whether a shape change was accepted.

## Checks accompanying this upgrade

```bash
python -m unittest -v test_shape_optimization.py test_optimization.py test_paneling.py
python validation/s702_mach_check.py --workers 2
python validation/incompressible_shape_check.py
```

The Mach study independently repeats the Mach-0.1 S702 5/10-degree fixed-lift
calculations and recomputes them at Mach zero. The latter are sensitivity
calculations and must not be relabeled as matched-condition experimental
validation. The single-element study checks three NACA0012 camber amplitudes
at fixed lift against the bundled native mfoil Newton/lift-constraint solver,
using identical actual surfaces. That verifies implementation consistency with
shared closures; it is not an independent physics model or MSES run.
Native Python mfoil's fixed-lift path needed two array-layout adapters, for its
edge-velocity basis and scalar lift residual, in addition to the existing sparse
reallocation adapter. The former was checked against native edge velocities and
their incidence derivative. These are confined to the comparison runner; the
vendored source is unchanged. Both failed unadapted attempts are preserved.

## Scientific basis and remaining limits

[Drela's LINDOP manual, section 14](https://web.mit.edu/drela/Public/web/mses/lindop.pdf)
motivates multiple operating points and limited shape freedom. It warns that
free-transition optimization can exploit fragile bubble behavior. The current
two-Ncrit calculation addresses sensitivity within the available model; forced
transition remains unimplemented.

[Fidkowski's mfoil paper](https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf)
documents the coupled potential-flow and integral-boundary-layer foundation.
Mach zero uses its incompressible relations directly, including wake momentum
drag; an incompressible lift calculation alone would not provide viscous drag.

The [S702 NASA report](https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf)
describes nearly incompressible tunnel conditions and differences in predicted
transition, bubbles and the low-drag range. Hinge-film disturbances may partly
explain pressure-distribution differences. These observations do not supply
a universal correction for drag or shape derivatives. Exact geometry, trips,
roughness, interacting/merging wakes, and a validated shape-change dataset
remain needed before using these predictions to certify a design improvement.
