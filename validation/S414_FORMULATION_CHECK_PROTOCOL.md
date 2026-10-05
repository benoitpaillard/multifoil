# S414 inviscid cross-formulation check (V3), 4 October 2026

Stage V3 of the verification plan in `DEVELOPMENT_PLAN.md`. Frozen before
execution; the script is `validation/report_s414_formulation_check.py`.
Linear analysis only: no viscous or nonlinear solve.

## Question

The model's excess nose suction already appears in its inviscid calculation.
Is that pressure specific to the production streamfunction/linear-vortex
formulation, or does an independently formulated potential-flow method give the
same pressure on the **identical final polygons**? Agreement weakens a
formulation-specific explanation. It does not validate the pressure against the
viscous experiment, and it says nothing about geometry.

## Formulations

- **Production**: the unchanged `CoupledSolver._build_inviscid` operator
  (linear vortex sheets, streamfunction boundary condition at nodes, finite
  trailing-edge panel with source and vortex, nodal Kutta condition), assembled
  on prescribed nodes without repaneling.
- **Hess–Smith**: the existing `multifoil.HessSmithSolver` (constant sources,
  one constant vortex per body, no-penetration at panel midpoints, base panel
  retained, equal tangential speed on the two trailing-edge panels).

The sharp-only linear-vortex backend is not used and no trailing edge is
sharpened. Neither method is assumed correct because its matrix solves.

## Qualification controls (before any S414 comparison)

1. Kármán–Trefftz airfoil (10° trailing-edge angle, cambered), exact surface
   speed and lift, at 120, 240 and 480 panels, α = 4°.
2. Symmetric finite-trailing-edge NACA 0012 at 0° and ±4°: zero lift, mirrored
   speeds.
3. Two finite-trailing-edge elements in both orders at 4°: identical speeds.
4. For every solve: linear residual, Kutta residual, closed-body source flux
   (Hess–Smith), and far-contour circulation and flux from each method's own
   velocity field (circle of radius 3 about (0.5, 0), 720 points) compared with
   its pressure-integrated lift.

Qualification screens, fixed now: linear and Kutta residuals below 1e−8;
symmetry, antisymmetry and ordering differences below 1e−8; Kármán–Trefftz
errors decreasing monotonically with resolution for both methods, with
finest-grid lift error below 0.5 % and RMS speed error (outside the last 2 % of
chord) below 0.005. A failed screen stops the S414 interpretation.

## S414 comparison

Geometries: `nominal` and `fore_source_cubic` final solver polygons at 340/240
panels, checked against archived saved-state nodes. Incidences, read from
archived records at full precision: the measured −3.07° and −1.04°; the
matched-lift CL 0.056, 0.175 and 0.284 nominal angles; both transition-boundary
bracket ends; and the cubic CL 0.056 angle (cubic geometry only).

Compared quantities: incompressible signed speed at **panel midpoints** for
both methods (the production nodal values are averaged on each panel; nodal and
midpoint samples are never mixed); then Cp with the same Kármán–Tsien
correction at Mach 0.10; stagnation position; fore-nose minimum Cp and its
position; Cp at the supported tap positions (printed-table coordinates, with
plotted coordinates as a sensitivity); the four surface loads; total lift from
pressure and from circulation.

## Three-resolution study

Schedule, fixed now: 340/240 (baseline), 480/340, 680/480 panels. The
continuous contour is held fixed: every grid samples the spline of the
unchanged baseline solver input with the production curvature spacing, so added
nodes carry no new geometry information. The baseline grid must reproduce the
production nodes exactly. There are no wakes in this linear analysis. Results
at each resolution are compared at the baseline midpoint arclength positions.

Reported: for each method the change of nose Cp, loads, lift, flux and
stagnation position with resolution; and the between-method difference at each
resolution. The formulations are said to agree within resolved discretization
spread only if their fore-nose Cp difference decreases under refinement and, at
the baseline, does not exceed the larger of the two methods' own
baseline-to-finest changes. Otherwise the gate stays open. An algebraic
residual is not spatial convergence. There is no fourth grid.
