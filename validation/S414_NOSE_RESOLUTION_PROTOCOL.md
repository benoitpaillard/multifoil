# S414 local nose-resolution pilot — 5 October 2026

Question: does the robust early fore-lower transition/drag excess survive local
streamwise refinement of the **same final polygon and legacy station metric**?
This addresses an open coupled-discretization question after V0–V5. It does not
calibrate geometry, change Ncrit, validate the manufactured contour, or replace
the separated-layer closure. No direct prescribed-pressure replay is admitted.

## Frozen six-case scope

Two reference states: the audited corrected-Jacobian jobs
`validation/results/s414_jacobian_20261003/cl0056_corrected_af` and
`cl028_corrected_af`. Use the exact angle in each `result.json`, approximately
−3.459015324304847° and −1.6724273335571296°. At each angle compare subdivision
factors 1, 2 and 4. The labels CL 0.056/0.284 identify the archived states;
**refined cases are at identical incidence, not necessarily identical lift**.
Do not report their drag difference as a matched-lift accuracy improvement.

Retain Re 1,000,000, Mach 0.10, Ncrit 9 on both elements, explicit corrected_af,
unit assembly reference chord (457.2 mm test model), nominal 340/240 input
panels, finite TE treatment, 1c wake, maximum 100 Newton iterations per attempted
point and residual tolerance 1e−7. Preserve the historical initialization:
replay the five recorded angles through selected index 4 with the archived
start angle and continuation controls. There are six final cases and at most
30 requested evaluations including 24 initialization evaluations; internal
continuation attempts remain separately recorded. No adaptive lift search,
alternative initialization or retry is authorized by resuming this study.

Both factor-1 cases run first and must match their archived selected state
arrays exactly under pinned Python 3.12.13, NumPy 2.3.5, SciPy 1.17.0 and
Matplotlib 3.10.8. Failure stops the four refined cases. Record mode, paths,
incidence, forces and all 18 state arrays. No historical manifest is rewritten.

## Geometry and wake intervention

Build the production nominal geometry once per case. In clockwise fore-node
order (lower TE to LE to upper TE), select parent panels before the geometric
LE whose **both endpoints** have total-chord x ≤ 0.04. The LE is the farthest
node from the TE midpoint, using the unchanged production definition. Freeze
the exact panel indices and generated node arrays in preparation artifacts.

Insert one or three collinear equally spaced points per selected parent panel.
Retain every original coordinate, `foil.s` and tangent exactly, with an explicit
old-to-new node map. Insert `foil.s` by linear interpolation of the original
values. At inserted nodes the tangent is the parent Δx/Δs; original tangents
stay unchanged. Do not renormalize old stations or repanel afterwards. The
original spline arclength differs slightly from polygon length; preserve and
quantify that convention rather than silently changing it with resolution.
This experiment tests the frozen polygon/metric, not smooth-contour convergence.

All original fore nodes, the complete aft polygon/sampling, element placement,
TE faces and geometry frame remain unchanged. Factor 1 follows the original
production operator path. Refined cases use the same influence kernels and
closures on the new node lists, with source/system slices rebuilt consistently.

Freeze each wake's original node count, first arclength spacing and nominal
extension; otherwise increasing fore nodes would also increase wake resolution.
Wakes still follow the recomputed inviscid velocity. Save and compare their
coordinates/tangents. A maximum displacement over 1e−5c relative to baseline is
a declared confounding flag: report the calculation, but do not claim a pure
nose-spacing attribution without a separately declared wake-isolation check.
No automatic follow-up or fourth resolution follows a failed screen.

## Verification and scientific screens

Before launch, test identity/subdivision, retained vertices/metric, nested
grids, unchanged aft/TE/wake controls, node/source indexing, archived control
selection, failures carrying null forces, manifest/resume rules and orchestration
recovery using fake commands. Independent review of the frozen protocol and
runner must pass before remote submission. Freeze all numerical source, test,
protocol, preparation and analysis dependencies in an isolated snapshot.
The linear preflight also builds actual wakes and source coupling at the
archived low-lift angle on all three grids. Check contiguous surface/system/source
blocks, complete dimensions, finite operators and the mass-to-source action
against explicit differences on every panel. Factor 1 must exactly match the
standard solver's source and mass operators. This does not initialize or solve
the nonlinear boundary layers. Save the coordinator recovery-test log as hashed
preparation evidence; the remote machine need only run the numerical tests.

After retrieval, verify transfer/source hashes and freshly rebuild every saved
converged state: geometry, wakes, station indexing, BL L2 residual, coupling
maximum residual, forces, Cp/Cf, transition and separation, and wake drag. Both
residual norms must remain ≤1e−7. Preserve failed attempts and missing final
cases; no accuracy screen can pass by omitting them.

Compare Cp, theta and displacement at all finest-grid fore-lower physical
stations within x≤0.04 supported by all three grids, using their shared station
parameter to interpolate without extrapolation. Relative thickness errors use
the finer value with a 1e−12 denominator floor. N is evaluated only where the
interpolation-bracketing stations in all grids are on the same lower branch and
laminar. Report coverage; never interpolate turbulent shear as amplification.
Also report full surface distributions, extrema and fore/aft forces.

| Metric, half-to-quarter spacing | Numerical screen | Noise floor for contraction |
|---|---:|---:|
| Maximum absolute Cp change | 0.01 | 1e−5 |
| Maximum relative theta change | 0.02 | 1e−4 |
| Maximum relative displacement change | 0.02 | 1e−4 |
| Maximum absolute N change on common laminar domain | 0.1 | 1e−3 |
| Absolute total CD change, counts | 1 | 0.01 |

For each metric the half-to-quarter difference must be no larger than the
baseline-to-half difference when the latter exceeds its noise floor. Otherwise
the half-to-quarter difference must remain below that floor. These are
predeclared operational screens relative to the ~55-count discrepancy, not
experimental uncertainties or proof of asymptotic convergence.

Require matching transition presence and ordered separation/TE-reaching topology
on all four surfaces across both adjacent grid comparisons (1→2 and 2→4).
Screen transition movement for the finest comparison (2→4): it must be smaller
than one enclosing finest-grid cell, using physical polygon arclength (report
that cell width). Report 1→2 transition movement for context; it does not have
to meet the finest-grid cell limit. Large coarse-grid movement alone does not
fail a comparison whose two finest grids agree, provided topology matches.
Missing/inadequately supported metrics are unknown, never a pass. A numerical
screen failure is a valid experimental finding if execution and evidence audit
pass. Record execution, evidence verification, scientific screens and physical
interpretation separately; a reviewer cannot turn one into another.

## Orchestration and resource authorization

The user authorized this pilot and its independent review workflow. Use isolated
source snapshots: the project is not currently a Git repository. Do not edit
existing Python modules or historical study outputs. New orchestration and
diagnostic files are allowed; speculative physical corrections are not.

At most two single-threaded workers on jony4, BLAS/OpenMP/MKL/NumExpr threads=1.
Maximum 1800 seconds per case and 5400 seconds for the entire detached remote
launch, including preflight and timeout cleanup. An exceeded
budget is an explicit execution failure, not a new convergence policy. One
authorized attempt per case; preserve interrupted attempts and reconcile the
remote job identity before considering any new submission. No duplicate launch
on a lost SSH response. Do not automatically escalate model or compute limits.

Workflow: frozen specification → isolated preparation → independent read-only
prelaunch review → detached remote submission → ordinary-process monitoring →
hash-verified collection → saved-state audit/report → independent read-only
interpretation review. The coordinator must persist state and release the
interactive agent after confirming launch. Reviewer time is bounded separately;
review failure or insufficient evidence prevents promoting a scientific conclusion.
Approval of this pilot is not approval of new physics, reference changes, extra
resolutions or sending the data-request email.

Output: immutable task/source manifest, preparation tests, reviewed specification,
remote launch/identity/logs, incremental requested evaluations and full states,
transfer receipt, fresh audit, tables/figures with missing data explicit and
independent review. The final decision is whether local resolution supports
moving to the independently coupled layer; it is not a declaration that the
experimental discrepancy is fixed.
