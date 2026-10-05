# Development priorities after the S702 transition study

**Usage hold, 5 October at 18:09 UTC:** the local r3 monitor was stopped after
the user reported excessive model usage. Remote CFD continues. No automatic
postanalysis AI review should resume without addressing its budget; see
[hold receipt](orchestration/runs/s414-nose-20261005-r3/monitor_hold.json).

Updated 5 October 2026. Priorities 1 and 2 are complete for the declared
S702 protocols. All seven combined-policy cases passed full selected-state
audit, including an exactly repeated difficult baseline. The original failed
integration matrix remains preserved. Priority 3's source audit records
unresolved author-coordinate/hinge geometry and transition-tape limitations.
The priority 4 wake-extent check passed all five cases and nine selected-state
audits; its largest effect on the flap drag difference is only 0.00821 count.
The reserved-condition 340/240 panel check is complete and audited. Paired
flap drag changes reach 1.05171 counts versus 260/180, exceeding the screening
threshold despite improved local resolution. The 440/320 third-grid check is
complete: all five jobs and thirteen selected-state audits pass, with an exact
control, but the numerical convergence screen fails. Latest individual/paired
changes reach 0.69746/0.75479 count. The CL-1.2 flap ranking remains reversed.
The automatic grid sequence is closed. The eleven-target contour/surface/wake
separation study is complete and all fifteen selected states pass audit:
surface sampling dominates the CL-0.8 paired change (-0.56128 count), followed
by first wake spacing (-0.19597). Numerical convergence remains unestablished.
Independent contour rankings remain unvalidated;
priorities 3–4 remain open and priorities 5–6 await further evidence.

## Verification plan for the S414 nose-pressure and transition discrepancy

Prepared 4 October 2026; execution started the same day at the user's request.
**Status: V0 complete, V1 gate passed, V2 complete and audited (both
tranches; allowance used), V3 complete with its gate open but bounded; V4
complete and audited (reversal rejected as a sufficient low-lift cause); V5
complete and numerically verified in two tranches (conditional growth comparison
passes; low-lift closure screen fails; coupled formation unresolved); V6 not
started because no remedy is justified.** The live stage
table, counts and decision table are in
[validation/S414_VERIFICATION_STATUS.md](validation/S414_VERIFICATION_STATUS.md).
The stage text below is the frozen plan, not a result. This section is the
current execution order; dated studies below retain their original findings
and historical next-step recommendations.

**Post-plan review, 5 October:** the
[evidence review](validation/results/s414_post_verification_review_20261005/FINDINGS.md)
matches 153 archived output hashes and re-derives the drag and V5 screen results.
The early-transition cause is still unresolved. V5 supports amplification on
the modeled bubble; it does not exclude closure/interaction effects on bubble
formation because archived displacement was imposed and the CL 0.056 momentum
thickness screen failed (9.3% versus 5%). The completed frozen protocol is not
a passed physical gate. The next selected work is a six-case local coupled
nose-resolution check, then an independently coupled layer if the discrepancy
survives; exact experimental/MSES inputs remain the parallel data priority.
See the [decision supplement](validation/S414_V5_NOTES.md#post-plan-decision-5-october-2026)
for bounds and controls. No remedy is promoted, and the original stages and
results remain preserved below.

**Orchestration pilot, 5 October:** the user authorized the six-case local
resolution study. [Revision 3 and live status](orchestration/runs/s414-nose-20261005-r3/README.md)
has passed isolated preparation (43 short checks, actual linear/source-coupling
checks on all three grids) and independent prelaunch review. The detached pilot
launched on jony4 at 18:06 UTC; its local monitor is running. Scientific results
remain pending. The [frozen protocol](validation/S414_NOSE_RESOLUTION_PROTOCOL.md)
replays five archived incidence requests per case, at most 30 requests, with
two single-threaded workers on jony4. Exact archived baseline controls precede
refinement. The 5400-second complete-launch limit includes preflight and cleanup.
Transfer verification, fresh saved-state audit, figures and an independent
postanalysis review are automated after submission; further simulations are not.
Revision 1 retains an unavailable CLI-model failure; revision 2 retains the
independent review's rejected timing/verification gaps. Revision 3 corrects
those gaps and clarifies topology across all grids versus transition movement
between the finest two. Scientific inputs/physics and historical receipts stay
unchanged. Read the live state before resuming; do not submit a second job.

### Objective and definition of resolution

Explain why the model predicts early fore-lower transition and a large drag
drop between CL 0.056 and 0.175 while the qualified experiment and published
MSES do not. Investigate three linked questions separately: the nose pressure
input, the reconstruction-induced turning reversal, and the boundary-layer/
amplification response. A smaller drag error alone does not resolve them.

Resolution requires a repeatable numerical result, an intervention that
distinguishes competing causes, and a pressure/transition/drag comparison under
the same declared conditions. If available geometry or measurements cannot
distinguish causes, the deliverable is a quantified limitation and a specific
data requirement. Do not promise a unique physical explanation from sparse taps.

| Established evidence to retain | What it leaves open |
|---|---|
| The corrected derivative leaves the CL 0.056/0.284 solutions effectively unchanged | Its effect on convergence/branch selection in the difficult intermediate interval |
| At CL 0.056, nominal/cubic drag errors are +55.35/+56.02 counts; both retain the nose event | This particular interpolator change does not isolate the local turning reversal |
| At −2.56888°, nominal/cubic fore transition is 0.75011/0.08282c and drag differs +37.53 counts | Whether the response is reproducible at identical incidence and under both continuation directions |
| Inviscid excess suction exists; measured incidence does not repair the full pressure profile | Geometry, formulation, discretization, boundary conditions and measurement reduction remain distinct candidates |
| Amplification balances and transition roots close; regularization supplies only 0.007 of 9 N units at the nominal low-lift state | Correct solution of these equations is not validation of the separated-layer/transition model |

Evidence: [fore-nose audit](validation/results/s414_fore_nose_20261003/FINDINGS.md),
[derivative comparison](validation/results/s414_jacobian_analysis_20261003/FINDINGS.md),
[incidence comparison](validation/results/s414_incidence_analysis_20261003/FINDINGS.md),
[load budget](validation/results/s414_pressure_budget_20261003/FINDINGS.md),
[inviscid decomposition](validation/results/s414_inviscid_budget_20261003/FINDINGS.md),
[failed pressure discrimination](validation/results/s414_prescribed_pressure_20261004/FINDINGS.md),
and [fore-only reconstruction](validation/results/s414_fore_reconstruction_analysis_20261004/FINDINGS.md).
Reuse their data and verification receipts; do not repeat completed sweeps.

### Common comparison rules and acceptance criteria

Freeze Re 1,000,000, Mach 0.10, Ncrit 9 on both elements, common reference
chord, finite trailing edges, nominal 340/240 panels and 1c wake, unless that
specific input is the declared intervention. New standard runs use explicit
`corrected_af`. Legacy controls retain explicit `legacy` and their own baseline.
Keep reported total-chord positions distinct from local element-chord positions.

Maintain three separate comparison bases: identical incidence for mechanism
isolation; identical lift for drag validation; identical prescribed edge velocity
for a conditional boundary-layer comparison. Never score drag improvement from
different lifts, equate inviscid Cp with the viscous experiment, or use a failed
target's nearby converged evaluation as its result. Preserve source, geometry,
code, runtime, initialization, intermediate failures and complete saved states.

Numerical admission retains both boundary-layer residual L2 and maximum coupling
residual below 1e−7, and absolute CL error at most 2e−5 for fixed-lift targets.
Independently rebuild residuals, forces, Cp/Cf, station coordinates, transition
roots and wake drag from each saved state. A convergence flag alone is insufficient.
Use the archived 2e−5 scaled finite-difference criterion for the affected Jacobian,
with unchanged residual values within its existing 1e−10 allowance.

For proposed same-angle path comparisons, predeclare equivalence screens of
2e−5 CL, 0.01 drag count, 1e−4 Cp at common stations and 1e−4c transition position,
alongside transition/separation topology. These are numerical decision thresholds,
not experimental uncertainties or universal grid-convergence criteria. Exceeding
them triggers investigation, not adjustment of the thresholds. Exact replay
controls under the original pinned runtime retain exact saved-array checks.

For each physical attribution, report the effect alongside numerical spread,
geometry sensitivity, pressure sampling and source/reduction uncertainty. Do not
combine unknown/systematic terms into a fictitious statistical confidence band.
An attribution is unresolved if its effect cannot be separated from those limits.

### V0 — Freeze the evidence and qualify the reference coverage

Reconcile existing manifests into one case/coverage table; acquire only missing
evidence. Preserve the caption-versus-appendix distinctions and existing tap
exclusions. The experimental normal-force/lift definitions and already applied
tunnel corrections remain explicit; do not apply another fitted angle correction.
The completed load budget already shows that force-definition differences are
too small to explain the measured-angle lift discrepancy.

| Reference | Admitted use | Missing evidence |
|---|---|---|
| 2010 free-transition experiment, CL −0.072/0.056/0.175/0.284 | Qualified pressure markers and matching drag records, with original masks and conditions | Exact case-specific transition/reattachment locations and sufficiently resolved nose pressure |
| Published 2010 envelope-MSES | Drag at 0.056/0.175/0.284; pressure at nominal 0.28 as the 0.284 control | Pressure and transition for 0.056/0.175, exact calculation geometry/deck and detailed solver settings |
| CL 1.041 pressure/drag and other completed S414 states | Regression controls after a proposed remedy | They are already inspected, not blind validation |
| 2014 full-e^N MSES and related measurements | Separate, condition-qualified supporting series | They must not be merged with the 2010 envelope series or counted as independent experiments where overlapping |

Use both printed-table and plotted tap coordinates. Distinguish extraction
repeatability, reported manufacturing tolerance, reconstruction spread and actual
as-built uncertainty. Public sharp-TE coordinate files are not substitutes for
the modified finite-TE test model. The model's unresolved between-tap nose peak
must remain visibly unmeasured. Finite tap size alone does not justify a fitted
pressure-averaging correction.

**Gate/output:** an evidence matrix with source page/file hashes, applicability,
measurement resolution, reported uncertainty and explicitly unknown quantities.
Strict matched-geometry validation stays conditional until geometry is qualified.
See the [source qualification](validation/results/s414_geometry_qualification_20261002/FINDINGS.md)
and [reference register](validation/REFERENCES.md); do not repeat a broad literature search.

### V1 — Verify the integrated default before using new flow results

When execution resumes, verify the standard constructor, CLI/config propagation,
copied model parameters, coexistence of corrected/legacy models, invalid-mode
rejection and mode provenance in successful and failed outputs. Check the
optimizer's correction-module hash and resume guard. Import/dispatch checks must
cover both standard-solver and compatibility-adapter entry points.

Reuse the previous derivative evidence and run the affected existing derivative/
assembled-Jacobian checks through the new standard entry point on the saved CL
0.056/0.284 states. Cover inactive/ramp/full amplification and Hk clipping;
finite differences must respect smooth branches, with branch boundaries checked
separately. Restore the archived states without a new nonlinear campaign and
compare standard corrected math to the frozen optional implementation. Legacy
must reproduce its original residual/Jacobian. The growth/residual equations
must remain unchanged. This verifies integration, not a new aerodynamic benefit.

**Gate/output:** a new source snapshot and implementation receipt. A failed
dispatch, derivative or state-replay check stops new nonlinear interpretation.
Do not weaken old manifests to accommodate the changed source. Preserve the known
baseline far-field test limitation; if touched behavior changes it, investigate.

### V2 — Determine the repeatability of the transition-boundary response

Keep the previously recommended bounded fixed-incidence diagnostic. First
reproduce nominal/cubic endpoints using their recorded initialization in explicit
legacy mode, then establish corrected endpoints with the same initialization.
The original full-precision angles come from the saved exports, approximately
−2.56888° and −2.46888°; do not replace them by rounded numbers in the run deck.
Endpoint exports lacking full Newton states receive a fresh saved-state audit;
do not describe their old export checks as equivalent to that audit.

The first tranche is four corrected paths (two geometries × increasing/decreasing
incidence), each on the same five equally spaced angles including both endpoints:
20 requested points, plus four legacy endpoint controls. Use independent endpoint
seeds for the two directions and save full states at every point. Retain 100
Newton iterations per attempt and the existing residual tolerance. Record all
internal continuation steps so a warm start is never mislabeled a cold solution.

Allow at most four additional incidence locations per geometry, evaluated in
both directions, and one independent cold check at a discriminating angle per
geometry: at most 42 requested points including the first tranche. Choose
additional locations by a predeclared bisection rule around path differences or
the transition change, not by lower drag. If nominal transition is already late
throughout the cubic bracket, locate its earlier transition bracket from the
archived CL 0.056-to-0.175 history and spend only the same bounded refinement
allowance there. Report achieved angular resolution; no automatic extension.

**Gate/output:** CL, CD, fore/aft drag, transition, minimum Cf, Hk, N and residuals
versus incidence, with each direction and failure visible. Agreement supports
repeatability only on sampled paths. Distinct residual-converged states at one
angle establish numerical multiplicity; they do not prove physical hysteresis.
If there is a fold, do not force the monotone fixed-lift solver through it.
Recover CL 0.175 only on an admitted bracket/branch; otherwise keep it missing.
Unresolved branching requires a separately bounded branch-continuation method
before unique geometry/closure effects can be claimed.

### V3 — Separate pressure-formulation error from geometry and sampling

After V1, this inexpensive linear-analysis track can be prepared locally while
V2 runs. Start with the identical archived final polygons and incidences,
including −3.07°/−1.04° and the low-lift/transition-boundary controls. Use the
existing Hess–Smith implementation as a cross-formulation diagnostic only after
qualifying boundary/Kutta residuals, closed-body source flux, symmetry, element
ordering and integrated force/circulation consistency on suitable controls.
Neither implementation is assumed correct merely because its matrix solves.

Compare incompressible edge velocities first, then apply the same admitted Mach
correction for Cp comparisons. Preserve each finite TE and document Kutta/base
treatment; unmatched boundary conditions cannot be blamed on a linear solver.
Compare common supported surface locations, signed velocity/stagnation location,
full nose profiles and all four surface loads. Do not confuse midpoint and nodal
samples or sharpen the TEs to access the sharp-only linear-vortex backend.

If formulations disagree materially, inspect self/near-panel influence,
normals, Kutta enforcement and force integration. A bounded three-resolution
linear study (baseline plus two declared refinements) must hold the underlying
contour fixed and preserve aft/wake inputs where applicable. Set its exact
sampling schedule before execution. Require residuals below the implementation's
1e−8 linear criterion and show refinement trends for flux/load diagnostics.
An algebraic residual is not spatial convergence. Persistent disagreement or
nonconvergent trends leave this gate open; there is no automatic fourth grid.

**Gate/output:** pressure/velocity differences and numerical sensitivity on
identical contours. Agreement within resolved discretization spread weakens a
formulation-specific explanation; it does not validate pressure against experiment.
Disagreement must be resolved or bounded before attributing the nose peak to geometry.

### V4 — Isolate the actual turning artefact, then test its aerodynamic effect

Trace the reported polyline through every transformation: source interpolation,
arclength resampling, curvature paneling, TE treatment and final solver nodes.
Inspect tangent, signed turn, spacing and smooth-curve curvature in an ordered
surface coordinate. Separate changes of the physical reconstructed curve from
changes of its numerical sampling. Adding nodes is not new geometry information.

Construct at most two fore-nose alternatives independently of pressure/drag
agreement: one removing an identified resampling artefact while retaining its
underlying curve, and, only if needed, one constrained local reconstruction.
Preserve the aft element, registration, reference chord, finite TEs and the fore
surface outside the declared nose support, with smooth joins. Do not impose
convexity if the primary geometry does not justify it. Admit an alternative
only after showing that the unwanted reversal is absent in the **final** solver
geometry and that its source departures are documented independently of drag.
The earlier cubic variant fails the reversal-removal requirement and remains
an informative sensitivity, not this controlled intervention.

Compare inviscid pressure cheaply first. If removing the artefact leaves the
relevant recovery unchanged within numerical spread, deprioritize it as the
cause; do not repeat a large viscous matrix. If pressure changes materially,
run only the admitted alternatives at CL 0.056/0.175 with the V2 branch protocol,
plus the CL 0.284 control. Also compare selected states at identical incidence.
Any panel-sampling check must use the same continuous geometry and fixed wake
controls; keep its bounded schedule separate from the geometry intervention.

**Gate/output:** before/after geometry, pressure recovery, Hk/N, transition and
drag, with a numerical sensitivity estimate. A causal sampling explanation needs
removal of the artefact, a repeatable pressure/transition response, and consistency
with supported observations. A smoother contour or lower drag alone is insufficient.
Without as-built geometry, retain the conditional nature of geometric attribution.

### V5 — Verify the boundary-layer/transition response independently

Proceed once pressure and branch uncertainty are bounded enough to discriminate
models. Reuse the already verified frozen-state amplification balances/root
checks. Inspect the momentum/shape closure and separated-laminar response, not
only Ncrit: compare theta, displacement thickness, Hk, Cf, edge-velocity gradient,
amplification growth and transition location through the nose recovery.

Qualify an independent boundary-layer/stability implementation against simple
attached and adverse-gradient controls before the S414 comparison. Use identical
edge velocity, Reynolds/Mach definitions, stagnation/upstream data, wall conditions
and disturbance assumptions. Shared inherited formulas are a code cross-check,
not independent physical validation. Full stability calculations require adequate
velocity profiles and disturbance specifications; an integral state alone does
not uniquely supply measured profiles. If separation defeats direct marching,
qualify an interactive/inverse method or a suitably resolved independent viscous
calculation; do not silently switch equations or call nonconvergence physical separation.

The previous experiment-tap replay is **inconclusive**. Do not simply relaunch it:
first make the model-at-tap sampling controls and seed-sensitivity checks pass,
then bracket plausible between-tap pressure reconstructions without fitting drag.
Only interpret a measured-pressure response when this sampling/seed dependence
is bounded. A result conditional on a model seed remains labeled conditional.
Near-threshold regularization may be checked as one explicit derivative-consistent
sensitivity if V2 shows it matters there; its small CL 0.056 contribution already
argues against repeating it as the main low-lift explanation. Keep Ncrit at 9.

**Gate/output:** a closure/pressure/initialization distinction supported by
admitted independent calculations. A closure limitation is supported if plausible
pressure inputs retain the faulty response while a separately qualified model
and applicable measurements consistently disagree. Agreement among codes sharing
the closure does not settle the physical question. If the evidence is insufficient,
state the missing measurement rather than introducing a fitted transition law.

### V6 — Validate a justified remedy and close the evidence chain

Change one demonstrated cause at a time. An implementation correction gets a
targeted regression and affected-case replay; a physical closure proposal stays
behind an explicit model option with its original baseline. Keep Ncrit, references
and scoring fixed. Evaluate every required target, including failed ones.

Primary acceptance returns to CL 0.056/0.175: absolute drag errors and the paired
change (archived MultiFoil −46.14 counts, experiment +5.12, published MSES +0.03),
fore transition and supported nose recovery, plus aft recovery and the full
pressure/load budget. Check the already available −0.072/0.284/1.041 states for
regressions and relevant S702 controls if a shared closure/core changes. These
are regression/consistency checks, not blind validation. Reserve an independently
qualified condition or S103 dataset before using it to choose a remedy; otherwise
explicitly report that independent validation remains outstanding.

Require complete audited coverage, reproducible branch selection and bounded
local discretization error before accepting an improvement. The correction must
explain the pressure-to-boundary-layer-to-drag chain, not just cancel errors in
total drag. No universal tolerance for experimental Cp/transition is invented:
use source-supported uncertainty and tap resolution, and retain unknowns. Small
design gains and global convergence remain uncertified unless separately established.

If pressure/transition agree while drag remains discrepant, then activate a
bounded wake/drag audit: momentum balance, reference normalization, wake-length/
sampling sensitivity, modeled fore/aft contributions and applicable experimental
wake-survey reduction. S702 wake insensitivity does not establish S414 insensitivity.
Do not assign the residual to experiment without specific inconsistent source
records, corrected reduction calculations or independent measurements supporting it.

Missing-data escalation is specific: modified as-built fore/aft coordinates and
rigging, low-lift nose pressure resolved around the peak, case-specific transition/
separation evidence, tunnel disturbance/reduction information, and the corresponding
MSES input/output decks. First test public access and existing archives. Preparing
a request does not authorize sending it. There is no installed MSES executable;
published curves remain references, and new MSES calculations are conditional on
obtaining legitimate access and a qualified matching setup.

### Execution, completion tracking and stop rules

V0 is an evidence-consolidation task; V1 is the next implementation gate; V2 is
the next bounded nonlinear run. V3's linear comparisons can proceed after V1,
independently of V2's execution. V4/V5 are conditional on those results; V6 follows
only a justified remedy. This order preserves the current question and prevents
another open-ended sweep. Stage states are tracked in the status file above.

Before each tranche, freeze a runnable protocol with exact inputs, cases, limits
and thresholds. Use a new directory under `validation/results/`, a pinned source
snapshot and environment, exclusive locking, incremental complete states, and
logged restart/resume instructions. On jony4, use at most two workers with
BLAS/OpenMP threads fixed to one. Inspect server availability before launching.
Launch long work detached, verify its start, provide exact status/log/resume
commands and end the turn; do not keep the agent waiting for completion.

Create a stage status table distinguishing planned, running, execution finished,
audit passed, interpreted and blocked. Track converged/failed/interrupted counts
separately from completed jobs. Each completed stage supplies source and result
manifests, verification receipt, machine-readable comparison tables, reviewed
figures and a concise finding linked here. Preserve failures. A failed numerical
gate stops its dependent interpretation; unresolved data can block physical
validation while independent numerical verification continues. No stage extends
its budget or changes acceptance criteria automatically.

The final decision table must classify each original issue as: demonstrated
implementation defect corrected; numerical sensitivity bounded; physical mechanism
supported by independent evidence; or unresolved with the exact missing evidence.
The project is not declared physically validated because all scripts finish.

## Current work after the literature review

**Verification plan execution, 4–5 October (V0–V5 complete):**

- **V0, evidence matrix.** Existing S414 manifests are reconciled into 95
  reference rows with 27 explicit unknowns; an independent verifier re-reads
  181 evidence files and 659 values. Experimental and MSES transition
  locations are missing at all five cases; MSES pressure exists only at nominal
  CL 0.28/1.04. Strict matched-geometry validation stays conditional. Two
  experimental drag conventions (Figure 13 markers, appendix table) coexist and
  are both retained. See
  [findings](validation/results/s414_evidence_matrix_20261004/FINDINGS.md).
- **V1, integrated default.** The standard constructor reproduces the frozen
  optional correction exactly at the archived CL 0.056/0.284 states (20 exact
  array comparisons; derivative-only residual change 2.75e−11), and explicit
  legacy reproduces the historical mathematics exactly. The archived 32 probes
  pass unchanged (maximum 4.47e−8). Two added ramp-branch probes exceeded 2e−5
  at the coarse step; refinement shows second-order truncation, and the first
  failed receipt is kept. Thirteen new dispatch/propagation/provenance tests
  pass. The complete suite has 38 new failures, all historical source guards
  refusing the five changed files; the far-field limitation is preserved. See
  [findings](validation/results/s414_default_integration_20261004/FINDINGS.md).
  This is implementation verification, not an aerodynamic result.
- **V2 tranche 1, complete and audited.** Six jobs under the frozen
  [protocol](validation/S414_TRANSITION_BOUNDARY_PROTOCOL.md). Both legacy
  endpoint controls reproduced the archived exports identically under the
  pinned runtime. Of 20 `corrected_af` path points, 19 converge; all 457
  retrieved files match and all 25 saved states pass a fresh residual/force/
  Cp/Cf/station/transition-root/wake-drag audit (maximum residual 4.65e−8).
  Increasing and decreasing paths agree at all nine comparable angles (CL to
  1.4e−9, CD to 1.1e−6 count); corrected and legacy endpoints agree; both round
  trips return to their seed. No second converged state was found. The nominal
  geometry is flat in this bracket (transition 0.750–0.794 of total chord, drag
  within 0.51 count). The cubic geometry drops 37.86 counts over 0.1° as
  fore-lower transition moves 0.083 → 0.767 through converged intermediate
  states, and its lift is non-monotone: CL 0.175 is bracketed in three
  intervals, so that matched-lift target stays **missing**. One decreasing-path
  point failed at −2.51888° where the increasing path converged; this is a path
  difference, not a fold or hysteresis. Repeatability holds on the sampled
  paths at 0.025° spacing only. See
  [findings and two graphs](validation/results/s414_transition_boundary_analysis_20261004/FINDINGS.md).
- **V2 tranche 2, complete and audited; V2 closed.** Eighteen evaluations
  under the frozen [protocol](validation/S414_TRANSITION_BOUNDARY_T2_PROTOCOL.md),
  using the full 42-point allowance. Continuation starts from restored audited
  states; both restore controls match the tranche-1 live results with zero
  difference. Of 16 new evaluations 12 converge; all 504 retrieved files match
  and all 14 saved states pass fresh audits (maximum residual 4.24e−8). All five
  comparable direction pairs agree, and both cold starts reach the path state
  (cubic at −2.54388°, a user-approved, recorded deviation from the declared
  −2.51888°; nominal at −2.81944°). No incidence has two distinct converged
  states. The nominal transition change is located between −2.694° and −2.569°
  (drag −38.87 counts, transition 0.060 → 0.750 of total chord, CL 0.155–0.169).
  The cubic change is confined to −2.556°…−2.519°; all four continuation
  attempts into −2.544°…−2.519° fail from both sides, although archived legacy
  states exist there. The interior of both intervals is unresolved and cubic
  CL 0.175 stays **missing**. See
  [findings and two graphs](validation/results/s414_transition_boundary_t2_analysis_20261004/FINDINGS.md).
- **V3, cross-formulation pressure.** Hess–Smith reproduces the production
  fore-nose suction peak within 0.017 Cp at the same between-tap position on
  identical polygons, and the nose profile within 0.047 Cp. The predeclared
  agreement criterion is met in 0 of 10 cases, and a lift offset of about 0.01
  does not vanish on refinement; a control diagnostic attributes it to
  Hess–Smith's solid blunt-base treatment. Gate open but bounded: a
  formulation-specific explanation of the nose suction is weakened, not
  excluded, and nothing is validated against experiment. See
  [findings](validation/results/s414_formulation_check_20261004/FINDINGS.md).

- **V4, geometry trace and inviscid pressure, complete.** The fore contour is
  known only as 83 figure-derived vertices, 15 of them at x/c ≤ 0.06. The
  reported polygon is convex in the nose but nearly straight at vertex 45
  (x/c 0.01449, +1.01° turn, where vertex spacing jumps fourfold). The
  production PCHIP curve overshoots into concavity there (x/c 0.0132–0.0187,
  −0.99°) and near the upper leading edge, and has a 12.1° tangent break at the
  branch-split vertex; the solver's spline and paneling only inherit this
  (−0.69° at x/c 0.0159 in the final polygon). The reversal is therefore a
  property of the interpolated curve, not a paneling artefact, and more nodes
  cannot remove it. Alternative A (single curve, no branch split) keeps the
  reversals and is not admitted. Alternative B (convex chain of circular arcs
  through the same vertices, minimum bending energy, departure ≤ 0.107 mm) has
  no negative turn in the final nose polygon and is admitted. B changes the
  inviscid lower-nose pressure by up to 0.056 Cp (two of six incidences above
  the 0.05 threshold frozen beforehand; Hess–Smith 0.052), moves the suction
  peak aft by 0.001c without lowering it, and gives 0.04–0.05 **more** suction
  at the first tap behind the peak with a roughly 30 % stronger deceleration.
  Removing the reversal does not explain the excess nose suction. See
  [findings and two graphs](validation/results/s414_nose_geometry_20261004/FINDINGS.md).
- **V4, viscous check, complete and audited; V4 closed.** Six jobs under the
  frozen [protocol](validation/S414_NOSE_VISCOUS_PROTOCOL.md). The nominal
  CL 0.056 control repeats the archived corrected state with 18 identical saved
  arrays; all 495 retrieved files match; all 19 saved states pass fresh audits
  (maximum residual 4.92e−8); 15 of 16 path points converge and both directions
  agree at all seven comparable incidences. On the convex-nose geometry the
  CL 0.056 error is +56.00 counts (nominal +55.35), the CL 0.284 control +6.19
  (unchanged), and the CL 0.175 error **+44.14** (nominal +4.09): fore-lower
  transition stays at x/c 0.025–0.070 over −2.945°…−2.469°, so its
  transition-driven drag change lies above −2.469° and below the laminar
  CL 0.284 state, and is not located. The CL 0.056 → 0.175 drag change is −6.74
  counts (nominal −46.14, experiment +5.12, published MSES +0.03); it is closer
  to the experiment only because both states now carry a large error. At the
  second fore-lower tap the model keeps 0.20–0.34 too much suction with either
  reconstruction and with early or late transition. **The reversal is not the
  cause of the discrepancy.** Three reconstructions within 0.11 mm of each
  other place the drag change at least 0.1° apart, so the nominal CL 0.175
  agreement is not robust and the paired drag change is not a reliable score on
  this figure-derived geometry. See
  [findings and two graphs](validation/results/s414_nose_viscous_analysis_20261004/FINDINGS.md).

- **V5 tranche 1, attached-laminar replay controls, complete and verified.**
  Frozen [protocol](validation/S414_LAMINAR_REPLAY_PROTOCOL.md); 114
  finite-difference and 121 closure marches, no coupled solve. The first
  execution showed a sawtooth in the finite-difference wall shear at kinks of
  the edge velocity; it is kept, the streamwise differencing was changed to a
  damped backward scheme and everything was rerun with unchanged acceptance
  values (recorded deviation, not approved beforehand). Results: (a) all six
  failures of the earlier replay lie at the direct-mode singular interval, and
  at that station the inherited equations have **two roots** for the same
  upstream state and edge velocity (shape factor 4.21 or 3.67 at CL 0.056, 4.83
  or 3.34 at CL 0.175); downstream N at x/c 0.0208 is 8.2 or 1.6, and 5.1 or
  2.1. The earlier native control passed only because it started from the
  archive. (b) A new independent finite-difference laminar solver passes 13
  gate controls (Howarth separation 0.1196, cylinder 104.3°) and agrees with
  the inherited closure under identical pressure in all three archived states
  (momentum thickness within 1 %, separation within 0.0007 of chord); the
  closure's shape factor is 0.1–0.2 higher just before the bubble. (c) 10 of 54
  non-trivial sampling controls pass; nothing is admitted at CL 0.056 or 0.175.
  (d) With the measured taps, taps-only reconstructions stay attached to x/c
  0.034 at all three lifts and model-shape reconstructions separate near 0.006.
  Outcomes by the frozen rule: not bounded at CL 0.056 and 0.175; "separation
  retained" at CL 0.284 on one family, conditional and not counted as
  discrimination. **The direct prescribed-pressure replay is closed as a
  method.** See
  [findings and two graphs](validation/results/s414_laminar_replay_20261004/FINDINGS.md)
  and the [branch diagnostic](validation/results/s414_direct_branch_diagnostic_20261004/summary.json).

- **V5 tranche 2, independent layer and linear stability through the bubble,
  complete and verified; V5 gate reached.** Frozen
  [protocol](validation/S414_BUBBLE_STABILITY_PROTOCOL.md), no deviation; 11
  inverse boundary-layer marches and 7 stability envelopes, no coupled solve.
  Three new modules share no formula with the inherited closure: an inverse
  finite-difference laminar solver that passes through separation (tabulated
  reversed Falkner–Skan wall shear reproduced to 3e−4), an Orr–Sommerfeld solver
  (Jordinson's Blasius eigenvalue 0.30859 − 0.00570i; critical Reynolds number
  519.4) and an e^N integrator (Blasius N = 9 within 3 % of the inherited
  correlation); 16 gate controls pass. On the archived bubble's displacement
  thickness the independent layer needs the model's pressure within 0.016 Cp at
  CL 0.056 and 0.009 at CL 0.175 (declared 0.05). Its momentum thickness is
  9.3 % lower at the end of the deep CL 0.056 bubble (outside the 5 % screen;
  inside at CL 0.175). Linear stability on its profiles gives N = 9.93 at x/c
  0.0208 at CL 0.056 against 8.45 from the inherited law on the same layer and
  8.23 archived; 5.98 against 5.04 and 5.13 at CL 0.175: ratios 1.18 and 1.19,
  inside the declared 0.80–1.25. **Neither the growth law nor the separated-layer
  closure produces the early transition; given the modeled bubble, linear
  stability also reaches N = 9 there.** On similar profiles the inherited
  correlation instead grows faster than linear stability (envelope slope ratio
  0.91 at H 4 down to 0.66 at H 11.3), a caution about the correlation and not
  a result for the bubble. Conditional on the archived displacement thickness;
  a coupled independent calculation was not made. See
  [findings and three graphs](validation/results/s414_bubble_stability_20261005/FINDINGS.md).

**V5 outcome.** Independent calculations give no support for a closure
limitation as the origin of the early fore-lower transition at CL ≤ 0.155. With
V1 (derivative), V2 (repeatability on sampled paths) and V4 (turning reversal)
also excluded, what remains is the pressure input that creates the nose bubble
(the model has 0.20–0.34 too much suction at the second tap) or assumptions
shared by all e^N calculations. Neither can be decided with the archived
measurements. No remedy is justified, so V6 has nothing to validate. The open
action is the missing-data request, drafted and not sent:
[draft](validation/results/s414_data_request_20261005/REQUEST_DRAFT.md).
Suggestions outside the plan are kept in
[validation/S414_V5_NOTES.md](validation/S414_V5_NOTES.md). No baseline, closure
or Ncrit change is justified. No simulation is running.

**Standard Jacobian correction integrated, 4 October:** new `CoupledSolver`
runs default to `corrected_af`, using the existing derivative-only correction.
The main CLI, validation runner and shape optimizer accept the explicit
`amplification_jacobian` setting; `legacy` remains available for controlled
comparisons. Solution diagnostics record the mode, and optimizer source hashes
include the correction module. The compatibility adapter now forwards to the
standard implementation. Vendored source, completed studies and their manifests
remain unchanged; replay historical studies using their frozen source snapshots.
This supersedes earlier descriptions of the correction as opt-in only.
No tests or simulations were launched for this integration, at the user's
request; verification was limited to static inspection at that time (runtime
verification followed in stage V1, above). The earlier optional-mode
comparison remains the numerical evidence: negligible drag/transition changes
at CL 0.056/0.284. Integration does not resolve the physical discrepancy or
establish solution uniqueness. The next diagnostic below remains unlaunched.

**Focus restored to S414's CL 0.056 → 0.175 transition/drag discrepancy,
4 October:** the priority is explaining the modeled −46.14-count drag drop
(−45.53 fore), versus experimental +5.12 counts and nearly flat published MSES.
The completed pressure/load investigations remain valid supporting evidence;
the proposed independent Hess–Smith comparison is deferred. Neither a lift
budget nor inviscid/viscous decomposition by itself diagnoses the drag error.

**Bounded prescribed-pressure diagnostic completed, inconclusive:** both
native fore-lower controls restore the saved laminar states. Of eight profiles,
four reach their target pressure, but only three remain below Ncrit through the
common pretransition domain. Both model-at-tap sampling controls and both
experimental-difference reconstructions fail continuation. Even completed
counterfactual roots have unresolved physical branch selection. All eight saved
states pass independent residual/pressure/seed checks and 352 input hashes
match. This verifies numerical records, not experimental-flow predictions.
See [findings and graph](validation/results/s414_prescribed_pressure_20261004/FINDINGS.md)
and [verification](validation/results/s414_prescribed_pressure_20261004/verification.json).
No drag or downstream transition prediction is admitted from this diagnostic.

**Fore-only reconstruction finished and audited, 4 October: 3/4 targets converge.**
Both nominal controls repeat 36 saved arrays exactly. All 430 retrieved files
and 278 frozen study inputs match hashes; all three available selected states
pass fresh residual/geometry/force/Cp/Cf/BL audits (maximum residual 3.66e−8).
Eleven converged request exports receive separate export/growth/root checks;
these are not fresh full-state audits. Thirteen targeted analysis/runner tests
pass; three graphs are reviewed. See the
[findings](validation/results/s414_fore_reconstruction_analysis_20261004/FINDINGS.md)
and [verification](validation/results/s414_fore_reconstruction_analysis_20261004/verification.json).

At CL 0.056, fore-only cubic reconstruction increases drag by 0.671 count:
error +55.35 → +56.02 counts, fore transition x/c 0.022466 → 0.022241. Both
retain the nose pressure/recovery event. The cubic CL 0.175 target fails all
three declared initialization paths: two hit Newton limits and one encounters
a converged but nonmonotone lift response. Its paired drag drop is **unknown**;
do not substitute a nearby converged request. All 28 evaluations and failures
are retained, with no finite force assigned to the failed target.

At the identical lower bracket angle −2.56888°, nominal/cubic transition is
x/c 0.75011/0.08282, with +37.53 counts cubic-minus-nominal drag. Their lifts
are 0.169392/0.170075: this is a same-angle sensitivity, not the missing
matched-lift comparison. At the upper angle both transition near the trailing
edge and drag differs by −0.03 count. Small reconstruction changes therefore
strongly affect this modeled transition boundary while leaving the original
low-lift error intact. Exact geometry and solution uniqueness remain open.

**Next recommended bounded diagnostic:** fixed-incidence continuation through
the narrow existing bracket in both directions, with nominal/cubic endpoint
controls and the corrected Jacobian explicitly identified. This examines the
local transition/drag response without the fixed-lift monotonicity assumption;
recover CL 0.175 only if that response supports it. This is proposed, not launched.
It is now stage V2 of the verification plan above, after the V0/V1 gates.
The four-job protocol is closed and no simulation is running. Preserve the
legacy-Jacobian comparison, Re 1m, M 0.1, Ncrit 9, 340/240 panels and 1c wake;
do not fit Ncrit or expand the grid campaign. The
[protocol](validation/S414_FORE_RECONSTRUCTION_PROTOCOL.md) and
[execution record](validation/results/s414_fore_reconstruction_launch_20261004/README.md)
retain the isolated/remote 35-test preflights and the initial endpoint-test
preparation failure. Previous studies and receipts remain unchanged.

**S414 same-angle inviscid/viscous load decomposition completed, 3 October:**
one unchanged linear operator and two measured-angle evaluations reproduce all
174 archived preflight taps. Geometry and saved viscous Cp restore exactly;
the linear matrix residual is 3.66e−15. Eighteen tests and independent pressure,
quadrature and budget checks pass. See
[findings and three graphs](validation/results/s414_inviscid_budget_20261003/FINDINGS.md)
and [verification](validation/results/s414_inviscid_budget_20261003/verification.json).

At −3.07°/−1.04°, inviscid CL 0.18283/0.42934 becomes viscous CL 0.10723/0.35918,
versus experimental captions 0.056/0.284. The net viscous response removes
62%/50% of the supported-tap normal-force difference (59%/48% on full contours)
and improves Cp RMSE on every surface, but fore-lower integrated errors become
negative while upper-surface errors remain positive. This is an uneven net
response, not evidence for a uniform stronger displacement correction. Because
experiment is viscous, the decomposition does not by itself prove an inviscid
solver error or establish which body/wake causes a surface response.

**Deferred candidate: independent potential-flow formulation on
the unchanged final polygons.** Qualify the existing Hess–Smith backend's
boundary/Kutta residuals and compare the same two angles, Mach correction,
fore/aft loads and tap masks. Preserve the finite trailing-edge faces; the
sharp-only linear-vortex backend cannot be substituted by sharpening geometry.
Agreement would reduce suspicion of a formulation-specific inviscid error;
disagreement would identify a numerical/representation question before changing
geometry or closure. No grid extension or fit is automatic. The current
decomposition is complete, with no new nonlinear solve or follow-up job running.

**S414 pressure-derived load budget completed, 3 October:** identical supported
tap quadrature accounts for 88%/92% of the measured-angle excess normal force.
The partial fore/aft errors are +0.01966/+0.02592 at −3.07° and
+0.04132/+0.02799 at −1.04°. At matched lift, fore errors −0.02508/−0.02558
largely cancel aft errors +0.01919/+0.02051. Published MSES at nominal CL 0.28
shows the same signs of surface-load error, with smaller aft-upper error.
See [findings and four plots](validation/results/s414_pressure_budget_20261003/FINDINGS.md)
and [independent verification](validation/results/s414_pressure_budget_20261003/verification.json).

The primary report defines CL using pressure normal force and wake drag;
applying that definition to the model changes CL by only +0.000539/+0.000123
at measured incidence. Table/plotted tap coordinates change the partial
normal-force gap by at most 0.000405. Conditional observed-Cp-range endpoint
scenarios cannot remove the excess, but are not measurement uncertainty bounds.
The experimental closure remainder lies inside those scenarios; this does
not establish an inconsistency in the experiment. Twelve analytic/admission
tests pass; independent checks cover 32 surface and eight MSES quadratures,
eight dense element forces and eight budgets. All previous outputs are intact.

The full-surface inviscid/viscous follow-up is now completed above. The earlier
pressure-budget task itself required no new flow solve. Neither task establishes
an exact geometry or justifies a geometry/closure fit.

**S414 measured-incidence comparison completed and audited, 3 October:** all
four jobs and 12 requests converge, with 36 corrected-control arrays repeated
exactly. All 1,977 retrieved files and 1,879 frozen study inputs pass hash checks.
Four fresh residual/force/pressure/transition audits pass (maximum residual
1.47e−10), 18 targeted analysis/runner tests pass, and three plots are reviewed.
See [findings and plots](validation/results/s414_incidence_analysis_20261003/FINDINGS.md)
and [verification](validation/results/s414_incidence_analysis_20261003/verification.json).

At experimental angles −3.07°/−1.04°, predicted CL is 0.10723/0.35918 versus
reported 0.056/0.284. Low-angle fore-lower transition remains near the nose
(x/c 0.02517). Overall supported-tap Cp RMSE worsens in both cases and the
strongest sampled aft recovery remains one tap interval downstream of experiment.
An incidence offset is therefore insufficient to explain the discrepancy in
these cases. Fixed-angle drag is context only; original matched-lift errors
remain +55.35/+6.19 counts versus published MSES +8.21/+5.33 counts. Exact
tested geometry, numerical convergence and solution uniqueness remain open.

The pressure-derived load follow-up is now completed above. Partial integrals
remain distinct from independently measured element forces. The four-job
incidence study is closed; no further campaign is automatic.

**S414 derivative comparison completed and audited, 3 October:** all four jobs
and 20 requests converge with no failed requests; both legacy repeats match
36 archived arrays exactly. All 987 retrieved files and 1,849 audit inputs
pass hash checks. All four fresh residual/force/Cp/Cf/transition/separation
audits pass, with maximum residual 6.40e−9. Fifteen targeted analysis tests
pass, and both new graphs are reviewed. See
[findings and plots](validation/results/s414_jacobian_analysis_20261003/FINDINGS.md)
and [verification](validation/results/s414_jacobian_analysis_20261003/verification.json).

The derivative correction leaves the aerodynamics unchanged in these cases:
CD changes −3.27e−9/+2.60e−9 counts at CL 0.056/0.284; maximum pressure change
is 7.38e−10 Cp and all transition flags match. Drag errors remain +55.35/+6.19
counts. Reported final-point Newton iterations fall 28→21 and 23→22, but wall
time increases; no speedup or aerodynamic improvement is established. Keep
the algebraic correction explicit and retain historical baselines. This study
is closed; the defect does not explain the observed drag discrepancy under
the declared protocol.

The completed measured-incidence follow-up retained geometry, Ncrit, panels
and physical inputs. Different solve histories retain the solution-branch
limitation; fixed-angle results cannot replace matched-lift drag validation.
See the
[frozen protocol](validation/S414_INCIDENCE_PROTOCOL.md) and
[completed execution/retrieval record](validation/results/s414_incidence_launch_20261003/README.md).
The isolated execution payload passed 58 tests before launch, separately from
the 18 analysis/runner tests reported above.

**S414 isolated derivative correction implementation, 3 October:**
`amplification_jacobian.AmplificationJacobianSolver` adds explicit
`legacy`/`corrected_af` modes without editing the historical solver or vendor
files. Growth values are unchanged; copied per-model parameters select the
derivative. All 32 assembled coupled-system directional checks pass at the
archived CL 0.056/0.284 states: maximum corrected discrepancy 4.47e−8 versus
0.24547 with the historical derivative. Legacy residuals/Jacobians and saved
arrays restore exactly. The isolated payload passes 50 tests, including two
existing viscous numerical checks. See
[implementation findings](validation/results/s414_jacobian_check_20261003/FINDINGS.md).

The four-job comparison finished on jony4 at 15:59:58 UTC under the
[frozen protocol](validation/S414_JACOBIAN_PROTOCOL.md). Its completed audit is
above; [execution and archive records](validation/results/s414_jacobian_launch_20261003/README.md)
remain preserved. No additional derivative target is automatic.

**S414 fore-nose diagnostic complete, 3 October:** archived CL 0.056/0.175
states and the 0.284 control pass amplification-balance, independently bracketed
transition-root and frozen-state march checks. Twenty-three targeted tests pass.
The baseline amplification Jacobian check **fails**: `get_damp` differentiates
`0.1*exp(-20*Hmi)` with coefficient −1 rather than −2. Eighteen of 72
component/step comparisons fail; a diagnostic-only correction passes all 72.
No production routine was changed. This may affect Newton paths, but does not
prove the cause of drag error or invalidate already converged residuals. See
the [findings and three graphs](validation/results/s414_fore_nose_20261003/FINDINGS.md),
[verification](validation/results/s414_fore_nose_20261003/verification.json) and
[bounded protocol](validation/S414_FORE_NOSE_PROTOCOL.md).

The excess fore-nose suction is already present in the unchanged inviscid
calculation at saved incidences, before viscous feedback. Geometry sampling
introduces a small lower-nose turning reversal absent from the reported
polyline; the archived cubic reconstruction retains a similar pressure peak.
Neither observation isolates the cause. Near-Ncrit regularization contributes
only 0.0070 of 9 N units at CL 0.056 on the saved state; base-envelope growth
through the computed high-Hk nose region dominates. No new nonlinear solve was
run, Ncrit remains 9, and all archived baselines remain unchanged.

The explicit derivative option and assembled-Jacobian check are now complete;
the bounded CL 0.056/0.284 comparison is completed above. Distinguish convergence/
branch effects from residual-equation changes; the growth value is unchanged.
Geometry/pressure-input work remains a separate question after this comparison.

**S414 free-transition lower-bucket comparison complete and audited, 3 October:**
all four nominal targets (CL −0.072, 0.056, 0.175, 0.284) pass full saved-state
and fresh residual/force/Cp/Cf/BL audits. All 22 requests converge on cold-local;
the control repeats all 18 arrays exactly with zero drag difference. Maximum
residual is 2.16875e−8 and maximum CL error 1.34774e−5. All 238 retrieved files
and 100 frozen inputs match hashes. Thirty-eight targeted analysis/runner tests
pass; independent table-token drag interpolation agrees exactly and 89 prior
audited outputs remain unchanged. The comparison scores 348 supported pressure
observations with explicit masks. Six graphs are reviewed. See the
[findings and plots](validation/results/s414_lower_bucket_analysis_20261003/FINDINGS.md),
[verification](validation/results/s414_lower_bucket_analysis_20261003/verification.json),
[frozen protocol](validation/S414_LOWER_BUCKET_PROTOCOL.md) and
[execution record](validation/results/s414_lower_bucket_launch_20261003/README.md).

The new discriminating result is **CL 0.056: MultiFoil drag error +55.35 counts,
published MSES +8.21**. MultiFoil retains fore-lower transition at x/c 0.02247,
then moves it to 0.76928 at CL 0.175. Of the computed −46.14-count drag change
between those targets, −45.53 comes from the fore element; experiment changes
by +5.12 counts. This locates a large model response at the fore element without
allocating experimental drag error by element. Supported leading-edge taps also
show pressure discrepancies. Aft recovery remains wrong, but its strongest
sampled rise is downstream of experiment at three cases and upstream at 0.175;
a fixed downstream-bubble diagnosis does not describe the full sequence.

The fore leading-edge geometry/sampling, pressure and amplification diagnostic
is now complete above, separating a confirmed derivative defect from the
unresolved pressure-input discrepancy. Do not fit a wake correction or retune Ncrit.
Exact geometry, numerical convergence, case-specific measured transition and
solution uniqueness remain unestablished. Missing MSES coverage is preserved:
drag at −0.072 and pressure at the three new cases. Re 1m, Mach 0.10, Ncrit 9,
340/240 panels, 1c wake and production equations remain unchanged. That four-job
study is closed; the separate derivative comparison is also completed above.

**S414 trip source qualification complete, 3 October:** six experimental
table curves (169 rows) and six 2010 MSES drag curves are recovered for Re
0.7m/Mach 0.07 and Re 1m/Mach 0.10, with free, fore-only-fixed and both-fixed
labels. All 108 independent plotted experimental markers agree with their
table rows within 0.002695 drag count. At Re 1m/CL 0.284, the fore-trip drag
increase is +37.61 counts experimentally versus +67.49 in MSES; at CL 1.041
it is +29.57 versus +37.36. These comparisons use matched lift, not matched
incidence. Nineteen targeted tests pass. See [findings and graphs](validation/results/s414_trip_qualification_20261003/FINDINGS.md)
and [qualification protocol](validation/S414_TRIP_QUALIFICATION_PROTOCOL.md).

The proposed clean intervention is **not admitted**: the author reports extra
trip displacement and ineffective fore-upper grit at low lift, and supplies
no paired tripped aft Cp, bubble endpoints or wake-turbulence profiles. The
total-force contrasts cannot isolate the missing wake-disturbance channel.
Retain them as qualified context, without a closure fit or forced-trip campaign.
The complementary free-transition Figure-7(a–d) pressure sequence is extracted
at CL −0.072/0.056/0.175/0.284. Its strongest sampled recovery location is not
monotonic, so it cannot substitute for measured bubble length. Independent
Figure-7/12 copies fail a strict precision check by up to 0.000338 Cp and
0.00008123c; offsets remain explicit, with no fitted alignment.

The bounded free-transition lower-boundary comparison is now fully audited
above, repeating the 0.284 control at Ncrit 9 and identifying a large fore
transition/drag discrepancy at 0.056. The nominal geometry
and numerical-convergence limits remain. Published MSES drag covers 0.056,
0.175 and 0.284 but not −0.072, and no additional MSES pressure curves match
the three new conditions. Report this coverage before scoring.
S103 qualification remains complementary; the newly found Twiss
2019 S414 thesis is an unqualified separate-campaign lead, not substituted data.

**S414 fore-wake information audit complete, 3 October:** restored both nominal
states and checked 16 finite-difference/Jacobian directions. All pass, with a
maximum receiving-block relative discrepancy of 1.2741e−10; saved states and
flags are restored exactly. The fore wake enters global mean-flow coupling
through `ue * delta_star`, but carries no direct disturbance/turbulence input
to aft transition. Its shear and momentum can still affect the aft element
indirectly by changing wake displacement after solving. All wakes use a fully
turbulent closure, including the low-lift case with two laminar fore surfaces.
This identifies model assumptions, not the cause or magnitude of the pressure
error. Eight reconstruction states supply geometry/integral-scale context;
displacement thickness does not establish wake width or disturbance overlap.
All 25 targeted tests pass. No production equations, Ncrit, archived results
or flow predictions changed. See the [findings and graphs](validation/results/s414_wake_transfer_20261003/FINDINGS.md)
and [bounded protocol](validation/S414_WAKE_TRANSFER_PROTOCOL.md).

The 2010 oil-flow/acoustic evidence puts CL 0.284 inside the laminar fore-surface
bucket. The authors' proposed wake-turbulence mechanism concerns its boundary
and observed aft-bubble shortening; it cannot be assigned to both current
pressure discrepancies without further evidence. The recommended free versus
fore-only-tripped qualification is now complete above: force data exist, but
trip effectiveness/displacement and absent paired aft pressure prevent clean
causal admission. Its complementary free-transition pressure sequence now
has a completed four-case model comparison above. No closure fit or trip simulation follows.

**S414 bounded reconstruction comparison complete and audited, 3 October:**
all nine jobs finished on jony4 at 21:49:27 UTC, 2 October. All 336 retrieved
files match remote hashes; all nine saved states pass full export and fresh
residual/force/Cp/Cf/BL checks. The repeat reproduces drag and all 18 arrays
exactly. All 46 requests converge without internal continuation failures, using
cold-local paths. Maximum residual is 2.8653e−8; 25 targeted tests pass.
No simulation is running. See the
[findings, pressure/BL/geometry/drag graphs](validation/results/s414_reconstruction_analysis_20261003/FINDINGS.md)
and [launch/retrieval record](validation/results/s414_reconstruction_launch_20261002/README.md).

At CL 0.284 / 1.041, nominal MultiFoil drag errors are +6.1883 / −10.0608
counts; published MSES gives +5.3319 / −13.8949. Both methods place the strongest
sampled aft pressure rise in x/c 0.89313–0.90370, while experiment places it in
0.88215–0.89313. All four reconstruction scenarios and both tap-coordinate
conventions retain this diagnosis. Nominal recovery-region Cp RMSE is
0.08360 / 0.16372, versus MSES 0.06884 / 0.11603. Maximum geometry-induced
drag changes are only 0.0465 / 0.1772 count. These are sensitivity cases, not
confidence bounds; exact geometry and numerical convergence remain unestablished.

This establishes a local pressure mismatch in addition to a drag discrepancy;
drag bookkeeping or wake-survey reduction alone cannot explain both observables.
It does not separate transition onset from reattachment/recovery error or prove
a wake-turbulence mechanism. Close the nine-job protocol without an automatic
extension or Ncrit tuning. The resulting fore-wake information/qualitative
oil-flow audit is now complete above; no new closure or flow run follows from
identifying the missing direct disturbance channel alone.
S103 tested-displacement qualification remains the complementary independent
case. Keep the 2010 envelope-MSES and 2014 full-e^N references distinct.

**S414 source qualification completed, 2 October:** recovered the original
2010 report from two public mirrors; archive supplied checksums pass and both
copies contain identical Figure-1 geometry. The modified nominal test shape
has 83 fore/68 aft vertices and matches four pressure-figure outlines within
7.04e−6c. Its lower fore modification reaches about 2.98 mm at the model scale.
Two public `.pts` candidates have undocumented registration and sharp trailing
edges in place of the report's finite edges; neither is adopted as exact input.
Three tap centers lack support on the reconstructed branch extents. See the
[qualification, geometry and published-pressure plots](validation/results/s414_geometry_qualification_20261002/FINDINGS.md).

The 2010 Figure-12 comparisons at Re 1 million/Mach 0.10 are extracted: 188
experimental markers (8 repeated endpoints) and 480 MSES vertices, nominal
CL 0.28/1.04 versus measured 0.284/1.041. Ordered Table-II tap associations
are unambiguous, but the original 2e−5c exact-coordinate check fails by up to
0.000260c; this failed check and the offsets are retained. The author-described
downstream aft-bubble/recovery prediction is a direct mechanism diagnostic.
These MSES curves use the **2010 envelope method**, distinct from the 2014
full-e^N reference. All seven targeted extraction/provenance tests pass.
No MultiFoil or MSES flow solve was run.

Nominal source reconstruction is qualified; strict matched-geometry benchmark
admission remains open pending original modified coordinates/reference placement
and edge treatment. The drag extraction and bounded protocol are now complete
above. The physical question is whether unretuned MultiFoil also delays aft
recovery at these two states. Preserve Ncrit 9 as primary independent setting,
the current core and the numerical-convergence limitation. Do not silently
promote public coordinates or this reconstruction to an exact tested deck.
S103 qualification remains complementary.

**Geometry cross-check and fore/aft transition verification complete,
2 October:** direct Figure-1 re-extraction reproduces all original inputs.
Four pressure-figure outlines match under their own x-axis calibrations to
1.4343e−5c; after one common fore-registered z translation, maximum z difference
is 4.7704e−6c. The foremost aft tap remains 0.000844–0.000853c ahead of those
contours. This corroborates the extraction but leaves a common source/tested-
geometry inconsistency. No justified correction or original coordinate deck
was recovered from the bounded public-source check. See the
[source reconciliation and plot](validation/results/s702_geometry_reconciliation_20261002/FINDINGS.md).

All ten transition-isolation jobs completed on jony4 at 13:11:39 UTC. The two
controls reproduce drag and all 36 saved-array comparisons exactly. All 18 new/
archived selected states pass audit with actual per-element Ncrit; 10 have fresh
residual rechecks and 8 reuse hash-verified earlier residual evidence. All 39 new
requested evaluations converge; one recovered internal continuation failure is
preserved. All factorial endpoints use the cold-local path. Twenty-three
targeted audit/report tests pass. No simulation is running. See the
[findings and pressure/BL/drag plots](validation/results/s702_transition_isolation_analysis_20261002/FINDINGS.md),
[protocol](validation/S702_TRANSITION_ISOLATION_PROTOCOL.md) and
[launch/retrieval record](validation/results/s702_transition_isolation_launch_20261002/README.md).

At flap 5 / CL 0.828, lowering fore Ncrit alone changes total drag by +23.710
counts; aft-only changes it by −0.410. The large calibration response is chiefly
fore-driven. Aft-only improves the higher-lift pre-hinge pressure rise from
0.03128 to 0.05823 (experiment 0.09031), but removes modeled lower-lift separation
and does not recover that measured pressure feature. At flap 10 / CL 0.774,
fore-only changes aft drag by +5.072 counts while fore drag barely changes;
the total-drag interaction is −3.394 counts. Component responses are coupled.
No tested threshold combination reproduces the flap-5 pressure evolution
consistently. Mixed Ncrit remains a counterfactual, not measured wake turbulence.

**Close this Ncrit diagnostic without further tuning.** The S414 source
qualification above advances the independent pressure/bubble comparison;
S103 remains the complementary MSES/experiment path. Exact
S702 tested geometry and hinge/transition evidence remain open prerequisites for
isolating its physical cause. Preserve the geometry, numerical core and frozen
calibration; no automatic simulation, closure correction or grid extension follows.

**Matched-pressure study complete and audited, 2 October:**
four S702 snapshots from Figures 7(b,c)/8(b,c) are qualified at Re 1 million,
Mach 0.10, free transition. The extraction recovers 386 plotted markers,
including 16 flagged repeated endpoints. Table-II tap positions, legend series,
axis calibration, original source hashes and experimental lift associations
were checked. See [source findings and plots](validation/results/s702_pressure_qualification_20261002/FINDINGS.md).
The two pre-hinge taps show the flap-5 pressure rise increasing from 0.00769 to
0.09031 as incidence rises from 2.03° to 6.11°; flap 10 retains a much smaller
rise at the selected low/high-lift snapshots. This is a pressure feature, not
an exact experimental transition or separation coordinate.

All eight bounded targets (four experimental CL values × Ncrit 9/6.25) completed
on jony4 under the
[matched-pressure protocol](validation/S702_MATCHED_PRESSURE_PROTOCOL.md).
All 272 retrieved files match remote hashes; all eight full saved-state audits
pass, including fresh coupled residual/force evaluation without a nonlinear
solve. All 29 requests converged; maximum residual is 8.3422e−8 and maximum CL
error is 1.6281e−5. The unchanged production solver uses 340/240 panels and the
existing combined initialization policy. Thirteen targeted reporting/data tests
pass. That eight-target study is closed. See the
[findings, pressure/BL plots and published-MSES drag comparison](validation/results/s702_matched_pressure_analysis_20261002/FINDINGS.md)
and [launch/retrieval record](validation/results/s702_matched_pressure_launch_20261002/README.md).

At flap 5, the measured pre-hinge pressure rise increases from 0.00769 to
0.09031 with lift; Ncrit 9 instead predicts 0.03525 to 0.03128 and retains a
substantial separated region. Ncrit 6.25 predicts 0.03667 to 0.05910, improving
the higher-lift recovery but failing to recover both profiles consistently.
Lowering Ncrit improves the low-lift flap-5 drag mainly through fore-element
transition (+24.66 fore counts, −1.70 aft), while the pre-hinge feature barely
changes. Conversely, high-lift flap-10 Ncrit-9 drag is within 1.29 counts of
experiment despite a hinge-region Cp RMSE of 0.0849. These local mismatches
strengthen the model/input investigation without isolating a unique cause.

**Geometry evidence motivating the source cross-check above:** Table II
places the foremost aft upper pressure tap at x/c 0.69941, ahead of our entire
Figure-1 reconstructed contour by 0.00084946c (0.3886 mm). The raw input already
has this inconsistency; paneling changes its nose x by only a few millionths
of chord. Pressure-figure outlines now corroborate the reconstruction but do not
resolve the inconsistency. Continue seeking exact tested coordinates and
slot/hinge/reference-frame evidence; preserve the baseline and do not fit a
translation to drag.
This finding neither quantifies a drag effect nor establishes an error in the
authors' MSES geometry. The unsupported tap is excluded explicitly from scoring.

Published MSES drag is retained at each snapshot's own CL; matching MSES Cp
curves and exact experimental transition locations are unavailable.
These snapshots are not same-CL flap pairs and do not replace Re-700,000 data.
Numerical convergence remains unestablished. No physical model, geometry or
Ncrit fit changed, and no automatic numerical extension follows.

The [reference inventory](validation/REFERENCES.md) adds stronger prospective
checks without invalidating the numerical gate above. Pomeroy 2014 supplies
three-element coordinate and placement tables with measured wakes; its 2015
companion varies placement and measures wake pressure/velocity. The completed
geometry extraction below finds an unresolved placement inconsistency, so
Pomeroy is not ready for quantitative benchmark admission. Its measured wake
cross-profiles also lack a directly matching output in the present solver.
These sources do not establish that wake bursting causes the S702 error.
S103/S414 are now the next priorities for published MSES/experiment comparison,
with tested geometry and transition formulation qualified before admission.
S103's tested element displacement
and S414's modified low-speed contour are material. The TMR near-wake experiment
has accessible data but requires forced-trip support; the TMR multielement case
has no experiment and cannot replace independent validation.

**Completed and audited, 1 October:** the eleven-job
[discretization protocol](validation/S702_DISCRETIZATION_PROTOCOL.md) separates
the existing 340/240-to-440/320 step at CL 0.8. All eleven targets and fifteen
full selected states pass. The four archived endpoints and fresh original-solver
control reproduce drag and shared arrays exactly. The production core remains
unchanged. The paired changes are +0.00730 count from source resampling,
-0.56128 from surface sampling, -0.00483 from wake count and -0.19597 from first
wake spacing, summing to -0.75479. The surface contribution is chiefly in the
aft drag term (-0.59039), offset by +0.02911 fore. Both flap-10 surface endpoints
use the same warm-anchor path. Two failed requests and three failed internal
continuation points remain preserved. Twelve reporting/audit regression tests
pass. See the
[findings and plots](validation/results/s702_discretization_analysis/FINDINGS.md).

These are conditional contributions along the declared path; they establish
neither a numerical error bound nor the physical cause of the experimental
discrepancy. The finest CL-0.8 flap delta remains 18.12595 counts above the
measurement. Small design gains and closure changes remain gated.

**Fore/aft sampling study completed and audited, 1 October:** all five new
targets and all nine selected states pass. The fresh C control reproduces drag
and all 24 saved array comparisons exactly. At fixed source/wake controls,
the B→F fore-sampling change is +0.030226 count in the paired flap response;
F→C aft sampling changes it by -0.591506, summing to -0.561280. Both increments
use matching selected paths within each flap. The alternate mixed A flap-10
state needs warm-entry, so the -0.048890 interaction and contrasts involving A
retain a path-change qualification. Four failed requests and eight failed
internal points among the 29 new requests remain preserved. Sixteen reporting
regression tests pass. That numerical study is closed. See
[the full report and experiment/MSES plots](validation/results/s702_surface_factorial_analysis/FINDINGS.md)
and the [launch record](validation/results/s702_surface_factorial_launch_20261001/README.md).

This closes the bounded numerical diagnosis: aft surface sampling dominates
the observed sub-count change, while the CL-0.8 experimental discrepancy
remains over 18 counts and the CL-1.2 ranking stays reversed. No error bound,
closure correction or validated design gain follows. Source qualification is
the next task; there is no automatic further numerical campaign.

**Ncrit influence on the flap error reanalysed, 1 October:** the user-requested
[focused analysis](validation/results/s702_ncrit_flap_error_20261001/FINDINGS.md)
uses all 168 archived calibration entries and both reserved Ncrit scenarios,
holding the original 260/180 discretization fixed. At Re=700,000/CL=0.8,
Ncrit 9→6.25 reduces MultiFoil's paired error from +35.89 to +19.38 counts;
91% of the response change comes from raising the 5° baseline drag. At Re=1
million the full sweep is non-monotonic: the closest complete sampled CL-0.8
pair is Ncrit 6, still +6.76 counts above experiment. Ncrit 6.25 retains the
wrong CL-1.2 sign at both Reynolds numbers. No new Ncrit was selected, no flow
solve was run, and no finer-grid values were mixed into the sensitivity.
MSES's own sensitivity is unknown because only published Ncrit-9 curves are
available. The S103/S414 qualification priority below remains unchanged.

**Physical diagnosis, 1 October:** the user prioritizes understanding over
perfect agreement. The [evidence assessment](validation/results/s702_physics_diagnosis_20261001/FINDINGS.md)
distinguishes a late predicted 5° low-drag-range exit from excessive 10° drag.
At reserved CL=0.8, measured incidence falls about 1.55° with the flap change;
the model captures fore unloading/laminar-flow benefit but nearly cancels it
with aft losses. Its aft bubbles reattach: neither aft stall nor wake bursting
has been established. The leading suspects are transition/bubble closure and
unrepresented wake disturbances or hinge-cover geometry. Mean-flow interaction
is present; external wake-induced transition is not. Bad experimental drag is
less supported, but raw S702 repeat surveys/uncertainties are insufficient to
exclude case-dependent bias. The observed numerical changes remain smaller
than the residual, without a rigorous error bound.

The recommended **5°/10° aft pressure and bubble comparison** in S702 Figures
7–8 at Re=1 million is now complete, as reported above, with each experimental
operating lift recovered from its own lift curve. These observations do not
replace Re=700,000 pressure measurements or interpolate moving bubble features.
S414's measured bubble shortening and drag dip, missed by MSES/OVERFLOW, provide
corroborating physics; they do not prove the S702 cause. S103/S414 geometry
qualification remains the broader independent-validation path. No new flow
campaign, Ncrit fit or closure correction was launched for this diagnosis.

**Pomeroy source qualification completed, admission unresolved, 1 October:**
all 243 coordinates (81 per element) were extracted with row-level provenance
and visually checked. Local and assembled contours are nonintersecting. Under
the declared clockwise transform, Table 1 placement yields gaps 0.023140c and
0.016441c, versus Table 2's 0.0100c for both. Eighteen checks of contour sampling
and trailing-edge endpoints do not reconcile them. This identifies a source
consistency problem under the stated frame, not which table matches the tunnel.
The 2015 companion has different absolute placement and an opposite overhang
sign convention; it cannot silently substitute for the 2014 geometry. Its
pressure maps are off-body probe measurements, not surface Cp. Six extraction
and geometry regression tests pass; zero flow solves were launched. See
[the source audit, coordinates and plot](validation/results/pomeroy_geometry_qualification_20261001/FINDINGS.md).

For independent validation, qualify S103's exact tested contours/displacement and S414's modified
low-speed contour with their corresponding experimental and published-MSES
curves. Pomeroy needs reconciliation against original conference/thesis data
or an author geometry deck before a flow benchmark; reviewed methods also leave
Mach and trip/Ncrit conditions unresolved. Do not fit placement to force table
agreement. New closure development and optimization remain gated by numerical
and independent experimental evidence.

### Deliverables toward experiment and MSES validation

The user reaffirmed experimental validation and MSES comparison on 1 October.
The bounded sampling audit and first Pomeroy qualification are complete. Continue
source qualification for independent cases; do not extend the global grid sequence. Numerical uncertainty remains
reported even if a useful experimental diagnostic can be qualified.

| Case | Next concrete deliverable | Comparison and admission limit |
|---|---|---|
| S702 | Transition factorial complete: retain the 18-state audit and pressure/BL/drag comparisons; seek exact tested geometry and hinge/transition evidence | Fore calibration and aft recovery responses are separated, with significant coupling at low-lift flap 10. Tested Ncrit changes fail to recover both flap-5 pressure profiles consistently. Geometry and numerical convergence remain unqualified. Published MSES N9 supplies drag context only. |
| Pomeroy 2014/2015 | Extraction complete; resolve the Table 1/Table 2 placement conflict from paper-specific original geometry before admission | No matched flow benchmark yet. Measured wake cross-profiles are not represented by the current solver's wake-edge speed or Squire–Young drag. The 2015 pressure maps are off-body, not surface Cp. No published MSES polar established here. |
| S103 | Qualify exact tested contours and the aft-element rigid placement; begin with the Re 3 million, Mach 0.10 free-transition comparison in Figures 10/11 | Three-way pressure and section-characteristic comparison. At this Reynolds number the aft element is rotated -1° about (0.8321,0.0002)c and translated (-0.0017,-0.0017)c. Table II gives pressure-orifice locations, not a complete contour deck. Published theory uses Ncrit 9; exact geometry and transition-method details remain gates. |
| S414 | Same-angle inviscid/viscous load decomposition complete; next recommended check is an independent potential-flow formulation on identical finite-TE polygons | Net viscosity reduces the supported load gap by 62%/50%, but fore-lower errors cross negative while upper errors remain positive. Viscous experiment is not an inviscid validation reference. Published MSES remains matched-lift context only. Preserve geometry and Mach convention; Ncrit 9 fixed, no closure fit or automatic extension. |

Freeze each admitted benchmark before scoring: geometry/source hashes,
reference chord, Re/M, transition/trips, tunnel corrections, operating range,
digitization method and uncertainty. Report requested/converged coverage first,
then signed bias, absolute/RMS and worst errors on common support, with plots
of pressure/transition/wake observables where available. Compare fixed-lift
drag and placement/flap changes using each configuration's own baseline.
Keep experimental uncertainty distinct from graphical reading allowance;
if uncertainty is unavailable, report raw errors without an uncertainty-based
pass claim. Declare case-specific accuracy criteria before inspecting new
scores. No per-case Ncrit tuning or omission of failed targets is permitted.

Source details above were checked against the preserved primary texts in
`validation/results/reference_review_20260924/sources/`; authoritative URLs
and hashes remain in that review's `source_manifest.json`. Published MSES
curves stay part of the deliverable. There is no installed MSES executable,
so a fresh MSES run is not claimed or required to use those curves honestly.

## 1. Locate the remaining drag error using existing flow states

**Completed 15 September 2026.** The runnable reporter and all outputs are in
`validation/report_s702_drag_mechanism.py` and
[`validation/results/s702_drag_mechanism/`](validation/results/s702_drag_mechanism/FINDINGS.md).
All 42 states were audited, including all 16 primary states; all six flap-delta
rows reproduce exactly. The calibrated CL=0.8 fore-element benefit is offset by
an aft-element penalty, almost completely at Re=700,000. Surface diagnostics tie
these modeled responses to fore upper transition and thicker aft upper layers
with aft lower separation. Experimental error cannot be allocated by element
from the available total-drag measurements. There were zero new solves or core
changes; 24 targeted tests pass. Continue priority 2 before mechanism-changing
experiments. The original acceptance specification follows for traceability.

First deliverable: `validation/report_s702_drag_mechanism.py` and a new results
directory containing an auditable table and paired plots. Start without new
solves. Use the 16 already computed states for flaps 5/10, CL 0.8/1.2, both
Reynolds numbers, and Ncrit 9/6.25. Retain all 21 operating conditions as the
broader comparison set; the CL 0.5 points are useful controls.

Read the authoritative converged target and its actual surface/wake exports.
Account for the old reused Ncrit 9 flow-directory alias using the existing
report helper. Plot Cp, Cf, edge speed, Hk, theta, displacement thickness,
transition positions and separated intervals for fore and aft elements, with
the achieved lift and incidence beside each comparison. Tabulate per-element
wake CD and friction diagnostics, and compute their contributions to the
5-to-10-degree drag change. Do not add friction to wake CD.

Check that total CD equals the sum of element wake terms within the exported
rounding precision, that the selected state satisfies the recorded residual
and lift tolerances, and that each plotted quantity comes from that state.
Reference pressure data at flap 0 / CL 1.04 can be shown as a separate check;
do not relabel them as measurements at the difficult flap settings.

Acceptance: reproduce all six published comparison deltas in
`s702_flap_5_to_10.csv`, identify which modeled element/surface changes dominate
each discrepancy, and distinguish observed numerical behavior from proposed
physical explanations. A diagnostic that rules out a hypothesis is useful;
no drag correction is required to call this task complete.

## 2. Make fixed-lift evaluation repeatable before using tiny differences

**Completed for the declared S702 protocols, 22 September 2026:** the original
11-job/19-target matrix is retrieved and reviewed in
[`s702_repeatability_analysis/FINDINGS.md`](validation/results/s702_repeatability_analysis/FINDINGS.md).
18 targets converged on the first attempt; the cold flap-10/CL=0.8 target failed.
Successful warm paths agree in force and state; the five fully covered conditions
have maximum CD span 2.77e-7 count. Maximum difference from older selected states
is 0.0027721 count. This supports addressing initialization before changing physics.

`fixed_lift_policy.py` now provides explicit fresh cold/warm paths, a bounded
computed-lift recenter, preserved failures and complete selected solver restoration
in memory. The seven-job follow-up passed 6/6 primary cases and 1/1 replay
control, including full selected-state review. A -1-degree warm anchor recovers
the difficult cold target without the archived seed prefix. See
[`s702_initialization_policy_analysis/FINDINGS.md`](validation/results/s702_initialization_policy_analysis/FINDINGS.md).
The original core and sweep worker remain unchanged.

`shape_optimization.py` now accepts an explicit per-point initialization policy,
records every attempt and saves restored selected states. Existing configurations
retain the legacy solve. All 34 targeted shape/policy/optimization/paneling checks
pass locally. The twelve-job matrix finished with 11/12 feasible cases and
5/6 complete comparisons. The largest complete-pair CD difference is 0.0033121
count. All eleven selected exports and five saved policy Newton states passed
audit, but the policy baseline at flap 10° failed. Its declared −1° warm start
also failed internally; the legacy control instead obtained a converged state
at center−0.4°. See
[`design_initialization_analysis/FINDINGS.md`](validation/results/design_initialization_analysis/FINDINGS.md).
That original matrix remains failed. The explicit-entry diagnostic subsequently
passed 3/3 cases and full selected-state audit, recovering the baseline without
changing physics or tolerances. Maximum CD difference from the controls is
0.0006445 count; two internal −1° starts still fail, but the declared entry point
converges and supplies a usable state. See
[`design_entry_seed_analysis/FINDINGS.md`](validation/results/design_entry_seed_analysis/FINDINGS.md).
The fresh seven-job combined check passed all six contour cases plus an
independent repeat of the difficult baseline. The sequence is cold-local,
fresh −1° start, then explicit entry, with first success retained. Full audit
covered thirteen selected exports (seven new and six controls), all seven new
Newton states, initialization histories and the independent repeat. Repeated
forces and every saved array field are exactly equal. Maximum CD difference
from legacy controls is 0.0033121 count; the largest difference in own-baseline
contour deltas is 0.0019178 count, below the 0.5-count screening threshold.
All five failed requested evaluations and thirteen failed internal continuation
points remain recorded. See
[`design_combined_policy_analysis/FINDINGS.md`](validation/results/design_combined_policy_analysis/FINDINGS.md)
and [`DESIGN_COMBINED_POLICY_PROTOCOL.md`](validation/DESIGN_COMBINED_POLICY_PROTOCOL.md).
This completes priority 2 within these conditions; it does not establish
universal path independence or independent contour-ranking validity.

The earlier remote preflight passed 45/46 tests but failed the unchanged 100-chord far-field
element-order test (CL difference 1.14e-6 against 5e-7). The untouched baseline
reproduced the identical failure. That known limitation remains; these changes
do not modify the flow core or claim that test passed.

Target `optimization.solve_fixed_lift`, `CoupledSolver.solve`, and the existing
local/warm-seed utilities. The reserved failures are useful regression cases.
The general fixed-lift routine requires admissible endpoints and can abort when
a trial incidence crosses an inadmissible region; separate drivers currently
provide the successful retries.

Extract a reusable, explicit initialization policy from those drivers. Use
computed lift to construct local brackets, retain branch checks and residual
tolerances, restore the selected solver state, and keep every failed attempt.
Record cold starts, continuation paths and warm replays separately. If full
state checkpoints are added, version them and verify exact physical/numerical
compatibility; existing resume utilities do not provide that facility.

Use a small matrix of cold, increasing-lift and decreasing-lift initializations
on the problematic states. Keep the original selected state as the reference.
Report CL, CD, CM, transition locations, residuals, work and path differences.
Do not select a path because its drag matches the experiment best.

Acceptance: a documented common policy handles the two known failures or
classifies the reason it cannot; repeatability is measured against the size of
the proposed design benefit. The previous TN808 one-count difference is a
case-specific observation, not a universal noise estimate. Record equivalent
solutions as such only when the relevant state and force comparisons support it.

## 3. Resolve geometry and transition-condition uncertainty

**Source audit completed; physical uncertainty remains, 22 September 2026.**
The report, existing geometry/pressure evidence and a bounded public-source
search are documented with snapshots and checksums in
[`s702_conditions_audit_20260922/FINDINGS.md`](validation/results/s702_conditions_audit_20260922/FINDINGS.md).
No original author coordinate deck or measured hinge-film shape was found.
The reported model tolerance does not specify film-step height or a statistical
geometry uncertainty. The report's tape measurements use different fore/aft
chord conventions and describe displacement effects and sometimes ineffective
tripping. They cannot cleanly validate an ideal zero-thickness trip condition.
No geometry or Ncrit was fitted and no pressure solve was repeated. Priority 3
remains open; the independent bounded wake-setting check below can proceed
without inventing missing geometry or trip dimensions.

Prioritize exact S702 author coordinates and a verified representation of the
hinge/cover region. The current figure reconstruction is adequate for a
qualified benchmark, but fine shape optimization requires tighter geometry
control. Preserve the reconstruction as a separate baseline if better geometry
is obtained. Check reference chord, element placement and flap definition
before comparing forces. If exact coordinates cannot be obtained, state that
limitation and use justified reconstruction variants for sensitivity only;
do not fit geometry to reproduce drag and then claim independent validation.

The report includes fixed-transition experiments, while the solver currently
supports natural transition only. Explicit per-surface forced-trip support could
reduce Ncrit ambiguity. If implemented, preserve free transition upstream and
apply a trip as a maximum laminar extent, with clear chord/surface conventions.
Verify the limit cases and compare with matched fixed-transition measurements.
Ncrit approaching zero is not a substitute for this boundary condition.

Keep the current S702 Ncrit 6.25 calibration frozen when comparing alternative
model implementations. Further tuning against the reserved cases would make
them development data; it would require a new untouched validation set. Values
near 6 and 9 remain useful scenarios, but the interval is not a probabilistic
description of the tunnel.

Acceptance: inputs are traceable and observables beyond total drag improve
under matching conditions. Per-case Ncrit or geometry tuning is not evidence
that a physical mechanism has been corrected.

## 4. Audit drag and wake interaction before extending the closures

**Wake-extent experiment complete and audited, 22 September 2026.**
[`S702_WAKE_EXTENT_PROTOCOL.md`](validation/S702_WAKE_EXTENT_PROTOCOL.md)
declares four 2c cases (flaps 5/10, CL 0.8/1.2, Re 700,000, Mach 0.07,
Ncrit 6.25) against four audited 1c references, plus one fresh difficult 1c
control. The unchanged solver derives wake node counts from extent, so this is
extent plus derived discretization sensitivity, not length at fixed resolution.
All five new cases converged and all nine selected states (including four
archived references) passed full audit. The fresh 1c control repeats drag and
all shared saved arrays exactly. Maximum individual drag sensitivity is
0.00553 count; maximum change of the matched flap delta is 0.00821 count.
This is small relative to the 0.5-count screening threshold and the ~19.38-count
CL-0.8 discrepancy. Longer modeled wakes do not explain that discrepancy in
these conditions. There are 21 requests, two failed requests and six failed
internal continuation points; all are retained. Source/runtime, nonzero-Mach
velocity conversion, Newton states, force sums, geometry, transition/separation
and wake clearance were checked. See
[`s702_wake_extent_analysis/FINDINGS.md`](validation/results/s702_wake_extent_analysis/FINDINGS.md)
for paired plots and tables. Positive centerline clearance does not prove
absent wake/BL interaction, and the Squire–Young station traces are not an
independent drag diagnostic. Launch/source verification and preflight logs are in
`validation/results/s702_wake_extent_launch_20260922/`; frozen manifest SHA-256
is `d79a1add8303d9f5b5de97d18396f7360a5dada3e44f7e88593ab3693c4d6f96`.

**Five-job panel-refinement check complete and audited, 22 September 2026.** The flap-10 aft
post-transition Hk drop reaches about 1.10 per panel at CL 0.8; velocity steps
exceed 0.05 at both lift targets. These are resolution diagnostics, not new
admissibility cutoffs. The previous 340/240 check covered Re=1 million, not
these Re=700,000 cases. The new
[`S702_HOLDOUT_PANELS_PROTOCOL.md`](validation/S702_HOLDOUT_PANELS_PROTOCOL.md)
declares four 340/240 cases at the same flaps/lifts against audited 260/180
references, plus a fresh difficult base-grid control. Wake extent remains 1c,
Ncrit 6.25, Mach 0.07, raw geometry and core equations unchanged. Surface count
also changes derived wake nodes and spacing, so this is operational panel
refinement with that coupling explicit. All five jobs and nine selected-state
audits passed. The base control repeats drag and all shared arrays exactly.
The largest individual CD change is 0.68448 count. Flap-delta changes are
−0.50432 and −1.05171 counts at CL 0.8 and 1.2; the remaining experimental
discrepancies are 18.88073 and 13.49186 counts. The CL-1.2 ranking stays
reversed. Post-transition Hk/velocity steps improve, but those observations
do not establish numerical resolution below the 0.5-count screening threshold.
All four fine cases select cold initialization. The base control still needs
the warm anchor; 19 requests include one failed request and four failed
internal points. See
[`s702_holdout_panels_analysis/FINDINGS.md`](validation/results/s702_holdout_panels_analysis/FINDINGS.md).
Actual solver-contour resampling changes reach 2.61e−5c by a sampled
vertex-to-other-polyline metric; TE endpoints agree to roundoff. No artificial
TE opening or monotone-paneling recovery was activated, so lowering the
inactive minimum-TE-gap threshold is not a useful explanation of these states.
The frozen manifest is
`cd861525d01a02d788bbaddabaacb0c81f9896e153ab2229743190b775fea6f1`.
All 27 targeted runner checks passed locally and in the pinned remote runtime.
The launch, source verification and restart instructions are saved in
`validation/results/s702_holdout_panels_launch_20260922/`.

**Five-job third-grid check complete; conclusion audited 24 September 2026.**
[`S702_PANEL_CONVERGENCE_PROTOCOL.md`](validation/S702_PANEL_CONVERGENCE_PROTOCOL.md)
declares four 440/320 cases at the same conditions, plus an independent
340/240 flap-10/CL-0.8 control. The immediate references are the audited fine
states, with original base references retained. Physics, raw geometry, 1c
extent and numerical incidence centers are fixed. All five jobs and thirteen
selected-state audits passed; the 340/240 control repeats CD and all eighteen
saved array fields exactly. All new states select cold initialization;
18 requests include no failed requests or internal continuation points.

The numerical screen fails: the largest individual step grows from 0.68448
to 0.69746 count; the largest paired step shrinks from 1.05171 to 0.75479 but
remains above 0.5 count. At CL 0.8, the paired step itself grows from −0.50432
to −0.75479. Flap-10/CL-1.2 CD reverses direction between the two steps.
The latest finest-grid flap deltas are −2.59757 and +4.52587 counts at CL 0.8
and 1.2, versus measured −20.72351 and −8.75843. Remaining discrepancies are
18.12595 and 13.28430 counts; the CL-1.2 ranking stays reversed. These observed
steps are not an error bound and do not isolate the discrepancy to a closure.

The largest second-step response (flap 10°, CL 0.8) is chiefly aft drag:
−0.73716 count, partially offset by +0.03971 fore. Local Hk/velocity steps
improve, but integrated drag does not consistently contract. Transition
presence and separation-region counts remain unchanged. Actual contour
resampling changes reach 1.54788e−5c in the second step, while derived wake
counts and first spacing also change. TE endpoints agree to roundoff, with
no artificial opening or paneling recovery. No wake/body intersection appears.
Full tables, plots and reproducible audit are in
[`s702_panel_convergence_analysis/FINDINGS.md`](validation/results/s702_panel_convergence_analysis/FINDINGS.md).

This closes the last automatic grid extension. No fourth grid or optimization
campaign was launched. Three nonuniform grids do not establish an asymptotic
error bound or independently validated design rankings. Thirty-six targeted
analysis/policy/optimization/paneling tests pass; the known baseline far-field
test limitation remains, and the full numerical suite was not rerun.
Frozen manifest:
`b522cead604443e3470eb7cb38df9efdd1863ea15a6bf8cc6440cbf6eeb43833`.
Twenty-six runner checks had passed locally and in the pinned remote runtime
before launch. The jony4 runner finished on 22 September. Launch verification,
logs and preserved restart instructions are in
`validation/results/s702_panel_convergence_launch_20260922/`.

The drag definition is the sum of separate Squire-Young wake estimates. The
wakes are traced in the global inviscid field; their merging and confluence
with other boundary layers are not modeled. These are concrete limitations,
but the existing evidence has not established them as the cause of the S702
ranking discrepancy.

**Discretization separation complete, 1 October:** explicit opt-in source,
surface-count and wake-count/first-step controls preserve production defaults.
The eleven-target study passed its full audit and identified surface sampling
as the largest conditional contribution; see the current-work section above
and [the report](validation/results/s702_discretization_analysis/FINDINGS.md).
Next isolate fore/aft surface sampling using mixed counts at fixed source and
wake controls. Keep the same fixed-lift policy and clearance checks; do not alter
the inactive minimum-TE-gap option as a supposed fix. If an independent
momentum/control-volume drag diagnostic is added, derive its assumptions and
normalization; re-evaluating the same Squire-Young formula is not independent
verification. Surface pressure plus friction requires its own discretization
and wake-consistency assessment before being adopted as a replacement CD.

Only implement displacement-aware wake relaxation, confluence or revised
transition/separation physics when the diagnostic evidence points to a need.
Keep the baseline available under an explicit option and preserve mfoil's
original source. Compare force differences, Cp/BL diagnostics and solver work.

Acceptance: the mechanism-specific test improves and a frozen regression matrix
does not acquire unexplained failures or errors. State separately whether the
change improves S702, other geometries, or only numerical consistency.

## 5. Validate geometry rankings and useful directional derivatives

Use exact, matched geometry variants with independent wind-tunnel or credible
RANS references. The S702 flap tests already supply a useful challenge but are
not independent contour perturbations. The TR664 candidate remains unscored
pending geometry/condition/uncertainty audits; TN808 is a coarse placement
challenge with strong cove-model limitations. Pomeroy's measured placement
variations and the X-57 unpowered alternatives are possible further sources,
subject to obtaining matched geometry and transition assumptions. Source links
and the existing evidence assessment are in `validation/OPTIMIZATION_READINESS.md`.

For a small set of smooth contour modes, compare plus/minus perturbations at
two step sizes and at least two panel densities, keeping physical assumptions
fixed. Inspect changes in actual repaneling, minimum TE opening, transition
topology and monotone-paneling recovery. A finite difference across a change
of numerical branch is not automatically a usable design derivative.

Acceptance should be declared before examining candidate reference results:
compare the sign and magnitude of delta-CD at fixed lift, worst-point error,
coverage and the optimum's location where experimentally resolved. Report
unresolved differences honestly when reference uncertainty is too large. Use
independent contour evidence before claiming validated contour gradients.

## 6. Expand optimization and improve performance after those checks

Retain the existing low-dimensional pattern search as a transparent baseline.
Use multipoint fixed-lift objectives, explicit geometric and moment limits,
separation/supersonic applicability checks, and paired transition scenarios.
Require predicted gains to exceed a justified case-specific combination of
numerical sensitivity and reference resolution. Verify finalists on finer
panels and independent calculations. Preserve baseline retention as a valid
outcome.

Profile before changing the dense solver. Cache only by a complete fingerprint:
actual geometry, Re, Mach, Ncrit/trips, panels, wake settings, TE regularization,
solver code/dependencies, tolerances and initialization policy. A path-sensitive
state must not be reused as if history were irrelevant. Use bounded worker
parallelism with one BLAS thread per worker and incremental durable results.

Optuna is a useful later option for budgeted study management and constrained
derivative-free exploration; it is not currently a dependency. Compare it with
the existing search on the same finite evaluation budget and scenario set.
Do not use a broad search to compensate for incorrect aerodynamic rankings.
An adjoint or automatic-differentiation effort should follow evidence that the
discrete response and model sensitivities are meaningful.

## How to judge progress efficiently

For each change, write a short hypothesis, a discriminating test matrix, the
expected observable, the source/input fingerprint, and a work budget before
execution. Report completed/failed counts, aggregate and worst errors, all six
S702 flap deltas, numerical work and any independent evidence. Preserve the
reference setting and evaluate one material change at a time.

A completed first iteration is the archived-state diagnostic in priority 1.
Priority 2's Mach-0.07 and Mach-zero policy/state audits have now passed.
The operational wake-extent experiment found negligible flap-delta sensitivity.
Both panel refinements improve local steps, but the three-grid screen fails:
latest individual/paired changes reach 0.69746/0.75479 count and some responses
do not contract. The automatic sweep is closed. The completed separation study
identifies surface sampling and first wake spacing as the main conditional
contributions at CL 0.8. The next numerical question separates fore/aft surface
sampling at fixed source/wake controls, as specified above. Preserve priority
3's unresolved geometry and trip limitations before deciding on any closure change.
A new full Ncrit sweep, a global drag
multiplier, a higher-dimensional optimizer or another Mach-zero-only check is
unlikely to resolve the specific gaps already identified without new evidence.
