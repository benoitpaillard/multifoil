# Working on MultiFoil

Usage hold (5 October, 18:09 UTC): the user reported excessive agent usage.
The r3 local monitor was stopped before any postanalysis model review; the
remote CFD job continues. See `orchestration/runs/s414-nose-20261005-r3/monitor_hold.json`.
Do not restart automatic AI reviews without addressing their budget. This
supersedes the monitor-running statements below.

Current solver default (4 October 2026): `CoupledSolver` now uses
`amplification_jacobian="corrected_af"`; explicit `legacy` mode remains available.
This supersedes historical opt-in-only descriptions below. See
`DEVELOPMENT_PLAN.md` for scope. Runtime verification of the integration is
complete (stage V1): see
`validation/results/s414_default_integration_20261004/FINDINGS.md`. It verifies
the implementation only. The complete test suite now shows 38 historical
source-guard failures plus the known far-field limitation; these guards are
correct and must not be relaxed. Keep historical source/hash guards intact and
use frozen snapshots for archived-study replay.

Latest decision (5 October): read
`validation/results/s414_post_verification_review_20261005/FINDINGS.md` and the
decision supplement in `validation/S414_V5_NOTES.md`. The review matches 153
archived output hashes. V5 prescribed the archived displacement and failed the
5% momentum-thickness screen at CL 0.056 (9.3% difference); it cannot exclude
closure/interaction effects on coupled bubble formation. Next selected is a
six-case local coupled nose-resolution diagnostic, then an independent coupled
layer if needed, with exact test/MSES inputs pursued in parallel. The user has
authorized the orchestration pilot. Its prepared revision and live status are
in `orchestration/runs/s414-nose-20261005-r3/README.md` and `state.json`.
Independent prelaunch review passed; the detached pilot launched on jony4 at
18:06 UTC on 5 October, with a local monitor for retrieval, audit and final
review. Read the live state before acting; do not submit another job. The
protocol is `validation/S414_NOSE_RESOLUTION_PROTOCOL.md`. Scientific results
remain pending; retain the model/Ncrit/geometry and all original study receipts.

The current staged S414 investigation is in `DEVELOPMENT_PLAN.md`, section
"Verification plan for the S414 nose-pressure and transition discrepancy".
It covers reference coverage, default integration, solution repeatability,
independent pressure/geometry checks, transition physics and remedy validation.
Execution started 4 October; the live stage table is
`validation/S414_VERIFICATION_STATUS.md`. V0 (evidence matrix) and V1 are
complete. V2 tranche 1 (fixed-incidence paths, frozen protocol
`validation/S414_TRANSITION_BOUNDARY_PROTOCOL.md`) is complete and audited; see
`validation/results/s414_transition_boundary_analysis_20261004/FINDINGS.md`.
Legacy endpoint controls matched the archive identically; 19/20 path points
converge and all 25 saved states pass fresh audits. Both incidence directions
and both derivative modes agree wherever comparable; no second converged state
was found. The cubic geometry's lift is non-monotone and brackets CL 0.175
three times: keep that matched-lift target missing. Repeatability holds on the
sampled paths only. Tranche 2 (frozen protocol
`validation/S414_TRANSITION_BOUNDARY_T2_PROTOCOL.md`) is also complete and
audited, and V2 is closed with its 42-point allowance used; see
`validation/results/s414_transition_boundary_t2_analysis_20261004/FINDINGS.md`.
Restored-state continuation reproduces live continuation exactly; 12/16 new
evaluations converge and all 14 saved states pass audit. No incidence has two
converged states, including two cold starts (the cubic one at −2.54388° is a
recorded, user-approved deviation from the declared rule). The nominal
transition change lies between −2.694° and −2.569° (−38.87 counts, CL
0.155–0.169); the cubic one between −2.556° and −2.519°, with four failed
continuation attempts inside. Interval interiors are unresolved; do not sample
them further without a separately bounded branch-continuation protocol.
V4's geometry and inviscid part is complete; see
`validation/results/s414_nose_geometry_20261004/FINDINGS.md`. The fore-nose
turning reversal is created by the PCHIP interpolant next to a nearly flat
figure-derived vertex (x/c 0.01449), not by paneling. A convex arc-chain
reconstruction through the same 83 vertices (`validation/s414_nose_geometry.py`,
frozen input `fore_input_convex_nose.dat`) removes it and is admitted; it is not
the tested contour. Its inviscid lower-nose pressure differs by up to 0.056 Cp
(declared threshold 0.05) and has more, not less, suction behind the peak.
V4's viscous check of that geometry (frozen protocol
`validation/S414_NOSE_VISCOUS_PROTOCOL.md`) is complete and audited; see
`validation/results/s414_nose_viscous_analysis_20261004/FINDINGS.md`. V4 is
closed: the reversal is **not** the cause. With it removed, the CL 0.056 error is
+56.00 counts (nominal +55.35), CL 0.284 is unchanged (+6.19), and CL 0.175
moves to early transition (+44.14 against +4.09). Three reconstructions within
0.11 mm place the transition-driven drag change at least 0.1° apart, so do not
treat the nominal CL 0.175 state as validated or use the CL 0.056 → 0.175 drag
change as a score on this figure-derived geometry. Robust across reconstructions:
early fore-lower transition at CL ≤ 0.155 and 0.20–0.34 excess suction at the
second fore-lower tap. Keep the nominal baseline; no geometry, closure or Ncrit
change is justified.
V5 tranche 1 (frozen protocol `validation/S414_LAMINAR_REPLAY_PROTOCOL.md`, with
one recorded, not pre-approved numerical amendment) is complete and verified; see
`validation/results/s414_laminar_replay_20261004/FINDINGS.md`. The
prescribed-pressure replay is **closed as a method**: with the model's own
pressure and the same upstream state the separation station has two roots
(archived bubble or attached layer; N at x/c 0.0208 is 8.2 or 1.6 at CL 0.056),
so its sampling and seed controls cannot pass (10 of 54 do). Do not relaunch a
direct replay in any form. A qualified independent finite-difference laminar
solver (`validation/s414_laminar_fd.py`) agrees with the inherited closure up to
separation under identical pressure; nothing is established about the separated
layer or amplification. Measured taps alone keep the layer attached to x/c
0.034; superposing the model's between-tap shape separates it near 0.006, so the
unmeasured pressure upstream of and between the first two taps decides. The
CL 0.284 "separation retained" outcome is conditional and is not discrimination.
V5 tranche 2 (frozen protocol `validation/S414_BUBBLE_STABILITY_PROTOCOL.md`, no
deviation) is complete and verified; see
`validation/results/s414_bubble_stability_20261005/FINDINGS.md`. On the archived
bubble's displacement thickness an independent inverse finite-difference layer
(`validation/s414_bubble_fd.py`) reproduces the pressure within 0.016 Cp, and
Orr–Sommerfeld stability (`validation/s414_stability.py`,
`validation/s414_amplification.py`) gives 18–19 % more amplification than the
inherited law: N 9.93 against 8.45 at x/c 0.0208 at CL 0.056, 5.98 against 5.04
at CL 0.175. This supports growth conditional on the archived bubble, without
excluding coupled-formation effects. One quantified closure difference: momentum thickness 9 %
high at the end of the deep CL 0.056 bubble. On similar profiles the inherited
correlation grows faster than linear stability (slope ratio 0.91–0.66 for H
4–11.3); do not use that as evidence about the bubble. **V5 execution is complete:
the growth screen passes, but the low-lift closure screen fails and coupled
bubble formation remains unresolved.** See the post-plan review above. V6
has nothing to validate. A missing-data request is drafted, not sent, in
`validation/results/s414_data_request_20261005/`; sending needs the user's
approval. Suggestions outside the plan are parked in
`validation/S414_V5_NOTES.md`; the user asked to continue the written plan and
keep suggestions and findings stored there. The authorized nose-resolution
pilot is now active; see the live orchestration state above.
Do not edit any existing
`*.py` in the repository root or `validation/`: the launched manifests hash them.
V3 is complete with its gate open
but bounded: Hess–Smith reproduces the nose suction peak within 0.017 Cp, the
predeclared agreement criterion is not met, and the residual offset follows
Hess–Smith's blunt-base treatment; see
`validation/results/s414_formulation_check_20261004/FINDINGS.md`. V6 is not
started. Follow the plan's gates and bounded runs rather than historical
next-step paragraphs.

Read `CODEX_HANDOVER.md` first for background, then `DEVELOPMENT_PLAN.md` for
the current state (4 October 2026). The latest archived-state diagnosis is in
`validation/results/s702_drag_mechanism/FINDINGS.md`; calibration and flap
comparisons remain in `validation/results/s702_ncrit_analysis/` and
`validation/results/s702_flap_effect/`. The fixed-lift repeatability matrix is
complete and audited in `validation/results/s702_repeatability_analysis/FINDINGS.md`.
The seven-job combined initialization policy is complete and audited in
`validation/results/design_combined_policy_analysis/FINDINGS.md`. The wake-extent
study is complete and audited in `validation/results/s702_wake_extent_analysis/FINDINGS.md`.
The first panel refinement is audited in
`validation/results/s702_holdout_panels_analysis/FINDINGS.md`. The third-grid
study is complete and all thirteen states passed audit, but the numerical
convergence screen failed; see
`validation/results/s702_panel_convergence_analysis/FINDINGS.md`.
The automatic grid-extension sequence is closed. The eleven-target diagnostic
separating contour/surface/wake controls is complete; all fifteen selected
states passed audit. Surface sampling dominates the conditional paired change,
followed by first wake spacing. See
`validation/results/s702_discretization_analysis/FINDINGS.md` and the bounded
fore/aft sampling follow-up in `DEVELOPMENT_PLAN.md`. That follow-up is complete:
all five targets and nine selected states pass, with an exactly repeated control.
Aft sampling dominates even in comparisons with matching selected paths;
one mixed case's different path qualifies the interaction estimate. See
`validation/results/s702_surface_factorial_analysis/FINDINGS.md` for the audit
and experiment/MSES plots. The Pomeroy source audit
recovered all 243 coordinates but found unresolved placement/gap inconsistency;
see `validation/results/pomeroy_geometry_qualification_20261001/FINDINGS.md`.
Do not admit that geometry as a matched benchmark or substitute the 2015 setup.
The earlier physical diagnosis is in
`validation/results/s702_physics_diagnosis_20261001/FINDINGS.md`: transition/bubble
behavior and physical hinge/wake inputs are leading suspects, not proven causes.
Four S702 flap-5/10 pressure snapshots are now qualified; see
`validation/results/s702_pressure_qualification_20261002/FINDINGS.md`.
All eight matched-lift pressure targets (Ncrit 9/6.25) completed on jony4 and
passed full saved-state/residual audit. See the latest completed pressure
diagnosis and plots in
`validation/results/s702_matched_pressure_analysis_20261002/FINDINGS.md`.
Flap-5 pressure/separation evolution remains mismatched; good total drag can
coexist with wrong local pressure. The foremost aft tap lies 0.00084946c ahead
of the raw reconstructed contour. The source cross-check in
`validation/results/s702_geometry_reconciliation_20261002/FINDINGS.md` reproduces
the contour in four pressure figures but leaves this inconsistency unresolved.
No geometry correction is justified. All ten fore/aft transition-isolation jobs
completed and all 18 new/archived states pass audit with actual per-element Ncrit;
both controls repeat exactly. Latest completed S702 findings:
`validation/results/s702_transition_isolation_analysis_20261002/FINDINGS.md`.
Fore transition drives the large low-lift flap-5 calibration response; aft
threshold changes still fail to recover both flap-5 pressure profiles. Strong
coupling appears at low-lift flap 10. Close the Ncrit diagnostic without tuning.
The adapted auditor is `validation/audit_s702_transition_state.py`; old receipts
remain intact. S414's nominal modified test outline and two 2010 published
MSES/pressure comparisons are now extracted; see
`validation/results/s414_geometry_qualification_20261002/FINDINGS.md`.
Repeated report outlines match, but exact geometry admission remains open:
public `.pts` files sharpen the trailing edges, three tap centers lack branch
support, and plotted/table tap coordinates differ. Keep 2010 envelope-MSES and
2014 full-e^N results distinct. Matching Figure-13 drag is now extracted; MSES
bias is +5.33 / −13.89 counts at CL 0.284 / 1.041. See
`validation/results/s414_drag_reference_20261002/FINDINGS.md`.
All nine bounded S414 reconstruction jobs completed and passed full saved-state
and fresh residual/force/Cp/Cf/BL audits, with an exact repeat. All 336 retrieved
files match remote hashes. See the latest findings and graphs in
`validation/results/s414_reconstruction_analysis_20261003/FINDINGS.md`.
At CL 0.284 / 1.041, nominal MultiFoil drag errors are +6.19 / −10.06 counts.
MultiFoil and MSES both delay the strongest sampled aft recovery by one tap
interval; all declared geometry variants and both tap-coordinate conventions
retain that mismatch. Small tested geometry changes do not explain it, but
exact tested geometry and grid convergence remain unestablished. Ncrit stays 9;
aft ±0.13 mm shifts are sensitivities, not measured rigging uncertainties.
The nine-job study is closed. The fore-wake information
audit is complete; see `validation/results/s414_wake_transfer_20261003/FINDINGS.md`.
Both nominal states pass 16 residual/Jacobian probes: fore-wake displacement
affects the aft mean flow, but there is no direct wake-disturbance transfer to
aft transition. Every wake uses a turbulent closure even after laminar fore
surfaces. Neither fact proves the cause of the pressure error; the low-lift
case lies inside the bucket, away from the authors' boundary-horn interpretation.
The free/fore-only/both-tripped source qualification is now complete; see
`validation/results/s414_trip_qualification_20261003/FINDINGS.md`. Six experimental
and six published-MSES drag curves are extracted at Re 0.7/1m. The author reports
excess trip displacement and ineffective fore-upper grit at low lift; no paired
tripped aft Cp or bubble endpoints are provided. Admit total-force context, not
a clean wake-disturbance benchmark. At Re 1m/CL 0.284, measured/MSES fore-trip
penalties are +37.61/+67.49 counts. Four free-transition boundary-pressure cases
are extracted (CL −0.072/0.056/0.175/0.284); pressure-copy offsets and a failed
strict precision check are retained. The four-case free-transition comparison
is complete and fully audited: all four saved states and fresh residual/force/
Cp/Cf/BL checks pass, with 18 exact control arrays, 22 converged requests and
38 passing targeted tests. All 238 retrieved files and 100 frozen inputs match
hashes. See `validation/results/s414_lower_bucket_analysis_20261003/FINDINGS.md`
for six graphs and the 348-supported-tap comparison. At CL 0.056 the drag error
is +55.35 counts versus published MSES +8.21; the modeled fore transition remains
at x/c 0.02247, moving to 0.76928 at CL 0.175. The −46.14-count computed drag
change between them is −45.53 fore/−0.61 aft, while experiment gains 5.12 counts.
Aft strongest-rise error reverses at 0.175; no uniform bubble shift is supported.
The fore-nose diagnostic is complete; see
`validation/results/s414_fore_nose_20261003/FINDINGS.md`. Saved amplification
balances/roots pass, but the inherited `get_damp` Jacobian misses a factor of
two in an exponential derivative (18/72 component checks fail). A diagnostic-only
correction passes; production code is unchanged. Excess nose suction already
appears in the inviscid calculation at saved incidences; geometry sampling adds
a local turning reversal. Neither finding proves the drag-error cause.
The separate `amplification_jacobian.AmplificationJacobianSolver` now supplies
explicit legacy/corrected modes, preserving solver/vendor files. All 32
assembled derivative checks and 50 isolated preflight tests pass. The four
CL 0.056/0.284 legacy/corrected jobs are completed and fully audited; see
`validation/results/s414_jacobian_analysis_20261003/FINDINGS.md`. All 987 files
match hashes; all four fresh state audits pass and 36 control arrays repeat
exactly. Drag changes by less than 4e−9 counts and Cp by less than 7.4e−10;
the defect does not explain the aerodynamic discrepancy in these cases.
The measured-incidence follow-up is complete and audited; see
`validation/results/s414_incidence_analysis_20261003/FINDINGS.md`. All four
states pass fresh audits, 36 control arrays repeat exactly and 18 tests pass.
At −3.07°/−1.04°, CL is 0.10723/0.35918 versus experiment 0.056/0.284;
overall Cp agreement worsens and early fore transition persists. Incidence
offset alone is insufficient. Fixed-angle drag does not replace matched-lift
validation. The identical-tap pressure budget is now complete; see
`validation/results/s414_pressure_budget_20261003/FINDINGS.md`. Supported
pressure differences account for 88%/92% of the excess normal force at measured
angle; fore/aft load errors largely cancel at matched lift, also in published
MSES at nominal CL 0.28. Twelve tests and independent quadrature checks pass.
Endpoint bounds are conditional and do not establish experimental uncertainty.
The same-angle inviscid/viscous decomposition is now complete; see
`validation/results/s414_inviscid_budget_20261003/FINDINGS.md`. Viscosity reduces
the supported-tap load gap by 62%/50% and improves all-surface Cp RMSE, but
fore-lower load errors cross negative while upper errors remain positive.
Eighteen tests pass; no nonlinear solve was needed. Comparing an inviscid model
with viscous experiment does not prove an inviscid defect. The independent
potential-flow/Hess–Smith candidate is deferred: focus is restored to the
CL 0.056 → 0.175 transition/drag discrepancy. The bounded prescribed-pressure
replay is complete but inconclusive: 4/8 profiles reach their requested
pressure; only 3/8 remain below Ncrit. Both sampling controls and both
peak-retaining experimental reconstructions fail continuation. Even completed
roots do not establish physical branch selection. See
`validation/results/s414_prescribed_pressure_20261004/FINDINGS.md`.
The four-job fore-only reconstruction check is finished and audited: 3/4 targets
converge, both controls repeat 36 arrays exactly, and all three fresh selected-state
audits pass. All 430 retrieved files match; 13 analysis/runner tests pass. See
`validation/results/s414_fore_reconstruction_analysis_20261004/FINDINGS.md`.
At CL 0.056 the cubic fore increases drag by 0.671 count and retains early
transition. Its CL 0.175 target fails; the paired drag change is unknown.
At the same −2.56888° angle, nominal/cubic transition differs 0.75011/0.08282
and drag differs +37.53 counts, but lifts differ: do not label this matched-lift.
Eleven additional request exports are checked separately, not full Newton audits.
Next recommendation is bounded fixed-angle continuation in both directions,
with nominal/cubic endpoint controls and the corrected Jacobian explicitly
identified, to diagnose the transition-sensitive interval. It is now V2
tranche 1 (launched 4 October; see above). Keep the legacy comparison, Re 1m, Mach 0.10, Ncrit 9,
340/240 panels and 1c wake; the aft input is unchanged. No extension is automatic.
The original matched-lift basis remains in the frozen protocol
`validation/S414_LOWER_BUCKET_PROTOCOL.md`. Exact geometry, grid convergence
and solution uniqueness remain unestablished.
Published MSES drag
does not cover −0.072; no new MSES pressure comparison exists for the three
new cases. No automatic extension follows these four jobs.
Do not retune Ncrit or extend the automatic grid campaign.
Preserve the baseline; no geometry/closure
correction or causal proof follows from the S702 counterfactual.
S103/S414 tested geometry remains the independent-validation path. Do not extend
the automatic numerical campaign. No MSES executable is installed.
Numerical convergence remains unestablished; do not certify small design gains.
A failed far-field preflight test was reproduced unchanged in the original
baseline and retained as a known limitation. The literature review and next
experimental candidates are in `validation/REFERENCES.md`.

The goal is a reproducible, open multi-element airfoil model useful for
incompressible shape optimization. It is currently a research/design-screening
prototype. Baseline drag calibration has improved; independent contour-design
rankings remain unvalidated. Published MSES curves are references, not runs of
an installed MSES executable. No proprietary MSES or Viiflow binary is bundled.

Keep scientific evidence traceable:

- Preserve original data, completed studies, failed targets and source snapshots.
  Put new calculations in new directories under `validation/results/`.
- Never turn a failed solve into a finite force value or a stall prediction.
  Report coverage before accuracy; incomplete candidates cannot win by omission.
- Compare designs at the same lift, reference chord, physical assumptions and
  numerical protocol. Use each scenario's own baseline.
- The default Ncrit remains 9. The frozen S702 fit is 6.25; 6 is nearly equivalent
  and better for absolute/flap-difference RMSE. Keep Ncrit fixed or use declared
  sensitivity scenarios during geometry optimization. Do not optimize it as a
  geometry variable or retune it per case using experimental drag.
- Preserve geometry, source/code hashes, solver settings, initialization history,
  wall time and convergence status with every new run. Check existing manifest
  guards before reusing results; do not bypass them after a source change.
- Keep vendored notices and historical source snapshots intact. Put proposed
  closure changes behind explicit model options and compare with the baseline.

Use targeted `rg` searches; avoid loading the entire results tree into context.
Start with archived outputs and the fast audit in the handover. Run flow solves
only for a concrete unresolved question. Set BLAS/OpenMP threads to 1 when using
multiple worker processes; the complete S702 sweep can take hours.

After core changes, run relevant existing numerical tests and the affected
benchmark cases. Distinguish implementation verification, calibration,
independent validation and model limitations in the resulting report.
The user prefers complete runnable scripts or a complete updated archive for
code changes, together with concise explanations and plots of scientific results.

`handover/verify_snapshot.py` checks the delivered snapshot before development.
Expected hash differences after intentional edits are not aerodynamic failures;
retain the original manifest and produce a new receipt for a new release.
