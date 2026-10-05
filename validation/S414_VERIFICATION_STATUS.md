# S414 verification plan: stage status

**Usage hold, 5 October at 18:09 UTC:** automatic local monitoring is stopped
to prevent further AI review usage; the remote simulation continues. See
[hold receipt](../orchestration/runs/s414-nose-20261005-r3/monitor_hold.json).
This supersedes the monitor-running statements below.

Status of the plan in `DEVELOPMENT_PLAN.md`, section "Verification plan for the
S414 nose-pressure and transition discrepancy". Last updated 5 October 2026,
after the post-plan evidence review. States: planned · running · execution finished · audit passed ·
interpreted · blocked. A finished script is not a physical validation.

**Current decision:** [post-plan review](results/s414_post_verification_review_20261005/FINDINGS.md).
V5's independent layer prescribed archived displacement, and its CL 0.056
momentum-thickness screen failed (9.3% versus 5%). It supports amplification
on that bubble; it does not exclude closure/interaction effects on coupled
bubble formation. Next selected: six local coupled-resolution cases, then an
independently coupled layer if needed, alongside exact test/MSES input recovery.
The local-resolution diagnostic is now an authorized orchestration pilot:
[revision 3 and live status](../orchestration/runs/s414-nose-20261005-r3/README.md),
[protocol](S414_NOSE_RESOLUTION_PROTOCOL.md). Preparation passed 43 short checks
and real linear/source-coupling checks; independent prelaunch review passed.
The detached pilot launched on jony4 at 18:06 UTC on 5 October. The local monitor
is running; baseline reproduction and scientific outcomes remain pending.
Original studies/receipts are unchanged.

| Stage | State | Evidence | What it does and does not establish |
|---|---|---|---|
| V0 evidence matrix | interpreted | [findings](results/s414_evidence_matrix_20261004/FINDINGS.md), [verification](results/s414_evidence_matrix_20261004/verification.json) | 95 reference rows, 27 explicit unknowns. Strict matched-geometry validation stays **conditional**. |
| V1 integrated default | gate passed | [findings](results/s414_default_integration_20261004/FINDINGS.md), [receipt](results/s414_default_integration_20261004/verification.json) | Implementation verified; no aerodynamic claim. |
| V2 tranche 1, fixed-incidence paths | **interpreted** (executed, audit passed) | [findings](results/s414_transition_boundary_analysis_20261004/FINDINGS.md), [verification](results/s414_transition_boundary_analysis_20261004/verification.json), [execution record](results/s414_transition_boundary_launch_20261004/README.md) | Both directions agree at 9 of 9 comparable angles; cubic lift is non-monotone and CL 0.175 is bracketed three times, so that target stays missing. Sampled paths only. |
| V2 tranche 2, refinement and cold checks | **interpreted** (executed, audit passed). **V2 closed** | [findings](results/s414_transition_boundary_t2_analysis_20261004/FINDINGS.md), [verification](results/s414_transition_boundary_t2_analysis_20261004/verification.json), [execution record](results/s414_transition_boundary_t2_launch_20261004/README.md) | No second state at any incidence, including two cold starts. Nominal change located between −2.694° and −2.569°. Interior of both intervals unresolved; cubic CL 0.175 stays missing. One recorded, user-approved deviation (cubic cold-check incidence). |
| V3 cross-formulation pressure | interpreted; **gate open, bounded** | [findings](results/s414_formulation_check_20261004/FINDINGS.md), [verification](results/s414_formulation_check_20261004/verification.json) | Nose peak reproduced by a second method within 0.017 Cp; predeclared agreement criterion not met (0/10). |
| V4 turning artefact: geometry trace and inviscid pressure | interpreted | [findings](results/s414_nose_geometry_20261004/FINDINGS.md), [verification](results/s414_nose_geometry_20261004/verification.json) | The reversal is made by the interpolant next to one nearly flat reported vertex, not by paneling. A convex reconstruction through the same vertices removes it (admitted). Inviscid lower-nose pressure changes by up to 0.056 Cp, just above the declared 0.05: material by the frozen rule. It does not reduce the nose suction. |
| V4 turning artefact: viscous check of the convex nose | **interpreted** (executed, audit passed). **V4 closed** | [findings](results/s414_nose_viscous_analysis_20261004/FINDINGS.md), [verification](results/s414_nose_viscous_analysis_20261004/verification.json), [execution record](results/s414_nose_viscous_launch_20261004/README.md) | Removing the reversal leaves the CL 0.056 error (+56.0 vs +55.4 counts) and the CL 0.284 control unchanged, and moves CL 0.175 to early transition (+44.1 vs +4.1). The reversal is not the cause. Transition position between CL 0.155 and 0.28 depends on sub-0.11 mm reconstruction choices. |
| V5 tranche 1, attached-laminar replay controls | **interpreted** (executed twice, verified). Direct replay **closed as a method** | [findings](results/s414_laminar_replay_20261004/FINDINGS.md), [verification](results/s414_laminar_replay_20261004/verification.json), [branch diagnostic](results/s414_direct_branch_diagnostic_20261004/summary.json), [protocol](S414_LAMINAR_REPLAY_PROTOCOL.md) | The prescribed-pressure station problem has two roots at separation, so the sampling and seed controls cannot pass (10 of 54 do). An independent laminar solver agrees with the closure up to separation under identical pressure. Measured taps do not decide separation. One recorded numerical deviation, not approved beforehand (first execution kept). |
| V5 tranche 2, independent layer and linear stability through the bubble | **interpreted** (executed, numerically verified); growth screen passes, low-lift closure screen fails | [findings](results/s414_bubble_stability_20261005/FINDINGS.md), [verification](results/s414_bubble_stability_20261005/verification.json), [protocol](S414_BUBBLE_STABILITY_PROTOCOL.md) | On archived displacement, independent pressure agrees within 0.016 Cp and stability gives 18–19% more amplification. The 9.3% momentum-thickness difference fails the 5% screen at CL 0.056. Conditional cross-check, not exclusion of coupled closure effects. |
| V5 outcome | amplification supported on the prescribed bubble; coupled formation and experimental input remain unresolved | [review](results/s414_post_verification_review_20261005/FINDINGS.md), [notes and decisions](S414_V5_NOTES.md) | No independent coupled calculation was made. Local coupled resolution is the next selected diagnostic. |
| V6 remedy validation | not started; **no remedy justified** | – | Await a demonstrated cause from the selected numerical checks or recovered experimental/MSES inputs. |

## Counts

| Stage | Jobs completed / required | Requested points | Converged | Failed | Interrupted |
|---|---|---|---|---|---|
| V2 tranche 1 (audited) | 6 / 6 | 20 path + 4 control + 2 seed | 19 path, 4 control, 2 seed; 25 saved states pass audit | 1 path point (cubic, decreasing, −2.51888°) | 0 |
| V2 tranche 2 (audited) | 18 / 18 | 14 continuation + 2 cold + 2 restore control | 10 continuation, 2 cold, 2 control; 14 saved states pass audit | 4 continuation (cubic, both directions, −2.53763° and −2.53138°) | 0 |
| V3 | 30 / 30 linear solves | – | all residuals below 4.5e−15 | – | – |
| V4 viscous (audited) | 6 / 6 | 16 path + 3 targets + 1 control | 15 path, 3 targets, 1 control; 19 saved states pass audit | 1 path point (convex, increasing, −2.51888°) | 0 |

| V5 tranche 1 (verified; second execution, local) | 1 / 1 | 114 finite-difference + 121 closure marches | 34 attached to x/c 0.034, 80 ended at laminar separation; closure agrees in topology in 112 of 114 | 0 paired marches; 1 archive-control march ends at the archived turbulent station | 0 |

| V5 tranche 2 (verified, local) | 1 / 1 | 11 inverse boundary-layer marches + 7 stability envelopes | 11 marches completed; 7 envelopes with no wave lost | 0 | 0 |

V2 and V4 counts are from the audited studies. All 457, 504 and 495 retrieved files match remote hashes. V2 used its 42-point allowance; nothing extends automatically.

## Decision table (original issues)

| Issue | Classification today | Missing evidence |
|---|---|---|
| Amplification-derivative defect | demonstrated implementation defect corrected; integration verified (V1); compared V2 endpoints agree between modes | None for this fix; global solution uniqueness is a separate question |
| Repeatability of the transition-boundary response | numerical sensitivity bounded on the sampled paths: one state per incidence from 14 direction pairs and 2 cold starts, both derivative modes at the endpoints | The interior of the two intervals (0.025° cubic, 0.125° nominal). Needs a separately bounded branch-continuation method. Uniqueness off the sampled incidences is not established |
| Matched-lift drag at CL 0.175, cubic geometry | unresolved: the target is met in three intervals | A declared branch rule, or the branch-continuation method above |
| Modeled −46-count drag drop, CL 0.056 → 0.175 (nominal) | located; not robust to reconstruction (V4) | Coupled formation of the low-lift bubble remains unresolved after V5 |
| Nose pressure: formulation | numerical sensitivity bounded in the nose (±0.05 Cp), not eliminated | A second method with a matched blunt-base model |
| Nose turning reversal | rejected as the cause: origin identified (interpolant overshoot at reported vertex 45), removed in an admitted alternative, low-lift error unchanged | None for this question |
| Transition position between CL 0.155 and 0.28 | numerical sensitivity shown, not bounded: three reconstructions within 0.11 mm place the transition-driven drag change at least 0.1° apart; the nominal CL 0.175 agreement is not robust | As-built fore-nose coordinates. Until then CL 0.175 is not a validation point and the CL 0.056 → 0.175 drag change is not a reliable score |
| Excess nose suction behind the first tap (0.20–0.34 Cp) | unresolved; present with every reconstruction and with early or late transition. The measured taps alone do not separate the laminar layer; adding the model's between-tap shape does | Nose pressure resolved upstream of and between the first two taps; as-built coordinates |
| Prescribed-pressure replay as a discriminating method | rejected for a demonstrated numerical reason: two roots at the separation station under identical pressure and upstream state; controls cannot pass | None; do not relaunch. An interactive method replaces it |
| Boundary-layer/transition closure, attached laminar part | numerical cross-check passed: an independent finite-difference solution agrees up to separation under identical pressure (not a validation against measurement) | None for the attached part |
| Boundary-layer/transition closure, separated layer | conditional cross-check: pressure within 0.016 Cp on archived displacement; independent momentum thickness 9.3% lower at CL 0.056, outside the 5% screen | Local coupled resolution and an independent calculation with displacement free; coupled effect of this difference is unknown |
| Amplification law in the bubble | numerical cross-check passed: linear stability gives 18–19 % more amplification than the inherited law on independent profiles; the law does not over-amplify | Assumptions shared by both (e^N with Ncrit 9, parallel two-dimensional linear theory) are untested; case-specific transition measurements do not exist |
| Cause of the early fore-lower transition at CL ≤ 0.155 | unresolved: derivative fix and reversal removal do not repair it; amplification is supported conditional on the bubble. Coupled formation/resolution, geometry/pressure input and shared e^N assumptions remain candidates | Local coupled-resolution check, independently coupled layer, and exact geometry/resolved nose-pressure/transition evidence |
| Experimental reference | not implicated | Uncertainty of Cp and of the 2010 force coefficients is not reported |

## Next actions, in order

1. Monitor the launched six-case local coupled nose-resolution pilot and review
   its audited outcome. Follow the live revision-3 state above: exact original
   controls gate baseline/half/quarter spacing at two archived incidences.
   Five archived requests per case preserve initialization; 30 requests maximum,
   two workers, 5400 seconds including preflight/cleanup. No automatic extension.
   Retrieval, audit, plots and independent postanalysis review follow automatically.
2. If the discrepancy survives that check, qualify an independently coupled
   nose layer with displacement and edge velocity free. Prescribed-displacement
   replay cannot close this question.
3. In parallel, seek exact experimental/MSES geometry and low-lift outputs.
   The missing-data request draft is in
   [results/s414_data_request_20261005](results/s414_data_request_20261005/REQUEST_DRAFT.md).
   It has not been sent; sending needs the user's approval.
4. Stage V6 has no admitted remedy to validate. Remaining unselected suggestions
   stay in [S414_V5_NOTES.md](S414_V5_NOTES.md); no fitting or direct replay.
5. Optional and separately bounded: a branch-continuation method for the
   interval interiors found in V2; locating the convex-nose transition change
   between −2.469° and −1.672°.

The detached nose-resolution pilot is active on jony4. Its ordinary local monitor
handles retrieval, audit, plots and independent review; inspect the live state
before resuming and never resubmit an unresolved remote identity.
