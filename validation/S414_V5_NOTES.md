# S414 stage V5: findings so far and suggestions kept for later

A running notebook for stage V5 of the verification plan in
`DEVELOPMENT_PLAN.md`. Findings are summaries with pointers; the dated
`FINDINGS.md` files and receipts remain the record. Suggestions are ideas that
are **not** part of the frozen plan. They are kept here so they are not lost;
none is scheduled, and none may start without its own bounded protocol.
Started 4 October 2026; last updated 5 October 2026, after the post-plan review.
The decision supplement below supersedes the earlier suggestion ordering and
overstrong exclusion of the separated-layer closure. Original study receipts
and their reported numerical results remain unchanged.

## Post-plan decision, 5 October 2026

The [review](results/s414_post_verification_review_20261005/FINDINGS.md) checks
153 archived output hashes and recomputes the V5 screens: momentum thickness
misses the declared 5% criterion by 9.3% at CL 0.056. The independent layer
prescribed archived displacement. Therefore closure/interaction effects on
**coupled bubble formation** remain unresolved, even though amplification on
that bubble is supported. Do not summarize V5 as ruling the closure out.

**Next selected diagnostic: six fixed-angle local-resolution cases.** This is a
decision/protocol outline, not an implemented runner or a launched campaign.
Use nominal geometry, corrected_af, Re 1m, Mach 0.10, Ncrit 9, existing residual
tolerance, and the archived initialization policies. Take the exact incidences
from the audited nominal CL 0.056 and CL 0.284 records (approximately −3.4590153°
and −1.6724273°). At each angle compare baseline, half-spacing and quarter-spacing
of fore-lower panels in the nose region through x/c 0.04. No CL 0.175 search,
geometry alternative or Ncrit scenario is included.

Before execution, implement the experiment in new files and freeze a complete
run deck with exact source hashes and station lists. A pure subdivision of the
archived final straight panels can preserve its polygon exactly; retain every
original vertex and bypass subsequent repaneling that would change that contour.
Confirm that the solver actually retains the prescribed nodes. Preserve all
other surface nodes, finite TE faces, element placement, reference chord, wake
node counts, first spacing and length. Record wake coordinates; if the operator
changes their traced positions materially, quantify that change or isolate it
before attributing the response solely to the nose spacing. The independent
source-to-velocity assembly must honor the refined indices. If these controls
cannot be met, revise the implementation before launch rather than treating
a changed geometry as refinement. This qualifies the frozen polygon only, not
the exact manufactured or an assumed smooth contour.

The two baseline controls must reproduce the selected archived states under
their recorded initialization/runtime before refined cases are admitted. Audit
each converged state freshly. Report geometry and inviscid Cp changes, the
coupled pressure profile, Cf, theta, displacement, Hk, N, transition, per-element
drag and all failures. Pressure comparisons use common locations; N is compared
only on the common pretransition laminar domain. Do not reinterpret turbulent
shear as N. These are fixed-incidence mechanism comparisons; any later accuracy
claim requires matching lift again in a separately declared tranche.

Proposed numerical screens, to freeze before execution: the half-to-quarter
change must contract relative to baseline-to-half and be below 0.01 Cp on the
common nose domain, 2% theta/displacement, 0.1 N and 1 drag count; transition
must have the same topology and shift by less than one finest local cell
(report its physical size). These are operational screens against a ~55-count
error, not experimental uncertainty or proof of asymptotic convergence. If a
quantity is already below numerical noise, record that explicitly rather than
requiring a meaningless contraction ratio. Failure leaves local resolution
unqualified; no automatic fourth level or wider sweep follows.

**Conditional next step:** if that bubble/drag response is stable, qualify an
independently coupled finite-difference nose layer with both displacement and
edge speed free, not supplied from the archive. Check the interaction mapping
and conservative interfaces to upstream/downstream and transition equations.
Use an identity-replacement control for the coupling machinery, then compare
the two same-angle states. Keep the integral baseline and shared assumptions
explicit. A local patch with the rest of the flow frozen is only a conditional
diagnostic; it cannot establish full coupled drag or exclude global interaction.
Prepare a separate bounded implementation/verification protocol before running.

**Parallel data priority:** exact test/MSES geometry, rigging and 2010 MSES input
and low-lift output, then as-built nose measurements and resolved pressure/
transition observations. The request is drafted, not sent. No pressure fit or
Ncrit adjustment is a remedy. No new simulation is running as part of this review.

## Findings

### Tranche 1, attached-laminar replay controls (4 October, verified)

Record: [findings](results/s414_laminar_replay_20261004/FINDINGS.md),
[verification](results/s414_laminar_replay_20261004/verification.json),
[branch diagnostic](results/s414_direct_branch_diagnostic_20261004/summary.json),
[protocol](S414_LAMINAR_REPLAY_PROTOCOL.md).

1. **The prescribed-pressure replay is closed as a method.** With the model's
   own edge velocity and the same upstream state, the separation station has
   two roots of the inherited equations (shape factor 4.21 or 3.67 at CL 0.056;
   4.83 or 3.34 at CL 0.175). Downstream, N at x/c 0.0208 is 8.2 or 1.6, and
   5.1 or 2.1. The earlier native control passed only because Newton started
   from the archive. Do not relaunch a direct replay in any form.
2. **Attached-layer closure cross-checked.** A new finite-difference laminar
   solver (`s414_laminar_fd.py`), qualified on Blasius, Hiemenz, Howarth and
   cylinder flows, agrees with the inherited momentum/shape closure up to
   separation under identical pressure: momentum thickness within 1 %,
   separation within 0.0007 of chord.
3. **The closure runs slightly ahead of the exact laminar solution near the
   nose separation.** Just before the bubble its shape factor is 0.1–0.2 higher
   (3.36 against 3.15 at CL 0.056). On Howarth and cylinder flows it separates
   later instead. Inside the declared screen; not a demonstrated defect.
4. **The model's own nose pressure is marginal for laminar separation at all
   three lifts.** Marched under the archived pressure, both methods stay
   attached where the coupled state has a bubble (CL 0.056, 0.175), and at
   CL 0.284 they separate or not depending on the interpolant between nodes.
   The archived pressure already contains the bubble's displacement relief, so
   this does not show that the bubble is wrong.
5. **Measured taps do not decide separation.** Taps-only reconstructions keep
   the laminar layer attached to x/c 0.034 at all three lifts; superposing the
   model's between-tap shape separates it near 0.006. 10 of 54 sampling
   controls pass; nothing is admitted at CL 0.056 or 0.175. The CL 0.284
   "separation retained" outcome is conditional and is not discrimination.
6. **Tap pattern.** The measured first tap (x/c 0.00181) has 0.27/0.20/0.17
   more suction than the model at CL 0.056/0.175/0.284, the second (0.00838)
   has 0.20/0.34/0.39 less. Read together they are compatible with a measured
   suction peak at or upstream of the first tap; this is an observation, not a
   result, because nothing is measured there.
7. **Numerics worth remembering.** The centred Keller box gives a sawtooth wall
   shear at kinks of the edge velocity (index up to 1.4); backward differences
   remove it (0.009). The earlier replay's seed treatment (an abrupt blend over
   0.00016 of chord) separates the layer on the spot when the target has less
   suction than the base. A 25 % thicker layer at x/c 0.00165 flips the direct
   outcome at both low lifts.

### Tranche 2, independent layer and linear stability through the bubble (5 October, verified)

Record: [findings](results/s414_bubble_stability_20261005/FINDINGS.md),
[verification](results/s414_bubble_stability_20261005/verification.json),
[protocol](S414_BUBBLE_STABILITY_PROTOCOL.md). No deviation.

1. **Given the modeled bubble, linear stability also reaches N ≈ 9.** On the
   archived displacement thickness, Orr–Sommerfeld N at x/c 0.0208 is 9.93 at
   CL 0.056 (inherited law on the same layer 8.45, archive 8.23) and 5.98 at
   CL 0.175 (5.04, 5.13). Ratios 1.18 and 1.19; sensitivities 1.17–1.19.
2. **The separated-layer closure reproduces the independent pressure** within
   0.016 Cp (CL 0.056) and 0.009 (CL 0.175).
3. **One quantified closure difference:** momentum thickness 9.3 % high at the
   end of the deep CL 0.056 bubble (shape factor 10.1 against 11.2). The
   reversed-flow approximation explains at most 2.2 % of it. It changes the
   inherited law's N by less than 0.25.
4. **Wall shear in reversed flow** is weaker in the closure (cf·Re_θ about
   −0.06) than in the independent layer (down to −0.10 and −0.13).
5. **On similar profiles the inherited correlation grows faster than linear
   stability** at Re_θ 100–150: envelope slope ratio 0.91 (H 4.0), 0.84 (6.4),
   0.75 (8.8), 0.66 (11.3). In the bubble the order is reversed. I did not
   separate the two likely reasons (non-similar profiles; momentum thickness
   that hardly grows, so one band of frequencies stays amplified).
6. **CL 0.175 sits close to the threshold under either law.** Its archived N
   peaks at 7.06 downstream of the tranche's domain; with the 1.19 ratio it
   would approach 8.4. Extrapolation, not a result.
7. **V5 outcome.** A closure limitation is not supported. The remaining
   candidates are the pressure input that creates the bubble and the
   assumptions shared by every e^N calculation (Ncrit 9, parallel
   two-dimensional linear theory).

Reusable tools from this tranche: `s414_bubble_fd.py` (inverse laminar layer
with reversed flow; `similar_profile` traces the Falkner–Skan family through
separation), `s414_stability.py` (Orr–Sommerfeld, temporal and spatial),
`s414_amplification.py` (e^N along a sequence of profiles).

## Suggestions kept for later (not scheduled)

1. **Tap-matched coupled counterfactual.** Perturb the inviscid fore-nose
   velocity until the *interacted* model pressure meets the measured taps, at
   fixed incidence, with model-to-model controls (taps of one archived state as
   the target for another). It needs coupled solves on jony4, a subclass of the
   solver in a new file and a bounded count. Caveat found in tranche 1: the
   answer will again depend on the assumed between-tap shape, so it brackets
   rather than decides unless nose pressure is resolved.
2. **Interactive finite-difference nose layer.** Replace the integral closure
   on the fore-lower nose by the finite-difference laminar equations, coupled
   through the solver's own source-influence matrix, everything else frozen at
   the archive. It would show whether exact laminar equations produce the same
   bubble under the same inviscid flow. More work than suggestion 1 and more
   decisive about the closure.
3. **Missing-data request, contents** (preparation only; sending needs the
   user's approval): as-built fore-nose coordinates; lower-nose pressure
   resolved upstream of the first tap and between the first two taps
   (x/c 0 to 0.01); case-specific transition, separation or reattachment
   evidence on the fore-lower surface at CL 0.056 and 0.175 (oil flow, hot
   film, stethoscope); the matching MSES input and output decks.
4. **Use the tap pattern as a geometry clue.** If the measured peak sits
   upstream of x/c 0.0018, a family of nose shapes or a small leading-edge
   droop could be screened inviscidly against the first three taps at all four
   lifts before any viscous run. This is a geometry hypothesis and must stay
   independent of drag.
5. **Resolve the marginality numerically.** The bubble spans six archived nodes.
   A local panel refinement of the fore nose on the same continuous contour
   would show whether the two-root station and the coupled bubble survive
   refinement. This belongs with the open grid-convergence question and would
   need a bounded schedule; the automatic grid campaign stays closed.
6. **Branch continuation for V2's unresolved intervals** (already listed as
   optional in the status file): pseudo-arclength continuation in incidence
   through the transition change, instead of fixed-incidence restarts.
7. **Finish the amplification comparison to the end of the CL 0.175 bubble.**
   Tranche 2 stopped at x/c 0.0208 by protocol. Extending the inverse layer and
   the stability envelope to reattachment (x/c about 0.03) would give linear
   stability's peak N for the state that sits near the threshold. Small, local,
   needs its own bounded protocol because inverse mode is ill-conditioned once
   the layer reattaches.
8. **Separate the two effects behind the similar-profile result.** Compare the
   local Orr–Sommerfeld growth rate per momentum thickness of bubble profiles
   and similar profiles at equal shape factor and Reynolds number. If the
   inherited correlation's slope is too high at large H and low Reynolds
   number, that matters for longer bubbles than this one (S702 flap cases).
   This would be a closure study with its own regression set, not an S414 fix.
9. **Do not tune Ncrit to this case.** If someone proposes it: linear stability
   on the model's own bubble supports the inherited law, so a lower growth or a
   higher Ncrit would be a fit to drag, which the project rules exclude.

## Open questions carried forward

- Does the modeled bubble exist with exact laminar equations under the same
  inviscid flow (suggestion 2), or is it a property of the integral closure?
- Is the amplification inside the bubble consistent with linear stability
  theory on independent velocity profiles? Answered in tranche 2: yes, linear
  stability gives 18–19 % more.
- Where is the measured suction peak, and how strong is it?
