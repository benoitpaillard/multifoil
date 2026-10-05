# S414 convex-nose viscous check (V4, viscous part), 4 October 2026

Frozen before execution; the runner is `validation/run_s414_nose_viscous.py`.
It follows the decision recorded in
`validation/results/s414_nose_geometry_20261004/FINDINGS.md`: the admitted
convex-nose reconstruction changes the inviscid lower-nose pressure by up to
0.056 Cp, above the 0.05 threshold declared beforehand.

## Question

With the turning reversal removed and nothing else changed, do the modeled
fore-lower transition and drag at CL 0.056 and near CL 0.175 change, and does
the transition-driven drag change move in incidence? This tests a
reconstruction sensitivity. It cannot show that either contour is the tested
one, and a lower drag would not by itself be evidence of anything.

## Fixed inputs

Re 1,000,000, Mach 0.10, Ncrit 9 on both elements, 340/240 panels, 1c wake,
100 Newton iterations per attempt, residual tolerance 1e−7, explicit
`corrected_af`, pinned Python 3.12.13 environment, two workers on jony4.

Geometry: the fore solver input is read from the frozen, hash-checked file
`fore_input_convex_nose.dat` (341 points; all 83 reported vertices kept; no
negative turn in the final nose polygon). The aft input is the nominal one.
The geometry is not regenerated on the server.

## Jobs

Phase 1, control gate: nominal geometry at CL 0.056 with the archived ordered
fixed-lift policy and center. It must reproduce the archived corrected state
(`s414_jacobian_20261003/cl0056_corrected_af`) with identical saved arrays and
the same selected path. Otherwise nothing else starts.

Phase 2, convex-nose geometry:

- **CL 0.056** and **CL 0.284** (control): the same ordered fixed-lift policy,
  centers, bounds, tolerances and budgets as the archived nominal jobs.
- **CL 0.175**: the same policy and the archived nominal center, with its
  monotone-lift guard unchanged. Its result is admitted in the analysis only
  if the fixed-incidence paths below show one monotone bracket of CL 0.175;
  otherwise the target is reported missing, whatever this job returns.
- **Two fixed-incidence paths** on the eight incidences already sampled for
  the other geometries (three in the nominal archived bracket, five in the
  original bracket, about −2.945° to −2.469°). Sixteen path points; full
  states saved at each. The increasing path starts cold at the lowest
  incidence. A local preparation check found that this geometry does not
  converge from a cold start at the highest incidence (100 iterations), so the
  decreasing path starts cold at −1.0°, the warm anchor of the archived
  policy, and reaches its first point by continuation; every internal step is
  recorded.

## Rules

- Ordered initialization, first success wins, lift alone selects the incidence;
  drag never selects a path, a seed or a branch.
- A failed target or path point has no force. A path whose first point fails is
  not started. After a mid-path failure the path continues from the last
  converged state in the same direction, and that base is recorded.
- Every evaluation is exported; selected fixed-lift states and all converged
  path points save full Newton states. Completion is execution only: each
  saved state needs the independent audit before interpretation.
- An interrupted job restarts from its beginning; the partial attempt is kept.

## Reading of outcomes

Same-incidence comparisons use the V2 screens (|ΔCL| ≤ 2e−5, |ΔCD| ≤ 0.01
count, |ΔCp| ≤ 1e−4, transition within 1e−4 c, same topology) between the two
directions. Differences between geometries are reported as sensitivities at
the same incidence and, where a target is admitted, at the same lift. They are
compared with the experiment and the published 2010 MSES drag at CL 0.056,
0.175 and 0.284 using the unchanged references.

A changed transition position or drag shows sensitivity to this reconstruction
choice. It supports a causal role of the reversal only if the response is
repeatable in both directions, the pressure change is the declared one, and it
is consistent with the supported taps. An unchanged response deprioritizes the
reversal. Ncrit stays 9. Nothing extends automatically after these six jobs.
