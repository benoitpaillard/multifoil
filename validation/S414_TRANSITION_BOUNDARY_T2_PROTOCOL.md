# S414 fixed-incidence transition-boundary paths (V2, tranche 2), 4 October 2026

Second and last tranche of stage V2 in `DEVELOPMENT_PLAN.md`. Frozen before
execution; the runner is `validation/run_s414_transition_boundary_t2.py`.
It applies the selection rule declared in
`validation/S414_TRANSITION_BOUNDARY_PROTOCOL.md` to the audited tranche-1
results, with one recorded deviation approved by the user.

## Question

Tranche 1 found the same state from both incidence directions at every
comparable angle, a 37.9-count drag change on the fore-only cubic geometry
through converged intermediate transition positions, non-monotone cubic lift,
and no transition change for the nominal geometry inside the bracket. Tranche 2
asks, within the remaining allowance: (i) does the cubic response stay
direction-independent when sampled at up to four times finer spacing; (ii) does
a cold start inside the transition change reach the same state; (iii) where in
its archived bracket does the nominal geometry's own transition change lie, and
is it direction-independent there.

This remains a numerical diagnostic at fixed incidence. It is not matched-lift
validation, cannot show physical hysteresis and does not prove uniqueness.

## Fixed inputs

Identical to tranche 1: Re 1,000,000, Mach 0.10, Ncrit 9 on both elements,
340/240 panels, 1c wake, 100 Newton iterations per attempt, residual tolerance
1e−7, solver continuation settings 1°/0.25°, explicit `corrected_af`, pinned
Python 3.12.13 environment, two workers on jony4. All incidences are read at
full precision from audited files; none is retyped.

## Requested evaluations (18; with tranche 1's 24 this uses the 42-point allowance)

**Restore controls, two.** Continuation here starts from saved states restored
into a new solver, which is not a bitwise copy of a live solver. For each
geometry, restore the tranche-1 increasing-path state at the lower bracket end
and continue to the next sampled angle. Gate: the result must match the
tranche-1 live result at that angle within the declared screens. If either
control fails, nothing else starts.

**Cubic geometry, four new incidences in both directions (eight).** From
`tranche2_selection.json`: about −2.55638°, −2.53763°, −2.53138° and −2.50638°.
Each is one continuation step from the nearest audited, converged tranche-1
state on that direction's own path: the state below for increasing, the state
above for decreasing. Where the adjacent decreasing-path point failed in
tranche 1 (−2.51888°), the next converged decreasing-path state above it is
used, as tranche 1 itself did.

**Nominal geometry, fallback bracket, three new incidences in both directions
(six).** No tranche-1 interval qualifies, so the declared rule moves to the
bracket from the archived history: the audited corrected state at the measured
incidence −3.07° (`s414_incidence_20261003/cl0056_measured_alpha`, fore-lower
transition 0.025 of total chord) and the tranche-1 state at the lower bracket
end (0.750). The single interval yields its midpoint and its two quarter
points; the fourth location of the allowance is not used. Increasing
evaluations continue from the restored −3.07° state, decreasing ones from the
restored tranche-1 decreasing-path state at the lower bracket end.

**Cold checks, two.** A new solver started directly at the incidence, one
attempt.

- Nominal: the midpoint of its fallback bracket, as the rule's second clause gives.
- Cubic: **deviation.** The declared rule names −2.51888°, the angle with a
  path difference. Tranche 1 already holds three failed cold starts there (one
  legacy, two corrected), so repeating it would add nothing. With the user's
  approval of 4 October the check is moved to the tranche-1 angle −2.54388°,
  inside the transition change, where both directions converged to the same
  state. The choice does not use drag level. It was made after seeing tranche
  1 and is reported as a deviation wherever this check is cited.

## Rules

- Every restored base state must reproduce its saved geometry, wake and station
  coordinates and a fresh residual below 1e−7 before continuation; otherwise
  that evaluation is not started and is reported as such.
- Each evaluation is independent: one process, one base state, one requested
  incidence. Internal continuation sub-steps, Newton history, wall time and the
  identity of the base state are recorded. A restored warm start is never
  described as cold.
- A nonconverged evaluation has no force and no saved state.
- Every converged evaluation saves both elements' full Newton states, exports,
  transition-cell audit and wake discretization, and needs the independent
  saved-state audit before interpretation.
- An interrupted job restarts from its beginning; the partial attempt is kept.

## Reading of outcomes

Same-angle comparisons use the tranche-1 screens: |ΔCL| ≤ 2e−5, |ΔCD| ≤ 0.01
count, |ΔCp| ≤ 1e−4, transition within 1e−4 c, identical topology. Exceeding
them triggers investigation, not adjustment.

Agreement supports repeatability on the sampled paths only. Two
residual-converged states at one incidence would establish numerical
multiplicity, not physical hysteresis. A cold check that converges to the path
state supports, and does not prove, uniqueness there; one that fails is a
convergence finding. CL 0.175 on the cubic geometry stays missing unless a
branch rule is separately declared and admitted. Achieved angular resolution
is reported.

**V2 ends with this tranche.** There is no further automatic refinement or
extension of the allowance.
