# S414 fixed-incidence transition-boundary paths (V2, tranche 1), 4 October 2026

Stage V2 of the verification plan in `DEVELOPMENT_PLAN.md`. Frozen before
execution; the runner is `validation/run_s414_transition_boundary.py`.

## Question

Near −2.57° the nominal and fore-only source-cubic reconstructions give fore
transition at x/c 0.750 and 0.083 with a +37.5-count drag difference, and the
cubic CL 0.175 target failed under the fixed-lift monotonicity guard. Is the
modeled transition/drag response through this narrow incidence interval
**repeatable on sampled paths**, independent of continuation direction, and
unchanged by the corrected amplification derivative? This is a numerical
diagnostic at fixed incidence. It is not matched-lift drag validation and
cannot demonstrate physical hysteresis.

## Fixed inputs

Re 1,000,000, Mach 0.10, Ncrit 9 on both elements, reference chord 1, finite
trailing edges, 340/240 panels, 1c wake, 100 Newton iterations per attempt,
residual tolerance 1e−7, solver continuation settings 1°/0.25°, minimum TE gap
1e−5. Geometries: `nominal`, and `fore_source_cubic` (fore contour interpolator
only; aft input identical). No fixed-lift search, no Ncrit change, no grid change.

The bracket ends and the start incidence are read at full precision from the
verified archived exports of `s414_fore_reconstruction_20261004`
(evaluations 0 and 1 of `cl0175_nominal` and `cl0175_fore_source_cubic`):
about −2.56888° and −2.46888°, started at about −2.51888°. They are not retyped.
The five path angles are the two archived ends and the three interior points
`lower + (upper − lower)·k/4`, k = 1, 2, 3.

## Jobs (two concurrent workers on jony4, pinned Python 3.12.13 environment)

Phase 1, explicit `legacy`, one job per geometry — **endpoint controls**:
fresh solver started at the archived start incidence; request the lower end,
then the upper end. This is the recorded initialization of the archived
evaluations. Gate: polar rows, Newton histories, continuation trails, transition
and force diagnostics, and all eight surface/wake export files per geometry must
be **identical** to the archive. If either control differs, no corrected path
starts.

Phase 2, explicit `corrected_af`, two jobs per geometry — **paths**:

- *increasing*: fresh solver, same start; request the five angles upward. Its
  first point is the corrected lower endpoint with the recorded initialization.
- *decreasing*: a separate fresh solver, same start; request the lower end as a
  labeled endpoint-seed evaluation, then the five angles downward from the upper
  end. Its first path point is the corrected upper endpoint with the recorded
  initialization (through the lower end in one step). The last point returns to
  the lower end after the round trip.

Requested points: 20 path points (4 paths × 5 angles) plus 4 legacy control
points. The two endpoint-seed evaluations are initialization, reported and
saved, not path points.

## Rules

- Every converged requested evaluation saves both elements' full Newton states,
  exports, transition-cell audit and wake discretization. Every evaluation
  records how it was reached (cold, continuation from which incidence), all
  internal continuation steps, the Newton history and wall time.
- A nonconverged point has no force. The path continues from the last converged
  state in the same direction, and that base is recorded.
- If a direction's own endpoint (or its seed) fails, that path is **not
  started**: remaining points stay unevaluated. No substitute cold start is used
  in this tranche.
- An interrupted job restarts from its beginning; the partial attempt is moved
  to `interrupted_attempts/`. A resumed job is not a restored Newton checkpoint.
- Completion is execution only. Each saved state then needs the independent
  residual/force/Cp/Cf/transition/wake audit before interpretation.

## Predeclared equivalence screens and reading of outcomes

Same-angle comparisons between directions, modes or controls use: |ΔCL| ≤ 2e−5,
|ΔCD| ≤ 0.01 count, |ΔCp| ≤ 1e−4 at common stations, transition position within
1e−4 c on every surface, and identical transition/separation topology. These
are numerical decision thresholds, not experimental uncertainties. Exceeding
them triggers investigation, not adjustment.

Agreement supports repeatability only on the sampled paths. Two
residual-converged states at one angle establish numerical multiplicity, not
physical hysteresis. A fold or failure is reported as such; the fixed-lift
solver is not forced through it, and CL 0.175 stays missing unless an admitted
branch brackets it. Same-angle results at different lifts are never scored as
matched-lift drag.

## Tranche 2 selection rule (declared now, frozen and launched separately)

Allowance: at most four additional incidence locations per geometry, each
evaluated in both directions, and one independent cold check per geometry
(42 requested points in total with tranche 1). Selection uses only audited
tranche-1 quantities and never drag level:

1. For each of the four intervals between adjacent sampled angles compute
   `d` = true if at either end the two directions both converge but differ
   beyond a screen above, or exactly one direction converges; and `t` = the
   largest change of fore-lower transition position (x/total chord) across the
   interval over directions with both ends converged, with `t` = ∞ when one end
   converges and the other fails on the same direction.
2. An interval is eligible if `d` is true or `t` > 0.05. Rank by `d` first,
   then `t` descending, then lower angle first.
3. New locations are midpoints in rank order; if allowance remains, the
   quarter points of the top-ranked interval (lower first).
4. If no interval is eligible for a geometry (for example nominal transition
   late throughout), apply the same rule to the bracket located from the
   archived nominal CL 0.056 → 0.175 request history, with the same allowance.
5. Each new location is reached by continuation from the audited saved state
   of the adjacent tranche-1 point on that direction. The cold check is a fresh
   solver started directly at the angle with the largest path difference, else
   the midpoint of the top-ranked interval.

Achieved angular resolution is reported. Nothing extends automatically.
