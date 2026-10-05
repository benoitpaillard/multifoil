# S414 isolated amplification-Jacobian comparison — 3 October 2026

Question: does correcting the missing factor of two in `get_damp` change
convergence, solution selection or the CL 0.056 drag discrepancy? The control
is CL 0.284. This is an implementation verification, not a fitted physical
closure or proof of improved experimental validation.

Use `amplification_jacobian.AmplificationJacobianSolver` with explicit modes
`legacy` and `corrected_af`. The latter adds only the missing derivative of
`0.1 exp(-20/(Hk-1))`; amplification values and physical equations are unchanged.
The standard solver and vendor files remain byte-identical to the archive.
A dispatcher reads the copied per-model parameter; untagged models remain
legacy. Record the mode on every flow evaluation and in every job/result.

Before launch, require finite-difference checks of the entire assembled coupled
Jacobian at the saved CL 0.056 and 0.284 states. Probe the largest-Hk active
laminar station and the last laminar station, all four state variables, at
relative steps 1e−4 and 1e−5. The maximum directional discrepancy divided by
the largest corrected derivative component (floor 1e−10) must be below 2e−5.
No residual rows may be omitted. Residual values may differ by at most 1e−10
due to transition-root roundoff; legacy residuals/Jacobians and all saved
state arrays must restore exactly. Preserve failed checks and stop on failure.
Also test growth values, inactive/ramp/full activation, Hk floor, copied
parameters, mode isolation and invalid settings.

Run exactly four jobs on jony4 with two single-threaded workers:

1. CL 0.056 legacy, repeating the archived `cl0056_nominal` job.
2. CL 0.284 legacy, repeating the archived `cl028_control` job.
3. CL 0.056 corrected_af, with the same declared numerical protocol.
4. CL 0.284 corrected_af, with the same declared numerical protocol.

Both fresh legacy controls must select the same initialization path, reproduce
all 18 saved arrays exactly and differ by no more than 0.01 drag count from
their respective archived states before either corrected case runs. Do not
bypass a failed control. Keep nominal S414 geometry, 340/240 panels, Re 1m,
Mach 0.10, Ncrit 9, 1c wake, 100 Newton iterations, 1e−7 residual and 2e−5
lift tolerances, archived incidence centers, fixed-lift bounds and ordered
cold-local/warm-anchor/warm-entry initialization paths. Use the archived
Python 3.12.13 and package versions. This matches initialization policies;
individual Newton trajectories and evaluated incidences may then diverge.

Preserve every attempted evaluation, convergence history, failure, selected
state, source hash and wall time. Completed failures count as completed work,
not successful comparisons. They cannot be silently retried or represented
by finite forces. Interrupted folders are preserved before restarting;
exclusive locks and immutable manifests guard resume.

After completion, retrieve and hash-check results. Audit selected states with
the correct per-model derivative setting, including forces, pressure, Cf,
transition and residuals. Compare coverage before accuracy; report matched-lift
CD and CL residuals, incidence, selected paths, Newton counts, failed requests,
fore/aft drag contributions, fore transition and pressure differences. A
same-state residual/derivative comparison separates implementation correctness
from solution changes. Available published MSES drag and CL 0.284 pressure
remain external references; no MSES executable is available.

Do not declare the aerodynamic discrepancy fixed merely because derivatives
pass, a solve converges or one drag error decreases. No Ncrit retuning,
geometry change, broader grid campaign or automatic additional target follows.
Numerical convergence, exact tested geometry and uniqueness remain unestablished.
