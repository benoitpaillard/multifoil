# S702 fore/aft surface sampling diagnostic

Declared 1 October 2026 before the new solves. The completed separated study
attributes -0.561280 count of its CL-0.8 paired change to surface sampling;
the aft drag term contributes -0.590387 and the fore term +0.029106. Force
allocation does not identify which element's sampling caused the change.
Test that question with a complete 2×2 surface-count design.

## Fixed conditions and five new targets

Re 700,000; Mach 0.07; Ncrit 6.25; CL 0.8; reference chord 1; flaps 5/10;
wake extension 1c. Keep the source contours at 440/320 panels, with the original
fine wake-count and requested-first-step recipes for each flap. Retain the
unchanged production core and opt-in DiscretizedSolver. The frozen previous
worker performs every new solve and exports source coordinates, controls,
initialization histories, failed requests/internal points and Newton states.

| Cell | Surface fore/aft panels | Flaps | Source of states |
|---|---|---|---|
| B | 340/240 | 5, 10 | Audited previous `contour` states |
| F | 440/240 | 5, 10 | Two new targets |
| A | 340/320 | 5, 10 | Two new targets |
| C | 440/320 | 5, 10 | Audited previous `surface` states |
| C repeat | 440/320 | 10 | One new explicit-control repeat |

The repeat uses the same explicit controls as C; it is not the earlier
automatic production-solver control. Counts above are panels, API counts
are nodes=panels+1. The regularized source and wake recipes are fixed; actual
wake location and curved spacing respond to each fixed-lift state.

Use the original incidence centers, bounds [-6,8] degrees, initial half-width
0.05 degree, at most two brackets and twenty evaluations per bracket, 100
Newton iterations, CL tolerance 2e-5 and residual tolerance 1e-7. Preserve
cold-local, fresh -1-degree warm anchor, then fresh anchor with explicit
center-0.4-degree entry; first success wins. Do not retune transition or select
initialization using drag agreement. Report selected paths. If they differ,
qualify interpretation and assess a matched-path check before claiming robust
factor effects; no automatic extra path or grid campaign is authorized here.

## Analysis fixed before execution

For each flap, report CD, fore/aft wake-drag contributions, CL, incidence, CM,
residual, source/actual geometry, transition/separation, local Hk/velocity
steps, wake count/spacing/clearance and cost. In drag counts, calculate:

* Fore effect at the low aft count: F−B; at the high aft count: C−A.
* Aft effect at the low fore count: A−B; at the high fore count: C−F.
* Interaction: C−F−A+B. Check `(F−B)+(A−B)+interaction = C−B`.
* Apply the same operations to paired flap-10-minus-5 CD, using each cell's
  own flap-5 baseline. Report individual and paired interactions; include
  per-element force contributions without treating them as causal locations.

Require 5/5 new feasible targets, all eight factorial cells, the repeat, full
selected-state audit and unchanged frozen parent/reference hashes before
complete attribution. Require repeat CD agreement with C within 0.01 count
and report exact agreement of every saved array field. Missing or failed cells
cannot be omitted from a passing score. Preserve failed force nulls and all
attempts. The runner reports provisional arithmetic only; full audit follows
completion. The 0.5-count screen is unchanged, and this study supplies neither
an asymptotic error bound nor independent physical validation.

## Experimental and MSES comparison

Keep the existing S702 measured and published-MSES flap deltas in the summary
as context, with explicit Re/M/CL and transition labels. At this condition the
archived experimental delta is -20.723512 counts and published MSES Ncrit-9
delta is +12.651491. MultiFoil remains the frozen Ncrit-6.25 calibration; that
is not a transition-matched code comparison. Published MSES is a third-party
curve, not a new executable run. Compare absolute CD and fixed-lift differences
in the subsequent report; no experimental value determines factor selection.

The broader validation deliverables remain: qualify Pomeroy's tabulated
three-element geometry and measured wake observables, then qualify S103/S414
tested geometries for experiment/MultiFoil/published-MSES comparisons. Record
reference chord, placement, Re/M, transition method/trips, tunnel corrections,
raw versus digitized values and uncertainty before admission. Use Ncrit 9 as
the unchanged default on new airfoils; declared sensitivity is separate from
calibration. Report coverage before accuracy, common-support errors and flap
or placement ranking; missing states never improve a score. Total-drag
agreement cannot substitute for pressure/transition/wake checks. Do not import
the S702 fitted Ncrit as an experimentally established condition on new cases.

## Execution and continuation

New isolated remote root:
`/home/jony/multifoil-runs/s702-surface-factorial-20261001`.
Use at most four worker processes, each with BLAS/OpenMP threads=1, in pinned
Python 3.12 / NumPy 2.3.5 / SciPy 1.17.0 / Matplotlib 3.10.8. Prepare and hash
the sources before running. Preserve interrupted attempts; resume only under
the identical manifest/runtime. No packages or production-core files change.

```bash
# Remote run/resume from its isolated root:
.venv/bin/python validation/run_s702_surface_factorial.py --workers 4
# Local progress (no numerical environment required):
python3 validation/monitor_s702_surface_factorial.py --watch
```

After completion retrieve and audit all nine selected states, assess the
conditional effects/interactions and return to experimental source qualification.
Do not launch a higher grid, closure change or optimizer automatically.
