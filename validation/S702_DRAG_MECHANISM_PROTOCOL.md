# Archived S702 drag mechanism diagnostic

Declared 14 September 2026, before computing the element comparison. This
implements priority 1 of `DEVELOPMENT_PLAN.md`; it changes no solver equations.

## Question and test

Determine whether the remaining 5-to-10-degree flap drag discrepancy is
associated primarily with the modeled fore or aft element, and which surface
transition, separation and trailing-edge thickness changes accompany it.
Association does not establish the missing physical mechanism: the measurements
resolve total drag, not an experimental decomposition by element or surface.

Use the frozen Ncrit 9 and 6.25 results at all 21 operating conditions (42 flow
states). The primary matrix has 16 states: flaps 5/10, CL 0.8/1.2, Re 1,000,000
at M 0.1 and Re 700,000 at M 0.07, both Ncrit values. CL 0.5 is a control.
Keep each Reynolds/Mach/Ncrit scenario paired at its own achieved target lift.
Retain the original cold/warm retry qualifications and initialization paths.

Budget: zero new flow solves; one archived-data report and targeted audit tests.
No drag, Ncrit, geometry, or reference-data fitting. No closure change.

## Audits and outputs

Resolve original/reused results through `flow_directory()` and existing study
manifest guards. Tie the selected result to its polar index, metadata and actual
surface/wake exports. Check convergence, recorded residual/lift tolerances,
physical settings, geometry and all consumed file hashes. Record the original
source hashes and solver/dependency versions separately from this reporting
environment. Preserve failed attempts; never replace their forces with a value.

Check total CD against the sum of element wake CD with 1e-12 absolute tolerance
(less than 1e-8 drag count), and reproduce every column of the six-row archived
`s702_flap_5_to_10.csv` with 1e-8-count tolerance. Neither is an aerodynamic
accuracy tolerance. Friction is a separate diagnostic, never added to wake CD.

Tabulate element forces, surface transition/separation and TE states, wake
endpoints, actual selected-state paths and initialization history. Plot paired
Cp, Cf, edge speed, Hk, theta and displacement thickness for fore/aft upper/lower
surfaces, marking transition and Cf-negative intervals. Include achieved CL and
incidence. Plot wake evolution and the element contributions to each flap delta.
Thickness uses the assembly reference chord; surface x uses local element chord;
wake x uses assembly chord. A null transition means laminar to the TE, not a
measured transition at x/c = 1. Do not mix in the flap-0/CL-1.04 pressure taps.

Acceptance is a reproducible diagnostic with full coverage, all six deltas,
explicit modeled associations and limitations. If no physical cause is resolved,
state that and nominate the smallest discriminating next experiment. This report
is neither new aerodynamic validation nor evidence for validated shape gradients.
