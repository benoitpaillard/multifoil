# S414 fore-wake information audit — 3 October 2026

## Question and scope

Identify which fore-wake quantities enter the aft mean-flow and boundary-layer
equations, and which proposed wake/transition mechanisms the current model can
represent. Use the completed S414 reconstruction study and its audited states.
Keep the numerical core, geometry, Ncrit 9, initialization and all saved states
unchanged. No nonlinear flow solve, fitted coefficient, closure option or
automatic simulation extension belongs to this check.

## Evidence and checks fixed before the probes

1. Trace the actual source code from the four station variables (momentum
   thickness, displacement thickness, amplification/root shear, edge velocity)
   through wake initialization, source coupling and aft transition. Distinguish
   a turbulent wake closure from measured turbulence intensity or intermittency.
2. Restore the nominal CL 0.284 and 1.041 states from their saved arrays. Rebuild
   the same wake and source operators and verify the coupled residual remains
   below 1e−7 before any diagnostic perturbation. Verify rebuilt foil/wake
   coordinates, branch indices and stagnation distances against the archive.
3. Inspect the Jacobian blocks linking fore-wake variables to aft BL equations
   and aft velocity-coupling equations. Independently perturb each of the four
   fore-wake variable fields by ±0.1% and ±0.01%, with all other state values,
   geometry, Ncrit and branch flags held fixed. Compare central differences
   with the Jacobian at both step sizes. A zero cross-element block must be
   exactly zero; nonzero relative directional-derivative discrepancies must
   be below 1e−7 (denominator at least 1e−12). Require a nonzero own-fore BL
   response as a positive control. Restore every saved state array exactly.
4. Resolve the unchanged mean-flow coupling into fore/aft surface/wake source
   columns. Their contributions must sum to the implemented displacement
   correction and reproduce the archived coupling residual. Report values
   around aft recovery without presenting them as independently physical drag
   components or a converged response to removing a source.
5. From all eight non-repeat archived states, measure shortest distances from
   aft upper-surface points to the fore-wake polyline. Interpolate saved wake
   displacement/momentum thickness and root shear along the nearest segment.
   Compare geometric distances with displacement thickness as a scale only:
   displacement thickness is not wake half-width, and the model does not supply
   a measured turbulence profile or a basis for an overlap/contact threshold.
6. Confront the trace with the 2010 report's oil-flow/acoustic observations and
   qualified 2014 interpretation. The low-lift CL-0.284 state is inside the drag
   bucket, not its lower-limit horn. Preserve that distinction. Do not extend
   a qualitative horn explanation to every pressure discrepancy. Keep 2010
   envelope and 2014 full-e^N MSES results distinct, and do not claim this code
   audit proves the undocumented internals of the published MSES setup.

## Interpretation and reproducibility

These are partial derivatives of the discrete residual at a fixed converged
state, not new equilibrated aerodynamic solutions. A zero direct shear-to-aft
block does not rule out indirect effects through fore thickness and the coupled
mean flow after solving. Perturbed states must never be scored as converged
forces, new drag predictions or measured turbulence conditions.

Save input/code hashes, source-function locations, runtime, baseline residuals,
each finite-difference check, source-column sums, geometry/thickness tables and
plots. The two operator checks are restartable with identical code/input/runtime
guards. Reuse the completed full state/force audit rather than rerunning all
cases. Preserve original results and receipts.

Success is a tested description of the available information channels and their
limitations, with a defensible next diagnostic. Missing physics is not thereby
proven to cause the measured discrepancy. No empirical wake-turbulence-to-Ncrit
mapping or closure correction is to be invented from the two lift conditions.
