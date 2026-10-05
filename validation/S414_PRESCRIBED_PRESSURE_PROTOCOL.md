# S414 fore-nose prescribed-pressure diagnostic, 4 October 2026

Question: does the resolved fore-lower pressure input account for the large
amplification difference between CL 0.056 and 0.175? The matched-lift drag
comparison remains the objective. The independent potential-flow comparison
is deferred; this diagnostic does not change or replace either archived state.

Use the audited nominal lower-bucket states, Re 1m, M 0.1, Ncrit 9. Retain the
saved geometry, streamwise distances and upstream boundary-layer seed. Reuse
the original momentum, shape and amplification residuals with the explicit
`corrected_af` Jacobian option (unchanged residual equations).

The domain starts at the last physical lower-surface node upstream of the first
pressure tap (x/c 0.00181). It ends at the last common native node before either
archived solution becomes turbulent, x/c about 0.02083. This contains most of
the CL 0.056 pretransition amplification but excludes its mixed transition cell.
This is a conditional laminar boundary-layer diagnostic: no total drag, downstream
transition location, or fully coupled aerodynamic solution will be inferred.

For each CL, compare four pressure inputs:

1. Native model pressure, with a perturbed initial guess, must restore the saved
   laminar states (theta/delta-star relative error <2e-5; N error <2e-5).
2. Model pressure linearly interpolated at the same experimental table taps.
   This controls for loss of the unresolved nose peak between taps.
3. Experimental pressure with the identical linear tap interpolation.
4. Native model pressure plus the linearly interpolated experiment-minus-model
   tap differences. This retains the model's unresolved curvature between taps
   and is an alternative conditional reconstruction, not measured pressure.

All modified profiles join the same saved pressure at the upstream seed. No
extrapolation downstream of the supported taps is permitted. Cp is inverted
using the solver's own Karman–Tsien convention; the boundary-layer equations
perform their normal velocity correction. The seed is model dependent and is
not experimental boundary-layer data. The pressure snapshots match CL, not
incidence; use that qualification throughout.

Continue pressure from native to target in steps at most 0.05; halve failed
steps down to 0.05/4096, with a maximum of 300 attempts. Solve each laminar
station with damped Newton and fixed outer velocity, never the initializer's
inverse fallback or silent extrapolation. Require maximum station residual
<1e-9. Allow separated laminar states; record nonlinear failures explicitly.
If continuation cannot reach the requested profile, it is incomplete, not a
prediction of physical separation or of the experimental flow. N>=9 is a
conditional threshold crossing; no post-crossing laminar result is interpreted.

Report coverage first, pressure/Hk/amplification profiles and endpoint N only
for complete below-threshold targets. Independently recompute station residuals,
pressure inversion, seed preservation and archive hashes. A large response to
model tap sampling alone limits what experimental taps can discriminate.
Agreement between two experimental reconstructions would strengthen the
pressure-input hypothesis; disagreement/failure cannot establish a closure
defect. No fitting, geometry edit, Ncrit sweep or automatic campaign follows.
