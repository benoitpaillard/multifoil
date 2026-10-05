# S414 fore-only reconstruction test, 4 October 2026

Primary question: does the fore-element coordinate reconstruction materially
control the erroneous CL 0.056 → 0.175 fore-transition/drag change? Existing
derivative and incidence checks did not remove that change. The bounded
prescribed-pressure replay could not complete both the sampling controls and
the alternative experimental reconstruction, so it cannot distinguish pressure
input from transition closure. Use the fully coupled solver for this test.

Four jobs, two concurrent workers, on jony4:

1. Exact nominal CL 0.056 repeat against the archived lower-bucket state.
2. Exact nominal CL 0.175 repeat against its archived lower-bucket state.
3. Fore-only source-cubic reconstruction at CL 0.056.
4. The same fore-only reconstruction at CL 0.175.

Both controls must pass the existing exact saved-array gate before either new
geometry job starts. A failed/incomplete target cannot contribute a finite force
or an accuracy score. Preserve all attempts and initialization histories.

The intervention uses the existing `s414_geometry.cubic_repanel` on the reported
fore coordinates. The aft element remains the nominal input exactly. Do not
change source coordinates, element placement, trailing-edge faces, panel counts
(340/240), 1c wake, Re 1m, M 0.1, Ncrit 9, tolerances, fixed-lift center/bounds,
or ordered initialization policy. Use the **legacy Jacobian in all four jobs**
to reproduce the original protocol exactly and isolate reconstruction alone.
Its identified derivative defect is retained as a known limitation; the completed
derivative-only study did not change the CL 0.056 aerodynamic solution.

This changes the fore contour interpolator, followed by the unchanged production
paneling procedure. It is not a pure local-nose edit and is not a measured model
uncertainty bound. The earlier cubic study changed both elements and evaluated
CL 0.284/1.041; it did not test this lower-bucket event.

After execution, independently audit all four selected states, including actual
geometry, Cp/Cf, residuals, transition cells, forces, hashes and protocol. Plot
total and element drag, fore-lower transition, nose pressure, Hk and N. Retain
the same experimental and published 2010 MSES drag references at both lifts;
no corresponding MSES nose-pressure/transition curves are available.

Primary diagnostic is the change in the **paired drag drop**, alongside changes
in fore-nose recovery, amplification and transition. Also report each absolute
experimental error and coverage. A changed transition/drag drop associated with
the reconstruction implicates sensitivity to that input, not accuracy of the
alternative geometry. An unchanged drag drop is informative only if the
intervention actually changes the relevant nose pressure/recovery. No material
pressure change means this alternative did not discriminate the mechanism.
Ncrit fitting, automatic grid expansion and a geometry optimum are excluded.

The run ends after these four jobs. Completion is execution only: full saved-state
audit and physical interpretation remain separate follow-up work.
