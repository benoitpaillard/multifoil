# S414 fore leading-edge diagnostic — 3 October 2026

Use the fully audited CL 0.056 and 0.175 nominal saved states, with CL 0.284
as control. The prior four-state audit, source hashes, Ncrit 9, geometry,
340/240 panels, Mach 0.10, Re 1m and initialization history remain unchanged.
Preserve every previous receipt. No nonlinear flow solves or model changes.

1. Trace geometry through the 83-vertex reported fore contour, PCHIP/cosine
   input sampling and mfoil's curvature paneling. Verify rebuilt geometry
   against saved arrays. Report near-nose spacing, turning and source-to-model
   distances over x/c <= 0.06. Distances to a published polyline are reconstruction
   diagnostics, not as-built errors. Retain the previously declared cubic-source
   sensitivity and its existing CL 0.284 state; do not substitute public .pts
   geometry or infer that a smoother contour is the true shape.
2. Rebuild the unchanged nominal inviscid linear operator once and evaluate
   pressure at the three archived incidences. Compare its geometric fore-lower
   branch with saved viscous pressure and qualified experimental taps under both
   x conventions. This is a same-incidence inviscid/viscous diagnostic, not a
   matched-lift inviscid validation. Do not change incidence or refit pressure.
3. Trace lower-flow-branch theta, delta*, Hk, Re_theta, edge velocity, pressure
   gradients, Cf, N and dN/ds. Interpret the third state as N only while laminar;
   after transition it is shear. Use distance from stagnation in the equations,
   preserving flow order even where global x is nonmonotone.
4. Check every laminar interval's amplification balance, separate the base
   envelope growth from the existing near-Ncrit regularization, independently
   bracket the transition root, and check representative rate derivatives.
   Re-march amplification only with frozen thickness/velocity states and restore
   all arrays exactly. A Gauss-quadrature comparison along linearly interpolated
   saved states is a discretization diagnostic, not a recomputed boundary layer
   or a valid prediction with the regularization removed.
5. Report which explanations the evidence supports or weakens. Keep geometry,
   discretization, closure and solution-branch uncertainties separate. Produce
   reproducible tables/graphs and recommend the next discriminating check.

No Ncrit fit, empirical wake correction, automatic grid campaign or new geometry
admission follows. A new coupled-flow comparison, if justified, needs a bounded
protocol after this diagnostic has been interpreted.
