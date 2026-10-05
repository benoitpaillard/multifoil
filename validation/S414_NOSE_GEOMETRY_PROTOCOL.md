# S414 fore-nose turning artefact (V4), 4 October 2026

Stage V4 of the verification plan in `DEVELOPMENT_PLAN.md`. The geometry trace
below was done first. The alternatives, their admission rule and the pressure
decision rule are frozen here **before any pressure is computed on them**. The
script is `validation/report_s414_nose_geometry.py`; the constructions are in
`validation/s414_nose_geometry.py`.

## Question

The final solver polygon has a small turning reversal on the fore lower nose.
Where in the geometry pipeline does it arise, can it be removed without adding
geometry information, and does removing it change the inviscid nose pressure
by more than the numerical sensitivity already measured?

## Trace (geometry only)

Stages: the 83 reported fore vertices; the continuous curve interpolated
through them; the 341 cosine-spaced input samples; the solver's spline and
curvature spacing; the final 341 solver nodes. At each stage: signed turn at
polygon vertices, signed curvature of continuous curves, spacing. Reported for
the nominal rule (PCHIP per branch), the archived fore-only cubic variant, and
the two alternatives below. The nose support is the run of reported vertices
with x / total chord ≤ 0.06 (15 vertices, both surfaces).

## Alternatives (at most two; both pass through all 83 reported vertices)

- **A, single curve.** The production rule interpolates the upper and lower
  branches separately and joins them at the leading-edge vertex without a
  tangent condition. A interpolates the whole contour with the same PCHIP rule
  as one curve. It changes only the two reported segments next to that vertex.
  It tests whether the branch split is the origin of a reversal.
- **B, convex nose.** Inside the nose support only, the interpolant is replaced
  by a tangent-continuous chain of circular arcs (two per reported segment,
  junction tangent parallel to the chord), which is convex whenever every
  vertex tangent lies between its two chords. The two join tangents are those
  of the production curve. Interior tangents are the minimum-bending-energy
  choice among convex-consistent ones, each kept at least 5 % of the local
  turn inside its chords. Convexity is imposed only because the reported
  polygon is itself convex throughout the support; the construction refuses
  otherwise. No aerodynamic quantity enters.

Both keep the production sampling parameters, and every input sample outside
the changed segments is the production sample bit for bit. The aft element,
registration, reference chord, panel counts and trailing edges are unchanged.

**Admission.** An alternative is admitted as a reversal-free geometry only if
its **final solver polygon** has no negative turn inside the nose support, it
reproduces all reported vertices, and its departures from the nominal curve are
tabulated (in mm of the 457.2 mm chord). An alternative that keeps a reversal
is reported and not carried into any viscous run. Neither alternative is
claimed to be the tested contour.

## Inviscid comparison (only after admission is decided)

Production operator and Hess–Smith on each variant's own final nodes with the
nominal aft polygon, 340/240 panels, at six archived incidences: measured
−3.07°; matched CL 0.056, 0.175 and 0.284 of the nominal geometry; and the two
ends of the nominal transition interval found in V2 (about −2.694° and
−2.569°). Cp uses the Kármán–Tsien correction at Mach 0.10.

Reported: fore-lower minimum Cp and its position; Cp at the supported nose
taps; the largest |ΔCp| against nominal on the fore-lower recovery
(0.008 ≤ x/c ≤ 0.06) and on the fore-upper nose (x/c ≤ 0.06), compared at
equal x on each surface; the strongest deceleration on the lower recovery; and
the same geometry differences evaluated with Hess–Smith.

**Decision rule, fixed now.** The pressure change is *material* if, at any of
the six incidences, the largest |ΔCp| on the fore-lower recovery or the change
of the minimum Cp exceeds 0.05, the formulation sensitivity found in V3.
Production grid sensitivity (0.013) is reported beside it.

- Not material: the reversal is deprioritized as the cause. No viscous run.
- Material, and B admitted: a separate, bounded viscous protocol is frozen for
  B at CL 0.056 and 0.175 with the V2 fixed-incidence rules, plus the CL 0.284
  control. Not launched from this protocol.

A smoother contour or a pressure change alone proves no cause. Without
as-built coordinates every geometric attribution stays conditional.
