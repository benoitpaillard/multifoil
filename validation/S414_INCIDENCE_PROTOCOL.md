# S414 matched-incidence diagnostic — 3 October 2026

The derivative correction leaves the audited CL 0.056/0.284 aerodynamics
unchanged. The next question is whether the fore-nose transition discrepancy
persists at the experimental pressure-figure incidence, while explicitly
reporting the resulting lift error. This is a complementary mechanism check;
the original matched-lift drag errors remain the primary validation result.

Launch four jobs, at most two single-threaded workers, on jony4:

1. Repeat corrected CL 0.056 with exactly the archived fixed-lift settings.
2. Repeat corrected CL 0.284 with exactly the archived fixed-lift settings.
3. At the Figure 7(b) measured incidence −3.07°, solve once from a local cold
   start with the same geometry and physical/numerical parameters.
4. At the Figure 7(d) measured incidence −1.04°, use the same fixed-angle rule.

Both fresh controls must reproduce all 18 arrays of their own corrected
reference exactly, with the same selected path and ≤0.01 count drag difference,
before the fixed-angle cases start. A failed control closes the gate.

Keep the nominal reported test geometry, 340/240 panels, Re 1m, Mach 0.10,
Ncrit 9, 1c wake, `corrected_af`, 100 Newton iterations and 1e−7 residual
tolerance. Fixed-lift controls retain 2e−5 lift tolerance, archived centers,
bounds and initialization policy. New fixed-angle jobs use `start_alpha`
equal to the measured caption angle; native solver safeguards remain in place,
but no additional warm path or retuned angle is allowed. A failed fixed-angle
job remains failed with no selected finite force. No lift target is imposed
on these two jobs; record their reference experimental CL separately.
These direct cold solves have different evaluation histories from the
five-evaluation fixed-lift controls. Any qualitative change remains conditional
on that initialization; it cannot by itself exclude branch selection or prove
that incidence is the sole cause. A continuation check would require a later
bounded decision if this ambiguity becomes material.

Use caption angles, not an optimized pressure match. The corresponding
appendix rows have angles −3.074° and −1.037° and CL 0.0556/0.2837; preserve
the caption/table rounding and data association. The source already states
that standard low-speed boundary corrections were applied. No extra tunnel
or incidence correction is fitted. The inviscid preflight evaluates both
reported angle conventions and all supported taps, including worsened fits.

Every job records its comparison basis, requested angle or lift, derivative
mode, geometry/source hashes, complete flow diagnostics, selected Newton
arrays, wall time and failures. Reuse the pinned Python 3.12.13 environment.
Immutable manifests, file locks and retained interrupted folders guard resume.
Completed failed targets may not be silently rerun to improve coverage.

After completion, retrieve and hash-check results, then audit the saved states
with the correct derivative mode and comparison basis. Score pressure on the
same supported taps under both coordinate conventions; show the excluded taps
and missing published MSES pressure at the low-lift case. Compare fore
transition, aft recovery, pressure shape and predicted CL at the measured
angle. Report CD as context, **not as a matched-lift improvement**. Published
MSES curves remain the existing matched-lift reference; do not relabel them as
matched-incidence data without source qualification.

If matching incidence improves the peak but not its neighboring taps or total
lift, an angle offset alone is inadequate. If early fore transition persists,
incidence selection does not remove it. Either result informs a later bounded
geometry/pressure-input investigation; no automatic additional case, shape
change, Ncrit fit or grid extension follows this four-job diagnostic.
