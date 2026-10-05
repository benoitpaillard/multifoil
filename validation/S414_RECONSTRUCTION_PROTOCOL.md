# Bounded S414 reconstruction comparison — 2 October 2026

## Question and evidence

Does unretuned MultiFoil reproduce the measured aft pressure recovery at two
independent S414 conditions, and does that diagnosis survive small, declared
changes to the reconstructed input? This is an exploratory comparison with
experiment and published MSES; exact tested-geometry admission remains open.
No result selects a preferred geometry by agreement with measured drag.

Use the modified nominal test outline extracted from Figure 1 of RDECOM
TR 10-D-112 (2010), with its finite trailing edges and common report coordinates.
The source qualification and its failed exact tap-coordinate check are retained
in `results/s414_geometry_qualification_20261002/`. Use Figure 12(a,b) Cp and
Figure 13(c) drag at Re 1,000,000, Mach 0.10, free transition. Published MSES
uses Ncrit 9 and the 2010 approximate-envelope method. The 2014 full-e^N curves
are a separate reference. No MSES executable or new MSES run is involved.

## Nine jobs, fixed before flow results

All jobs use Ncrit 9 on both elements, unit reference chord (physical model
457.2 mm), 340 fore/240 aft panels, and a 1c wake. No geometry renormalization
or Ncrit fit is allowed. Numerical convergence remains unestablished.

| Scenario | Source contour treatment | Aft z displacement |
| --- | --- | --- |
| nominal | Production PCHIP per surface branch in arc length | 0 |
| source_cubic | Natural cubic spline instead of PCHIP at the same cosine arc-length samples | 0 |
| aft_z_minus | Nominal, then rigid aft displacement | −0.13 mm = −0.0002843394576c |
| aft_z_plus | Nominal, then rigid aft displacement | +0.13 mm = +0.0002843394576c |

Run every scenario at measured CL 0.284 and 1.041 (eight jobs), plus one fresh
exact nominal CL-0.284 repeat. Preserve source TE endpoints and the LE branch
split. The cubic alternative changes source regularization only; the production
coupled solver and its internal boundary-layer paneling remain unchanged.
Initial geometry checks must verify nonintersection and finite TE thickness.

The reported ±0.13 mm manufacturing contour tolerance motivates the displacement
scale; it is **not measured rigging uncertainty, a probability distribution, or
a confidence bound**. These scenarios do not exhaust fore contour, nose,
x-registration, edge or as-built uncertainty. Cubic interpolation is an
alternative reconstruction, not a correction. Report all responses separately.

The two nominal jobs and the repeat run first. Require the low-lift repeat to
select the same initialization path and reproduce all saved arrays exactly
(drag difference at most 0.01 count), before launching six sensitivity jobs.
Stop the campaign if that gate fails. Failed targets retain null selected
forces; completed failures are not silently retried on resume.

## Numerical and execution protocol

Use the existing combined initialization policy: cold local, fresh −1° anchor,
then fresh −1° anchor with an explicit entry 0.4° below the center. Centers are
−1.5717960133708968° / 5.063277366434972°, obtained from the increasing-lift branch
of the published MSES lift curve at each target. They are only initialization
heuristics. Preserve numerical CL-based recentering and all failed attempts.
Allow alpha in [−6, 8]°, initial half-width 0.05°, at most two brackets with
20 requested evaluations each, 100 Newton iterations, residual tolerance 1e−7,
and absolute CL tolerance 2e−5. Measured drag never enters the flow evaluator.

Freeze source/input hashes, Python 3.12.13, NumPy 2.3.5, SciPy 1.17.0 and
Matplotlib 3.10.8. Run at most four processes, each with BLAS/OpenMP threads 1.
Save full selected states, surface/wake exports, geometry, path history, residuals
and wall time. Resume only with the identical manifest/runtime; preserve partial
attempts in `interrupted_attempts/`. A lock prevents duplicate execution.

## Interpretation after completion

First report coverage, exact-repeat results, selected paths and a full saved-state
audit with freshly recomputed residuals/forces. Completion alone is not validation.
Do not launch an automatic grid extension or rescue a failed scenario by fitting.

Overlay measured Cp, published MSES Cp, and all MultiFoil cases in the common
source frame. Tabulated tap coordinates are the primary scoring coordinates;
also report sensitivity to the plotted coordinates without fitting translations.
Exclude the eight repeated lower-curve endpoints and unsupported tap centers
from numeric scores, and report the excluded locations. Recheck branch support
for each reconstruction and use common supported taps for scenario comparisons.
Preserve every original marker in source plots. Partition aft upper-surface Cp
errors into x/c < 0.80, [0.80,0.86), [0.86,0.94), and >= 0.94 regions, and show
the signed pressure rises between successive supported taps in 0.80–0.98.
Compare the location of recovery with predicted Cf, separation and transition;
do not label a sparse pressure feature as an exact measured transition endpoint.

Experimental captions round CL to three decimals; retain the ±0.0005 rounding
interval as a reference-coordinate limitation. MSES Cp is published at nominal
CL 0.28/1.04, whereas MultiFoil targets measured 0.284/1.041. Disclose this small
mismatch. The drag comparison interpolates both original curves at the same
measured CL, within their support, with no extrapolation. Retain caption drag,
nearest original marker, and their plotting offsets separately. Reference
interpolation precision is not measurement uncertainty.

Compare each sensitivity to its own nominal at the same lift; qualify differences
when selected paths differ. If the pressure diagnosis persists across all tested
inputs, it is robust only to these declared perturbations. If it changes, geometry
uncertainty prevents that diagnosis. Neither outcome isolates a unique closure
error or certifies small design gains. Close this bounded study after its audit.
