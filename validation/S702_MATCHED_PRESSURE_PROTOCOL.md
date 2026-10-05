# S702 matched pressure diagnostic

Declared 2 October 2026 before new flow solves. The question is whether the
model reproduces the measured fore/slot loading, aft pressure recovery and
hinge-region features at the operating conditions actually plotted. This is a
bounded physical-observable check, not a new automatic refinement sequence.

## Four experimental snapshots, eight numerical targets

All measurements are free transition, Re 1,000,000, M 0.10, c = 457.48 mm.
The common report reference frame and reconstructed flap geometry are retained.

| Flap | Figure / PDF page | Caption incidence | Lift marker CL | Purpose |
|---|---|---:|---:|---|
| 5° | 7(b) / 43 | 2.03° | 0.8280869042 | Reported indication of an aft bubble near the low-drag comparison |
| 5° | 7(c) / 44 | 6.11° | 1.2343679194 | Reported disappearance of that pressure feature |
| 10° | 8(b) / 48 | −0.01° | 0.7739518681 | Low-drag-region aft/hinge recovery |
| 10° | 8(c) / 49 | 4.07° | 1.1763141482 | Higher-lift aft/hinge recovery |

These are digitized lift values, not raw-data precision. Associate each caption
with its nearest experimental lift marker within a declared 0.02° limit. Actual
offsets are at most 0.0073°. Save offsets and drawing IDs. Do not interpolate Cp
between incidences or force the snapshots to CL 0.8/1.2. Each model run matches
its own snapshot CL; the flap-5 and flap-10 snapshots are not a same-lift pair.
This does not directly test the larger Re-700,000 discrepancy.

Run Ncrit 9 and frozen 6.25 for every snapshot, common to both elements. No Cp
or drag observation chooses Ncrit, geometry, initialization or a winning run.
Use the unchanged production CoupledSolver with 340/240 panels, 1c wake
extension, 100 Newton iterations, residual tolerance 1e−7 and CL tolerance
2e−5. This is the existing intermediate resolution, not a claim of numerical
convergence. Source regularization, boundary-layer sampling and default wake
recipe remain coupled as in the production solver. No explicit discretization
override, closure change, trip, or hinge modification is introduced.

Incidence centers come only from the nearest two archived converged numerical
CL/incidence points for the same flap, Re/M and Ncrit; extrapolation is limited
to 0.1 in CL. Cold-local, fresh −1° anchor, then fresh −1° anchor with explicit
center−0.4° entry are tried in that order. First success wins. Bounds [−6,8]°,
half-width 0.05°, at most two brackets and 20 requested evaluations per bracket.
Retain failed requests/internal continuation and null failed target forces.

## Source qualification

`results/s702_pressure_qualification_20261002/source.json` records the original
NASA PDF checksum/URL, axis transforms, marker selectors, central dots and
original-page rasters. Four selected polylines per snapshot contain 386 plotted
markers. Lower polylines duplicate the upper LE/TE endpoints; preserve and flag
these 16 duplicates, exclude them from lower-surface scores. These are not 386
independent physical taps. Verify all 36 fore upper positions against Table II,
and available aft upper positions ahead of the hinge against its zero-flap table.
Do not use the zero-flap table unchanged for deflected post-hinge positions.
Use marker central dots, because triangle bounding-box centers are biased.

## Analysis after retrieval and audit

Require all eight requested targets for complete coverage. Recompute selected
state residual/forces, confirm saved surface and Newton arrays correspond to the
selected history, inspect transition/separation, initialization and geometry,
and verify source/runtime hashes. Completion alone is not a physics audit.

Overlay measured taps and calculated Cp in the common x/c frame, separately
for upper/lower and fore/aft. Interpolate the model to measured tap positions
only within model support; disclose any excluded taps. Report per-surface
tap RMSE/bias and aft-upper regional residuals for x/c 0.80–0.86, 0.86–0.94 and
0.94–TE. Show sparse tap spacing and avoid fitting a smooth experimental bubble
or assigning an exact transition location. Plot calculated Cf, transition flag,
amplification and separation/reattachment with Cp using consistent coordinates.
Compare changes between low/high lift within each flap, especially the reported
flap-5 bubble disappearance. Neither absent pressure plateau nor matching Cp
alone proves a boundary layer is correct: pressure can agree while transition,
momentum thickness and drag differ. Independent transition/wake measurements
would still be needed for definitive causal attribution.

Show the source geometry and hinge location; the actual film/steps are unknown.
Include total CD and published MSES Ncrit-9 CD interpolated on common support
at each snapshot's CL as context. The report provides no matched MSES Cp curves
for these four flap-5/10 snapshots. Do not draw invented MSES pressure curves;
the existing zero-flap Figure-22 pressure comparison remains separate. No MSES
executable is installed, so this study cannot evaluate MSES Ncrit sensitivity.

## Execution

Isolated jony4 root `/home/jony/multifoil-runs/s702-matched-pressure-20261002`;
reuse the pinned Python 3.12.13 environment with NumPy 2.3.5, SciPy 1.17.0 and
Matplotlib 3.10.8. Maximum four workers, BLAS/OpenMP threads one. Freeze and
verify sources, geometry, extracted data, protocol and runtime before launch.
Keep source snapshots, incremental progress and interrupted attempts. Resuming
with modified sources/runtime is rejected; completed failures are retained.

```bash
# From the isolated remote root (also the restart command):
bash run_matched_pressure.sh
# From the local repository:
python3 validation/monitor_s702_matched_pressure.py --watch
```

After completion retrieve, audit and interpret these eight states. No automatic
extra grid, tuning sweep, geometry adjustment or optimization follows.
