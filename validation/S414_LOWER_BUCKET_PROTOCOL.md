# S414 free-transition lower-bucket diagnostic — 3 October 2026

## Question and bounded scope

Does the unchanged model reproduce the joint evolution of fore-lower loading/
transition and aft-upper pressure recovery as lift crosses the lower drag-bucket
boundary? Use the qualified Figure 7(a–d) free-transition measurements, avoiding
the physical-trip effectiveness/displacement confounds. Four nominal-geometry
jobs only: CL −0.072, 0.056, 0.175 and an exact fresh CL 0.284 control against
the audited 2 October reconstruction study. This is an exploratory mechanism
comparison, not a new geometry calibration, grid campaign or causal proof.

Source: `results/s414_trip_qualification_20261003/`. Use the rounded pressure
caption lifts, retaining exact appendix rows and their offsets. Re 1,000,000,
Mach 0.10, Ncrit 9 on both elements; unchanged reconstructed modified tested
outline, finite trailing edges, 340/240 panels, unit reference chord (457.2 mm
model) and 1c wake. Source geometry and numerical convergence remain unqualified
for small design gains. No changes to transition/wake closures or physical trips.

## Frozen execution and gate

Run the fresh CL 0.284 control first, using its original center
−1.5717960133708968 degrees and original initialization policy/settings. Require
the same selected path, all 18 saved arrays exactly equal to the audited control,
and at most 0.01 drag count difference. A failed control prevents the three new
jobs from launching, including on resume. This fresh run validates continuity;
the archived state remains the reference rather than being overwritten.

New CL 0.056 and 0.175 centers come from the previously qualified increasing-lift
branch of the published MSES lift curve. At CL −0.072 that curve lacks coverage;
use the measured incidence −4.09 degrees only as an initialization heuristic.
No extrapolated MSES value is a validation reference. Centers do not select a
physical result or fit experimental drag; the same numerical lift-based
recentering and ordered paths apply to every job.

Retain cold-local, fresh −1-degree anchor, and fresh −1-degree anchor with entry
0.4 degrees below the center. Incidence bounds [−6,8] degrees, initial half-width
0.05 degrees, at most two brackets/20 evaluations per bracket/path, 100 Newton
iterations, residual tolerance 1e−7 and absolute CL tolerance 2e−5. Failed
requests and internal continuation failures remain recorded; failed final targets
have null forces and are not retried merely by resuming.

Freeze all numerical/source hashes and use the existing jony4 Python 3.12.13,
NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 environment. Reuse the unchanged
S414 worker. At most three independent workers after the control, one BLAS/
OpenMP thread each. Save source snapshots, geometry, runtime, all evaluations,
selected arrays, surface/wake exports, wall times and status incrementally.
An exclusive lock prevents duplicate execution; preserve interrupted attempts
before restarting with identical code/settings/runtime. No automatic extension.

## Analysis after completion

Report four-target coverage, control reproducibility, requested/internal failures
and selected paths first. Audit full exports and freshly recompute residuals,
forces, Cp/Cf and boundary-layer summaries before interpreting flow results.

Compare the fore lower surface and aft upper recovery jointly. Score all
commonly supported independent taps, with Table-II coordinates primary and
original plotted coordinates as a sensitivity. Exclude repeated endpoints and
the existing unsupported tap centers; disclose masks and preserve every source
marker. Retain the failed strict Figure-7/12 copy check (up to 0.000338 Cp and
0.00008123c), with no fitted source alignment.

Plot measured/model Cp, computed separation/transition/reattachment and fore/wake
states. Use fore-lower x/c <0.60 and ≥0.60, plus aft-upper x/c <0.80,
[0.80,0.86), [0.86,0.94) and ≥0.94 regional pressure diagnostics. In
0.80≤x/c≤0.98 retain all adjacent-tap pressure rises and slopes; strongest-rise
locations are sampling statistics, not exact bubble endpoints. The experiment's
largest-rise interval is nonmonotonic across these snapshots; do not impose a
monotonic bubble-shortening target or infer numerical transition from a Cp kink.
The report qualitatively describes turbulent fore-lower flow at −0.072 and
laminar flow surviving the pressure peak near the lower boundary at 0.056;
no measured case-specific transition coordinates or wake turbulence are available.

Compare total drag at the same lift using appendix free-transition data. Published
2010 envelope-MSES drag covers 0.056, 0.175 and 0.284, not −0.072. Only the 0.284
case has a matching published MSES pressure reference (nominal CL 0.28); mark the
other three pressure references missing. Keep the 2014 full-e^N results separate.
Decompose model fore/aft drag for interpretation only, never as measured element
drag. Distinguish altered mean-flow coupling from the absent direct turbulence
channel: matching or mismatching this sequence alone does not prove a cause.

Success is a fully audited, coverage-qualified joint pressure/transition/drag
comparison. Any new closure or investigation requires a further evidence-based
decision. Stop after these four jobs; do not tune Ncrit or rescue missing cases
by changing the frozen protocol.
