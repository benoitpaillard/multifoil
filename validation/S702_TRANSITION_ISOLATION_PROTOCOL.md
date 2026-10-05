# S702 element-by-element transition verification

Declared 2 October 2026 before launching the new calculations. The user asked
for the most informative next approach and its verification. The preceding
matched-pressure study shows a flap-5 pressure-recovery discrepancy, while
global Ncrit reduction improves low-lift total drag mainly through the fore
element. Separate these responses at the four measured pressure conditions.

## Source gate and choice of experiment

The new geometry reconciliation directly re-extracts Figure 1 and compares
the outlines on PDF pages 43, 44, 48 and 49 under their own x-axis calibrations.
Maximum x difference is 1.44e−5c; after one common vertical registration to the
fore element, maximum z difference is 4.78e−6c. The aft-tap inconsistency remains
0.000844–0.000853c. These repeated drawings corroborate our reconstruction;
they do not establish the as-built geometry. No correction or justified rigid
aft displacement follows. Exact tested coordinates and film shape remain open.
See `results/s702_geometry_reconciliation_20261002/FINDINGS.md`.

A geometry perturbation sized from this one unsupported tap would be arbitrary.
A further global Ncrit sweep would mix fore and aft changes again. Instead,
complete the two-by-two Ncrit factorial using the existing 9 and 6.25 values,
without choosing a new threshold from experimental agreement. This tests the
model's sensitivity to element-specific transition assumptions. It does not
model transported wake turbulence, prescribed transition, film steps, or a
measured element-specific Ncrit.

## Eight new mixed states, two fresh controls

Use the four qualified pressure snapshots from the completed matched-pressure
study, each at its own measured CL. Re = 1 million, M = 0.10, 340/240 panels,
1c wakes, unmodified Figure-1 geometry, original solver and closure equations.

| Fore Ncrit | Aft Ncrit | States |
|---:|---:|---|
| 9 | 9 | Four audited archived states |
| 6.25 | 9 | Four new fore-only changes |
| 9 | 6.25 | Four new aft-only changes |
| 6.25 | 6.25 | Four audited archived states |

Fresh controls reproduce flap-5 low-lift N9/N9 and high-lift N6.25/N6.25, using
their exact previous initialization centers and policies. Run these first.
Both must converge, repeat their selected paths, meet 0.01-count CD agreement
and reproduce every saved Newton/geometry/branch/wake array exactly before
the eight mixed states start. Failure stops the campaign for investigation.
No failed control or mixed target may be omitted from coverage.

Mixed-state incidence centers are the mean of the two archived **declared
numerical centers** for that snapshot, with no experimental drag/Cp selection.
Use unchanged cold-local, fresh −1° anchor, then fresh −1° anchor with
center−0.4° entry. First success wins. Domain [−6,8]°, initial half-width 0.05°,
at most two brackets and 20 requested evaluations per bracket. All fresh
solver instances receive the declared fore/aft Ncrit values before any solve.
CL tolerance 2e−5, coupled residual tolerance 1e−7, 100 Newton iterations.
The scalar legacy Ncrit diagnostic describes the fore element; the explicit
`ncrit_by_element` field is authoritative and saved at every evaluation.

## Predeclared analysis

Retrieve and audit all 10 new selected states and the 8 archived endpoints.
Recompute residuals/forces with the actual per-element thresholds. Verify
runtime/source/geometry hashes, failed attempts, selected histories, transition
cells and surface/wake/Newton exports. Compare initialization paths explicitly;
if they change, qualify causal attribution to the threshold contrast.

For each snapshot let Q00 = Q(9,9), Q10 = Q(6.25,9), Q01 = Q(9,6.25),
Q11 = Q(6.25,6.25). Report fore effect Q10−Q00, aft effect Q01−Q00,
and interaction Q11−Q10−Q01+Q00. Also show the conditional aft effect
Q11−Q10 and conditional fore effect Q11−Q01. Apply these to total/element
drag, incidence, pre-hinge ΔCp and regional Cp error. Report transition and
Cf-negative intervals as topology/location observations rather than smooth
derivatives. Keep the same measured taps, regions, duplicate exclusions and
unsupported foremost aft tap as the original analysis. No geometry changes
mean support should be unchanged; reject unexpected support changes.

Key discrimination: does reducing **aft** Ncrit with fore Ncrit held fixed
recover the observed flap-5 pressure evolution, or does the better total drag
come mainly from fore transition while aft Cp remains wrong? Require improvement
of the lower/higher-lift evolution together before interpreting the mechanism
as supported. Even a successful sensitivity is not causal proof of wake-induced
transition. If pressure remains mismatched, prioritize source geometry/film
and local transition/closure limitations; do not extrapolate to a fitted Ncrit.

Retain experimental and published MSES N9 drag at each snapshot's CL as context.
No matching MSES Cp or MSES executable is available. The snapshots have unequal
CL across flap settings; they are not fixed-lift flap differences and do not
replace the Re-700,000 comparison. Numerical convergence and independent
contour-design validation remain unestablished. The automatic grid campaign
stays closed. No further parameter sweep is authorized by this runner.

## Execution and continuation

Isolated jony4 root: `/home/jony/multifoil-runs/s702-transition-isolation-20261002`.
Reuse pinned Python 3.12.13 / NumPy 2.3.5 / SciPy 1.17.0 / Matplotlib 3.10.8;
four workers maximum, one BLAS/OpenMP thread each. Resume only with identical
code, data, settings and runtime; retain interrupted attempts and completed
failures. Each evaluation and selected full state is saved incrementally.

```bash
# In the isolated remote directory; also the guarded restart command:
bash run_transition_isolation.sh
# From the local repository:
python3 validation/monitor_s702_transition_isolation.py --watch
```

After completion, retrieval and full selected-state audit/interpretation remain
required. The automatic summary checks execution and controls only.
