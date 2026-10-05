# First optimization-readiness pilot

This pilot implements and exercises a fixed-lift design-screening workflow.
It checks numerical sensitivity of predicted geometry changes. It does not
establish agreement with independent MSES/RANS or experimental geometry rankings,
and does not perform continuous shape optimization.

[Executed results and plots](results/optimization_pilot/REPORT.md).
[Literature and remaining validation needs](OPTIMIZATION_READINESS.md).

## Declared study

Nominal NLR7301 ECARP coordinates; Re = 2.51 million; Mach = 0.185;
target CL = 2.484; reference chord = 1; moment reference = (0.25, 0).
Ncrit = 9 except the Ncrit = 6 sensitivity cases. The target lift is the printed
CL of [MIT's official MSES example](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf).
Its exact geometry is unconfirmed; it supplies context, not a matched-input score.

| Variable | Tested change | Definition |
|---|---|---|
| Main-element camber | +/-0.001 local chord; +/-0.0005 checks | Add `16*a*c*x^2*(1-x)^2` along the local chord normal to both surfaces |
| Flap placement | +/-0.002 reference chord | Translate the entire flap vertically in assembly coordinates |

Here `x` is normalized local chordwise position. The camber mode fixes the
leading- and trailing-edge endpoints and end slopes, and preserves thickness
at fixed chordwise position in the input shape. It is a smooth contour change.
The flap change is a placement change, not a prescribed minimum-gap increment.
The assembly is not renormalized after moving the flap. Geometry is modified
before repaneling, so every panel density samples the same underlying variant.

| Numerical setting | Main/flap panels | Wake extension beyond aftmost TE | Ncrit | Initialization |
|---|---:|---:|---:|---:|
| coarse | 180/120 | 1 reference chord | 9 | 4 degrees |
| base | 260/180 | 1 reference chord | 9 | 4 degrees |
| fine | 340/240 | 1 reference chord | 9 | 4 degrees |
| long_wake | 260/180 | 2 reference chords | 9 | 4 degrees |
| ncrit6 | 260/180 | 1 reference chord | 6 | 4 degrees |
| start0 | 260/180 | 1 reference chord | 9 | 0 degrees |

Baseline and four full-sized variants run at the first four settings: 20 jobs.
Two half-sized camber changes run at base settings. Baseline and full-sized
camber variants run at Ncrit 6. Baseline runs once with the alternate start:
**26 geometry/settings combinations** in total. Each job requires several
incidence evaluations; the number of jobs is not the number of flow solves.
Surface panel counts are nominal input counts; actual surface and wake nodes
are retained in each job's exported distributions.

The initial pilot found that changing Ncrit from 9 to 6 reversed the preferred
camber direction. Three **additional adaptive checks** repeat baseline and
full-sized camber variants at Ncrit 6 on 340/240 panels. They are selected after
seeing the reversal, not independent holdouts, and are stored separately under
`results/transition_followup/`. The original 26-job study remains intact.

Ncrit sensitivity is not a roughness, forced-trip or fully turbulent model.
Only one alternate initialization path and one smaller camber step are checked;
passing them does not prove global branch uniqueness or gradient convergence.
Lengthening a wake also changes its node count under the existing solver rule.
This checks the combined wake truncation/discretization choice, not either in isolation.

## Fixed-lift objective and constraints

`optimization.py:solve_fixed_lift` solves incidence within [4, 7] degrees on an
increasing-lift branch, using a safeguarded secant procedure. Both endpoints
must converge and bracket the target; no extrapolation is allowed. CL tolerance
is 0.00002. Each flow solve uses residual tolerance 1e-7 and at most 100 Newton
iterations per continuation point. Every geometry starts with a new solver;
flow-state reuse occurs only within its incidence solve.

Intersections, minimum clearance below 0.015 reference chords, nonconvergence,
nonpositive drag, local supersonic flags and unreachable target lift make a
job inadmissible. Clearance is checked on both input and final repaneled
contours. These clearance and incidence limits are declared pilot restrictions,
not universal physical safety limits. Limited modeled separation is allowed
and its arclength fraction is reported; massive separation is not validated.
CM is reported without imposing an aircraft-specific trim constraint.

`CD` remains total modeled two-dimensional profile drag from the independent
Squire-Young wake estimates. No friction term is added to it. There is no
induced drag or shock drag. The flow solver and boundary-layer closures were
not changed or fitted for this pilot.

The analysis compares each variant with its baseline under identical numerical
settings. The spread across settings is an observed sensitivity range, not a
rigorous error bound. A stable sign alone is insufficient when changes are tiny.
Central differences use independently re-solved fixed-lift states; they are
not an adjoint or proof of an exact geometry derivative.

## Run or resume

From the extracted project root, after installing `requirements.txt`:

```bash
python -m unittest test_optimization test_validation
OPENBLAS_NUM_THREADS=1 python validation/optimization_pilot.py --workers 3 --out pilot_rerun
OPENBLAS_NUM_THREADS=1 python validation/transition_followup.py --workers 3 --out transition_rerun
python validation/retry_pilot_job.py --study transition_rerun/study.json --geometry baseline --setting ncrit6_fine --start-alpha 0 --out transition_retry
python validation/rebracket_transition.py --study transition_rerun/study.json --out transition_rebracket
python validation/report_optimization_pilot.py --out pilot_rerun --transition-followup transition_rerun --transition-retry transition_retry --transition-retry transition_rebracket
python validation/plot_three_way.py
```

The study can take several minutes per group of geometries. `--workers` controls
independent operating-system processes; the default is 1. Each flow evaluation
is checkpointed. Rerunning the same command resumes completed jobs only if the
study settings and code/input hashes match. Failed completed jobs remain failed;
use a new folder for a changed study. A changed code or study cannot silently
reuse cached results. `--summarize-only` refreshes tables without solving flow.

The initial fine-grid Ncrit 6 baseline exhausted its Newton budget when started
directly at 4 degrees. Its follow-up command returns exit code 2 for this recorded
failure. The explicit retry above changes initialization only, starts at 0 degrees,
and preserves the original attempt in its own folder. The report includes both;
the retry is not a change to Reynolds number, Mach number, geometry or Ncrit.
That initialization retry also failed at the lower bracketing point (4 degrees).
A further adaptive attempt narrows the bracket to [5,6] inside the original
[4,7] domain and starts at 5 degrees. This distinguishes bracketing difficulty
from failure at the requested lift. The report keeps the initial attempt and
both retries; no failed CD is replaced with a numerical penalty.
The [5,6] retry reached the target. Its baseline-relative comparisons are
qualified by the changed procedure; the two modified geometries retain their
original matched procedures and still reverse their drag ordering between
Ncrit 9 and Ncrit 6 on both panel densities.

Every job retains original modified input coordinates, actual surface/wake
coordinates, all requested incidence solutions, transition/separation diagnostics,
software versions and final status. Internal continuation points are listed in
solver diagnostics. Results never substitute a penalty value for physical CD.

The seven unit tests exercise the fixed-lift contract and geometry operations;
they are not aerodynamic validation. Existing comparison integrity tests are
retained. The report records actual completion and failures separately.

## What would justify accepting an optimum

Obtain exact baseline and changed coordinates with independent MSES/RANS or
experimental results at matched CL, Reynolds number, Mach and transition
conditions. Check rankings and changes in drag across multiple operating points,
and verify finalists independently. The published S702 report offers coordinates
through its provider but does not tabulate them. An outline drawing or the single
MIT coefficient point cannot substitute for the missing geometry-ranking data.

This pilot prepares the numerical machinery and exposes potential sensitivity;
it does not close that independent-validation gap. The existing high-incidence
NLR drag error and missing wake-confluence physics remain unresolved.

The measured transition sensitivity in this code is a concrete reason to include
multiple transition and operating conditions in a future optimization objective.
A geometry that improves drag under one Ncrit assumption can degrade it under
another. The numerical pilot does not identify which assumption matches an
intended installation; that requires matched experimental or independent data.
