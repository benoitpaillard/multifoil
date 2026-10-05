# S702 operational wake-extent sensitivity

Declared 22 September 2026 before new flow solves. Priority 2 has passed the
seven-case combined policy/state audit. Priority 3 source limitations remain
in `results/s702_conditions_audit_20260922/FINDINGS.md`.

Question: does extending the modeled wake materially change the aft penalty
that cancels the fore benefit in the 5-to-10-degree flap comparison? This is a
sensitivity experiment, not a drag fit or independent validation.

Four fresh 2c cases: flaps 5 and 10 degrees, each at CL 0.8 and 1.2.
Re = 700,000, Mach = 0.07, Ncrit = 6.25, reference chord 1, reconstructed
geometry and 260/180 surface panels are fixed. Four corresponding audited 1c
states from `s702_initialization_policy_20260921` are the archived references.
A fifth fresh job repeats the 1c flap-10/CL-0.8 case as a driver control.
The numerical core and selected archived reference files must match their
previous manifest/audit hashes; no cache guard is bypassed.

Every new job uses the same declared ordered initialization sequence:
fresh cold at the previous numerical incidence center, fresh −1° anchor,
then fresh −1° anchor followed by an explicit center−0.4° seed. The first
successful path wins. Bounds [−6°, 8°], initial half-width 0.05°, at most
two brackets and twenty target evaluations per bracket, 100 Newton iterations,
CL tolerance 2e−5, residual tolerance 1e−7. Archived references used the first
two paths, so the third path extends the available fallback but does not
change their successful path. Record all requested and internal failures.
No experimental drag enters incidence or path selection.

The solver's wake-length parameter controls downstream extent in reference
chords; the fore wake is longer to approach the common downstream plane.
Its existing node-count formula also changes with extent. Therefore this
test changes extent plus its derived discretization, not length at fixed
resolution. Save actual wake coordinates, arclengths, tangents, node counts
and spacing with selected Newton arrays. Do not attribute a response to
length alone. No wake merging/confluence model is added.

Before interpretation require 5/5 feasible jobs, full selected-state/export
audit, and fresh 1c control agreement within 0.01 CD count of its archived
reference. A failed target remains failed and excludes a complete paired
conclusion. A change in the 2c result is an observable, not an accuracy failure
or improvement. For each CL compare flap-10 minus flap-5 at each wake setting,
then their difference; compare individual cases with their own 1c reference.
Report fore/aft wake drag, CL/CM/incidence, residual, transition and separation,
wake/body clearance and downstream momentum evolution. Relate sensitivity
to the 0.5-count design screening threshold and the previously reported
CL-0.8 flap-delta discrepancy (~19.38 counts), without selecting an extent
to minimize experimental error. Further independent resolution/length tests
require an informative observed response; no broad sweep is preauthorized
by this protocol.

Budget: five independent workers, one BLAS/OpenMP thread each, one attempt per
declared path. A source/configuration change requires a new output directory.
Interrupted unfinished job directories are retained before fresh restart;
completed results are reused only under the exact manifest/runtime. No
portable Newton restart is claimed.

Run/resume from the isolated remote root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python validation/run_s702_wake_extent.py --out validation/results/s702_wake_extent_20260922 --workers 5
```

Monitor locally:

```bash
python3 validation/run_s702_wake_extent.py --remote-status --watch
```

Remote root: `/home/jony/multifoil-runs/s702-wake-extent-20260922` on
`jony4-anywhere`. `launch.log`, study `status.json`, per-job `run.log` and
`progress.json` remain available independently of the monitor. Ctrl-C stops
only the monitor. Retrieve with:

```bash
rsync -azc -e 'ssh -o BatchMode=yes -o StrictHostKeyChecking=yes' jony4-anywhere:/home/jony/multifoil-runs/s702-wake-extent-20260922/validation/results/s702_wake_extent_20260922/ validation/results/s702_wake_extent_20260922/
```

The runner summary checks coverage and the control only. Selected-state audit,
paired plots and physical interpretation are the next interactive step after
completion; a finished process alone does not complete priority 4.
