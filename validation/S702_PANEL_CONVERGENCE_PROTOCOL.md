# S702 third-grid comparison

Declared 22 September 2026 before new solves. The completed 260/180 to 340/240
comparison passed all five jobs and nine selected-state audits, but changed
the flap drag difference by −0.50432 count at CL 0.8 and −1.05171 at CL 1.2.
Those changes exceed the existing 0.5-count screening threshold. Local aft
boundary-layer steps improve with refinement, while the remaining model/data
discrepancies remain much larger. This motivates one additional grid.

Four fresh cases use 440 fore / 320 aft panels, flaps 5/10 at CL 0.8/1.2.
A fifth independent case repeats the 340/240 flap-10/CL-0.8 state as a control.
All immediate references are the newly audited 340/240 states; their original
260/180 references remain pinned for the three-grid comparison. Re=700,000,
Mach=0.07, Ncrit=6.25, reference chord/axes, raw source geometry, 1c wake extent,
TE regularization, tolerances and ordered initialization remain fixed.
No measured force or retuned Ncrit chooses the grid, incidence or path.

The mesh sequence is 260/180 → 340/240 → 440/320. These are nonuniform
curvature-based meshes of a reconstructed source geometry. Resampling changes
the solver contour representation and the derived wake count/first spacing.
Record these effects explicitly; do not call this pure surface spacing with
fixed geometry representation or a formal asymptotic convergence study.

Every job uses fresh cold-local, fresh −1° anchor, then fresh −1° anchor with
explicit center−0.4° entry, first success wins. Retain the same original
numerical incidence centers, [−6°,8°] bounds, initial half-width 0.05°, at most
two brackets and twenty target evaluations per bracket, 100 Newton iterations,
CL tolerance 2e−5 and residual tolerance 1e−7. Every requested/internal failure
remains saved. No additional retry path is introduced.

`study.json` contains the complete frozen matrix. Each job declares `panels`;
its exported `metadata.case.study` equals the manifest with top-level `panels`
set to that job's list. The unchanged worker from `run_s702_wake_extent.py`
executes that explicit configuration. New source/input hashes are checked
against the completed fine-grid study; no existing manifest guard is bypassed.

Before interpretation require 5/5 feasible jobs, full selected-state/metadata
audit and fresh 340/240 control agreement within 0.01 CD count. A failed target
prevents complete paired conclusions. For each CL and each grid, compute
flap-10 minus flap-5 using that grid's own baseline. Report per-case drag,
CL/CM/incidence, fore/aft contributions, transition/separation, actual geometry,
TE endpoints/opening, local Hk/velocity steps, wake discretization/clearance,
initialization and cost. Retain experimental discrepancies as context only.

The numerical decision is whether the maximum absolute individual and paired
changes from 340/240 to 440/320 both shrink relative to the preceding step
and fall below 0.5 count. Also show every nonmonotone per-condition response;
do not hide it in the maxima. A passing screen qualifies only these states
and does not supply an error bound, independent shape validation, or a new
universal minimum design benefit. Three nonuniform grids alone do not justify
Richardson extrapolation or an asserted convergence order.

This is the last automatic grid extension in this diagnostic sequence. If
coverage fails or changes do not shrink enough, inspect the panel/geometry
representation and state changes before another run. Keep priority 3's source
limitations and priority 4's wake/BL interaction questions open; do not change
closures solely to match total drag.

Budget: five workers with one BLAS/OpenMP thread each, detached and logged.
Interrupted attempts are archived before restart. Completed jobs resume only
with exact source, inputs, settings and runtime; any change requires a new
study directory.

Remote root: `/home/jony/multifoil-runs/s702-panel-convergence-20260922` on
`jony4-anywhere`. Run/resume from that root:

```bash
.venv/bin/python validation/run_s702_panel_convergence.py --out validation/results/s702_panel_convergence_20260922 --workers 5
```

Local monitor (closing it does not stop the run):

```bash
python3 validation/run_s702_panel_convergence.py --remote-status --watch
```

Retrieve after completion:

```bash
rsync -azc -e 'ssh -o BatchMode=yes -o StrictHostKeyChecking=yes' jony4-anywhere:/home/jony/multifoil-runs/s702-panel-convergence-20260922/validation/results/s702_panel_convergence_20260922/ validation/results/s702_panel_convergence_20260922/
```

The runner reports coverage/control checks only. State audit, three-grid plots
and the contraction decision follow after completion.
