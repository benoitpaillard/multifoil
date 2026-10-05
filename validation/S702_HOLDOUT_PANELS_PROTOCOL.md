# S702 reserved-condition surface-panel sensitivity

Declared 22 September 2026 before new solves. The completed wake-extent audit
(`results/s702_wake_extent_analysis/FINDINGS.md`) found a maximum change of
0.00821 count in the flap drag difference, but identified post-transition aft
Hk/velocity steps that warrant a panel check. The existing 340/240 calibration
study used Re=1 million; this experiment covers the reserved Re=700,000 cases.

Question: does finer surface discretization change the aft drag penalty or
flap ranking at these conditions, and reduce the observed local panel steps?
No transition threshold, contour or closure is fitted or modified.

Four fresh jobs use 340 fore / 240 aft panels, flaps 5/10 and CL 0.8/1.2.
Their four archived references use 260/180 panels. A fifth fresh 260/180
flap-10/CL-0.8 control checks the driver against its audited reference.
Keep Re=700,000, Mach=0.07, Ncrit=6.25, wake extent 1c, chord/axes, raw input
geometry, minimum TE opening 1e−5 local chord, tolerances and initialization
policy fixed. Default Ncrit remains 9 outside the declared S702 scenario.

Each job declares its own `panels` list in `study.json`. Exported
`metadata.case.study` equals that manifest with its top-level `panels` replaced
by the job's list. This explicit per-job configuration is passed to the
unchanged, previously audited wake-study worker. Preserve the complete frozen
manifest and verify this single declared substitution during the next audit.
The current automatic wake-node count and first wake spacing also depend on
surface panels: this is operational panel refinement, not surface-only
refinement at fixed wake discretization. Save actual solver geometry, wake
coordinates/spacing, states and all requested/internal failures.

All five jobs use fresh cold-local, fresh −1° warm anchor, then fresh −1°
anchor plus explicit center−0.4° seed, in that order; first success wins.
Centers come from archived numerical incidences only. Bounds [−6°,8°],
0.05° initial half-width, up to two brackets and twenty target evaluations per
bracket; 100 Newton iterations, CL tolerance 2e−5, residual tolerance 1e−7.
No measured drag chooses a seed, path, grid or target.

Before paired interpretation require 5/5 feasible jobs, selected-state/export
audit, and fresh base-control CD agreement within 0.01 count. Incomplete
coverage cannot win by omission. Report fine-minus-base CD/CL/CM/incidence,
fore/aft contributions and flap-10-minus-5 changes at each CL; each grid uses
its own flap baseline. Preserve all failed targets and initialization work.
Audit Cp/BL/transition/separation, local post-transition Hk drops and corrected
velocity steps, actual repaneling, TE openings, and wake/body clearance.

Relate drag sensitivity to the 0.5-count screening threshold and the remaining
experimental discrepancy, without fitting either. The XFOIL documentation's
approximate Hk-drop/velocity-step guides are diagnostics, not MultiFoil validity
tolerances. A two-grid result is not an asymptotic error estimate. An apparent
improvement in total drag alone does not validate geometry rankings or identify
the correct closure. Priority 3 source/hinge/trip limitations remain open.

Budget: five workers with one BLAS/OpenMP thread each. Source or runtime changes
require a fresh study directory; completed results resume only under an exact
manifest, and interrupted attempts are archived. The core remains unchanged.

Remote root: `/home/jony/multifoil-runs/s702-holdout-panels-20260922`, accessed
through `jony4-anywhere`. Run/resume from that root:

```bash
.venv/bin/python validation/run_s702_holdout_panels.py --out validation/results/s702_holdout_panels_20260922 --workers 5
```

Monitor locally (Ctrl-C affects only the monitor):

```bash
python3 validation/run_s702_holdout_panels.py --remote-status --watch
```

Retrieve after completion:

```bash
rsync -azc -e 'ssh -o BatchMode=yes -o StrictHostKeyChecking=yes' jony4-anywhere:/home/jony/multifoil-runs/s702-holdout-panels-20260922/validation/results/s702_holdout_panels_20260922/ validation/results/s702_holdout_panels_20260922/
```

The automatic summary checks coverage and control forces only. Full state
auditing, paired plots, resolution diagnostics and interpretation follow in
the next interactive session.
