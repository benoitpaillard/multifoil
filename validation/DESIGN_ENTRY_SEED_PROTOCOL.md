# Explicit entry seed: bounded initialization diagnostic

The twelve-job Mach-zero study finished 11/12 feasible, with 5/6 complete
comparisons. Its failed flap-10° baseline is preserved. The nominal warm anchor
at −1° failed internally, so that path never acquired a converged starting state.
The legacy control converged first at center−0.4° (−0.5662261150448021°).

Hypothesis: requesting that same entry point explicitly before the local lift
solve supplies a usable starting state. Run three independent flap-10° contours:
fore-camber amplitude 0, −0.001 and +0.001. Each uses one declared path with
start −1°, explicit seed center−0.4°, then the existing local bracket centered
at −0.16622611504480211° with half-width 0.05°. The center and seed are identical
for all three shapes and come from the previous numerical protocol, not drag.
This diagnostic does not repeat the two failed baseline paths or yet add the
entry path to the common configuration.

Keep Mach=0, Re=700,000, Ncrit=6.25, CL=0.8, 260/180 panels, unit reference
chord/wake length, existing design restrictions, Newton limit 100, residual
tolerance 1e-7, CL tolerance 2e-5 and incidence domain [−6°, 8°]. One lift-based
recenter is allowed; each of the two brackets permits 20 requested evaluations.
Including the explicit seed gives at most 41 requests/job, 123 total; internal
continuation is additional. No flow-core or fixed-lift-helper changes are made.

Require 3/3 feasible results and absolute CD difference ≤0.01 count from each
contour's already computed legacy control. Preserve all failed requests. Review
the selected surfaces, transition, force sums, Newton arrays and entry histories
before proposing a common-policy regression. Agreement would support this
initialization path, not independently validate contour rankings.

The new study snapshots and hashes source/configuration/reference files, pins
the existing Python 3.12 numerical environment, uses three single-threaded
workers, and resumes complete jobs only. Partial attempts are archived before
restarting. Changed sources or inputs require a fresh output directory.

Remote root: `/home/jony/multifoil-runs/design-entry-seed-20260921` on jony4.
Launch/resume there with:

```bash
.venv/bin/python -u validation/run_design_entry_seed.py \
  --out validation/results/design_entry_seed_20260921 --workers 3
```

Follow from the local project with:

```bash
python3 validation/run_design_entry_seed.py --remote-status --watch
```

Results are under `validation/results/design_entry_seed_20260921/`; launch and
preflight logs are in the remote root. The earlier failed integration study is
unchanged. Its audit is in `validation/results/design_initialization_analysis/`.
