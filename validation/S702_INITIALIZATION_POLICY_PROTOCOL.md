# S702 common initialization policy: bounded follow-up

Declared 21 September 2026 after auditing the original 19-target matrix and
before executing this follow-up. The original results are retained unchanged.

## Question and acceptance

The original cold flap-10 / CL=0.8 run failed all three cold attempts. Increasing,
decreasing and archived replay paths converged to matching states. Test whether a
single prescribed warm anchor at -1 degree, following one cold-path failure,
suffices without the two historical intermediate replay incidences.

Primary matrix: six independent fixed-lift jobs, flaps 5/10 and CL 0.5/0.8/1.2.
Each uses `fixed_lift_policy.solve_with_initialization`: cold at the existing
numerical incidence estimate, then a fresh solver initialized at -1 degree if
needed. No candidate ranking or seed uses experimental drag. A seventh control
uses only the original two-incidence warm replay at flap 10 / CL=0.8 through the
same new helper. This separates wrapper regression from anchor sufficiency.

Keep Re=700,000, M=0.07, Ncrit=6.25, panels 260/180, chord 1, wake length 1,
100 Newton iterations, residual tolerance 1e-7 and lift tolerance 2e-5. The local
half-width remains 0.05 degree. Declare an overall incidence domain [-6, 8]
degrees, including initialization. An unbracketed increasing pair may recenter
once using computed CL, with slope 0.02–0.3 per degree and at most one degree of
extrapolation beyond its samples. Inadmissible solves do not supply coefficients
for recentering. A failed cold solve immediately advances to the warm path;
the original repeated cold retries are not silently rerun.

Acceptance: 6/6 primary targets and the replay control converge; every selected
state meets lift/residual/admissibility checks; maximum absolute CD difference
from its archived selected reference is <=0.01 count. This last threshold is a
declared implementation-regression screen, chosen above the audited 0.00278-count
archive difference and below the historical 0.5-count design screen. It is not
an aerodynamic uncertainty bound. Audit full selected-state/topology differences
before promoting the policy to shape optimization. A failure remains a failure;
do not claim universal reliability from these seven jobs.

## Work budget and state handling

At most two brackets per initialization path, 20 requested solves per bracket,
and two paths per primary target: <=480 primary requested evaluations. The warm
control has two explicit seed evaluations and <=40 bracket evaluations. Total
nominal bound: 522 requested solves; internal continuation is additional and
preserved in the flow metadata. Seven independent workers, one BLAS thread each.

The first successful path wins. Every failed flow, initialization phase, bracket,
wall time and selected index is retained. The helper copies the complete solver
in memory at each candidate matching the lift tolerance; when the lower endpoint
is returned after evaluating the upper endpoint, continuation therefore starts
from the selected state, including its coupling arrays. These copies are not
disk-restorable checkpoints and cannot be transferred across geometries.

The wrapper and benchmark are explicit additions. The original
`optimization.solve_fixed_lift`, viscous equations and sweep worker remain
unchanged. The existing shape optimizer does not adopt this policy until the
follow-up passes and its selected states have been reviewed.

## Run and resume

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  validation/run_s702_initialization_policy.py \
  --out validation/results/s702_initialization_policy_20260921 --workers 7
```

Add `--prepare-only` to record sources/settings without solving. The same command
resumes completed whole jobs only. Partial job folders are archived before being
restarted. Manifest/source/runtime changes stop execution and require a new study
directory. A lock prevents duplicate execution. Status, logs, result JSON, polar,
surface, wake and selected state NPZ files are saved incrementally. The detached
launcher and exact remote paths are in `handover/S702_REMOTE_RUN.md`.
