# S702 initialization repeatability: first remote matrix

Declared 15 September 2026, before the new runs. This is the first experiment in
development priority 2. It measures the existing fixed-lift/retry policy before
changing that policy or the numerical core.

## Question and matrix

How much do selected forces and boundary-layer states depend on initialization,
relative to the 0.19/0.47-count historical design benefits and their 0.5-count
screening threshold? These numbers are comparisons, not uncertainty limits.

Keep Re=700,000, M=0.07, common Ncrit=6.25, panels 260/180, reference chord 1,
wake length 1, Newton tolerance 1e-7, lift tolerance 2e-5 and 100 Newton iterations.
Use flaps 5 and 10 degrees and target CL 0.5, 0.8 and 1.2. Compare:

- Cold local starts, independently for each target: six targets in six jobs.
- Increasing target-lift order 0.5, 0.8, 1.2: two jobs, six targets.
- Decreasing target-lift order 1.2, 0.8, 0.5: two jobs, six targets.
- Exact archived warm-replay initialization for flap 10 / CL 0.8: one control.

Total: 11 independent jobs, 19 requested fixed-lift targets. The three repeated
paths cover the same six operating points, including both original failures and
the paired flap-5 baseline. The warm control checks the successful historical
replay. It is not selected preferentially using drag.

For the two originally failed flap-10 targets, use `local_seed()` on their
original computed lift records. For other targets use the archived selected
incidence. These are numerical seeds, not fits to measured drag. No flow state
or force is reused as a newly solved result. The original selected states and
all original failures stay unchanged.

## Existing policy and interpretation

Execute the existing audited `s702_ncrit_sweep.run_group()` worker unchanged,
with initial half-width 0.05 degrees and explicit start-at-seed for local cold
starts. The worker retains its three-attempt retry logic: bracket widening when
lift is not bracketed, otherwise cold restart with a narrower local bracket.
All endpoints, failures, retries, selected-state indices, internal continuation,
BL/wake exports and source snapshots are retained.

Increasing/decreasing means target order. Secant/bracket evaluations within a
target may reverse incidence. An automatic cold fallback is explicitly flagged;
a success after fallback is not evidence that the original warm path converged.
Compare first-attempt outcomes separately from final-policy outcomes. The warm
control replays the two originally admissible requested incidences, using the
original start alpha -1, before entering the same local bracket.

The source worker restores selected model/solution state before proceeding, as
in the archived driver. This does not add a disk-restorable Newton checkpoint or
claim general exact-checkpoint support. That behavior remains a subject for the
subsequent reusable-policy implementation.

Each target has at most three fixed-lift attempts, each with up to 20 requested
flow evaluations; warm replay adds two. The nominal outer bound is 1,142 requested
evaluations; internal continuation steps are additional and are logged. No
automatic extra sweep or closure change is launched. Use 1–8 worker processes,
one BLAS/OpenMP thread each, leaving capacity for existing server work.

## Restart, reporting and acceptance

Launch with `run_s702_repeatability.py --out PATH --workers N`. `--prepare-only`
writes the immutable manifest without flow solves. Resume with the same command;
only completed whole jobs are reused. An incomplete job, including its finished
targets, is moved to `interrupted_attempts/` and restarted from its original
initialization, preserving the intended increasing/decreasing path. Source or
runtime changes cause a hard stop; use a new output directory instead.

Produce incremental `status.json`, per-job logs/exports and a final comparison
with archived CL/CD/CM/incidence, transition and separation. Report all coverage
before accuracy, per-point path spans and paired flap drag changes. Failed targets
retain null forces. Missing paths prevent a complete repeatability score. No path
is chosen to minimize experimental drag error. Exact physical equivalence is not
asserted from scalar forces alone; inspect the selected states and path histories.

After the job finishes, audit selected-state and source compatibility, inspect
the comparison plots, then decide the smallest common-policy change. The run
itself is a numerical repeatability experiment, not new aerodynamic validation.
