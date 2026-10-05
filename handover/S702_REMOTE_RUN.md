# S702 repeatability run on jony4

## Current follow-up: common initialization policy, 21 September

The original matrix below is complete, retrieved and audited. Its findings are
in `validation/results/s702_repeatability_analysis/FINDINGS.md`. A separate
seven-job follow-up is deployed at:

```text
/home/jony/multifoil-runs/s702-policy-20260921
```

**Continuation launched after the baseline failure was reproduced.** The original
preflight passed 45/46 tests. The unchanged
`test_multiple_elements_order_and_far_limit` failed with CL order difference
1.1415574480699853e-6 versus 5e-7. The separate test in the untouched original
deployment reproduced exactly the same value and failure. Its logs are retained:

```bash
ssh jony4-anywhere 'cat /home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime/original_far_field.json'
ssh jony4-anywhere 'tail -n 30 /home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime/original_far_field.log'
```

This is classified as an existing limitation of the 100-chord far-field test,
not evidence of a policy regression. Its threshold and failed outcome remain
unchanged. The continuation permits only the declared S702 experiment after
checking both failure logs, all 92 original payload hashes, baseline core/test
hashes, the shared environment and its pinned package versions. It reuses the
seven unaffected viscous checks and reruns the 38 fast checks. No full-suite pass
is claimed. Classification evidence is in
`handover/policy_runtime/baseline_classification.json` on jony4.

### Follow completion from the laptop

```bash
cd /home/jony/multifoil
python3 handover/policy_status.py --watch
```

This refreshes every 10 seconds and exits when the launcher finishes or fails.
It displays completed jobs separately from converged targets, requested solves,
failed requests and initialization paths. Closing it or pressing Ctrl-C stops
only the monitor; the detached remote run continues. For a single check, omit
`--watch`; add `--json` for structured output. Exit 0 means a completed force
regression passed (or a one-shot check observed an active run); 2 means completed
with regression failures, 1 means execution stopped, and 3 means a one-shot
connection failure. Completion still requires scientific review of selected states.

This new directory contains an isolated copy of the baseline project and the new
policy/helper/tests. Its `.venv` links to the unchanged pinned environment of the
original run; do not remove that original environment while this follow-up is
in use. jony4 was idle with about 29 GiB available RAM at launch. Seven workers,
one BLAS/OpenMP thread each. The original launcher and its failed preflight are
retained; `handover/resume_s702_policy.py` performs the classified continuation.

Protocol: `validation/S702_INITIALIZATION_POLICY_PROTOCOL.md`. Latest launch
observation: `handover/policy_resume_launch_20260921.json`; the initial failed
attempt remains in `handover/policy_launch_20260921.json`. The local new-study folder
contains preparation files until results are retrieved; do not mistake it for
remote progress.

```bash
ssh jony4-anywhere 'cat /home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime/launch_status.json'
ssh jony4-anywhere 'cat /home/jony/multifoil-runs/s702-policy-20260921/validation/results/s702_initialization_policy_20260921/status.json'
ssh jony4-anywhere 'tail -n 40 /home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime/continuation.log'
```

Retrieve completed outputs without deleting original studies:

```bash
rsync -az --checksum jony4-anywhere:/home/jony/multifoil-runs/s702-policy-20260921/validation/results/s702_initialization_policy_20260921/ validation/results/s702_initialization_policy_20260921/
rsync -az jony4-anywhere:/home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime/ handover/policy_runtime_20260921/
```

Inspect `summary.json` and `comparison.csv`, distinguishing the six primary
conditions from the one replay control. Inspect each selected polar/state export
and all failed attempts before adopting the helper in shape optimization. A
passed force regression is not independent aerodynamic validation.

Only after an actual interruption, resume on the remote host with:

```bash
cd /home/jony/multifoil-runs/s702-policy-20260921
nohup setsid .venv/bin/python handover/resume_s702_policy.py --workers 7 </dev/null >>handover/policy_runtime/continuation.log 2>&1 &
```

Launcher/study locks prevent duplicates. Completed whole jobs are retained;
partial jobs are archived and restarted. Source/runtime changes require a fresh
study directory. Changed baseline evidence also stops continuation; the failed
physics test cannot be silently reclassified as passed.

## Original matrix: completed, retrieved and reviewed

Launch date: 16 September 2026. User authorized selecting the best available
server and launching the prepared study. See `repeatability_launch.json` for
the latest confirmed launch state and process IDs.

## Machine and deployment

Selected **jony4**, AMD Threadripper 3990X, 64 physical cores, idle at inspection,
29 GiB available memory and 386 GiB free disk. Eight worker processes, one
BLAS/OpenMP thread each. `server_capacity_20260916.json` records all three machines;
the selection is based on hardware/load, not a measured flow benchmark.

Launch used the existing key-based `jony4` LAN SSH alias. On 17 September, the
shared Tailscale policy was changed to accept non-root SSH to owned devices;
direct root SSH retains check mode. All three `*-anywhere` aliases passed fresh
unattended connections with strict host-key verification. Use `jony4-anywhere`
for status and retrieval, including when away from the LAN. Policy backup,
verification and rollback are in `/home/jony/remoteConnection/README.md`.

Latest observation: the matrix finished at **16 September 12:48:40 UTC** with
11/11 jobs completed and 18/19 targets converged. The launcher reports
`finished_with_numerical_failures` (exit 2). See
`remote_status_observed_20260917.json`. Retrieval and numerical repeatability
review were completed on 21 September; do not relaunch this matrix. This review
does not establish independent aerodynamic accuracy.

Remote project root:

```text
/home/jony/multifoil-runs/s702-repeatability-20260916
```

The complete project was copied into this isolated directory. The launcher checks
87 source/baseline/manifest files against `remote_payload_20260916.json`. The
prepared immutable study lives at `validation/results/s702_repeatability_20260915/`;
its name retains the protocol preparation date. The runner checks the study and
source hashes again before any solve.

The environment is isolated under the remote root: `.tools/uv`, `.python`,
`.venv`, and `.cache/uv`. Python 3.12.13 is used because the available uv 0.12.1
catalogue has no 3.12.14 build; 3.12.13 also appears in the historical S702 metadata.
The study requires Python 3.12 and exact NumPy 2.3.5, SciPy 1.17.0 and Matplotlib
3.10.8. The failed first setup attempt at 3.12.14 launched no flow solves and is
preserved in `handover/remote_runtime/bootstrap_attempt_001/` and the append-only
launch log. All 25 preflight tests run in the new environment before the matrix.

## Inspect progress

From the laptop, inspect the launcher stage and log:

```bash
ssh jony4-anywhere 'cat /home/jony/multifoil-runs/s702-repeatability-20260916/handover/remote_runtime/launch_status.json'
ssh jony4-anywhere 'tail -n 40 /home/jony/multifoil-runs/s702-repeatability-20260916/handover/remote_runtime/launch.log'
```

Once the launcher says `running_study`, inspect matrix progress:

```bash
ssh jony4-anywhere 'cat /home/jony/multifoil-runs/s702-repeatability-20260916/validation/results/s702_repeatability_20260915/status.json'
```

Each of the 11 job directories has a `run.log`, incremental target `progress.json`
and flow exports. Completed targets have `result.json`; completed jobs have their
own group `result.json`. The matrix contains 19 targets, not 19 Newton evaluations.
The launcher was detached with `start_new_session=True` and all standard streams
redirected, so it does not depend on the interactive SSH session.

## Collect and continue after completion

From `/home/jony/multifoil`, copy results and runtime logs back without deleting
or replacing unrelated archived studies:

```bash
rsync -az jony4-anywhere:/home/jony/multifoil-runs/s702-repeatability-20260916/validation/results/s702_repeatability_20260915/ validation/results/s702_repeatability_20260915/
rsync -az jony4-anywhere:/home/jony/multifoil-runs/s702-repeatability-20260916/handover/remote_runtime/ handover/remote_runtime_20260916/
```

Inspect coverage first, then `comparison.csv`, `path_spans.csv`,
`paired_flap_deltas.csv`, `boundary_layer_comparison.csv` and `summary.json`.
Separate first-attempt convergence from cold fallback success. Audit indexed
selected states and plot comparisons before asserting force/state equivalence or
changing the common policy. No post-run aerodynamic conclusion is part of launch.

## Restart after interruption

The launcher holds a lock, and the study runner holds a separate study lock.
Do not start a second copy while either is running. On the remote machine, the
same launcher reuses its isolated environment, reruns preflight tests, checks
runtime/source compatibility and resumes the study:

```bash
cd /home/jony/multifoil-runs/s702-repeatability-20260916
nohup setsid bash handover/launch_s702_repeatability.sh 8 </dev/null >>handover/remote_runtime/bootstrap_startup.log 2>&1 &
```

Completed whole jobs are retained. An incomplete job, including any completed
targets within it, is archived under `interrupted_attempts/` and restarted from
its original initialization. This preserves the target-order experiment; it is
not a restored Newton checkpoint. A changed source or runtime requires a fresh
study directory. Numerical failures remain explicit and do not trigger a new
sweep or a change of model.
