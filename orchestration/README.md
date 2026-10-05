# S414 orchestration pilot

`pilot.py` connects a frozen six-case protocol to the existing scientific
runners. It is a linear workflow, with a fresh independent Codex review before
remote execution and after verification. It does not change solver physics or
retry failed scientific cases. No Git repository is required: every revision
has an isolated source snapshot and explicit input hashes.

The user authorizes the six-case budget once. The coordinator then enforces two
reserved remote workers across runs in this project. The remote adapter must
also enforce two actual workers and the six-case task list; metadata alone
cannot constrain an arbitrary program.

## Commands

For this pilot, first generate its explicit source/input list with
`python orchestration/build_s414_plan.py`. This writes
`orchestration/s414_nose_resolution.yaml`; it reads archived control receipts
but runs no tests, simulations or model sessions. Review the plan before freezing.

```sh
python orchestration/pilot.py freeze orchestration/PLAN.yaml --run RUN_DIRECTORY
python orchestration/pilot.py prepare --run RUN_DIRECTORY
python orchestration/pilot.py review --run RUN_DIRECTORY --phase prelaunch
python orchestration/pilot.py submit --run RUN_DIRECTORY
python orchestration/pilot.py status --run RUN_DIRECTORY
python orchestration/pilot.py poll --run RUN_DIRECTORY
python orchestration/pilot.py collect --run RUN_DIRECTORY
python orchestration/pilot.py audit --run RUN_DIRECTORY
python orchestration/pilot.py review --run RUN_DIRECTORY --phase postanalysis
```

`advance --run RUN_DIRECTORY` performs exactly one next operation. After
submission, `watch --run RUN_DIRECTORY --interval 60 --max-hours 3` can run
detached with stdout/stderr saved to a monitor log. It polls, collects, audits
and starts the final independent review. It cannot submit or extend the
campaign. It stops on command failure, evidence mismatch or a terminal review
decision. Remote execution remains independent if this local monitor stops;
restart it to reconcile the recorded job. Do not keep a Codex session waiting.

`state.json` is the concise status view. Each operation retains an immutable
`invocation.json`, prompt when applicable, stdout/stderr logs, receipt and
completion record under `attempts/STEP/NNN/`. The frozen plan and source
manifest, invocation hashes, result hashes and evidence hashes are checked
before subsequent operations. The original working tree is never executed by
the coordinator's commands; commands use the snapshot as their working directory.
The coordinator itself must match its frozen source hash. If the working copy
later changes, use `RUN_DIRECTORY/snapshot/orchestration/pilot.py` to inspect or
resume the old run; never relax its source guard.

Execution outcome, artifact verification, scientific criteria, and independent
review are separate state fields. A verified experiment that fails its numerical
convergence criterion still proceeds to interpretation and review. It must not
be retried until a preferred physical answer appears.

## Frozen plan format

The plan may be JSON (also valid YAML) or YAML when PyYAML is installed.

```json
{
  "schema": 1,
  "id": "s414-nose-resolution-r1",
  "project_root": "/home/jony/multifoil",
  "resources": {"remote_workers": 2, "max_final_cases": 6},
  "files": ["orchestration/pilot.py", "validation/runner.py"],
  "criteria": [{"id": "resolution", "description": "Frozen criterion and tolerance"}],
  "reviewer": {"effort": "high", "timeout_seconds": 900},
  "commands": {
    "prepare": {"argv": ["python", "validation/runner.py", "prepare", "--result", "{result}"], "timeout_seconds": 120},
    "submit": {"argv": ["python", "validation/runner.py", "submit", "--job-id", "{job_id}", "--workers", "{workers}", "--result", "{result}"], "timeout_seconds": 120},
    "poll": {"argv": ["python", "validation/runner.py", "poll", "--job-id", "{job_id}", "--result", "{result}"], "timeout_seconds": 60},
    "collect": {"argv": ["python", "validation/runner.py", "collect", "--job-id", "{job_id}", "--result", "{result}"], "timeout_seconds": 300},
    "audit": {"argv": ["python", "validation/runner.py", "audit", "--result", "{result}"], "timeout_seconds": 300}
  }
}
```

The filenames above illustrate the contract; they are not a runnable study.
Snapshot entries can also be objects with `source` and relative `destination`.
Use explicit files, including all transitive source and reference inputs. No
mutable-directory copying or implicit current-source dependency is supported.
Commands are argv lists, never shell strings; available placeholders are
`{snapshot}`, `{run}`, `{result}`, `{job_id}` and `{workers}`. Outputs belong
under `{run}`. For scientific tools with existing relative-path conventions,
the frozen plan may declare exact `writable_snapshot_paths` (new result-directory
paths relative to the snapshot). New files outside those directories are
rejected; all frozen input files remain protected even inside a writable tree.
Commands must write their JSON
receipt to `{result}`. Subprocess timeouts bound each operation. BLAS/OpenMP
threads are set to one. This is not an OS sandbox for command adapters: they
must be reviewed trusted programs. Codex reviewers use the read-only sandbox.

## Adapter receipt contract

- Prepare: `status: prepared`, plus nonempty `evidence` array.
- Submit/poll: `job_id` and `status: running|complete|failed|unknown`.
- Collect: `job_id`, `status: collected`, plus nonempty `evidence` array.
- Audit: `evidence_verified: true|false`, nonempty `evidence`, and `criteria`.

Every artifact in `evidence` is an object with `path` (absolute or relative to
the run directory) and `sha256`. Every criterion has `id`,
`outcome: pass|fail|inconclusive|not_evaluated`, and an `evidence` list of paths or
measurement references. The audit and reviewers must report each frozen ID
exactly once. Failed solves must remain failed/missing measurements.

Submission first writes a persistent `SUBMISSION_PENDING` state and reserves
resources. If the reply is lost or the coordinator crashes, **poll the same job
ID**. A second submit is rejected. The remote adapter must likewise reconcile
the same ID, payload hash and scheduler/process identity before starting work.
An unknown state retains its worker reservation. A confirmed terminal remote
receipt releases resources; failed execution still proceeds to collection and
audit when usable evidence exists. No automatic retries or refinement follow.

If a protocol/source correction is required, preserve the rejected run and
freeze a new revision. Do not edit its plan, receipts, source hashes or acceptance
criteria in place. Do not delete the resource ledger while a remote job is
unresolved. Command failure details and the resume point remain in `state.json`.

## Coordinator verification

```sh
python -m unittest discover -s orchestration -p test_pilot.py -v
```

These checks use temporary directories and fake local commands only. They cover
submission ambiguity, budget enforcement, immutable inputs/receipts, evidence
hashes, independent review gating, timeouts and the distinction between a failed
scientific criterion and failed evidence verification. They do not validate the
aerodynamic runner or impersonate an actual independent Codex review.
