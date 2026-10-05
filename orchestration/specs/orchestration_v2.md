# Codex Orchestration v2 — Bounded Task Capsules, Hard Usage Controls, and Manual Launch Approval

## Purpose

This document supersedes the earlier orchestration brief.

The first orchestration pilot demonstrated that a multi-agent approach can consume extremely large amounts of model usage even when most input tokens are cached. The main failure was not the engineering task itself, but the architecture of the control plane:

- the main conversation remained an expensive model-driven coordinator;
- several agents inherited the full parent conversation;
- long-lived agents accumulated growing tool and shell histories;
- corrective follow-ups resumed those large contexts;
- broad AI reviews repeated substantial repository inspection;
- an unavailable model triggered an upward fallback to a much more expensive model;
- timeouts existed, but no hard token/tool/model budget governed model work.

The revised design must make **model calls scarce, explicit, isolated, bounded, and auditable**.

The core principle is:

> **Automate everything between AI calls. Do not automate away visibility of the AI calls themselves until the system has proved that its budgets and isolation work correctly.**

During the initial debug phase, **every Codex command launch must require explicit human approval before execution**.

---

# 1. Target workflow

The desired workflow is:

```text
Human objective
      │
      ▼
Planner
one bounded planning session
      │
      ▼
tasks.yaml
      │
      ▼
Human reviews plan
      │
      ▼
PURE PYTHON ORCHESTRATOR
      │
      ├─ resolve task dependencies
      ├─ build minimal task capsule
      ├─ estimate/check capsule size
      ├─ select configured model + effort
      ├─ prepare exact Codex command
      ├─ show expected scope and budget
      │
      ▼
HUMAN APPROVES THIS CODEX LAUNCH
      │
      ▼
fresh Codex process
small isolated context
hard limits
      │
      ▼
deterministic validation
      │
      ▼
result / patch / evidence
      │
      ▼
human reviews or approves next launch
```

The orchestrator itself must not use an LLM.

---

# 2. Fundamental design rules

## 2.1 The control plane is deterministic

The following must be implemented in ordinary Python/shell/git logic rather than by an LLM:

- task scheduling;
- dependency resolution;
- capsule construction;
- copying/symlinking files;
- process launching;
- process supervision;
- timeout handling;
- token accounting;
- tool-call accounting where observable;
- command-count accounting;
- validation command execution;
- status transitions;
- artifact collection;
- log collection;
- process completion checks;
- CFD/process monitoring;
- branch/worktree management;
- diff generation;
- result packaging.

An LLM may **design** a task graph or solve a bounded task.

An LLM must not remain active as the scheduler, monitor, or runtime coordinator.

---

## 2.2 No inherited conversation context

Every executor or reviewer starts as a **fresh Codex process**.

Forbidden by default:

```text
fork_turns="all"
spawn_agent inheriting parent conversation
follow-up continuation of a large executor context
```

Each call receives only its explicit task capsule.

If a correction is needed, prefer:

```text
large failed run
      ↓
small deterministic result/failure summary
      ↓
new small correction capsule
      ↓
fresh Codex run
```

Do not resume the large previous context.

---

## 2.3 Task capsules are the default workspace

A task should normally not receive the full repository.

Instead, build a small directory containing only what the task needs.

Example:

```text
.task_capsules/T17/
├── TASK.yaml
├── TASK.md
├── AGENTS.md
├── work/
│   ├── scripts/
│   │   └── postprocess_cp.py
│   └── tests/
│       └── test_postprocessing.py
├── reference/
│   ├── refinement_2.csv
│   ├── refinement_3.csv
│   ├── refinement_4.csv
│   ├── refinement_5.csv
│   └── T12_result.md
└── metadata/
    └── capsule_manifest.json
```

The capsule is deliberately incomplete.

If required information is absent, Codex must report `BLOCKED` or `MISSING_CONTEXT`.

It must not be allowed to wander through the full repository to discover additional context.

---

# 3. Workspace levels

Support three workspace levels.

## Level A — Tiny capsule

Default.

Typical contents:

- task specification;
- 2–20 relevant source files;
- relevant tests;
- small reference data;
- concise dependency summaries;
- minimal build/configuration files if needed.

Use for:

- bug fixes;
- plotting;
- scripts;
- local algorithms;
- bounded engineering implementation;
- numerical analysis based on preselected datasets.

## Level B — Sparse project

Use when normal package/import/build structure is necessary.

Expose selected subtrees, for example:

```text
src/aero/
tests/aero/
validation/naca0012/
pyproject.toml
```

This may be implemented with deterministic copying, symlinks, or a git sparse checkout.

Use only when Level A is insufficient.

## Level C — Full repository

Exceptional.

Use only for tasks such as:

- architecture-wide analysis;
- large refactoring;
- global dependency inspection;
- integration review that genuinely requires full source access.

Level C must require explicit human approval and should have tighter model and token controls.

---

# 4. Read/write separation

Prefer a capsule structure that makes write scope obvious:

```text
capsule/
├── work/        # editable
├── reference/   # read-only source/reference material
├── TASK.yaml
└── AGENTS.md
```

Task metadata should define:

```yaml
workspace:
  level: A

  writable:
    - work/scripts/postprocess_cp.py
    - work/tests/

  readonly:
    - reference/refinement_data/
    - reference/T12_result.md
```

After execution:

1. collect modifications under `work/`;
2. generate a patch against the real repository;
3. run validation;
4. present the patch/result;
5. apply/merge only after the appropriate approval.

Do not let Codex directly edit the main working copy by default.

---

# 5. Planner contract

The planner should produce a task graph, not execute it.

Each task should contain enough information for deterministic capsule construction.

Example:

```yaml
tasks:
  - id: T17
    title: Diagnose mesh convergence

    depends_on:
      - T12

    task_class: investigation
    model_class: reasoning
    effort: high

    objective: >
      Determine why Cp stops converging after refinement level 4.

    context:
      files:
        - scripts/postprocess_cp.py
        - tests/postprocessing/**
        - cases/refinement/summary/*.csv

      dependency_artifacts:
        - task: T12
          artifact: result.md

      exclude:
        - conversation_history
        - unrelated_docs
        - unrelated_source
        - full_git_history

    workspace:
      level: A

    actions:
      - compare refinement levels 2-5
      - inspect y+, skewness and non-orthogonality
      - distinguish discretization error from solver convergence
      - propose the next discriminating experiment

    acceptance:
      - convergence plot generated
      - numerical evidence supports the diagnosis
      - recommendation written to docs/validation.md

    invariants:
      - baseline numerical data must remain unchanged
      - solver source code must not be modified

    may_modify:
      - scripts/postprocess_cp.py
      - tests/postprocessing/
      - docs/validation.md

    validation:
      - python scripts/postprocess_cp.py
      - pytest tests/postprocessing

    stop_and_escalate_if:
      - required solver modification is discovered
      - expected physical interpretation conflicts with evidence
      - necessary input data is absent

    expected_artifacts:
      - convergence.png
      - refinement_metrics.csv
      - result.json
```

The planner should prefer tasks that are independently executable and independently verifiable.

---

# 6. Human approval gate — mandatory during debug phase

This is a hard requirement.

## 6.1 No Codex process may launch automatically

Before **every** model-executing Codex invocation, the orchestrator must stop and show a launch preview.

For example:

```text
============================================================
CODEX LAUNCH REQUEST
============================================================

Task:           T17 — Diagnose mesh convergence
Purpose:        Executor
Workspace:      Level A capsule
Capsule path:   .task_capsules/T17/
Files exposed:  11
Approx size:    184 KB
Approx LOC:     3,120

Model:          configured-balanced-model
Effort:         high

Hard limits:
  Wall time:       300 s
  Input budget:    250,000 tokens
  Output budget:    10,000 tokens
  Shell commands:  12

Writable:
  scripts/postprocess_cp.py
  tests/postprocessing/
  docs/validation.md

Command:
  codex exec ...

Launch this Codex process? [y/N]
```

Default answer must be **No**.

Only explicit `y` / `yes` launches the process.

## 6.2 Approval applies to every type of model call

The gate applies to:

- executors;
- reviewers;
- retries;
- corrective runs;
- planners launched by the orchestrator;
- integration reviews;
- model escalation;
- fallback to another model.

No hidden/background AI launch is permitted in debug mode.

## 6.3 No chained launches after one approval

Approval is per process.

If T17 succeeds and makes T18 ready, the orchestrator may prepare T18 but must stop again:

```text
T17: PASS
T18: READY

Codex launch required for T18.
Approve? [y/N]
```

One approval must never authorize a chain of subsequent model calls.

---

# 7. Model policy

The orchestrator should use an explicit allowlist.

Example:

```yaml
model_policy:
  default_model_class: balanced
  default_effort: medium

  automatic_upward_fallback: false
  automatic_model_escalation: false

  require_explicit_approval_for_every_launch: true

  classes:
    light:
      model: <configured light model>
      max_effort: medium

    balanced:
      model: <configured balanced model>
      max_effort: high

    reasoning:
      model: <configured reasoning model>
      max_effort: high

    exceptional:
      model: <expensive frontier model>
      max_effort: xhigh
      explicit_warning: true
```

The exact model names must be discovered from the local Codex configuration/version rather than hard-coded from assumptions.

## 7.1 Upward fallback is forbidden

If a requested model is unavailable:

```text
MODEL_UNAVAILABLE
       ↓
STOP
```

Do not silently substitute a more expensive model.

The orchestrator may present alternatives to the human, but launching one requires new approval.

## 7.2 `xhigh` and expensive frontier models are exceptional

During the debug phase, any launch involving:

- `xhigh`;
- the most expensive available model;
- a Level C full-repository workspace;

must display a conspicuous warning.

Then require explicit approval.

---

# 8. Hard usage controls

Timeouts alone are insufficient.

Every Codex process must have multiple independent limits.

Example initial configuration:

```yaml
budgets:

  executor:
    max_wall_seconds: 300
    max_input_tokens: 250000
    max_output_tokens: 10000
    max_shell_commands: 12

  focused_reviewer:
    max_wall_seconds: 240
    max_input_tokens: 150000
    max_output_tokens: 6000
    max_shell_commands: 8

  exceptional_task:
    max_wall_seconds: 600
    max_input_tokens: 500000
    max_output_tokens: 15000
    max_shell_commands: 20
```

These values are initial debug limits, not permanent truths.

Tune them from measured runs.

---

# 9. Live usage governor

Where the installed Codex CLI exposes usage events through JSON/JSONL output, the Python launcher should parse them while the process is running.

Track, when available:

- total input tokens;
- cached input tokens;
- uncached input tokens;
- output tokens;
- reasoning output tokens;
- number of model generations/turns;
- shell/tool command count;
- wall time.

Pseudo-logic:

```python
while codex_process_running:
    event = read_next_json_event()

    update_usage(event)

    if input_tokens > max_input_tokens:
        terminate("INPUT_TOKEN_BUDGET_EXCEEDED")

    if output_tokens > max_output_tokens:
        terminate("OUTPUT_TOKEN_BUDGET_EXCEEDED")

    if shell_commands > max_shell_commands:
        terminate("COMMAND_BUDGET_EXCEEDED")

    if elapsed > max_wall_seconds:
        terminate("TIME_BUDGET_EXCEEDED")
```

If the installed CLI version does not expose reliable live token usage, record that limitation explicitly and rely temporarily on:

- wall time;
- shell-command count;
- capsule size;
- model restriction;
- per-launch human approval.

Do not pretend a token ceiling is enforced if it is not technically observable.

---

# 10. Capsule size guard

Before presenting the launch approval, calculate deterministic capsule metrics.

At minimum:

```text
files
total bytes
source LOC
largest files
reference-data size
```

Optionally estimate prompt/token size.

Refuse or warn on oversized capsules.

Example:

```yaml
capsule_limits:
  level_A:
    max_files: 30
    max_source_lines: 10000
    max_total_bytes: 2000000

  level_B:
    max_files: 150
    max_source_lines: 50000
    max_total_bytes: 15000000
```

Do not automatically promote an oversized Level A capsule to Level B or Level C.

Report it and require a decision.

---

# 11. Missing context behavior

A capsule is intentionally restrictive.

If Codex discovers that it needs a missing file, dependency, dataset, or API contract, it should not search outside the capsule.

Standard result:

```json
{
  "task_id": "T17",
  "status": "MISSING_CONTEXT",
  "requested_context": [
    "src/common/coefficients.py"
  ],
  "reason": "postprocess_cp.py imports this module and its behavior is required to interpret Cp."
}
```

The orchestrator then:

1. shows the request to the human;
2. optionally adds the approved file;
3. builds a new capsule;
4. presents a **new Codex launch approval**.

Do not resume the previous large conversation.

---

# 12. Codex invocation

Prefer a fresh non-interactive `codex exec` invocation.

Conceptually:

```bash
codex exec   --model <approved-model>   --sandbox workspace-write   --ask-for-approval never   --json   -
```

The outer Python orchestrator provides the human approval gate.

`--ask-for-approval never` refers only to Codex-internal tool approval inside the already-approved capsule. It does **not** mean that the orchestrator can launch Codex without human approval.

Adapt flags to the installed Codex CLI version.

Do not assume unsupported flags or model names.

---

# 13. Executor prompt

The executor receives only:

- `TASK.yaml` / `TASK.md`;
- capsule-local `AGENTS.md`;
- allowed work files;
- reference files;
- concise dependency artifacts.

Suggested execution rules:

```text
You are executing one bounded engineering task.

1. Work only inside this capsule.
2. Do not seek files outside the capsule.
3. Modify only paths explicitly marked writable.
4. Do not weaken acceptance criteria.
5. Do not weaken or remove tests merely to obtain a pass.
6. Do not silently change numerical or physical assumptions.
7. Run the specified validation commands.
8. Base conclusions on evidence.
9. If required context is missing, return MISSING_CONTEXT.
10. If the task cannot be completed within scope, return BLOCKED.
11. Do not claim PASS unless every acceptance criterion has evidence.
12. Keep exploration minimal and relevant.
```

---

# 14. Standard result

Prefer structured output.

Example:

```json
{
  "task_id": "T17",
  "status": "PASS",
  "summary": "Solver convergence rather than discretization dominates Cp variation after refinement level 4.",
  "changed_files": [
    "scripts/postprocess_cp.py",
    "docs/validation.md"
  ],
  "validation": [
    {
      "command": "python scripts/postprocess_cp.py",
      "status": "PASS"
    },
    {
      "command": "pytest tests/postprocessing",
      "status": "PASS"
    }
  ],
  "evidence": [
    "artifacts/convergence.png",
    "artifacts/refinement_metrics.csv"
  ],
  "unresolved": [],
  "missing_context": []
}
```

Store per-run outputs under:

```text
runs/T17/<attempt-id>/
├── task.yaml
├── capsule_manifest.json
├── prompt.txt
├── codex_command.txt
├── stdout.jsonl
├── stderr.log
├── usage.json
├── result.json
├── patch.diff
├── validation.log
└── artifacts/
```

---

# 15. Deterministic validation first

After Codex exits, the Python orchestrator—not another LLM—should run the declared validation commands.

Examples:

```text
pytest
type checking
lint
comparison scripts
numerical tolerances
file existence
schema validation
reference-case comparison
mass/conservation checks
grid-convergence checks
```

A model's statement that tests passed is not sufficient.

The orchestrator must execute or independently verify them.

---

# 16. Three levels of correctness

For scientific/engineering tasks, acceptance may require:

## Level 1 — Software correctness

- builds/imports/runs;
- unit/integration tests;
- no obvious regression.

## Level 2 — Numerical correctness

- convergence;
- conditioning;
- conservation;
- tolerance/grid/time-step sensitivity;
- robustness to initial conditions when relevant.

## Level 3 — Physical/model correctness

- agreement with theory;
- experiment;
- trusted benchmark;
- valid physical assumptions.

The planner must specify which levels apply.

Do not use an AI reviewer as a substitute for deterministic Level 1/2 checks where executable checks are possible.

---

# 17. AI review policy

AI review is **not automatic by default**.

Use it only when deterministic validation leaves meaningful uncertainty, for example:

- subtle numerical algorithms;
- physics/model changes;
- architecture changes;
- surprising result despite passing tests;
- security-sensitive changes;
- high-impact engineering decisions.

Before a reviewer launches, present the same human approval gate as any other Codex process.

## 17.1 Reviewer receives a small review packet

Do not give a reviewer the full repository by default.

Prepare:

```text
review_packet/
├── TASK.yaml
├── patch.diff
├── validation.json
├── relevant_source/
├── relevant_tests/
└── evidence/
```

Reviewer instruction:

```text
Review only this bounded packet.

Try to falsify the claimed result.

Check:
- acceptance criteria;
- unintended scope changes;
- weakened tests;
- numerical errors;
- unit/sign/reference-frame errors;
- hidden assumptions;
- unsupported engineering conclusions.

Do not perform a general repository review unless a specific finding requires additional context.
```

If more context is required, return `MISSING_CONTEXT`.

## 17.2 Corrective review is narrow

If the first review reports three findings and they are corrected, the second review must default to:

> Verify whether findings 1–3 are fixed.

Do not automatically repeat a comprehensive review.

A second reviewer launch requires separate human approval.

---

# 18. Retries and correction runs

Automatic retries are disabled during debug.

If a task fails:

```text
FAIL / BLOCKED / MISSING_CONTEXT
             │
             ▼
Python prepares failure summary
             │
             ▼
Human decides
```

If a correction run is desired:

- construct a new compact task;
- include only the failure evidence and relevant files;
- start a fresh Codex process;
- require new approval.

Do not use a long-lived follow-up on a large existing agent context.

---

# 19. Dependency graph

The orchestrator may automatically determine which tasks are `READY`, but readiness does not imply launch.

Example:

```text
T12 PASS
   │
   ├──── T17 READY
   └──── T18 READY
```

The orchestrator can show:

```text
READY:
  T17 — estimated capsule 180 KB — balanced/high
  T18 — estimated capsule 95 KB  — light/medium
```

The human chooses what to launch.

During debug mode, even parallel execution must require approval for each individual Codex process.

---

# 20. Parallelism

Initial default:

```yaml
max_parallel_codex_processes: 1
```

This is intentional.

First establish:

- capsule quality;
- actual token usage;
- process accounting;
- budget enforcement;
- correctness of launch gating.

Only then consider `2`.

Never start a large fan-out of AI workers merely because several DAG nodes are ready.

---

# 21. Monitoring long-running simulations

Do not use an LLM to repeatedly monitor CFD or other long-running jobs.

Use deterministic process logic:

```text
PID
scheduler status
log tail
output file timestamps
completion marker
exit code
residual parser
```

Only invoke a model when a human intentionally requests interpretation of a completed or problematic state.

---

# 22. Usage ledger

Every run must record usage in a central machine-readable ledger.

Example:

```json
{
  "timestamp": "2026-10-05T20:00:00+02:00",
  "task_id": "T17",
  "attempt": 1,
  "purpose": "executor",
  "model": "configured-model",
  "effort": "medium",
  "workspace_level": "A",
  "capsule_files": 11,
  "capsule_bytes": 184000,
  "capsule_loc": 3120,
  "wall_seconds": 84.2,
  "input_tokens": 84000,
  "cached_input_tokens": 51000,
  "output_tokens": 4200,
  "shell_commands": 6,
  "status": "PASS"
}
```

Aggregate reports should show:

- usage by task;
- usage by model;
- usage by effort;
- usage by purpose;
- retries;
- reviewer overhead;
- cached vs uncached input;
- capsule size vs usage;
- pass/fail/block rates.

The objective is to learn which task categories genuinely need expensive reasoning.

---

# 23. Debug-mode configuration

Initial mode should be deliberately conservative.

Example:

```yaml
debug_mode: true

launch_policy:
  require_human_approval_for_every_codex_process: true
  default_answer: no
  allow_chained_launches: false
  automatic_retry: false
  automatic_review: false
  automatic_model_escalation: false
  automatic_upward_fallback: false

context_policy:
  inherit_parent_conversation: false
  continue_previous_agent_context: false
  full_repository_default: false

execution_policy:
  max_parallel_codex_processes: 1

model_policy:
  xhigh_requires_explicit_warning: true
  frontier_model_requires_explicit_warning: true

workspace_policy:
  default_level: A
  full_repository_requires_explicit_warning: true
```

---

# 24. Suggested CLI

Keep the interface simple.

```bash
# inspect graph and readiness
python tools/orchestrate.py status

# build capsule without launching anything
python tools/orchestrate.py prepare T17

# inspect capsule
python tools/orchestrate.py inspect T17

# request launch; command is shown and confirmation requested
python tools/orchestrate.py run T17

# deterministic validation only
python tools/orchestrate.py validate T17

# prepare a review packet
python tools/orchestrate.py prepare-review T17

# launch reviewer, again with explicit approval
python tools/orchestrate.py review T17

# show usage ledger
python tools/orchestrate.py usage
```

`prepare` must never execute Codex.

`run` must never skip confirmation while debug mode is enabled.

---

# 25. Suggested state machine

```text
PENDING
   │ dependencies satisfied
   ▼
READY
   │
   ▼
CAPSULE_PREPARED
   │
   ▼
AWAITING_LAUNCH_APPROVAL
   │ human yes
   ▼
RUNNING
   │
   ▼
VALIDATING
   │
   ├── PASS
   ├── FAIL
   ├── BLOCKED
   ├── MISSING_CONTEXT
   └── BUDGET_EXCEEDED
```

Optional review states:

```text
PASS
  │
  ▼
REVIEW_AVAILABLE
  │
  ▼
AWAITING_REVIEW_APPROVAL
  │ human yes
  ▼
REVIEWING
```

---

# 26. Safety against another usage explosion

The system must specifically prevent recurrence of the failed pilot pattern.

The following should be impossible by default:

```text
main frontier/xhigh conversation acts as coordinator
        +
3 agents inherit entire parent conversation
        +
follow-up tasks resume giant contexts
        +
reviewers inspect full repository
        +
second broad reviewer repeats everything
        +
model failure silently escalates to frontier model
        +
no token/tool budget
```

Instead enforce:

```text
one bounded planner
        ↓
planner stops
        ↓
small deterministic capsule
        ↓
human approval
        ↓
fresh Codex process
        ↓
hard budget
        ↓
deterministic validation
        ↓
compact result
```

---

# 27. Implementation order

Do not build the full system at once.

## Phase 0 — Audit installed Codex capabilities

Before implementation, inspect the local Codex CLI and determine:

- exact version;
- supported `exec` flags;
- available model-selection mechanism;
- JSON/JSONL event format;
- whether live token usage is available;
- how final usage is reported;
- available reasoning-effort values;
- output-schema support.

Document findings.

Do not launch expensive test workloads merely to discover this.

## Phase 1 — Manual-gated single-task capsule runner

Implement only:

```bash
python tools/orchestrate.py prepare T17
python tools/orchestrate.py run T17
```

Required behavior:

- parse `tasks.yaml`;
- build Level A capsule;
- compute capsule metrics;
- construct exact Codex command;
- show launch preview;
- require `y/N`;
- launch a fresh Codex process;
- record logs and usage;
- enforce available hard limits;
- collect result;
- run deterministic validation.

This is the most important milestone.

## Phase 2 — Usage governor

Before DAG automation, ensure:

- usage events are captured;
- token limits work if technically available;
- shell-command limits work;
- wall-time limits work;
- budget termination is tested;
- ledger is correct.

Create unit/integration tests for all budget exits.

## Phase 3 — Missing-context workflow

Implement:

- `MISSING_CONTEXT`;
- requested-file reporting;
- human-approved capsule expansion;
- fresh re-run.

## Phase 4 — DAG readiness

Add:

- dependency parsing;
- READY state;
- downstream blocking.

Still require manual approval for every launch.

## Phase 5 — Optional review packets

Only after executor usage is under control.

Add focused reviewer support with the same per-launch approval.

## Phase 6 — Consider reducing manual gating

Only after enough real runs have demonstrated that:

- capsule sizes are predictable;
- model choice is predictable;
- usage governor works;
- no silent escalation exists;
- typical usage is acceptable.

Any future automatic-launch mode should be opt-in.

Debug mode must remain available.

---

# 28. Tests for the orchestrator itself

The orchestration layer must have ordinary automated tests.

At minimum test:

- task parsing;
- dependency resolution;
- capsule file selection;
- writable/read-only separation;
- excluded files are absent;
- path traversal is rejected;
- oversized capsule warning/refusal;
- default approval is `No`;
- no process starts after `No`;
- one `Yes` launches exactly one Codex process;
- second task requires a second approval;
- unavailable model does not upward-fallback;
- retry does not happen automatically;
- reviewer does not happen automatically;
- timeout termination;
- command-count termination;
- token-budget termination if supported;
- result parsing;
- usage ledger writing;
- missing-context workflow;
- deterministic validation execution;
- patch generation;
- main repository remains unchanged before explicit application.

Use mocks/fake subprocesses for most tests.

Do not consume real model usage to test ordinary scheduler logic.

---

# 29. Success criteria

The revised orchestration is successful when a normal interaction looks like:

```text
1. Planner writes tasks.yaml.
2. I inspect/edit it.
3. Python builds a small capsule for T17.
4. Python shows me:
   - files exposed;
   - workspace level;
   - model;
   - effort;
   - hard limits;
   - exact Codex command.
5. I explicitly approve or reject that launch.
6. A fresh bounded Codex process runs.
7. Python validates the result.
8. I see:
   - outcome;
   - evidence;
   - patch;
   - actual usage.
9. No further Codex process starts until I explicitly approve it.
```

And critically:

> A task should not be able to consume millions of accumulated model-input tokens merely because an agent repeatedly explores a growing context.

---

# 30. Instructions to Codex implementing this brief

When you receive this document:

1. Inspect the existing orchestration code and repository structure.
2. Do not launch subagents to implement the orchestrator.
3. Do not use `fork_turns="all"` or equivalent inherited-context workers.
4. Do not launch additional Codex processes while designing the implementation unless the human explicitly approves each launch.
5. Prefer ordinary Python, YAML, subprocess, pathlib, git, and small utilities.
6. Reuse existing code only where it is simpler and consistent with this design.
7. Remove or disable model-driven scheduling/monitoring from the critical path.
8. Make task capsules the default workspace.
9. Make the manual per-launch approval gate impossible to bypass while `debug_mode: true`.
10. Implement usage limits before adding automatic DAG execution.
11. Implement tests using mocks rather than real Codex calls wherever possible.
12. Verify assumptions against the locally installed Codex CLI.
13. If a requested hard token limit cannot be enforced with the local CLI, state that explicitly and implement the strongest observable fallback rather than pretending it is enforced.
14. Do not silently choose a more expensive model if the requested model is unavailable.
15. Do not implement autonomous parallel multi-agent behavior as part of the first version.

Before modifying code, produce a concise implementation plan containing:

- files to add/change;
- current functionality to remove or disable;
- capsule representation;
- approval-gate implementation;
- budget enforcement strategy;
- tests;
- any Codex CLI limitations discovered.

The target is a **small, inspectable, deterministic orchestration layer around deliberately bounded Codex calls**, not an autonomous-agent framework.
