Continue development of this MultiFoil project for reliable incompressible
multi-element shape optimization. Read AGENTS.md, CODEX_HANDOVER.md and
DEVELOPMENT_PLAN.md, then verify the delivered snapshot and run the documented
fast S702 audits.

First implement the archived-state drag diagnostic described in priority 1.
Use the existing completed results for flaps 5/10, CL 0.8/1.2, both Reynolds
numbers and Ncrit 9/6.25. Resolve the correct selected flow states and compare
per-element wake drag, friction diagnostics, transition, Cp, Cf, theta,
displacement thickness and Hk. Reproduce the six recorded flap drag differences
and explain which modeled changes dominate the remaining discrepancies.

Preserve the original results and the frozen calibration. Distinguish the two
successful numerical retries from the original seven-of-nine reserved run.
Do not rerun the full sweep or alter the physics just to improve a fitted score.
After inspecting the diagnostic, choose the smallest discriminating experiment
for the next development step and record its rationale and evaluation budget.

Deliver complete runnable code, a concise findings report with plots, exact
commands and validation results. Keep the handover and development status current.
