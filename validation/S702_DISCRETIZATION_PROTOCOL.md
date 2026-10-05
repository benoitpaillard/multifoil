# S702 separation of contour, surface and wake controls

Declared 1 October 2026, before new study solves. The previous three-grid
screen failed; another coupled grid increase cannot identify the cause.
The new literature identifies Pomeroy 2014/2015 as a promising independent
wake-mechanism benchmark with coordinate tables, but does not establish that
wake bursting causes the S702 discrepancy. Qualify those measured observables
after this numerical diagnostic; keep S103/S414 queued for geometry audits.

## Hypothesis and scope

The prior 340/240 to 440/320 step changed the input spline representation,
surface station count, automatically derived wake count and first wake step
together. Isolate their conditional contributions along a declared path.
Keep Re 700000, Mach 0.07, Ncrit 6.25, CL 0.8, flaps 5/10, reference chord 1,
wake extension 1c and the original fixed-lift initialization policy unchanged.
No new closure, transition fit or additional surface resolution is proposed.

The opt-in `DiscretizedSolver` in `diagnostic_discretization.py` subclasses the
unchanged production solver. It samples the same regularized source spline at
an explicit surface count and uses explicit wake counts/requested first steps.
All boundary-layer equations, force estimates and wake tracing remain unchanged.
The copied wake-tracing path is verified against the production path, including
an exact converged single-element state comparison and an S702 repeated control.

## Eleven-job matrix

Run both flap settings for each of the five rows below, then repeat flap 10
using the original production solver with automatic controls (eleven jobs).
Counts in the first two columns are panels; the API records nodes = panels + 1.

| Stage | Source resampling | Surface panels | Wake node recipe | Requested first-step recipe |
|---|---|---|---|---|
| A: fine | 340/240 | 340/240 | fine | fine |
| B: contour | 440/320 | 340/240 | fine | fine |
| C: surface | 440/320 | 440/320 | fine | fine |
| D: wake count | 440/320 | 440/320 | finer | fine |
| E: wake spacing | 440/320 | 440/320 | finer | finer |

Recipes are evaluated from the corresponding original 340/240 or 440/320
geometry using the unchanged production `new_model`, count formula and mean
terminal surface arclength step. Compute them once per job and save the
numerical values, source coordinates and source hashes. They remain fixed
through incidence continuation. Extent follows the same original rule; the
traced wake location still responds to geometry and fixed-lift incidence.
Fixing the requested first tracing step does not fix its curved chord length.

Stage B varies the spline source at fixed surface count; it cannot hold the
actual surface node coordinates fixed. Stage C holds that source fixed while
changing its surface sampling. This distinction is part of the diagnostic.
The path measures conditional contributions, not unique independent physical
causes: discretization factors may interact. The four increments telescope
to the original fine-to-finer change if both endpoints reproduce.

Use original model-derived incidence centers, bounds [-6,8] degrees, 0.05-degree
initial half-width, two brackets, at most twenty target evaluations per bracket,
100 Newton iterations, CL tolerance 2e-5 and residual tolerance 1e-7.
Paths remain cold-local, fresh -1-degree anchor, then fresh -1-degree anchor
with explicit center-0.4-degree entry; first success wins. Preserve every failed
request/internal continuation point. Do not select a path based on drag error.

## Gates and interpretation

Before any complete attribution require 11/11 feasible jobs, selected-state
audit, unchanged frozen parent source/input hashes and matching runtime.
Both A states must reproduce archived fine values and both E states archived
finer values within 0.01 CD count. The fresh original-solver flap-10 control
must reproduce A within that tolerance; report exact saved-array agreement too.
Missing or failed cases cannot be omitted from a passing result.

Report each individual increment and each change of flap-10-minus-flap-5 drag,
using the same stage's own baseline. Include fore/aft contributions, CL/CM/
incidence, transition/separation, actual contours and source arrays, surface
Hk/velocity steps, requested/realized wake spacing, clearance, initialization
and cost. Verify telescoping numerically. Retain experimental deltas as context.
A large contribution identifies a numerical sensitivity to investigate; it
does not validate the model or justify a drag correction. The existing 0.5-count
screen remains unchanged; this experiment alone supplies no error bound.

No automatic further grid or optimization campaign follows. If gates fail,
audit the failed target or compatibility discrepancy. If a factor dominates,
design a follow-up around that factor before changing physics. Otherwise use
the Pomeroy measurements to formulate a separately qualified mechanism test.

## Execution and monitoring

At most four workers, one BLAS/OpenMP thread per worker, detached and logged on
jony4. Eleven targets; no unbounded retries. The frozen runtime is Python 3.12
with NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Hash all code, parent
references and settings. Interrupted attempts are archived before restart;
completed results are reused only within the same manifest and runtime.

Remote root: `/home/jony/multifoil-runs/s702-discretization-20261001`.

```bash
# Remote run or resume, from that isolated root:
.venv/bin/python validation/run_s702_discretization.py --workers 4

# Local monitoring; closing the monitor leaves the remote job running:
python3 validation/run_s702_discretization.py --remote-status --watch
```

The runner reports coverage, endpoint/control checks and provisional increments.
Full selected-state audit and scientific interpretation follow completion.
