"""Explicit, bounded initialization paths around the unchanged fixed-lift solver.

This module changes initialization, not flow equations or admissibility. Each
path owns a fresh solver. The first successful path wins; drag never selects a
seed or a path. Returned solvers are complete in-memory copies at the selected
evaluation, not portable or persistent Newton checkpoints.
"""
from dataclasses import dataclass, field
import copy
import math

from optimization import admissible, solve_fixed_lift, LiftTargetError


@dataclass(frozen=True)
class InitializationPath:
    name: str
    start_alpha: float
    seed_alphas: tuple = ()


@dataclass
class FixedLiftOutcome:
    solution: object = None
    solver: object = None
    selected_flow_index: int | None = None
    selected_path: str | None = None
    evaluations: list = field(default_factory=list)
    attempts: list = field(default_factory=list)
    error: str = ''

    @property
    def feasible(self):
        return self.solution is not None


def lift_seed(samples, target_cl, alpha_bounds, slope_bounds=(.02, .3), max_extrapolation=1.):
    """Use two nearest admissible numerical lift values; never inspect drag error.

    Slope limits are declared per-degree subcritical assumptions. An unsuitable
    slope or seed outside the incidence domain is a failure, not clipping.
    """
    points = {float(s.alpha): float(s.cl) for s in samples
              if admissible(s) and math.isfinite(s.alpha)}
    nearest = sorted(points, key=lambda a: abs(points[a]-target_cl))[:2]
    if len(nearest) != 2:
        raise LiftTargetError('Local seed requires two distinct admissible incidences')
    a, b = sorted(nearest)
    slope = (points[b]-points[a])/(b-a)
    if not slope_bounds[0] <= slope <= slope_bounds[1]:
        raise LiftTargetError('Computed lift slope is unsuitable for a local retry')
    center = a+(target_cl-points[a])/slope
    if not (a-max_extrapolation <= center <= b+max_extrapolation and
            alpha_bounds[0] <= center <= alpha_bounds[1]):
        raise LiftTargetError('Computed lift seed exceeds the declared incidence domain')
    return center


def solve_with_initialization(solver_factory, target_cl, center, alpha_bounds, paths,
                              half_width=.05, cl_tolerance=2e-5,
                              residual_tolerance=1e-7, max_evaluations=20,
                              max_brackets=2, on_evaluation=None):
    """Try ordered fresh initialization paths with one optional lift-based recenter.

    ``solver_factory(start_alpha)`` must return a new solver with the same
    geometry/physics/settings for every call. All start/seed/requested incidences
    lie inside ``alpha_bounds``. Each bracket is center +/- half_width intersected
    with that domain. Recenter only after two admissible increasing-lift endpoints
    fail to bracket the target; an inadmissible solve advances to the next path.

    Callback signature: ``callback(solution, context)``. It receives every
    requested evaluation, including failed flows and explicit warm seeds. Internal
    continuation remains in each solution's diagnostics. Exceptions in caller
    callbacks propagate rather than being misclassified as numerical failure.
    """
    lo, hi = map(float, alpha_bounds)
    paths = tuple(paths)
    if (not all(math.isfinite(v) for v in (lo, hi, target_cl, center, half_width,
                                          cl_tolerance, residual_tolerance)) or
            not lo < hi or not lo <= center <= hi or
            min(half_width, cl_tolerance, residual_tolerance) <= 0 or
            not isinstance(max_evaluations, int) or max_evaluations < 3 or
            not isinstance(max_brackets, int) or not 1 <= max_brackets <= 2 or
            not paths or len({p.name for p in paths}) != len(paths)):
        raise ValueError('Invalid initialization policy, bounds, tolerances or budgets')
    for path in paths:
        if (not path.name or not all(math.isfinite(a) and lo <= a <= hi
                                    for a in (path.start_alpha, *path.seed_alphas))):
            raise ValueError('Initialization path exceeds the declared incidence domain')
    outcome = FixedLiftOutcome()
    # A private exception prevents a caller's RuntimeError from being swallowed
    # by the numerical failure handling below.
    class CallbackError(Exception):
        pass

    for path in paths:
        solver = solver_factory(path.start_alpha)
        local_center = float(center)
        for bracket in range(max_brackets):
            bounds = (max(lo, local_center-half_width), min(hi, local_center+half_width))
            attempt = dict(path=path.name, start_alpha=path.start_alpha,
                           seed_alphas=list(path.seed_alphas) if bracket == 0 else [],
                           bracket=bracket, center=local_center, bounds=list(bounds),
                           first_record=len(outcome.evaluations), feasible=False)
            outcome.attempts.append(attempt)
            selected = {}
            phase = 'warm_seed'

            def checkpoint(sol):
                index = len(outcome.evaluations)
                context = dict(index=index, path=path.name, bracket=bracket, phase=phase)
                outcome.evaluations.append(dict(solution=sol, **context))
                if on_evaluation is not None:
                    try:
                        on_evaluation(sol, dict(context))
                    except Exception as exc:
                        raise CallbackError() from exc
                if (not admissible(sol) or not math.isfinite(sol.residual) or
                        sol.residual > residual_tolerance):
                    raise LiftTargetError('Flow fails admissibility or residual tolerance')
                if phase == 'fixed_lift' and abs(sol.cl-target_cl) <= cl_tolerance:
                    selected[id(sol)] = (index, copy.deepcopy(solver))

            try:
                if bracket == 0:
                    for alpha in path.seed_alphas:
                        sol = solver.solve(alpha)
                        checkpoint(sol)
                        if not math.isfinite(sol.alpha) or abs(sol.alpha-alpha) > 1e-8:
                            raise LiftTargetError('Warm seed returned a different incidence')
                phase = 'fixed_lift'
                sol, samples = solve_fixed_lift(solver, target_cl, bounds, cl_tolerance,
                                               max_evaluations, checkpoint)
                if not math.isfinite(sol.alpha) or id(sol) not in selected:
                    raise LiftTargetError('Selected state is unavailable or invalid')
                index, restored = selected[id(sol)]
                attempt.update(feasible=True, last_record=len(outcome.evaluations)-1)
                outcome.solution, outcome.solver = copy.deepcopy(sol), restored
                outcome.selected_flow_index, outcome.selected_path = index, path.name
                outcome.error = ''
                return outcome
            except CallbackError as exc:
                raise exc.__cause__
            except (LiftTargetError, ValueError, FloatingPointError, AssertionError, RuntimeError) as exc:
                outcome.error = f'{type(exc).__name__}: {exc}'
                attempt.update(error=outcome.error, last_record=len(outcome.evaluations)-1)
                samples = [e['solution'] for e in outcome.evaluations[attempt['first_record']:]
                           if e['phase'] == 'fixed_lift']
                # Restrict the only adaptive step to a genuine unbracketed,
                # increasing pair. Do not infer a usable branch across a failure.
                if (bracket+1 < max_brackets and len(samples) == 2 and
                        all(admissible(s) and math.isfinite(s.residual) and
                            s.residual <= residual_tolerance for s in samples) and
                        samples[0].alpha < samples[1].alpha and samples[0].cl < samples[1].cl and
                        not samples[0].cl-cl_tolerance <= target_cl <= samples[1].cl+cl_tolerance):
                    try:
                        local_center = lift_seed(samples, target_cl, (lo, hi))
                        continue
                    except LiftTargetError as seed_error:
                        attempt['recenter_error'] = str(seed_error)
                break
    return outcome
