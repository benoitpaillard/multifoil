"""Small, auditable helpers for fixed-lift design screening.

No shape adjoint or validated optimizer is claimed. Each geometry needs its own
solver; reuse only the flow state while finding incidence for that geometry.
"""
from __future__ import annotations

import math
import numpy as np


class LiftTargetError(RuntimeError):
    """No admissible, converged solution at the required lift."""


def admissible(solution):
    return (solution.converged and
            all(v is not None and math.isfinite(v) for v in (solution.cl, solution.cd, solution.cm)) and
            solution.cd > 0 and
            not any(solution.diagnostics.get('local_supersonic_detected', {}).values()))


def solve_fixed_lift(solver, target_cl, bounds=(4., 7.), cl_tolerance=2e-5,
                     max_evaluations=16, on_evaluation=None):
    """Safeguarded secant solve on a prescribed increasing-lift branch.

    Bounds are a design-domain restriction, not a stall detector. Both endpoints
    must be admissible and bracket the target. No extrapolation or use of failed
    coefficients. The callback receives every requested flow evaluation.
    """
    lo, hi = map(float, bounds)
    if not (np.isfinite([lo, hi, target_cl, cl_tolerance]).all() and lo < hi and
            cl_tolerance > 0 and isinstance(max_evaluations, int) and max_evaluations >= 3):
        raise ValueError('Invalid fixed-lift bounds, target, tolerance or budget')
    trace = []

    def evaluate(alpha):
        result = solver.solve(alpha)
        trace.append(result)
        if on_evaluation is not None:
            on_evaluation(result)
        if abs(result.alpha-alpha) > 1e-8 or not admissible(result):
            raise LiftTargetError(f'Inadmissible flow evaluation at alpha={alpha:g}')
        return result

    left, right = evaluate(lo), evaluate(hi)
    if right.cl <= left.cl:
        raise LiftTargetError('Bounds do not identify an increasing-lift branch')
    for result in (left, right):
        if abs(result.cl-target_cl) <= cl_tolerance:
            return result, trace
    if not left.cl < target_cl < right.cl:
        raise LiftTargetError('Required lift is not bracketed inside the allowed incidence range')
    for _ in range(max_evaluations-2):
        alpha = left.alpha + (target_cl-left.cl)*(right.alpha-left.alpha)/(right.cl-left.cl)
        if not left.alpha+1e-9 < alpha < right.alpha-1e-9:
            alpha = .5*(left.alpha+right.alpha)
        result = evaluate(alpha)
        if not left.cl-cl_tolerance <= result.cl <= right.cl+cl_tolerance:
            raise LiftTargetError('Nonmonotone lift response inside the selected branch')
        if abs(result.cl-target_cl) <= cl_tolerance:
            return result, trace
        if result.cl < target_cl:
            left = result
        else:
            right = result
        # Prevent endpoint locking on curved responses.
        if len(trace) % 4 == 0 and len(trace) < max_evaluations:
            result = evaluate(.5*(left.alpha+right.alpha))
            if not left.cl-cl_tolerance <= result.cl <= right.cl+cl_tolerance:
                raise LiftTargetError('Nonmonotone lift response inside the selected branch')
            if abs(result.cl-target_cl) <= cl_tolerance:
                return result, trace
            if result.cl < target_cl:
                left = result
            else:
                right = result
        if len(trace) >= max_evaluations:
            break
    raise LiftTargetError('Fixed-lift evaluation budget exhausted')


def camber_bump(xy, amplitude):
    """Shift both surfaces by 16*a*c*x^2*(1-x)^2 along the chord normal.

    Amplitude is a fraction of local chord. Endpoints and end slopes are fixed;
    thickness at fixed chordwise position is preserved in the input geometry.
    """
    xy = np.array(xy, dtype=float, copy=True)
    if not np.isfinite(amplitude):
        raise ValueError('Camber amplitude must be finite')
    te = .5*(xy[0]+xy[-1])
    le = xy[np.argmax(np.linalg.norm(xy-te, axis=1))].copy()
    chord = np.linalg.norm(te-le)
    axis = (te-le)/chord
    normal = np.array([-axis[1], axis[0]])
    x = np.clip((xy-le)@axis/chord, 0, 1)
    return xy + (16*float(amplitude)*chord*x*x*(1-x)**2)[:, None]*normal


def minimum_clearance(first, second):
    """Exact segment-distance minimum for two nonintersecting closed contours.

    Call validate_geometry first: intersecting segments invalidate this formula.
    """
    def distance(points, contour):
        start, end = contour, np.roll(contour, -1, axis=0)
        delta = end-start
        squared = np.sum(delta*delta, axis=1)
        valid = squared > 1e-24
        start, delta, squared = start[valid], delta[valid], squared[valid]
        offset = points[:, None, :]-start[None, :, :]
        t = np.clip(np.sum(offset*delta, axis=2)/squared, 0, 1)
        return np.min(np.linalg.norm(offset-t[:, :, None]*delta, axis=2))
    return float(min(distance(first, second), distance(second, first)))
