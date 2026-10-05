"""Calibration metrics with explicit coverage and grouped holdout selection.

No estimated experimental standard deviations are invented. MAPE is the stated
engineering loss; absolute drag errors and configuration differences are also
reported. Incomplete candidates cannot win by omitting difficult points.
"""
import itertools
import math
import numpy as np


def score(rows, required):
    """Score one candidate on exactly the required (flap, CL) keys."""
    lookup = {(r['flap'], r['target_cl']): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError('Duplicate operating points')
    subset = [lookup.get(key) for key in required]
    good = [r for r in subset if r and r['feasible'] and r.get('CD') is not None
            and math.isfinite(r['CD']) and r['CD'] > 0]
    complete = len(good) == len(required)
    lower_bound = float(sum(abs(r['CD']/r['CD_experiment']-1) for r in good)*100/len(required))
    if not complete:
        return dict(complete=False, evaluated=len(good), required=len(required),
                    mape_percent=None, rmse_counts=None, maximum_error_percent=None,
                    minimum_possible_mape_percent=lower_bound)
    error = np.array([r['CD']-r['CD_experiment'] for r in good])
    relative = error/np.array([r['CD_experiment'] for r in good])
    return dict(complete=True, evaluated=len(good), required=len(required),
        mape_percent=float(np.mean(abs(relative))*100),
        rmse_counts=float(np.sqrt(np.mean(error**2))*1e4),
        maximum_error_percent=float(np.max(abs(relative))*100),
        minimum_possible_mape_percent=lower_bound)


def select(rows, required, candidates=None):
    if candidates is None:
        candidates = sorted({r['ncrit'] for r in rows})
    scored = []
    for n in candidates:
        result = score([r for r in rows if r['ncrit'] == n], required)
        if result['complete']:
            scored.append((result['mape_percent'], abs(n-9.), n, result))
    if not scored:
        raise ValueError('No candidate has complete calibration coverage')
    _, _, n, result = min(scored)
    return dict(ncrit=n, **result)


def leave_one_flap_out(rows, flaps, targets, candidates):
    """Select using only training flap errors; then expose every held-out point.

    Use the predeclared coarse grid. Refinement chosen from all flap errors is
    deliberately excluded from this check to avoid adaptive selection leakage.
    """
    folds = []
    heldout = []
    for flap in flaps:
        train = [(f, c) for f in flaps if f != flap for c in targets]
        test = [(flap, c) for c in targets]
        winner = select(rows, train, candidates)
        chosen = [r for r in rows if r['ncrit'] == winner['ncrit']]
        test_score = score(chosen, test)
        baseline = score([r for r in rows if r['ncrit'] == 9.], test)
        folds.append(dict(heldout_flap=flap, selected_ncrit=winner['ncrit'],
            training_mape_percent=winner['mape_percent'],
            test_mape_percent=test_score['mape_percent'], test_complete=test_score['complete'],
            baseline_n9_mape_percent=baseline['mape_percent']))
        heldout.extend(r for r in chosen if r['flap'] == flap)
    overall = score(heldout, [(f, c) for f in flaps for c in targets])
    return folds, overall


def pairwise(rows, flaps, targets):
    """Drag changes are reported separately; no hidden ranking penalty in loss."""
    lookup = {(r['flap'], r['target_cl']): r for r in rows}
    result = []
    for cl in targets:
        for first, second in itertools.combinations(flaps, 2):
            a, b = lookup.get((first, cl)), lookup.get((second, cl))
            if a is None or b is None or not a['feasible'] or not b['feasible']:
                continue
            predicted = (b['CD']-a['CD'])*1e4
            measured = (b['CD_experiment']-a['CD_experiment'])*1e4
            result.append(dict(first_flap=first, second_flap=second, target_cl=cl,
                predicted_delta_counts=predicted, measured_delta_counts=measured,
                delta_error_counts=predicted-measured,
                same_sign=bool(np.sign(predicted) == np.sign(measured))))
    return result
