#!/usr/bin/env python3
"""Post-execution diagnostic for V5 tranche 1: is the prescribed-pressure station problem unique?

Added after the frozen laminar-replay run, to explain why its closure marches
and the archive disagree downstream of separation. It changes no outcome of
that run. Same equations, same native nodes, same archived Mach number and edge
velocity as the archived states; only the Newton starting point differs.
No coupled solve, no forces.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import json
import warnings
from pathlib import Path
import numpy as np
import run_s414_laminar_replay as replay
from run_s414_laminar_replay import ROOT, RESULTS, CASES, STATES, X_END, MACH, mf, station, digest, write, write_csv

OUT = RESULTS/'s414_direct_branch_diagnostic_20261004'
STUDY = RESULTS/'s414_laminar_replay_20261004'
ENDPOINT = .02082853867960402  # last station of the earlier replay


def solve(left, guess, speed, xi, param):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            state, error, _ = station(left, guess.copy(), speed, xi, param)
        return state, error
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, RuntimeWarning):
        return None, None


def main():
    if (OUT/'summary.json').exists():
        raise SystemExit('A completed diagnostic exists here')
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((STUDY/'input_manifest.json').read_text())['files']
    receipt = json.loads(replay.STATE_RECEIPT.read_text())
    for directory, _ in CASES.values():
        assert digest(STATES/directory/'fore_state.npz') == receipt[directory+'/fore_state.npz']
    inputs = {name: digest(ROOT/name) for name in ('validation/diagnose_s414_direct_branches.py',
        'validation/run_s414_laminar_replay.py', 'validation/run_s414_prescribed_pressure.py', 'vendor/mfoil.py',
        'viscous.py', 'amplification_jacobian.py')}
    for name in ('validation/run_s414_laminar_replay.py', 'vendor/mfoil.py', 'viscous.py'):
        assert inputs[name] == manifest[name], 'Source differs from the frozen run: '+name
    for directory, _ in CASES.values():
        name = str((STATES/directory/'fore_state.npz').relative_to(ROOT)); inputs[name] = manifest[name]
    write(OUT/'input_manifest.json', dict(files=inputs))

    _, param = replay.parameters(MACH)
    pressure = replay.Pressure(param)
    stations, roots, summary = [], [], {}
    for case, (directory, lift) in CASES.items():
        s = replay.Surface(np.load(STATES/directory/'fore_state.npz'), pressure)
        end = int(np.sum(s.xi < s.xi_of_x(X_END)))
        hk = lambda u: float(mf.get_Hk(u, param)[0])
        # 1. Direct march from the stagnation start, each station started from its upstream neighbour.
        march = replay.closure_march(s.xi[:end], s.U[3, :end], param)
        n = len(march['xi'])
        archive_residual = 0.
        for j in range(n):
            laminar = not s.turbulent[j]
            if j and laminar:
                r = mf.residual_station(param, s.xi[j-1:j+1], s.U[:, j-1:j+1], np.zeros(2))[0]
                archive_residual = max(archive_residual, float(np.max(abs(r))))
            stations.append(dict(case=case, station=j, x=float(s.x[j]), archive_laminar=bool(laminar),
                archive_theta=float(s.U[0, j]), archive_Hk=hk(s.U[:, j]), archive_N=float(s.U[2, j]) if laminar else None,
                march_theta=float(march['states'][0, j]), march_Hk=float(march['Hk'][j]), march_N=float(march['states'][2, j]),
                march_cf_Re_theta=float(march['friction'][j])))
        # 2. Identical upstream state (the archive) and edge velocity; two Newton starting points.
        doubles = []
        for j in range(1, end):
            if s.turbulent[j]:
                break
            left, right = s.U[:, j-1], s.U[:, j]
            a, a_error = solve(left, right, right[3], s.xi[j-1:j+1], param)
            b, b_error = solve(left, left, right[3], s.xi[j-1:j+1], param)
            two = a is not None and b is not None and abs(hk(a)-hk(b)) > 1e-6
            if two:
                doubles.append(j)
            roots.append(dict(case=case, station=j, x=float(s.x[j]), archive_Hk=hk(right),
                root_from_archive_Hk=None if a is None else hk(a), residual_from_archive=a_error,
                root_from_upstream_Hk=None if b is None else hk(b), residual_from_upstream=b_error,
                root_from_upstream_delta_star_ratio=None if b is None else float(b[1]/right[1]),
                two_distinct_roots=bool(two)))
        # Agreement with the archive from the earlier replay's seed node up to the first two-root station.
        seed = s.le+int(np.where(s.x[s.le:] >= .00165)[0][0])
        split = min(doubles[0] if doubles else n, n)
        k = int(np.argmin(abs(s.x[:n]-ENDPOINT))) if n else None
        reached = k is not None and abs(s.x[k]-ENDPOINT) < 1e-9
        summary[case] = dict(target_cl=lift, march_status=march['status'], march_message=march['message'],
            march_max_residual=march['residual'], archive_laminar_max_residual=archive_residual,
            compared_from_x=float(s.x[seed]), compared_stations=int(split-seed),
            upstream_max_Hk_difference=float(max(abs(march['Hk'][j]-hk(s.U[:, j])) for j in range(seed, split))),
            upstream_max_theta_relative=float(np.max(abs(march['states'][0, seed:split]/s.U[0, seed:split]-1))),
            stagnation_node_Hk_difference=float(march['Hk'][0]-hk(s.U[:, 0])),
            stations_with_two_roots=[dict(station=j, x=float(s.x[j])) for j in doubles],
            archive_max_Hk=max(hk(s.U[:, j]) for j in range(n) if not s.turbulent[j]),
            march_max_Hk=float(np.max(march['Hk'])), march_min_cf_Re_theta=float(np.min(march['friction'])),
            endpoint_x=ENDPOINT if reached else None,
            archive_N_at_endpoint=float(s.U[2, k]) if reached and not s.turbulent[k] else None,
            march_N_at_endpoint=float(march['states'][2, k]) if reached else None)
        print(case, json.dumps(summary[case]))
    write_csv(OUT/'stations.csv', stations); write_csv(OUT/'station_roots.csv', roots)
    write(OUT/'summary.json', dict(cases=summary, coupled_flow_solves=0, forces_computed=False,
        output_sha256={name: digest(OUT/name) for name in ('stations.csv', 'station_roots.csv', 'input_manifest.json')}))


if __name__ == '__main__':
    main()
