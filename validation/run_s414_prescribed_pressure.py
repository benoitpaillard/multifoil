#!/usr/bin/env python3
"""Strict, bounded laminar BL replay; never computes forces or changes archives."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import hashlib
import json
import platform
import time
from pathlib import Path
import sys
import warnings
import importlib.metadata
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from viscous import new_model, mf
from amplification_jacobian import configure
from s414_geometry import elements
from audit_s414_lower_bucket import Inputs, load_states

OUT = ROOT/'validation/results/s414_prescribed_pressure_20261004'
PRIOR = ROOT/'validation/results/s414_fore_nose_20261003'
KINDS = ('native', 'model_taps', 'experiment_taps', 'experiment_delta')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def speed(cp, param):
    cp = np.asarray(cp)
    b = param.KTb if param.Minf > 0 else 1.
    k = .5*param.KTl*(1+b) if param.Minf > 0 else 0.
    raw_cp = b*cp/(1-k*cp)
    if np.any(1-raw_cp <= 0):
        raise ValueError('Pressure does not define positive outer speed')
    return param.Vinf*np.sqrt(1-raw_cp)


def station(left, guess, ue, xi, param):
    """Solve only theta, delta-star and N; fixed ue, strict residual/line search."""
    u = guess.copy(); u[3] = ue
    for iteration in range(60):
        with np.errstate(all='raise'):
            r, jac, _ = mf.residual_station(param, xi, np.column_stack((left, u)), np.zeros(2))
        norm = float(np.max(np.abs(r)))
        if norm < 1e-9:
            return u, norm, iteration
        scale = np.array([u[0], u[1], 10.])
        step = np.linalg.solve(jac[:, 4:7]*scale, -r)
        cost = np.linalg.norm(r/np.array([1., 1., 10.]))
        accepted = False
        for power in range(20):
            alpha = 2.**(-power)
            trial = u.copy()
            trial[:2] *= np.exp(np.clip(alpha*step[:2], -10, 10))
            trial[2] += alpha*10*step[2]
            try:
                with np.errstate(all='raise'):
                    if trial[2] < -1e-7 or mf.get_Hk(trial, param)[0] <= 1.000001:
                        continue
                    rr = mf.residual_station(param, xi, np.column_stack((left, trial)), np.zeros(2))[0]
                    newcost = np.linalg.norm(rr/np.array([1., 1., 10.]))
                if newcost < cost:
                    u = trial; accepted = True; break
            except (FloatingPointError, ValueError, OverflowError):
                continue
        if not accepted:
            raise ValueError(f'Newton line search failed; residual={norm:.9g}')
    raise ValueError(f'Newton iteration limit; residual={norm:.9g}')


def march(initial, seed, ue, xi, param):
    state = initial.copy(); state[:, 0] = seed
    maximum, iterations = 0., 0
    for j in range(1, len(xi)):
        try:
            state[:, j], error, nit = station(state[:, j-1], state[:, j], ue[j], xi[j-1:j+1], param)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ValueError(f'station {j}: {exc}') from exc
        maximum = max(maximum, error); iterations += nit
    return state, maximum, iterations


def profiles(x, native, tx, model, experiment):
    assert x[0] < tx[0] and x[-1] <= tx[-1] and np.all(np.diff(x) > 0)
    qx = np.r_[x[0], tx]
    return {
        'native': native.copy(),
        'model_taps': np.interp(x, qx, np.r_[native[0], model]),
        'experiment_taps': np.interp(x, qx, np.r_[native[0], experiment]),
        'experiment_delta': native+np.interp(x, qx, np.r_[0., experiment-model]),
    }


def run_profile(original, cp, xi, param, native):
    seed = original[:, 0].copy()
    initial = original.copy()
    initial[0, 1:] *= 1.0001; initial[1, 1:] *= .9999; initial[2, 1:] += .0001
    # Every target first independently restores the unchanged prescribed-pressure control.
    current, error, nit = march(initial, seed, original[3], xi, param)
    relative = float(np.max(abs(current[:2]/original[:2]-1)))
    nerror = float(np.max(abs(current[2]-original[2])))
    if relative >= 2e-5 or nerror >= 2e-5:
        raise RuntimeError('Native laminar restoration failed')
    control = dict(relative_thickness_error=relative, N_error=nerror, maximum_residual=error)
    history = []
    amount, step = (1., .05) if native else (0., .05)
    baseline_cp = mf.get_cp(original[3], param)[0]
    failure = None
    for attempt in range(300):
        if amount >= 1: break
        target = min(1., amount+step)
        try:
            candidate, error, nit = march(current, seed, speed((1-target)*baseline_cp+target*cp, param), xi, param)
            history.append(dict(fraction=target, accepted=True, residual=error, iterations=nit,
                                endpoint_N=float(candidate[2, -1])))
            current = candidate; amount = target; step = min(.05, step*1.5)
        except ValueError as exc:
            failure = str(exc)
            history.append(dict(fraction=target, accepted=False, error=failure))
            step *= .5
            if step < .05/4096: break
    complete = amount >= 1
    crossing = bool(np.any(current[2] >= param.ncrit))
    return current, dict(complete=complete, reached_fraction=amount, control=control,
        last_failure=failure, continuation=history, conditional_threshold_crossing=crossing,
        endpoint_N=float(current[2, -1]) if complete and not crossing else None,
        endpoint_Hk=float(mf.get_Hk(current[:, -1], param)[0]) if complete and not crossing else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(); out = args.output
    out.mkdir(parents=True, exist_ok=True)
    inputs = Inputs()
    _, states, _ = load_states(inputs)  # Retain all existing study/hash admission guards.
    prior = json.loads((PRIOR/'verification.json').read_text())
    pressure_file = PRIOR/'pressure_at_taps.csv'
    assert digest(pressure_file) == prior['output_sha256']['pressure_at_taps.csv']
    for path in (pressure_file, PRIOR/'verification.json', Path(__file__),
                 ROOT/'amplification_jacobian.py', ROOT/'validation/S414_PRESCRIBED_PRESSURE_PROTOCOL.md'):
        inputs.add(path)
    environment = dict(python=platform.python_version(),
        packages={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')},
        threads={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    manifest = dict(files=inputs.files, environment=environment)
    manifest_path = out/'input_manifest.json'
    if manifest_path.exists():
        assert json.loads(manifest_path.read_text()) == manifest, 'Changed inputs: choose a new output directory'
    else:
        write(manifest_path, manifest)
    rows = list(csv.DictReader(pressure_file.open()))
    model = new_model('fore', elements('nominal')[0][1], 1e6, .1, 9, (.25, 0), 1e-5)
    configure([model], 'corrected_af'); param = mf.build_param(model, 0)
    selected = {case: states[case+'_nominal']['arrays']['fore'] for case in ('cl0056', 'cl0175')}
    # Common physical node subset; exclude the mixed transition cell in either state.
    common_end = min(float(a['xy'][0, a['lower'][a['turbulent'][a['lower']] == 0][-1]])
                     for a in selected.values())
    summaries = []
    for case, a in selected.items():
        assert np.max(abs(model.foil.x-a['xy'])) < 1e-13
        taps = [r for r in rows if r['case'] == case and r['coordinate'] == 'table']
        tx = np.array([float(r['x']) for r in taps])
        tcp = np.array([float(r['viscous_Cp']) for r in taps])
        ecp = np.array([float(r['experiment_Cp']) for r in taps])
        branch = a['lower'][a['lower'] <= np.argmin(a['xy'][0])]
        xx = a['xy'][0, branch]
        first = np.where(xx < tx[0])[0][-1]
        ids = branch[first:][a['xy'][0, branch[first:]] <= common_end+1e-12]
        assert np.all(a['turbulent'][ids] == 0)
        x, xi, original = a['xy'][0, ids], a['xi'][ids], a['U'][:, ids]
        native = mf.get_cp(original[3], param)[0]
        pressure_profiles = profiles(x, native, tx, tcp, ecp)
        for kind in KINDS:
            path = out/(case+'_'+kind+'.json')
            if path.exists():
                saved = json.loads(path.read_text())
                assert digest(out/(case+'_'+kind+'.npz')) == saved['state_sha256']
                summaries.append(saved); continue
            started = time.monotonic()
            with warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                u, summary = run_profile(original, pressure_profiles[kind], xi, param, kind == 'native')
            summary.update(case=case, profile=kind, seed_x=float(x[0]), endpoint_x=float(x[-1]),
                number_of_stations=len(x), wall_seconds=time.monotonic()-started)
            statepath = out/(case+'_'+kind+'.npz')
            np.savez(statepath, U=u, original=original, x=x, xi=xi, ids=ids,
                requested_cp=pressure_profiles[kind], actual_cp=mf.get_cp(u[3], param)[0],
                taps_x=tx, taps_model=tcp, taps_experiment=ecp)
            summary['state_sha256'] = digest(statepath)
            write(path, summary); summaries.append(summary)
            print(f'{case} {kind}: complete={summary["complete"]}, '
                  f'fraction={summary["reached_fraction"]:.6f}, N={summary["endpoint_N"]}', flush=True)
            write(out/'status.json', dict(state='running', completed=len(summaries), total=8))
    write(out/'summary.json', dict(cases=summaries, completed_profiles=sum(s['complete'] for s in summaries),
        required_profiles=8, coupled_flow_solves=0, forces_computed=False))
    write(out/'status.json', dict(state='finished', completed=8, total=8))


if __name__ == '__main__':
    main()
