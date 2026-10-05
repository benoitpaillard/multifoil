#!/usr/bin/env python3
"""V5 tranche 1: attached-laminar replay controls. No coupled solve, no forces.

This is the amended runner (see the protocol's amendment). The first execution
and the sources that produced it are kept in initial_run_box_scheme/.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import re
import sys
import time
import warnings
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'validation'))
from viscous import new_model, mf
from amplification_jacobian import configure
from s414_geometry import elements
from run_s414_prescribed_pressure import station  # the earlier station solver, unchanged
import s414_laminar_fd as fd

RESULTS = ROOT/'validation/results'
OUT = RESULTS/'s414_laminar_replay_20261004'
STATES = RESULTS/'s414_lower_bucket_20261003'
STATE_RECEIPT = RESULTS/'s414_lower_bucket_launch_20261003/retrieved_files_sha256.json'
TAPS = RESULTS/'s414_fore_nose_20261003/pressure_at_taps.csv'
TAP_RECEIPT = RESULTS/'s414_fore_nose_20261003/verification.json'
EARLIER = RESULTS/'s414_prescribed_pressure_20261004'
PROTOCOL = ROOT/'validation/S414_LAMINAR_REPLAY_PROTOCOL.md'
CASES = dict(cl0056=('cl0056_nominal', .056), cl0175=('cl0175_nominal', .175), cl028=('cl028_control', .284))
CROSS = (('cl0056', 'cl0175'), ('cl0175', 'cl0056'), ('cl0175', 'cl028'), ('cl028', 'cl0175'))
FAMILIES = ('linear', 'monotone', 'additive', 'ratio')
UPSTREAM = ('hold', 'ramp', 'late')
X_END, TAP_COUNT, SUB, NU, REYNOLDS, MACH = .03424, 5, 8, 1e-6, 1e6, .1
TOL_X, TOL_THETA = .0012, .05


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def parameters(mach):
    model = new_model('fore', elements('nominal')[0][1], REYNOLDS, mach, 9, (.25, 0), 1e-5)
    configure([model], 'corrected_af')
    return model, mf.build_param(model, 0)


class Pressure:
    """The model's own pressure/speed relations at the archived Mach number."""
    def __init__(self, param):
        self.l, self.b = param.KTl, param.KTb

    def physical(self, u):
        return u*(1-self.l)/(1-self.l*u**2)

    def incompressible(self, q):
        q = np.asarray(q, dtype=float)
        return (-(1-self.l)+np.sqrt((1-self.l)**2+4*self.l*q**2))/(2*self.l*q)

    def cp(self, q):
        raw = 1-self.incompressible(q)**2
        return raw/(self.b+.5*self.l*(1+self.b)*raw)

    def speed(self, cp):
        cp = np.asarray(cp, dtype=float)
        raw = self.b*cp/(1-.5*self.l*(1+self.b)*cp)
        if np.any(1-raw <= 0):
            raise ValueError('Pressure does not define a positive speed')
        return self.physical(np.sqrt(1-raw))


class Surface:
    """Fore-lower surface of one archived state, ordered from the stagnation node."""
    def __init__(self, arrays, pressure, interpolant='pchip'):
        low = arrays['lower']
        self.node, self.x, self.xi = low, arrays['xy'][0, low], arrays['xi'][low]
        self.U, self.turbulent = arrays['U'][:, low], arrays['turbulent'][low]
        self.le = int(np.argmin(self.x))
        self.q_nodes = pressure.physical(self.U[3])
        count = int(np.where(self.x[self.le:] > .06)[0][0])+self.le+1
        build = PchipInterpolator if interpolant == 'pchip' else CubicSpline
        self.q = build(self.xi[:count], self.q_nodes[:count])
        self.pressure = pressure

    def xi_of_x(self, x):
        return np.interp(x, self.x[self.le:], self.xi[self.le:])

    def x_of_xi(self, xi):
        return float(np.interp(xi, self.xi, self.x))

    def cp(self, xi):
        return self.pressure.cp(self.q(xi))

    def nodes(self, extra=(), sub=SUB):
        end = float(self.xi_of_x(X_END))
        base = np.r_[self.xi[self.xi < end-1e-12], end]
        fine = np.concatenate([np.linspace(a, b, sub+1)[:-1] for a, b in zip(base[:-1], base[1:])]+[[end]])
        extra = np.asarray(extra, dtype=float)
        merged = np.unique(np.r_[fine, extra[(extra > base[0]) & (extra < end)]])
        return merged[np.r_[True, np.diff(merged) > 1e-10]]


def reconstruction(base, tap_xi, target_cp, family, upstream):
    """Speed along the base surface that reproduces the target tap pressures."""
    pressure = base.pressure
    target_q = pressure.speed(target_cp)
    base_q, base_cp = base.q(tap_xi), base.cp(tap_xi)
    first, ratio = tap_xi[0], target_q[0]/base_q[0]
    seed = base.xi[base.xi < first][-1]
    monotone = PchipInterpolator(tap_xi, target_cp)

    def speed(xi):
        xi = np.asarray(xi, dtype=float)
        if upstream == 'hold':
            factor = np.full(xi.shape, ratio)
        elif upstream == 'ramp':
            factor = 1+(ratio-1)*xi/first
        else:
            factor = 1+(ratio-1)*np.clip((xi-seed)/(first-seed), 0, 1)
        inside = np.clip(xi, first, tap_xi[-1])
        if family == 'linear':
            down = pressure.speed(np.interp(inside, tap_xi, target_cp))
        elif family == 'monotone':
            down = pressure.speed(monotone(inside))
        elif family == 'additive':
            down = pressure.speed(base.cp(inside)+np.interp(inside, tap_xi, target_cp-base_cp))
        else:
            down = base.q(inside)*np.interp(inside, tap_xi, target_q/base_q)
        return np.where(xi < first, base.q(np.minimum(xi, first))*factor, down)
    return speed


def stagnation(speed, xi, param):
    """Similarity start of the inherited initializer."""
    nu = param.mu0/param.rho0
    K = speed/xi
    theta, delta = mf.thwaites_init(K, nu)
    state = np.array([theta, delta, 0., K*1e-6])
    for _ in range(30):
        param.simi = True
        r, jac, _ = mf.residual_station(param, np.r_[1e-6, 1e-6], np.stack((state, state), axis=-1), np.zeros(2))
        param.simi = False
        if np.linalg.norm(r) < 1e-10:
            break
        step = np.append(np.linalg.solve(jac[:, 4:7]+jac[:, 0:3], -r), 0)
        size = max(abs(step[0]/state[0]), abs(step[1]/state[1]))
        state = state+step*(1 if size < .2 else .2/size)
    else:
        raise ValueError('Stagnation state did not converge')
    state[3] = speed
    return state


def friction(state, param):
    return float(mf.get_cf(state, param)[0]*mf.get_Ret(state, param)[0])


def closure_march(xi, speed, param, start=None, function=None, depth=12):
    """Direct march of the inherited integral equations; stops at cf <= 0 or at a failure.

    With a speed function a failed station solve is halved (depth levels)
    before the failure is recorded; without one every node gets one attempt.
    """
    position = [float(xi[0])]
    states = [stagnation(speed[0], xi[0], param) if start is None else np.array(start, dtype=float)]
    values = [friction(states[0], param)]
    out = dict(status='completed', message='', xi_separation=None, residual=0.)
    for j in range(1, len(xi)):
        target = float(xi[j]); smallest, failed = (target-position[-1])*2.**(-depth), None
        while position[-1] < target:
            final = failed is not None and failed-position[-1] < smallest
            trial = target if failed is None else (failed if final else .5*(position[-1]+failed))
            value = float(speed[j]) if trial == target else float(function(trial))
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('error', RuntimeWarning)
                    state, error, _ = station(states[-1], states[-1].copy(), value, np.array([position[-1], trial]), param)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError, RuntimeWarning) as exc:
                if function is None or final:
                    out.update(status='failed', message=f'station {j}: {exc}'); break
                failed = trial
                continue
            out['residual'] = max(out['residual'], error)
            position.append(float(trial)); states.append(state); values.append(friction(state, param))
            if values[-1] <= 0:
                weight = values[-2]/(values[-2]-values[-1])
                out.update(status='separated', xi_separation=float(position[-2]+weight*(position[-1]-position[-2])),
                           theta_reference=float(states[-2][0]+weight*(states[-1][0]-states[-2][0])))
                break
            if final:
                failed = None
        if out['status'] != 'completed':
            break
    states = np.array(states).T
    out.update(xi=np.array(position), states=states, Hk=np.array([mf.get_Hk(s, param)[0] for s in states.T]),
               friction=np.array(values))
    if out['status'] != 'separated':
        out['theta_reference'] = float(states[0, -1])
    return out


def response(surface, nodes, speed_function, param, with_thwaites=False):
    """FD, closure (and optionally Thwaites) responses to one edge-velocity distribution."""
    speed = speed_function(nodes)
    started = time.monotonic()
    finite = fd.march(nodes, lambda s: float(speed_function(s)), NU, m_start=1.)
    fd_seconds = time.monotonic()-started
    closure = closure_march(nodes, speed, param, function=speed_function)
    result = dict(
        fd=dict(status=finite['status'], reason=finite['stop_reason'],
                x_separation=None if finite['xi_separation'] is None else surface.x_of_xi(finite['xi_separation']),
                xi_separation=finite['xi_separation'], theta_reference=float(finite['theta'][-1]),
                H_last=float(finite['H'][-1]), stations=int(len(finite['xi'])), edge_resolved=finite['edge_resolved'],
                oscillation=fd.oscillation(finite['wall_shear']), wall_seconds=fd_seconds),
        closure=dict(status=closure['status'], reason=closure['message'],
                     x_separation=None if closure['xi_separation'] is None else surface.x_of_xi(closure['xi_separation']),
                     xi_separation=closure['xi_separation'], theta_reference=closure['theta_reference'],
                     Hk_last=float(closure['Hk'][-1]), stations=int(len(closure['xi'])), residual=closure['residual'],
                     N_last=float(closure['states'][2, -1])))
    arrays = dict(nodes=nodes, speed=speed, cp=surface.pressure.cp(speed), x=np.array([surface.x_of_xi(s) for s in nodes]),
                  fd_xi=finite['xi'], fd_x=np.array([surface.x_of_xi(s) for s in finite['xi']]), fd_theta=finite['theta'],
                  fd_delta_star=finite['delta_star'], fd_H=finite['H'], fd_cf=finite['cf'], fd_speed=finite['ue'],
                  closure_xi=closure['xi'], closure_x=np.array([surface.x_of_xi(s) for s in closure['xi']]),
                  closure_states=closure['states'], closure_Hk=closure['Hk'], closure_friction=closure['friction'])
    if with_thwaites:
        simple = fd.thwaites(nodes, speed, np.gradient(speed, nodes), NU, k_stagnation=True)
        result['thwaites'] = dict(x_separation=None if simple['xi_separation'] is None
                                  else surface.x_of_xi(simple['xi_separation']))
        arrays.update(thwaites_theta=simple['theta'], thwaites_lambda=simple['lam'], thwaites_H=simple['H'])
    return result, arrays


def compare(result, truth):
    """Declared control screen: topology, separation position and reference theta."""
    same = result['status'] == truth['status'] and result['status'] in ('completed', 'separated')
    shift = None
    if same and result['status'] == 'separated':
        shift = result['x_separation']-truth['x_separation']
    theta = result['theta_reference']/truth['theta_reference']-1
    passed = bool(same and (shift is None or abs(shift) <= TOL_X) and abs(theta) <= TOL_THETA)
    return dict(same_topology=bool(same), separation_shift=shift, theta_relative=float(theta), passed=passed)


def qualification(param):
    """Part B. Returns the control table and the gate decision."""
    eta = fd.eta_grid(); rows = []

    def add(name, value, reference, tolerance, gate=True):
        rows.append(dict(control=name, value=float(value), reference=reference, tolerance=tolerance,
                         passed=bool(abs(value-reference) <= tolerance) if gate else None, gate=gate))
    d1, d2, shear = fd.measures(eta, *fd.similarity(0., eta))
    add('Blasius wall shear', shear, .332057, 1e-4); add('Blasius shape factor', d1/d2, 2.5911, 1e-3)
    add('Hiemenz wall shear', fd.measures(eta, *fd.similarity(1., eta))[2], 1.232588, 2e-4)
    for m in (1., 0., 1/3, -.05):
        r = fd.march(np.geomspace(1e-3, 1., 41), lambda x, m=m: x**m, NU, m_start=m)
        add(f'similarity preserved, m={m:.4g}', np.max(abs(r['wall_shear']-r['wall_shear'][0])), 0., 1e-8)
    howarth = lambda count: np.r_[np.geomspace(1e-4, 2e-3, 8)[:-1], np.linspace(2e-3, .13, count)]
    for count in (26, 51, 101, 201):
        r = fd.march(howarth(count), lambda x: 1-x, NU, m_start=0.)
        add(f'Howarth separation, {count} nodes', r['xi_separation'] if r['status'] == 'separated' else np.nan,
            .1198, 5e-4, gate=count >= 101)
    for count in (46, 91, 181):
        r = fd.march(np.linspace(.02, 2., count), lambda x: 2*np.sin(x), NU, m_start=1.)
        add(f'cylinder separation [deg], {count} nodes',
            np.degrees(r['xi_separation']) if r['status'] == 'separated' else np.nan, 104.5, .3, gate=count >= 91)
    # Amendment: a kink in the edge velocity must not produce a sawtooth in the wall shear.
    kinked = lambda x: 1. if x <= .5 else 1-.1*(x-.5)
    for count in (101, 201):
        r = fd.march(np.r_[np.linspace(1e-3, .5, count)[:-1], np.linspace(.5, 1., count)], kinked, NU, m_start=0.)
        add(f'kinked edge velocity, sawtooth index, {count} nodes per segment',
            fd.oscillation(r['wall_shear'][r['xi'] > .5]) if r['status'] == 'completed' else np.nan, 0., 1e-3)
    # The inherited closure and Thwaites on the same adverse-gradient controls: reported, never a gate.
    for name, nodes, speed, slope, start, reference, scale in (
            ('Howarth', howarth(201), lambda x: 1-x, lambda x: -np.ones_like(x), 'flat', .1198, 1.),
            ('cylinder [deg]', np.linspace(.02, 2., 181), lambda x: 2*np.sin(x), lambda x: 2*np.cos(x), 'stagnation', 104.5,
             180/np.pi)):
        if start == 'flat':
            # Blasius start: theta = 0.664 sqrt(nu x / ue), H = 2.591.
            theta = .66411*np.sqrt(NU*nodes[0]/speed(nodes[0]))
            first = np.array([theta, 2.5911*theta, 0., speed(nodes[0])])
        else:
            first = None
        r = closure_march(nodes, speed(nodes), param, start=first, function=speed)
        add(f'{name} separation, inherited closure', r['xi_separation']*scale if r['status'] == 'separated' else np.nan,
            reference, 0., gate=False)
        t = fd.thwaites(nodes, speed(nodes), slope(nodes), NU, k_stagnation=True if start == 'stagnation' else None)
        add(f'{name} separation, Thwaites', t['xi_separation']*scale if t['xi_separation'] else np.nan, reference, 0., gate=False)
    return rows, all(r['passed'] for r in rows if r['gate'])


def classify_earlier(surfaces, param):
    """Part A: earlier failures against the direct-mode singular interval."""
    earlier = {(r['case'], r['profile']): r for r in json.loads((EARLIER/'summary.json').read_text())['cases']}
    rows, folds, determinants = [], [], []

    def determinant(xi, states, j):
        jac = mf.residual_station(param, xi[j-1:j+1], states[:, j-1:j+1], np.zeros(2))[1]
        return float(np.linalg.det(jac[:2, 4:6]*states[:2, j]))

    def first_change(values):
        for k in range(1, len(values)):
            if np.isfinite(values[k-1]) and values[k-1] < 0 < values[k]:
                return k
        return None
    for case in ('cl0056', 'cl0175'):
        s = surfaces[case]
        laminar = int(np.where(s.turbulent > 0)[0][0])
        native = [np.nan]+[determinant(s.xi, s.U, j) for j in range(1, laminar)]
        change = first_change(native)
        for j, value in enumerate(native):
            determinants.append(dict(case=case, station=j, node=int(s.node[j]), x=float(s.x[j]),
                                     Hk=float(mf.get_Hk(s.U[:, j], param)[0]), determinant=None if j == 0 else value))
        for kind in ('native', 'model_taps', 'experiment_taps', 'experiment_delta'):
            record = earlier[(case, kind)]
            saved = np.load(EARLIER/f'{case}_{kind}.npz')
            ids = saved['ids']
            interval = int(np.where(ids == s.node[change])[0][0])  # downstream end, earlier numbering
            own = [np.nan]+[determinant(saved['xi'], saved['U'], j) for j in range(1, len(ids))]
            own_change = first_change(own)
            match = re.match(r'station (\d+):', record['last_failure'] or '')
            failed = int(match.group(1)) if match else None
            rows.append(dict(case=case, profile=kind, complete=record['complete'], reached_fraction=record['reached_fraction'],
                last_failed_station=failed, native_singular_interval=f'{interval-1}-{interval}',
                x_interval=f'{s.x[change-1]:.5f}-{s.x[change]:.5f}',
                at_or_after_interval=None if failed is None else bool(failed >= interval-1),
                saved_state_first_sign_change=own_change))
        left, old = s.U[:, change-1], s.U[:, change]
        for eps in (-1e-3, -1e-4, -1e-5, -1e-6, 1e-6, 1e-5, 1e-4, 1e-3):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('error', RuntimeWarning)
                    new, error, _ = station(left, old.copy(), old[3]*(1+eps), s.xi[change-1:change+1], param)
                folds.append(dict(case=case, x=float(s.x[change]), relative_speed_change=eps, root_found=True,
                    Hk=float(mf.get_Hk(new, param)[0]), delta_star_relative_change=float(new[1]/old[1]-1), message=''))
            except (ValueError, FloatingPointError, np.linalg.LinAlgError, RuntimeWarning) as exc:
                folds.append(dict(case=case, x=float(s.x[change]), relative_speed_change=eps, root_found=False,
                    Hk=None, delta_star_relative_change=None, message=str(exc)))
    return rows, folds, determinants


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    out = parser.parse_args().output
    if (out/'summary.json').exists():
        raise SystemExit('A completed study exists here; choose a new output directory')
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    log = (out/'run.log').open('w')

    def say(text):
        print(text, flush=True); log.write(text+'\n'); log.flush()

    # Frozen inputs.
    receipt = json.loads(STATE_RECEIPT.read_text()); files = {}
    for case, (directory, _) in CASES.items():
        path = STATES/directory/'fore_state.npz'
        assert digest(path) == receipt[directory+'/fore_state.npz'], 'Archived state changed: '+directory
    assert digest(TAPS) == json.loads(TAP_RECEIPT.read_text())['output_sha256']['pressure_at_taps.csv'], 'Tap table changed'
    earlier_receipt = json.loads((EARLIER/'verification.json').read_text())['output_sha256']
    for name, expected in earlier_receipt.items():
        assert digest(EARLIER/name) == expected, 'Earlier replay output changed: '+name
    tracked = [STATES/d/'fore_state.npz' for d, _ in CASES.values()]+[STATE_RECEIPT, TAPS, TAP_RECEIPT, PROTOCOL,
        EARLIER/'summary.json', EARLIER/'verification.json', Path(__file__), ROOT/'validation/s414_laminar_fd.py',
        ROOT/'validation/run_s414_prescribed_pressure.py', ROOT/'validation/s414_geometry.py', ROOT/'viscous.py',
        ROOT/'amplification_jacobian.py', ROOT/'vendor/mfoil.py', ROOT/'multifoil.py', ROOT/'paneling.py']
    tracked += sorted(EARLIER.glob('*.npz'))
    tracked += [RESULTS/'s414_geometry_qualification_20261002'/f'report_test_{n}.dat' for n in ('fore', 'aft')]
    for path in tracked:
        files[str(Path(path).relative_to(ROOT))] = digest(path)
    environment = dict(python=platform.python_version(),
        packages={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')},
        threads={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    write(out/'input_manifest.json', dict(files=files, environment=environment))

    model10, param10 = parameters(MACH)
    _, param0 = parameters(0.)
    assert abs(param0.mu0/param0.rho0-NU) < 1e-18
    pressure = Pressure(param10)
    arrays = {case: np.load(STATES/directory/'fore_state.npz') for case, (directory, _) in CASES.items()}
    for a in arrays.values():
        assert np.max(abs(model10.foil.x-a['xy'])) < 1e-13, 'Rebuilt geometry differs from the archive'
    surfaces = {case: Surface(a, pressure) for case, a in arrays.items()}
    taps = {}
    for row in csv.DictReader(TAPS.open()):
        taps.setdefault((row['case'], row['coordinate']), []).append((float(row['x']), float(row['experiment_Cp']), float(row['viscous_Cp'])))
    taps = {key: np.array(value[:TAP_COUNT]) for key, value in taps.items()}

    # Part A.
    rows_a, folds, determinants = classify_earlier(surfaces, param10)
    write_csv(out/'earlier_failures.csv', rows_a); write_csv(out/'fold_demonstration.csv', folds)
    write_csv(out/'direct_determinant.csv', determinants)
    classified = all(r['at_or_after_interval'] for r in rows_a if r['last_failed_station'] is not None)
    say(f'Part A: {sum(r["last_failed_station"] is not None for r in rows_a)} recorded failures, '
        f'all at or after the singular interval: {classified}')

    # Part B (gate).
    rows_b, gate = qualification(param0)
    write_csv(out/'qualification.csv', rows_b)
    say(f'Part B: independent-solver gate {"passed" if gate else "FAILED"}')
    if not gate:
        write(out/'summary.json', dict(state='stopped', reason='Part B gate failed',
                                       failed_controls=[r['control'] for r in rows_b if r['gate'] and not r['passed']]))
        raise SystemExit('Qualification gate failed; Parts C-F not run')

    marches, store = [], {}
    counts = dict(fd=0, closure=0)

    def run(identifier, part, surface, nodes, function, **labels):
        result, data = response(surface, nodes, function, param0, with_thwaites=part == 'C')
        counts['fd'] += 1; counts['closure'] += 1
        for key, value in data.items():
            store[identifier+'/'+key] = value
        record = dict(id=identifier, part=part, **labels, **{f'fd_{k}': v for k, v in result['fd'].items()},
                      **{f'closure_{k}': v for k, v in result['closure'].items()})
        if 'thwaites' in result:
            record['thwaites_x_separation'] = result['thwaites']['x_separation']
        marches.append(record)
        say(f'{identifier}: FD {result["fd"]["status"]} x={result["fd"]["x_separation"]}, '
            f'closure {result["closure"]["status"]} x={result["closure"]["x_separation"]}')
        return result

    # Part C: the model's own pressure.
    native, part_c, sensitivity, archive_control = {}, [], [], []
    for case, s in surfaces.items():
        tap_xi = s.xi_of_x(taps[(case, 'table')][:, 0])
        native[case] = run('native_'+case, 'C', s, s.nodes(tap_xi), s.q, base=case, target=case)
        screen = compare(native[case]['closure'], native[case]['fd'])
        part_c.append(dict(case=case, fd_status=native[case]['fd']['status'], closure_status=native[case]['closure']['status'],
            fd_x_separation=native[case]['fd']['x_separation'], closure_x_separation=native[case]['closure']['x_separation'],
            thwaites_x_separation=native[case]['thwaites']['x_separation'], **screen))
        for label, surface, sub in (('sub4', s, 4), ('sub16', s, 16), ('spline', Surface(arrays[case], pressure, 'spline'), SUB)):
            r = run(f'native_{case}_{label}', 'C-sensitivity', surface, surface.nodes(tap_xi, sub), surface.q, base=case, target=case)
            sensitivity.append(dict(case=case, variant=label,
                fd_shift=None if r['fd']['x_separation'] is None or native[case]['fd']['x_separation'] is None
                else r['fd']['x_separation']-native[case]['fd']['x_separation'], fd_status=r['fd']['status'],
                closure_shift=None if r['closure']['x_separation'] is None or native[case]['closure']['x_separation'] is None
                else r['closure']['x_separation']-native[case]['closure']['x_separation'], closure_status=r['closure']['status']))
        # Closure at the archived Mach number on native nodes, against the archive.
        end = int(np.sum(s.xi < s.xi_of_x(X_END)))
        r = closure_march(s.xi[:end], s.U[3, :end], param10); counts['closure'] += 1
        n = len(r['xi'])-(1 if r['status'] == 'separated' else 0)
        archive_control.append(dict(case=case, status=r['status'], stations_compared=n,
            theta_max_relative=float(np.max(abs(r['states'][0, :n]/s.U[0, :n]-1))),
            delta_star_max_relative=float(np.max(abs(r['states'][1, :n]/s.U[1, :n]-1))),
            theta_last_relative=float(r['states'][0, n-1]/s.U[0, n-1]-1)))
        store[f'archive_{case}/x'] = s.x; store[f'archive_{case}/xi'] = s.xi; store[f'archive_{case}/U'] = s.U
        store[f'archive_{case}/turbulent'] = s.turbulent
        store[f'archive_{case}/Hk'] = np.array([mf.get_Hk(u, param10)[0] for u in s.U.T])
        store[f'archive_{case}/cp'] = np.array([mf.get_cp(u, param10)[0] for u in s.U[3]])
    write_csv(out/'native_comparison.csv', part_c); write_csv(out/'native_sensitivity.csv', sensitivity)
    write_csv(out/'archive_control.csv', archive_control)

    # Part D: sampling controls.
    controls = []
    for case, s in surfaces.items():
        tap_xi = s.xi_of_x(taps[(case, 'table')][:, 0]); own = s.cp(tap_xi); nodes = s.nodes(tap_xi)
        table_error = float(np.max(abs(own-taps[(case, 'table')][:, 2])))
        for family in FAMILIES:
            function = reconstruction(s, tap_xi, own, family, 'hold')
            if family in ('additive', 'ratio'):
                identity = float(np.max(abs(function(nodes)-s.q(nodes))))
                assert identity < 1e-12
                controls.append(dict(kind='self', base=case, target=case, family=family, upstream='all', method='fd',
                    trivial=True, same_topology=True, separation_shift=0., theta_relative=0., passed=True,
                    identity_error=identity, tap_table_difference=table_error))
                continue
            for other in ('ramp', 'late'):
                assert np.max(abs(reconstruction(s, tap_xi, own, family, other)(nodes)-function(nodes))) < 1e-12
            r = run(f'self_{case}_{family}', 'D-self', s, nodes, function, base=case, target=case, family=family, upstream='all')
            for method in ('fd', 'closure'):
                controls.append(dict(kind='self', base=case, target=case, family=family, upstream='all', method=method,
                    trivial=False, **compare(r[method], native[case][method]), identity_error=None, tap_table_difference=table_error))
    for base, target in CROSS:
        a, b = surfaces[base], surfaces[target]
        x_taps = taps[(base, 'table')][:, 0]
        tap_xi = a.xi_of_x(x_taps); values = b.cp(b.xi_of_x(x_taps)); nodes = a.nodes(tap_xi)
        for family in FAMILIES:
            for upstream in UPSTREAM:
                r = run(f'cross_{base}_to_{target}_{family}_{upstream}', 'D-cross', a, nodes,
                        reconstruction(a, tap_xi, values, family, upstream), base=base, target=target, family=family, upstream=upstream)
                for method in ('fd', 'closure'):
                    controls.append(dict(kind='cross', base=base, target=target, family=family, upstream=upstream, method=method,
                        trivial=False, **compare(r[method], native[target][method]), identity_error=None, tap_table_difference=None))
    write_csv(out/'controls.csv', controls)

    # Part E, closure only: restart at the earlier seed node with scaled thicknesses.
    seeds = []
    for case in ('cl0056', 'cl0175'):
        s = surfaces[case]
        tap_xi = s.xi_of_x(taps[(case, 'table')][:, 0]); nodes = s.nodes(tap_xi)
        seed_xi = s.xi[s.xi < tap_xi[0]][-1]; k = int(np.argmin(abs(nodes-seed_xi)))
        assert abs(nodes[k]-seed_xi) < 1e-12
        reference = closure_march(nodes, s.q(nodes), param0, function=s.q)  # the native march above; not counted again
        for factor in (.8, 1.25):
            at = int(np.argmin(abs(reference['xi']-seed_xi))); assert abs(reference['xi'][at]-seed_xi) < 1e-12
            start = reference['states'][:, at].copy(); start[:2] *= factor
            r = closure_march(nodes[k:], s.q(nodes[k:]), param0, start=start, function=s.q); counts['closure'] += 1
            position = lambda m: None if m['xi_separation'] is None else s.x_of_xi(m['xi_separation'])
            screen = compare(dict(status=r['status'], x_separation=position(r), theta_reference=r['theta_reference']),
                dict(status=reference['status'], x_separation=position(reference), theta_reference=reference['theta_reference']))
            seeds.append(dict(case=case, seed_x=s.x_of_xi(seed_xi), thickness_factor=factor, status=r['status'],
                reference_status=reference['status'], separation_shift=screen['separation_shift'],
                theta_relative=screen['theta_relative'], bounded=screen['passed']))
    write_csv(out/'seed_sensitivity.csv', seeds)

    # Part F: measured taps.
    experiment = []
    for case, s in surfaces.items():
        for coordinate, upstreams in (('table', UPSTREAM), ('plotted', ('hold',))):
            data = taps[(case, coordinate)]; tap_xi = s.xi_of_x(data[:, 0]); nodes = s.nodes(tap_xi)
            for family in FAMILIES:
                for upstream in upstreams:
                    identifier = f'exp_{case}_{family}_{upstream}'+('' if coordinate == 'table' else '_plotted')
                    r = run(identifier, 'F', s, nodes, reconstruction(s, tap_xi, data[:, 1], family, upstream),
                            base=case, target='experiment', family=family, upstream=upstream, coordinate=coordinate)
                    relevant = [c for c in controls if c['method'] == 'fd' and c['base'] == case and c['family'] == family
                                and c['upstream'] in (upstream, 'all')]
                    experiment.append(dict(case=case, coordinate=coordinate, family=family, upstream=upstream,
                        admitted=bool(all(c['passed'] for c in relevant)), controls_passed=sum(c['passed'] for c in relevant),
                        controls=len(relevant), fd_status=r['fd']['status'], fd_x_separation=r['fd']['x_separation'],
                        fd_theta=r['fd']['theta_reference'], closure_status=r['closure']['status'],
                        closure_x_separation=r['closure']['x_separation'], closure_theta=r['closure']['theta_reference'],
                        methods_agree=bool(r['fd']['status'] == r['closure']['status'])))
    write_csv(out/'experiment.csv', experiment)

    # Part E spreads over the upstream treatments, and the declared outcome per lift.
    spreads = []
    groups = [('cross', b, t) for b, t in CROSS]+[('experiment', c, 'experiment') for c in CASES]
    for kind, base, target in groups:
        for family in FAMILIES:
            prefix = f'cross_{base}_to_{target}_{family}_' if kind == 'cross' else f'exp_{base}_{family}_'
            rows = [m for m in marches if m['id'] in [prefix+u for u in UPSTREAM]]
            for method in ('fd', 'closure'):
                statuses = {m[f'{method}_status'] for m in rows}
                positions = [m[f'{method}_x_separation'] for m in rows if m[f'{method}_x_separation'] is not None]
                thetas = [m[f'{method}_theta_reference'] for m in rows]
                spread = max(positions)-min(positions) if len(positions) == len(rows) else None
                theta_spread = max(thetas)/min(thetas)-1
                spreads.append(dict(kind=kind, base=base, target=target, family=family, method=method,
                    one_topology=len(statuses) == 1, separation_spread=spread, theta_spread=float(theta_spread),
                    bounded=bool(len(statuses) == 1 and (spread is None or spread <= TOL_X) and theta_spread <= TOL_THETA)))
    write_csv(out/'upstream_spread.csv', spreads)
    outcomes = {}
    for case in CASES:
        admitted = [e for e in experiment if e['case'] == case and e['coordinate'] == 'table' and e['admitted']]
        if not admitted or not all(e['methods_agree'] for e in admitted):
            outcome = 'not bounded'
        elif all(e['fd_status'] == 'separated' for e in admitted):
            outcome = 'separation retained'
        elif all(e['fd_status'] == 'completed' for e in admitted):
            outcome = 'separation removed'
        else:
            outcome = 'not bounded'
        everything = [e for e in experiment if e['case'] == case and e['coordinate'] == 'table']
        outcomes[case] = dict(outcome=outcome, admitted=len(admitted), reconstructions=len(everything),
            fd_separated_all=sum(e['fd_status'] == 'separated' for e in everything),
            fd_separated_admitted=sum(e['fd_status'] == 'separated' for e in admitted))
    write_csv(out/'marches.csv', [{k: m.get(k) for k in sorted({k for m in marches for k in m})} for m in marches])
    np.savez_compressed(out/'profiles.npz', **store)
    summary = dict(state='finished', part_a=dict(failures_classified_as_direct_singularity=classified, rows=rows_a),
        part_b_gate_passed=gate, part_c=part_c, archive_control=archive_control,
        controls=dict(total=len(controls), passed=sum(c['passed'] for c in controls),
                      fd_nontrivial_passed=sum(c['passed'] for c in controls if c['method'] == 'fd' and not c['trivial']),
                      fd_nontrivial=sum(1 for c in controls if c['method'] == 'fd' and not c['trivial'])),
        seed_sensitivity=seeds, outcomes=outcomes, counts=dict(fd_marches=counts['fd'], closure_marches=counts['closure'],
        coupled_flow_solves=0, forces_computed=False), tolerances=dict(x=TOL_X, theta=TOL_THETA),
        fd_scheme='bdf2', fd_max_sawtooth_index=max(m['fd_oscillation'] for m in marches),
        wall_seconds=time.monotonic()-started)
    assert counts['fd'] <= 120 and counts['closure'] <= 130
    write(out/'summary.json', summary)
    say(json.dumps(dict(outcomes=outcomes, counts=summary['counts'], controls=summary['controls']), indent=2))
    log.close()


if __name__ == '__main__':
    main()
