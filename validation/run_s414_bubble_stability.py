#!/usr/bin/env python3
"""V5 tranche 2: independent boundary layer and linear-stability amplification through the archived nose bubble.

No coupled solve, no forces. See validation/S414_BUBBLE_STABILITY_PROTOCOL.md.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_s414_laminar_replay as replay
from run_s414_laminar_replay import ROOT, RESULTS, CASES, STATES, STATE_RECEIPT, NU, MACH, mf, digest, write, write_csv
import s414_laminar_fd as fd
import s414_bubble_fd as bf
import s414_stability as stability
import s414_amplification as amp

OUT = RESULTS/'s414_bubble_stability_20261005'
PROTOCOL = ROOT/'validation/S414_BUBBLE_STABILITY_PROTOCOL.md'
X_END, SWITCH_HK, SUB, FREQUENCIES, BAND = .02082853867960402, 2.9, 8, 52, (.004, .6)
TOL_CP, TOL_THETA, RATIO = .05, .05, (.80, 1.25)
BUBBLES = ('cl0056', 'cl0175')
SIMILAR = ((3.5, np.linspace(8, 152, 37)), (5., np.linspace(8, 152, 37)), (6., np.linspace(8, 152, 37)),
           (6.75, np.linspace(8, 152, 37)))


def collocation():
    return stability.Grid(140, 60., 3.)


def similar_stations(m, profile, re_theta, eta, grid):
    """Stations of a similar flow ue = x**m at the requested momentum-thickness Reynolds numbers."""
    d1, d2 = bf.measures(eta, *profile)[:2]
    out = []
    for value in re_theta:
        x = (NU*(value/d2)**2)**(1/(1+m)); ue = x**m; scale = np.sqrt(NU*x/ue)
        item = amp.station(x, d1*scale, ue, NU, eta, profile[1], np.gradient(profile[2], eta), d1, grid)
        item['theta'] = d2*scale; out.append(item)
    return out, d1/d2


def frequencies(stations, count=FREQUENCIES):
    return np.geomspace(*BAND, count)*np.median([s['ue']/s['theta'] for s in stations])


def stability_envelope(stations, grid, count=FREQUENCIES):
    omegas = frequencies(stations, count)
    rate, wavenumber = amp.growth(stations, omegas, grid)
    xi = np.array([s['xi'] for s in stations])
    each, total = amp.envelope(xi, rate)
    lost = int(np.sum(np.isnan(each[:, -1]) & np.any(each > 0, axis=1)))
    return dict(xi=xi, omegas=omegas, rate=rate, each=each, N=total, lost_after_start=lost,
                tracked=int(np.isfinite(rate).sum()))


def correlation(xi, theta, delta_star, ue, param, start=0.):
    """The inherited growth function integrated with the inherited trapezoid rule."""
    damp = lambda j, n: mf.get_damp(np.array([theta[j], delta_star[j], n, ue[j]]), param)[0]
    values = [start]
    for j in range(1, len(xi)):
        n, step = values[-1], xi[j]-xi[j-1]
        predicted = n+damp(j-1, n)*step
        values.append(n+.5*(damp(j-1, n)+damp(j, predicted))*step)
    return np.array(values)


def qualification(param0, say):
    """Part A. Returns the control rows, the similar-profile table and the gate decision."""
    rows, eta, grid = [], bf.eta_grid(), collocation()

    def add(name, value, reference, tolerance, gate=True):
        rows.append(dict(control=name, value=float(value), reference=reference, tolerance=tolerance,
                         passed=bool(abs(value-reference) <= tolerance) if gate else None, gate=gate))
    blasius = fd.similarity(0., eta)
    d1, d2 = bf.measures(eta, *blasius)[:2]
    profile = lambda g: stability.on_grid(g, eta/d1, blasius[1], np.gradient(blasius[2], eta)*d1**2)
    u, upp = profile(grid)
    frequency, shape = stability.temporal(.308, 998., u, upp, grid)
    add('Blasius temporal phase speed, real part', (frequency/.308).real, .3641, 5e-4)
    add('Blasius temporal phase speed, imaginary part', (frequency/.308).imag, .0080, 5e-4)
    start, shape = stability.seed(.1122, stability.scan(998., u, upp, grid, np.linspace(.2, .4, 11)))
    alpha, _, ok = stability.spatial(.1122, 998., u, upp, grid, start, shape)
    add('Blasius spatial wavenumber, real part', alpha.real if ok else np.nan, .3086, 5e-4)
    add('Blasius spatial wavenumber, imaginary part', alpha.imag if ok else np.nan, -.0057, 5e-4)
    fine = stability.Grid(200, 60., 3.); uf, uppf = profile(fine)
    startf, shapef = stability.seed(.1122, stability.scan(998., uf, uppf, fine, np.linspace(.2, .4, 11)))
    alphaf, _, okf = stability.spatial(.1122, 998., uf, uppf, fine, startf, shapef)
    add('collocation 141 against 201 points', abs(alpha-alphaf) if ok and okf else np.nan, 0., 1e-6)
    for reynolds, sign in ((515., -1), (525., 1)):
        found = stability.scan(reynolds, u, upp, grid, np.linspace(.2, .4, 11))
        best = -np.inf
        for omega in np.linspace(.10, .14, 21):
            start, shape = stability.seed(omega, found)
            if start is not None:
                value, _, ok = stability.spatial(omega, reynolds, u, upp, grid, start, shape)
                best = max(best, -value.imag) if ok else best
        rows.append(dict(control=f'Blasius largest spatial growth at Re_delta* {reynolds:.0f}', value=float(best),
                         reference=0., tolerance=0., passed=bool(sign*best > 0), gate=True))
    # Falkner-Skan family through separation.
    members = bf.family(np.round(np.arange(1.75, 10.01, .25), 2), eta)
    table = np.array([(m, 2*m/(m+1), bf.measures(eta, *p)[2]/np.sqrt((m+1)/2), bf.measures(eta, *p)[0]/bf.measures(eta, *p)[1])
                      for _, m, p in members])
    k = int(np.argmin(table[:, 1])); lower = table[k:]
    add('Falkner-Skan separation, beta', table[k, 1], -.1988, 5e-4)
    add('Falkner-Skan separation, H', table[k, 3], 4.03, .01)
    for beta, reference in ((-.18, -.0976), (-.15, -.1334), (-.10, -.1405), (-.05, -.1082)):
        add(f'Falkner-Skan reversed-flow wall shear, beta {beta}', np.interp(beta, lower[:, 1], lower[:, 2]), reference, 1e-3)
    nodes = np.linspace(.02, 1.7, 85); cylinder = lambda x: 2*np.sin(x)
    new, old = bf.march(nodes, cylinder, NU), fd.march(nodes, cylinder, NU, m_start=1.)
    add('direct mode against the tranche-1 solver, theta', np.max(abs(new['theta']/old['theta']-1)), 0., 1e-8)
    back = bf.march(nodes, cylinder, NU, switch=.6, displacement=PchipInterpolator(new['xi'], new['delta_star']))
    add('direct solution fed back as displacement, edge-speed ratio', np.max(abs(back['ratio']-1)), 0., 1e-8)
    say('Part A: eigenvalue, similar-profile and inverse-mode controls done')
    # e^N envelope of similar flows against the inherited correlation.
    similar = []

    def compare(name, m, profile, re_theta, count=FREQUENCIES):
        stations, shape_factor = similar_stations(m, profile, re_theta, eta, grid)
        result = stability_envelope(stations, grid, count)
        theta = np.array([s['theta'] for s in stations]); delta = np.array([s['delta_star'] for s in stations])
        inherited = correlation(result['xi'], theta, delta, np.array([s['ue'] for s in stations]), param0)
        return shape_factor, result, inherited
    re_blasius = np.linspace(150, 1300, 47)
    h, result, inherited = compare('Blasius', 0., blasius, re_blasius)
    at6 = (np.interp(6, result['N'], re_blasius), np.interp(6, inherited, re_blasius))
    add('Blasius e^N: Re_theta at N = 6, Orr-Sommerfeld over inherited correlation', at6[0]/at6[1], 1., .10)
    add('Blasius e^N: Re_theta at N = 9, Orr-Sommerfeld', np.interp(9, result['N'], re_blasius), np.interp(9, inherited, re_blasius),
        0., gate=False)
    coarse = compare('Blasius', 0., blasius, re_blasius, 26)[1]
    add('Blasius N at Re_theta 1050, 52 minus 26 frequencies', np.interp(1050, re_blasius, result['N'])
        - np.interp(1050, re_blasius, coarse['N']), 0., 0., gate=False)
    for value in (300, 600, 900, 1200):
        similar.append(dict(profile='Blasius', H=h, re_theta=value, N_orr_sommerfeld=float(np.interp(value, re_blasius, result['N'])),
                            N_inherited=float(np.interp(value, re_blasius, inherited)), lost_waves=result['lost_after_start']))
    say('Part A: Blasius envelope done')
    wanted = [d for d, _ in SIMILAR]
    exact = {d: (m, p) for d, m, p in bf.family(wanted, eta)}
    approximate = {d: (m, p) for d, m, p in bf.family(wanted, eta, flare=True)}
    for displacement, re_theta in SIMILAR:
        m, profile = exact[displacement]
        h, result, inherited = compare(f'similar {displacement}', m, profile, re_theta)
        for value in (32, 56, 80, 104, 128, 152):
            similar.append(dict(profile=f'reversed similar, scaled displacement {displacement}', H=h, re_theta=value,
                N_orr_sommerfeld=float(np.interp(value, re_theta, result['N'])),
                N_inherited=float(np.interp(value, re_theta, inherited)), lost_waves=result['lost_after_start']))
        slope = lambda n: (np.interp(152, re_theta, n)-np.interp(104, re_theta, n))/48
        add(f'similar H {h:.2f}: envelope slope dN/dRe_theta over 104-152, Orr-Sommerfeld over inherited',
            slope(result['N'])/slope(inherited), 1., 0., gate=False)
        mf_, pf = approximate[displacement]
        a, b = bf.measures(eta, *profile), bf.measures(eta, *pf)
        add(f'similar H {h:.2f}: FLARE change of wall shear, relative', b[2]/a[2]-1 if a[2] else 0., 0., 0., gate=False)
        add(f'similar H {h:.2f}: FLARE change of H, relative', (b[0]/b[1])/(a[0]/a[1])-1, 0., 0., gate=False)
        add(f'similar H {h:.2f}: strongest reversed velocity over edge speed', a[3], 0., 0., gate=False)
        say(f'Part A: similar profile H {h:.2f} done')
    return rows, similar, all(r['passed'] for r in rows if r['gate'])


class Displacement:
    """Prescribed displacement thickness: the solver's own value at the switch node, archived values beyond."""
    def __init__(self, xi, archived, kind='pchip', scaled=False):
        self.xi, self.archived, self.kind, self.scaled = np.asarray(xi), np.asarray(archived, dtype=float), kind, scaled
        self.own = None

    def join(self, xi, own):
        assert abs(xi-self.xi[0]) < 1e-12, 'The march did not stop at the switch node'
        values = self.archived*(own/self.archived[0] if self.scaled else 1.)
        values[0] = own; self.own = float(own); self.values = values
        self.function = PchipInterpolator(self.xi, values) if self.kind == 'pchip' else None

    def __call__(self, xi):
        return float(self.function(xi)) if self.function is not None else float(np.interp(xi, self.xi, self.values))


def refine(xi, sub):
    return np.concatenate([np.linspace(a, b, sub+1)[:-1] for a, b in zip(xi[:-1], xi[1:])]+[[xi[-1]]])


def bubble(s, param10, sub=SUB, switch_shift=0, kind='pchip', scaled=False):
    """Part B for one archived state. Returns the march, the node table and bookkeeping."""
    hk = np.array([mf.get_Hk(u, param10)[0] for u in s.U.T])
    end = int(np.argmin(abs(s.x-X_END))); assert abs(s.x[end]-X_END) < 1e-9 and not s.turbulent[end]
    switch = s.le+int(np.where(hk[s.le:end+1] >= SWITCH_HK)[0][0])-switch_shift
    nodes = refine(s.xi[:end+1], sub)
    target = Displacement(s.xi[switch:end+1], s.U[1, switch:end+1], kind, scaled)
    started = time.monotonic()
    result = bf.march(nodes, s.q, NU, switch=float(s.xi[switch]), displacement=target, keep_profiles=True)
    seconds = time.monotonic()-started
    peak = s.le+int(np.argmax(s.q_nodes[s.le:end+1]))
    rows = []
    if result['status'] == 'completed':
        for j in range(peak, end+1):
            k = int(np.argmin(abs(result['xi']-s.xi[j]))); assert abs(result['xi'][k]-s.xi[j]) < 1e-12
            param10.turb = False
            rows.append(dict(node=int(s.node[j]), x=float(s.x[j]), mode='inverse' if j > switch else 'direct',
                Cp_archive=float(s.pressure.cp(s.q_nodes[j])), Cp_independent=float(s.pressure.cp(result['ue'][k])),
                theta_archive=float(s.U[0, j]), theta_independent=float(result['theta'][k]),
                delta_star_archive=float(s.U[1, j]), delta_star_independent=float(result['delta_star'][k]),
                H_archive=float(hk[j]), H_independent=float(result['H'][k]),
                cf_Re_theta_archive=replay.friction(s.U[:, j], param10), cf_Re_theta_independent=float(result['cf_re_theta'][k]),
                reversed_velocity=float(result['reverse'][k]), N_archive=float(s.U[2, j])))
    return dict(result=result, rows=rows, switch=switch, end=end, peak=peak, nodes=nodes, seconds=seconds, own=target.own)


def amplification(s, run, param0, grid):
    """Part C for one march: Orr-Sommerfeld envelope and the inherited correlation on the independent state."""
    result = run['result']
    marching = np.array([np.min(abs(run['nodes']-x)) < 1e-12 for x in result['xi']])
    index = np.where(marching)[0]
    xi, theta, delta, ue = (result[k][index] for k in ('xi', 'theta', 'delta_star', 'ue'))
    inherited = correlation(xi, theta, delta, ue, param0)
    first = int(np.argmin(abs(xi-s.xi[run['peak']])))
    stations = []
    for k in index[first:]:
        p = result['profiles'][k]
        item = amp.station(p['xi'], result['delta_star'][k], result['ue'][k], NU, p['eta'], p['ratio'], p['curvature'],
                           p['scaled_displacement'], grid)
        item['theta'] = float(result['theta'][k]); stations.append(item)
    started = time.monotonic()
    waves = stability_envelope(stations, grid)
    return dict(xi=xi, inherited=inherited, waves=waves, seconds=time.monotonic()-started,
                N_at=lambda x: float(np.interp(x, waves['xi'], waves['N'], left=0.)),
                inherited_at=lambda x: float(np.interp(x, xi, inherited)))


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

    receipt = json.loads(STATE_RECEIPT.read_text()); files = {}
    for directory, _ in CASES.values():
        assert digest(STATES/directory/'fore_state.npz') == receipt[directory+'/fore_state.npz'], 'Archived state changed'
    tracked = [STATES/d/'fore_state.npz' for d, _ in CASES.values()]+[STATE_RECEIPT, PROTOCOL, Path(__file__)]
    tracked += [ROOT/'validation'/n for n in ('s414_bubble_fd.py', 's414_stability.py', 's414_amplification.py', 's414_laminar_fd.py',
        'run_s414_laminar_replay.py', 'run_s414_prescribed_pressure.py', 's414_geometry.py')]
    tracked += [ROOT/n for n in ('viscous.py', 'amplification_jacobian.py', 'vendor/mfoil.py', 'multifoil.py', 'paneling.py')]
    for path in tracked:
        files[str(Path(path).relative_to(ROOT))] = digest(path)
    environment = dict(python=platform.python_version(), packages={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')},
        threads={k: os.environ[k] for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    write(out/'input_manifest.json', dict(files=files, environment=environment))

    _, param10 = replay.parameters(MACH)
    _, param0 = replay.parameters(0.)
    pressure = replay.Pressure(param10)

    rows_a, similar, gate = qualification(param0, say)
    write_csv(out/'qualification.csv', rows_a); write_csv(out/'similar_profiles.csv', similar)
    say(f'Part A: gate {"passed" if gate else "FAILED"}')
    if not gate:
        write(out/'summary.json', dict(state='stopped', reason='Part A gate failed',
                                       failed_controls=[r['control'] for r in rows_a if r['gate'] and not r['passed']]))
        raise SystemExit('Qualification gate failed; Parts B and C not run')

    grid = collocation()
    surfaces = {case: replay.Surface(np.load(STATES/directory/'fore_state.npz'), pressure) for case, (directory, _) in CASES.items()}
    variants = [('primary', {})]+[(name, options) for name, options in (('switch_upstream', dict(switch_shift=1)),
        ('linear_displacement', dict(kind='linear')), ('scaled_displacement', dict(scaled=True)), ('sub16', dict(sub=16)))]
    node_rows, amp_rows, screens, store = [], [], [], {}
    counts = dict(marches=0, envelopes=0, failed_marches=0, lost_waves=0)
    for case, s in surfaces.items():
        for name, options in variants:
            if name != 'primary' and case not in BUBBLES:
                continue
            run = bubble(s, param10, **options); counts['marches'] += 1
            result = run['result']; label = f'{case}/{name}'
            say(f'{label}: march {result["status"]} {result["stop_reason"]} ({len(result["xi"])} stations, {run["seconds"]:.0f} s), '
                f'switch x/c {s.x[run["switch"]]:.5f}')
            if result['status'] != 'completed':
                counts['failed_marches'] += 1
                screens.append(dict(case=case, variant=name, march=result['status'], reason=result['stop_reason']))
                continue
            for row in run['rows']:
                node_rows.append(dict(case=case, variant=name, **row))
            inverse = [r for r in run['rows'] if r['mode'] == 'inverse']
            cp_error = max(abs(r['Cp_independent']-r['Cp_archive']) for r in inverse)
            theta_error = run['rows'][-1]['theta_independent']/run['rows'][-1]['theta_archive']-1
            record = dict(case=case, variant=name, march='completed', reason='', switch_x=float(s.x[run['switch']]),
                own_over_archived_displacement_at_switch=run['own']/float(s.U[1, run['switch']]),
                largest_Cp_difference=float(cp_error), theta_endpoint_relative=float(theta_error),
                strongest_reversed_velocity=float(min(r['reversed_velocity'] for r in run['rows'])),
                H_endpoint_archive=run['rows'][-1]['H_archive'], H_endpoint_independent=run['rows'][-1]['H_independent'],
                closure_screen=bool(cp_error <= TOL_CP and abs(theta_error) <= TOL_THETA) if case in BUBBLES else None)
            x = np.array([s.x_of_xi(v) for v in result['xi']])
            for key in ('xi', 'ue', 'ratio', 'delta_star', 'theta', 'H', 'cf_re_theta', 'reverse'):
                store[f'{label}/{key}'] = result[key]
            store[f'{label}/x'] = x; store[f'{label}/cp'] = pressure.cp(result['ue'])
            if name in ('primary', 'switch_upstream', 'linear_displacement'):
                waves = amplification(s, run, param0, grid); counts['envelopes'] += 1
                counts['lost_waves'] += waves['waves']['lost_after_start']
                end_xi = float(s.xi[run['end']])
                n_os, n_inherited = waves['N_at'](end_xi), waves['inherited_at'](end_xi)
                crossing = None
                if waves['waves']['N'].max() >= 9:
                    crossing = s.x_of_xi(float(np.interp(9, waves['waves']['N'], waves['waves']['xi'])))
                record.update(N_endpoint_archive=float(s.U[2, run['end']]), N_endpoint_inherited_on_independent_state=n_inherited,
                    N_endpoint_orr_sommerfeld=n_os, ratio_orr_sommerfeld_over_inherited=n_os/n_inherited if n_inherited > 0 else None,
                    growth_screen=bool(RATIO[0] <= n_os/n_inherited <= RATIO[1]) if case in BUBBLES and n_inherited > 0 else None,
                    orr_sommerfeld_reaches_9_at_x=crossing, lost_waves=waves['waves']['lost_after_start'],
                    waves_tracked=waves['waves']['tracked'], envelope_seconds=waves['seconds'])
                for row in run['rows']:
                    amp_rows.append(dict(case=case, variant=name, x=row['x'], mode=row['mode'], N_archive=row['N_archive'],
                        N_inherited_on_independent_state=waves['inherited_at'](float(s.xi[np.argmin(abs(s.x-row['x']))])),
                        N_orr_sommerfeld=waves['N_at'](float(s.xi[np.argmin(abs(s.x-row['x']))]))))
                store[f'{label}/stability_xi'] = waves['waves']['xi']; store[f'{label}/stability_N'] = waves['waves']['N']
                store[f'{label}/stability_x'] = np.array([s.x_of_xi(v) for v in waves['waves']['xi']])
                store[f'{label}/stability_rate'] = waves['waves']['rate']; store[f'{label}/stability_omega'] = waves['waves']['omegas']
                store[f'{label}/stability_each'] = waves['waves']['each']
                store[f'{label}/inherited_xi'] = waves['xi']; store[f'{label}/inherited_N'] = waves['inherited']
                store[f'{label}/inherited_x'] = np.array([s.x_of_xi(v) for v in waves['xi']])
                say(f'{label}: N at x/c 0.0208: archive {record["N_endpoint_archive"]:.2f}, inherited law on independent state '
                    f'{n_inherited:.2f}, Orr-Sommerfeld {n_os:.2f} ({waves["seconds"]:.0f} s, lost waves {record["lost_waves"]})')
            screens.append(record)
        store[f'{case}/archive_x'] = s.x; store[f'{case}/archive_U'] = s.U; store[f'{case}/archive_turbulent'] = s.turbulent
        store[f'{case}/archive_cp'] = pressure.cp(s.q_nodes)
        store[f'{case}/archive_Hk'] = np.array([mf.get_Hk(u, param10)[0] for u in s.U.T])
    write_csv(out/'nodes.csv', node_rows); write_csv(out/'amplification.csv', amp_rows)
    keys = sorted({k for r in screens for k in r})
    write_csv(out/'screens.csv', [{k: r.get(k) for k in keys} for r in screens])
    np.savez_compressed(out/'profiles.npz', **store)
    primary = {r['case']: r for r in screens if r['variant'] == 'primary'}
    assert counts['marches'] <= 15 and counts['envelopes'] <= 8
    summary = dict(state='finished', part_a_gate_passed=gate, primary=primary, counts=dict(**counts, coupled_flow_solves=0,
        forces_computed=False), tolerances=dict(Cp=TOL_CP, theta=TOL_THETA, growth_ratio=RATIO),
        closure_screen={c: primary[c].get('closure_screen') for c in BUBBLES},
        growth_screen={c: primary[c].get('growth_screen') for c in BUBBLES}, wall_seconds=time.monotonic()-started)
    write(out/'summary.json', json.loads(json.dumps(summary, default=float)))
    say(json.dumps(dict(closure_screen=summary['closure_screen'], growth_screen=summary['growth_screen'], counts=summary['counts']),
                   indent=2))
    log.close()


if __name__ == '__main__':
    main()
