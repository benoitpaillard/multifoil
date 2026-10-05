#!/usr/bin/env python3
"""Independent re-check of V5 tranche 2 from its saved outputs; writes the receipt.

Repeats the primary marches and compares them with the stored arrays, repeats
them on a doubled cross-layer grid, recomputes stored growth rates from fresh
temporal scans on a finer collocation grid, evaluates the momentum-integral
balance of the stored boundary layers, and re-derives both declared screens
from the node tables. No coupled solve, no forces.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import csv
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_s414_bubble_stability as study
import s414_bubble_fd as bf
import s414_stability as stability
import s414_amplification as amp
from run_s414_bubble_stability import ROOT, OUT, CASES, STATES, NU, MACH, BUBBLES, TOL_CP, TOL_THETA, RATIO, replay, digest, write

FIGURES = ('bubble_closure', 'bubble_amplification', 'similar_profiles')


def table(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def number(text):
    return None if text in ('', 'None') else float(text)


def main():
    manifest = json.loads((OUT/'input_manifest.json').read_text())
    for name, expected in manifest['files'].items():
        assert digest(ROOT/name) == expected, 'Changed input: '+name
    summary = json.loads((OUT/'summary.json').read_text())
    assert summary['state'] == 'finished' and summary['part_a_gate_passed']
    qualification = table(OUT/'qualification.csv')
    assert all(r['passed'] == 'True' for r in qualification if r['gate'] == 'True')
    z = np.load(OUT/'profiles.npz')
    nodes, growth, screens = table(OUT/'nodes.csv'), table(OUT/'amplification.csv'), table(OUT/'screens.csv')
    _, param10 = replay.parameters(MACH)
    pressure = replay.Pressure(param10)
    surfaces = {case: replay.Surface(np.load(STATES/directory/'fore_state.npz'), pressure) for case, (directory, _) in CASES.items()}

    # 1. Repeat the primary marches; compare with the stored arrays and with a doubled cross-layer grid.
    repeated, balance, rates = [], [], []
    fine_eta = bf.eta_grid(30., 2401)
    for case, s in surfaces.items():
        run = study.bubble(s, param10)
        result = run['result']; label = f'{case}/primary'
        assert result['status'] == 'completed'
        assert np.array_equal(result['xi'], z[label+'/xi']) and np.array_equal(result['ue'], z[label+'/ue']), label
        target = study.Displacement(s.xi[run['switch']:run['end']+1], s.U[1, run['switch']:run['end']+1])
        fine = bf.march(run['nodes'], s.q, NU, switch=float(s.xi[run['switch']]), displacement=target, eta=fine_eta)
        assert fine['status'] == 'completed'
        common = np.array([np.min(abs(run['nodes']-x)) < 1e-12 for x in result['xi']])
        common_fine = np.array([np.min(abs(run['nodes']-x)) < 1e-12 for x in fine['xi']])
        speed = float(np.max(abs(fine['ue'][common_fine]/result['ue'][common]-1)))
        theta = float(np.max(abs(fine['theta'][common_fine]/result['theta'][common]-1)))
        assert speed < 1e-3 and theta < 2e-3, (case, speed, theta)
        repeated.append(dict(case=case, stored_arrays_identical=True, edge_speed_change_on_doubled_grid=speed,
                             theta_change_on_doubled_grid=theta))
        # Momentum-integral balance of the stored layer (FLARE drops part of it in reversed flow).
        xi, th, h, ue = (result[k][common] for k in ('xi', 'theta', 'H', 'ue'))
        cf = result['cf'][common]
        left = np.log(th[1:]/th[:-1])
        pressure_term = (2+.5*(h[1:]+h[:-1]))*np.log(ue[1:]/ue[:-1])
        friction = .25*(cf[1:]+cf[:-1])*(xi[1:]-xi[:-1])/(.5*(th[1:]+th[:-1]))
        residual = np.abs(left+pressure_term-friction)/np.maximum.reduce([np.abs(left), np.abs(pressure_term), np.abs(friction)])
        attached = result['reverse'][common][1:] >= 0
        # Summed over the reversed stations, the unbalanced part is the change of ln(theta) that FLARE leaves out.
        balance.append(dict(case=case, attached_median=float(np.median(residual[attached])),
                            attached_largest=float(np.max(residual[attached])),
                            reversed_largest=float(np.max(residual[~attached])) if np.any(~attached) else None,
                            reversed_cumulative_log_theta=float(np.sum((left+pressure_term-friction)[~attached]))))
        assert np.median(residual[attached]) < .02, (case, np.median(residual[attached]))
        # 2. Stored growth rates recomputed from fresh scans on a finer collocation grid.
        grid = stability.Grid(200, 60., 3.)
        stored_xi, stored_rate, omegas = z[label+'/stability_xi'], z[label+'/stability_rate'], z[label+'/stability_omega']
        index = np.where(common)[0]
        for position in (len(stored_xi)-1, int(.75*len(stored_xi)), len(stored_xi)//2):
            if not np.any(np.isfinite(stored_rate[:, position])):
                continue
            k = int(np.nanargmax(stored_rate[:, position]))
            j = index[int(np.argmin(abs(result['xi'][index]-stored_xi[position])))]
            p = result['profiles'][j]
            item = amp.station(p['xi'], result['delta_star'][j], result['ue'][j], NU, p['eta'], p['ratio'], p['curvature'],
                               p['scaled_displacement'], grid)
            reduced = omegas[k]*item['delta_star']/item['ue']
            start, shape = stability.seed(reduced, stability.scan(item['re'], item['U'], item['Upp'], grid, amp.WAVENUMBERS))
            alpha, _, ok = stability.spatial(reduced, item['re'], item['U'], item['Upp'], grid, start, shape)
            fresh = -alpha.imag/item['delta_star']
            assert ok and abs(fresh/stored_rate[k, position]-1) < 5e-3, (case, position, fresh, stored_rate[k, position])
            rates.append(dict(case=case, x=float(s.x_of_xi(stored_xi[position])), stored=float(stored_rate[k, position]),
                              recomputed=float(fresh)))

    # 3. Screens re-derived from the node tables.
    rederived = {}
    for case in CASES:
        rows = [r for r in nodes if r['case'] == case and r['variant'] == 'primary']
        inverse = [r for r in rows if r['mode'] == 'inverse']
        cp = max(abs(float(r['Cp_independent'])-float(r['Cp_archive'])) for r in inverse)
        theta = float(rows[-1]['theta_independent'])/float(rows[-1]['theta_archive'])-1
        last = [r for r in growth if r['case'] == case and r['variant'] == 'primary'][-1]
        ratio = float(last['N_orr_sommerfeld'])/float(last['N_inherited_on_independent_state'])
        stored = next(r for r in screens if r['case'] == case and r['variant'] == 'primary')
        assert abs(cp-float(stored['largest_Cp_difference'])) < 1e-12 and abs(theta-float(stored['theta_endpoint_relative'])) < 1e-12
        assert abs(ratio-float(stored['ratio_orr_sommerfeld_over_inherited'])) < 1e-9
        rederived[case] = dict(largest_Cp_difference=cp, theta_endpoint_relative=theta, growth_ratio=ratio,
            closure_screen=bool(cp <= TOL_CP and abs(theta) <= TOL_THETA) if case in BUBBLES else None,
            growth_screen=bool(RATIO[0] <= ratio <= RATIO[1]) if case in BUBBLES else None)
        if case in BUBBLES:
            assert rederived[case]['closure_screen'] == summary['closure_screen'][case]
            assert rederived[case]['growth_screen'] == summary['growth_screen'][case]
    counts = summary['counts']
    assert counts['marches'] <= 15 and counts['envelopes'] <= 8 and counts['coupled_flow_solves'] == 0
    figures = {f'{name}.{ext}': digest(OUT/f'{name}.{ext}') for name in FIGURES for ext in ('png', 'pdf', 'svg')}
    receipt = dict(complete=True, inputs_verified=len(manifest['files']), gate_controls_passed=sum(r['gate'] == 'True' for r in qualification),
        primary_marches_repeated=repeated, momentum_integral_balance=balance, growth_rates_recomputed=rates, screens=rederived,
        counts=counts, numerical_verification_pass=True,
        physical_interpretation='Conditional on the archived displacement thickness, parallel-flow stability theory and, in '
        'reversed flow, the FLARE approximation. No transition location, drag or closure change follows.',
        figures_reviewed=figures, source_sha256={name: digest(ROOT/name) for name in (
            'validation/verify_s414_bubble_stability.py', 'validation/report_s414_bubble_stability.py')})
    receipt['output_sha256'] = {p.name: digest(p) for p in sorted(OUT.iterdir())
                                if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md', 'tests.log')}
    write(OUT/'verification.json', json.loads(json.dumps(receipt, default=float)))
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('output_sha256', 'figures_reviewed', 'source_sha256')},
                     indent=1, default=float))


if __name__ == '__main__':
    main()
