#!/usr/bin/env python3
"""Independent re-check of V5 tranche 1 from its saved outputs; writes the receipt.

Recomputes every stored closure residual, repeats a sample of finite-difference
marches on a refined cross-layer grid, re-derives the control decisions and the
declared outcomes from the march table, and checks hashes, counts and the
preserved first execution. No coupled solve, no forces.
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
import run_s414_laminar_replay as replay
import s414_laminar_fd as fd
from run_s414_laminar_replay import ROOT, OUT, CASES, CROSS, FAMILIES, UPSTREAM, TOL_X, TOL_THETA, NU, digest, write, mf

BRANCH = ROOT/'validation/results/s414_direct_branch_diagnostic_20261004'
FIRST = OUT/'initial_run_box_scheme'
FIGURES = ('two_branches', 'measured_taps')
SAMPLE = ('native_cl0056', 'native_cl0175', 'native_cl028', 'exp_cl0056_linear_hold', 'exp_cl0056_additive_hold',
          'exp_cl0175_ratio_ramp', 'cross_cl0175_to_cl0056_monotone_hold', 'self_cl028_monotone')


def table(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def number(text):
    return None if text in ('', 'None') else float(text)


def screen(result, truth):
    same = result[0] == truth[0] and result[0] in ('completed', 'separated')
    shift = result[1]-truth[1] if same and result[0] == 'separated' else None
    theta = result[2]/truth[2]-1
    return bool(same and (shift is None or abs(shift) <= TOL_X) and abs(theta) <= TOL_THETA)


def main():
    manifest = json.loads((OUT/'input_manifest.json').read_text())
    for name, expected in manifest['files'].items():
        assert digest(ROOT/name) == expected, 'Changed input: '+name
    summary = json.loads((OUT/'summary.json').read_text())
    assert summary['state'] == 'finished' and summary['part_b_gate_passed']
    marches = {m['id']: m for m in table(OUT/'marches.csv')}
    z = np.load(OUT/'profiles.npz')
    _, param0 = replay.parameters(0.)

    # 1. Every stored closure station satisfies the inherited equations.
    worst = 0.
    for name in marches:
        xi, states = z[name+'/closure_xi'], z[name+'/closure_states']
        for j in range(1, len(xi)):
            r = mf.residual_station(param0, xi[j-1:j+1], states[:, j-1:j+1], np.zeros(2))[0]
            worst = max(worst, float(np.max(abs(r))))
    assert worst < 1e-9, worst

    # 2. Finite-difference records: resolved outer edge, no sawtooth, and a refined cross-layer grid.
    assert all(m['fd_edge_resolved'] == 'True' for m in marches.values())
    sawtooth = max(float(m['fd_oscillation']) for m in marches.values())
    assert sawtooth < 1e-2, sawtooth
    pressure = replay.Pressure(replay.parameters(replay.MACH)[1])
    surfaces = {case: replay.Surface(np.load(replay.STATES/directory/'fore_state.npz'), pressure)
                for case, (directory, _) in CASES.items()}
    taps = {}
    for row in table(replay.TAPS):
        if row['coordinate'] == 'table':
            taps.setdefault(row['case'], []).append((float(row['x']), float(row['experiment_Cp'])))
    taps = {case: np.array(rows[:replay.TAP_COUNT]) for case, rows in taps.items()}

    def rebuilt(name):
        """The edge-velocity function of one march, rebuilt from its definition."""
        part = name.split('_')
        base = surfaces[part[1]]; tap_xi = base.xi_of_x(taps[part[1]][:, 0])
        if part[0] == 'native':
            return base.q
        if part[0] == 'self':
            return replay.reconstruction(base, tap_xi, base.cp(tap_xi), part[2], 'hold')
        if part[0] == 'cross':
            target = surfaces[part[3]]
            return replay.reconstruction(base, tap_xi, target.cp(target.xi_of_x(taps[part[1]][:, 0])), part[4], part[5])
        return replay.reconstruction(base, tap_xi, taps[part[1]][:, 1], part[2], part[3])
    refined = []
    for name in SAMPLE:
        nodes = z[name+'/nodes']; function = rebuilt(name)
        assert np.max(abs(function(nodes)-z[name+'/speed'])) < 1e-13, name
        base = fd.march(nodes, lambda s: float(function(s)), NU, m_start=1.)
        fine = fd.march(nodes, lambda s: float(function(s)), NU, m_start=1., eta=fd.eta_grid(14., 1121))
        assert base['status'] == fine['status'] == marches[name]['fd_status'], name
        assert np.array_equal(base['xi'], z[name+'/fd_xi']) and np.array_equal(base['theta'], z[name+'/fd_theta']), name
        shift = None if base['xi_separation'] is None else float(fine['xi_separation']-base['xi_separation'])
        theta = float(fine['theta'][-1]/base['theta'][-1]-1) if base['status'] == 'completed' else None
        assert shift is None or abs(shift) < 1e-4, (name, shift)
        assert theta is None or abs(theta) < 1e-3, (name, theta)
        refined.append(dict(march=name, status=fine['status'], xi_separation_shift=shift, theta_relative=theta))

    # 3. Control decisions and outcomes re-derived from the march table alone.
    def response(name, method):
        m = marches[name]
        return m[f'{method}_status'], number(m[f'{method}_x_separation']), float(m[f'{method}_theta_reference'])
    controls = table(OUT/'controls.csv')
    decided = {}
    for c in controls:
        if c['trivial'] == 'True':
            decided[(c['base'], c['target'], c['family'], c['upstream'], c['method'])] = True
            continue
        name = f"self_{c['base']}_{c['family']}" if c['kind'] == 'self' else \
            f"cross_{c['base']}_to_{c['target']}_{c['family']}_{c['upstream']}"
        passed = screen(response(name, c['method']), response('native_'+c['target'], c['method']))
        assert passed == (c['passed'] == 'True'), (name, c['method'])
        decided[(c['base'], c['target'], c['family'], c['upstream'], c['method'])] = passed
    assert len(controls) == summary['controls']['total'] and sum(decided.values()) == summary['controls']['passed']
    outcomes = {}
    for case in CASES:
        admitted = []
        for family in FAMILIES:
            for upstream in UPSTREAM:
                keys = [k for k in decided if k[0] == case and k[2] == family and k[3] in (upstream, 'all') and k[4] == 'fd']
                if all(decided[k] for k in keys):
                    admitted.append(f'exp_{case}_{family}_{upstream}')
        agree = all(marches[n]['fd_status'] == marches[n]['closure_status'] for n in admitted)
        states = {marches[n]['fd_status'] for n in admitted}
        outcome = ('not bounded' if not admitted or not agree or len(states) != 1 else
                   'separation retained' if states == {'separated'} else 'separation removed')
        assert outcome == summary['outcomes'][case]['outcome'] and len(admitted) == summary['outcomes'][case]['admitted'], case
        outcomes[case] = dict(outcome=outcome, admitted=admitted)
    experiment = table(OUT/'experiment.csv')
    for e in experiment:
        name = f"exp_{e['case']}_{e['family']}_{e['upstream']}"+('' if e['coordinate'] == 'table' else '_plotted')
        assert e['fd_status'] == marches[name]['fd_status'] and e['closure_status'] == marches[name]['closure_status']

    # 4. Limits, the Part A classification and the branch diagnostic.
    counts = summary['counts']
    assert counts['fd_marches'] == len(marches) <= 120 and counts['closure_marches'] <= 130
    assert counts['coupled_flow_solves'] == 0 and not counts['forces_computed']
    assert summary['part_a']['failures_classified_as_direct_singularity']
    branch = json.loads((BRANCH/'summary.json').read_text())
    for name, expected in branch['output_sha256'].items():
        assert digest(BRANCH/name) == expected, 'Branch diagnostic output changed: '+name
    roots = [r for r in table(BRANCH/'station_roots.csv') if r['two_distinct_roots'] == 'True']
    assert all(float(r['residual_from_archive']) < 1e-9 and float(r['residual_from_upstream']) < 1e-9 for r in roots)
    two_roots = [dict(case=r['case'], x=float(r['x']), Hk_archive=float(r['root_from_archive_Hk']),
                      Hk_other=float(r['root_from_upstream_Hk'])) for r in roots]

    # 5. The first execution is intact and its outcomes are recorded for comparison.
    first_manifest = json.loads((FIRST/'input_manifest.json').read_text())['files']
    for name in ('run_s414_laminar_replay.py', 's414_laminar_fd.py'):
        assert digest(FIRST/'source_snapshot'/name) == first_manifest['validation/'+name], 'First-execution source changed'
    first = json.loads((FIRST/'summary.json').read_text())
    first_marches = {m['id']: m for m in table(FIRST/'marches.csv')}
    changed = [name for name in marches if marches[name]['fd_status'] != first_marches[name]['fd_status']]
    moved = [abs(number(marches[n]['fd_x_separation'])-number(first_marches[n]['fd_x_separation'])) for n in marches
             if marches[n]['fd_status'] == first_marches[n]['fd_status'] == 'separated']

    figures = {f'{name}.{ext}': digest(OUT/f'{name}.{ext}') for name in FIGURES for ext in ('png', 'pdf', 'svg')}
    receipt = dict(complete=True, inputs_verified=len(manifest['files']), closure_max_residual=worst,
        fd_max_sawtooth_index=sawtooth, fd_refined_grid_sample=refined, controls_rederived=len(controls),
        controls_passed=summary['controls']['passed'], fd_nontrivial_controls_passed=summary['controls']['fd_nontrivial_passed'],
        fd_nontrivial_controls=summary['controls']['fd_nontrivial'], outcomes=outcomes, counts=counts,
        stations_with_two_roots=two_roots,
        first_execution=dict(outcomes={k: v['outcome'] for k, v in first['outcomes'].items()},
            fd_nontrivial_controls_passed=first['controls']['fd_nontrivial_passed'], fd_topology_changed=changed,
            largest_separation_shift=max(moved), sources_match_manifest=True),
        numerical_verification_pass=True,
        physical_interpretation='Outcomes are those of the frozen rule. An admitted outcome remains conditional on the '
        'between-tap pressure shape of the admitted reconstruction; no drag, transition or closure conclusion follows.',
        figures_reviewed=figures,
        source_sha256={name: digest(ROOT/name) for name in ('validation/verify_s414_laminar_replay.py',
            'validation/report_s414_laminar_replay.py', 'validation/diagnose_s414_direct_branches.py')})
    receipt['output_sha256'] = {p.name: digest(p) for p in sorted(OUT.iterdir())
                                if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md', 'tests.log')}
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('output_sha256', 'figures_reviewed', 'source_sha256')}, indent=1))


if __name__ == '__main__':
    main()
