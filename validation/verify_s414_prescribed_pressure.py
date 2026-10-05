#!/usr/bin/env python3
"""Re-evaluate every saved prescribed-pressure station and plot admitted results."""
import csv
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_s414_prescribed_pressure import ROOT, OUT, digest, write, new_model, elements, mf, configure


def main():
    manifest = json.loads((OUT/'input_manifest.json').read_text())
    for name, expected in manifest['files'].items():
        assert digest(ROOT/name) == expected, 'Changed input: '+name
    model = new_model('fore', elements('nominal')[0][1], 1e6, .1, 9, (.25, 0), 1e-5)
    configure([model], 'corrected_af'); param = mf.build_param(model, 0)
    summary = json.loads((OUT/'summary.json').read_text())
    checks, rows = [], []
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex='col')
    colors = dict(native='#222222', model_taps='#7d3c98', experiment_taps='#1f77b4', experiment_delta='#d55e00')
    labels = dict(native='Native model', model_taps='Model at taps',
                  experiment_taps='Linear experimental taps', experiment_delta='Model + tap differences')
    for r in summary['cases']:
        name = r['case']+'_'+r['profile']; path = OUT/(name+'.npz')
        assert digest(path) == r['state_sha256']
        a = np.load(path); u, old, x, xi = (a[k] for k in ('U', 'original', 'x', 'xi'))
        assert np.array_equal(u[:, 0], old[:, 0])
        residuals = [np.max(abs(mf.residual_station(param, xi[j-1:j+1], u[:, j-1:j+1], np.zeros(2))[0]))
                     for j in range(1, len(x))]
        assert max(residuals) < 1e-9
        cp = mf.get_cp(u[3], param)[0]
        expected_cp = (1-r['reached_fraction'])*mf.get_cp(old[3], param)[0]+r['reached_fraction']*a['requested_cp']
        cp_error = float(np.max(abs(cp-expected_cp)))
        assert cp_error < 5e-14
        valid = r['complete'] and not r['conditional_threshold_crossing']
        assert (r['endpoint_N'] is not None) == valid
        if valid: assert abs(r['endpoint_N']-u[2, -1]) < 1e-13
        checks.append(dict(profile=name, maximum_residual=float(max(residuals)), Cp_error=cp_error,
            seed_exact=True, complete=r['complete'], admitted_endpoint=valid))
        column = 0 if r['case'] == 'cl0056' else 1
        color = colors[r['profile']]
        axes[0, column].plot(x, a['requested_cp'], color=color, label=labels[r['profile']])
        # Incomplete target BL profiles are deliberately absent. A threshold crossing
        # terminates the displayed laminar profile at the first N>=9 station.
        if r['complete']:
            hk = np.array([mf.get_Hk(v, param)[0] for v in u.T])
            end = int(np.where(u[2] >= 9)[0][0])+1 if r['conditional_threshold_crossing'] else len(x)
            axes[1, column].plot(x[:end], hk[:end], color=color)
            axes[2, column].plot(x[:end], u[2, :end], color=color)
            for j in range(end):
                rows.append(dict(case=r['case'], profile=r['profile'], x=float(x[j]),
                    Cp=float(cp[j]), Hk=float(hk[j]), N=float(u[2, j]), theta=float(u[0, j]),
                    delta_star=float(u[1, j]), threshold_reached=bool(u[2, j] >= 9)))
        if r['profile'] == 'native':
            mask = a['taps_x'] <= .03424
            axes[0, column].scatter(a['taps_x'][mask], a['taps_experiment'][mask], color='#1f77b4', marker='o', label='Experimental taps', zorder=5)
            axes[0, column].scatter(a['taps_x'][mask], a['taps_model'][mask], color='#7d3c98', marker='x', label='Model at tap centers', zorder=5)
        axes[0, column].set_title('CL '+('0.056' if column == 0 else '0.175'))
    for column in range(2):
        axes[0, column].invert_yaxis()
        axes[0, column].legend(fontsize=8, loc='best')
        axes[0, column].set_ylabel('Requested Cp (all profiles)')
        axes[1, column].set_ylabel('Hk (complete profiles only)')
        axes[2, column].set_ylabel('N (laminar diagnostic)')
        axes[2, column].axhline(9, color='gray', ls='--', label='Ncrit = 9')
        axes[2, column].set_xlabel('x / total chord')
        axes[2, column].set_xlim(.001, .022)
        for row in range(3): axes[row, column].grid(alpha=.2)
    fig.suptitle('Prescribed-pressure replay is inconclusive\nBoth model-tap and peak-preserving experimental targets fail continuation', fontsize=13)
    fig.text(.5, .012, 'Fixed model-dependent inlet; matched-lift pressure snapshots. No downstream transition or drag prediction.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .94))
    for ext in ('png', 'pdf', 'svg'): fig.savefig(OUT/('pressure_replay.'+ext), dpi=180)
    plt.close(fig)
    with (OUT/'admitted_stations.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    receipt = dict(complete=True, inputs_verified=len(manifest['files']), profiles_checked=8,
        requested_profiles_reached=sum(r['complete'] for r in summary['cases']),
        complete_below_threshold_profiles=sum(r['admitted_endpoint'] for r in checks),
        checks=checks, numerical_verification_pass=True, causal_diagnosis='inconclusive',
        source_sha256={str(Path(__file__).relative_to(ROOT)): digest(__file__)})
    receipt['output_sha256'] = {p.name: digest(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'verification.json'}
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('checks', 'output_sha256', 'source_sha256')}, indent=2))


if __name__ == '__main__': main()
