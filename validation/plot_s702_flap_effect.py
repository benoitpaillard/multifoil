#!/usr/bin/env python3
"""Plot audited S702 flap-angle comparisons; this script does not run a solver.

In the project:
    python validation/plot_s702_flap_effect.py
In the standalone bundle (requires numpy and matplotlib):
    python plot_s702_flap_effect.py --data s702_flap_drag_data.csv --out figures
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


SOURCE = 'https://ntrs.nasa.gov/citations/20220014964'
TARGETS = (0.5, 0.8, 1.2)
CASES = (
    ('calibration', 1000000, 0.1, (-2, 0, 5, 10), 'points.csv'),
    ('reserved', 700000, 0.07, (0, 5, 10), 'holdout_after_retries.csv'),
)
SERIES = {
    'experiment': dict(label='Experiment', color='#212529', marker='o', ls='-', lw=1.7),
    'mses_n9': dict(label=r'Published MSES, $N_{crit}=9$', color='#8554A3', marker='s', ls='--', lw=1.8),
    'multifoil_n9': dict(label=r'MultiFoil, $N_{crit}=9$', color='#277DA8', marker='^', ls=':', lw=2),
    'multifoil_n6p25': dict(label=r'MultiFoil, $N_{crit}=6.25$ (calibrated)', color='#D26919', marker='D', ls='-', lw=2.1),
}
FIGURES = {
    ('calibration', -2): '23(c)', ('calibration', 0): '24(a)',
    ('calibration', 5): '25(b)', ('calibration', 10): '26(c)',
    ('reserved', 0): '23(b)', ('reserved', 5): '25(a)',
    ('reserved', 10): '26(b)',
}


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def from_analysis(analysis):
    """Keep the frozen transition selection, matching references and retry flags."""
    frozen = json.loads((analysis / 'frozen_selection.json').read_text())
    if frozen['selected_ncrit'] != 6.25:
        raise ValueError('Expected the frozen Ncrit=6.25 selection')
    if sha(analysis / 'points.csv') != frozen['calibration_points_sha256']:
        raise ValueError('Calibration points differ from the frozen selection')
    rows = []
    for dataset, re, mach, flaps, filename in CASES:
        source = read_csv(analysis / filename)
        selected = {}
        for row in source:
            n = float(row['ncrit'])
            if n not in (9, 6.25):
                continue
            key = (int(row['flap']), float(row['target_cl']), n)
            if key in selected:
                raise ValueError(f'Duplicate point: {dataset}, {key}')
            if row['feasible'] != 'True' or row['reference_conditions_matched'] != 'True':
                raise ValueError(f'Infeasible or unmatched point: {dataset}, {key}')
            if abs(float(row['CL']) - key[1]) > 2e-5:
                raise ValueError(f'Unmatched lift: {dataset}, {key}')
            selected[key] = row
        expected = {(f, cl, n) for f in flaps for cl in TARGETS for n in (9, 6.25)}
        if set(selected) != expected:
            raise ValueError(f'Incomplete or unexpected {dataset} coverage')
        for cl in TARGETS:
            for flap in flaps:
                original, calibrated = (selected[(flap, cl, n)] for n in (9, 6.25))
                for key in ('CD_experiment', 'CD_mses_n9'):
                    if float(original[key]) != float(calibrated[key]):
                        raise ValueError('Reference differs between transition settings')
                for series in SERIES:
                    point = calibrated if series == 'multifoil_n6p25' else original
                    cd_key = {'experiment': 'CD_experiment', 'mses_n9': 'CD_mses_n9'}.get(series, 'CD')
                    numerical = series.startswith('multifoil')
                    rows.append(dict(
                        dataset=dataset, reynolds=re, mach=mach, target_cl=cl,
                        flap_deg=flap, series=series,
                        ncrit='' if series == 'experiment' else (6.25 if series == 'multifoil_n6p25' else 9),
                        CD=float(point[cd_key]),
                        numerical_retry=bool(point.get('local_retry_of')) if numerical else False,
                        achieved_CL=float(point['CL']) if numerical else '',
                        result_path=point['result_path'] if numerical else '',
                        retry_of=point.get('local_retry_of', '') if numerical else '',
                        source_table=filename,
                        reference_figure=FIGURES[(dataset, flap)],
                        reference_url=SOURCE,
                    ))
    return validate(rows)


def validate(rows):
    """Reject incomplete groups, invalid numbers or changed convergence labels."""
    for row in rows:
        for field in ('reynolds', 'mach', 'target_cl', 'flap_deg', 'CD'):
            row[field] = float(row[field])
            if not np.isfinite(row[field]):
                raise ValueError(f'Nonfinite {field}')
        row['numerical_retry'] = str(row['numerical_retry']).lower() == 'true'
        if row['CD'] <= 0:
            raise ValueError('Nonpositive drag')
        if row['series'].startswith('multifoil'):
            if abs(float(row['achieved_CL']) - row['target_cl']) > 2e-5:
                raise ValueError('Lift tolerance is not met')
        expected_n = {'experiment': '', 'mses_n9': 9., 'multifoil_n9': 9., 'multifoil_n6p25': 6.25}
        n = '' if row['ncrit'] == '' else float(row['ncrit'])
        if n != expected_n[row['series']]:
            raise ValueError('Transition label mismatch')
    lookup = {(r['dataset'], r['target_cl'], r['flap_deg'], r['series']): r for r in rows}
    expected = {(d, cl, f, s) for d, _, _, flaps, _ in CASES for cl in TARGETS for f in flaps for s in SERIES}
    if set(lookup) != expected or len(rows) != len(expected):
        raise ValueError('Expected 84 unique source/operating-point rows')
    for dataset, re, mach, _, _ in CASES:
        if any(r['reynolds'] != re or r['mach'] != mach for r in rows if r['dataset'] == dataset):
            raise ValueError('Reynolds or Mach mismatch')
    retry_keys = {key for key, row in lookup.items() if row['numerical_retry']}
    if retry_keys != {('reserved', cl, 10, 'multifoil_n6p25') for cl in (0.5, 0.8)}:
        raise ValueError('Numerical retry provenance differs from the audited study')
    for key, row in lookup.items():
        baseline = lookup[(key[0], key[1], 0, key[3])]['CD']
        row['CD_counts'] = row['CD'] * 1e4
        row['delta_CD_from_zero_counts'] = (row['CD'] - baseline) * 1e4
    return rows


def make_figure(rows, out, delta=False):
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.labelcolor': '#28333D', 'text.color': '#202B35',
        'axes.edgecolor': '#A9B2BA', 'xtick.color': '#46515B',
        'ytick.color': '#46515B', 'svg.fonttype': 'none',
    })
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5), sharex=True, sharey='col')
    fig.subplots_adjust(left=.076, right=.975, top=.756, bottom=.19, hspace=.62, wspace=.27)
    title = 'S702 | Drag change from 0° flap' if delta else 'S702 | Drag versus flap angle'
    fig.text(.076, .966, title, fontsize=20, weight='bold', va='top')
    subtitle = (r'$\Delta C_D=C_D(\delta)-C_D(0°)$ at equal lift; negative means lower drag.' if delta
                else r'Equal-lift comparisons at $C_L=0.5$, $0.8$ and $1.2$; one drag count = $10^{-4}$.')
    fig.text(.076, .918, subtitle, fontsize=11.5, va='top')
    handles = {}
    ykey = 'delta_CD_from_zero_counts' if delta else 'CD_counts'
    row_y = (.817, .457)
    for i, (dataset, re, mach, flaps, _) in enumerate(CASES):
        label = 'Calibration conditions' if dataset == 'calibration' else 'Reserved conditions · same frozen calibration'
        fig.text(.076, row_y[i], f'Re = {re:,}  |  M = {mach:g}     {label}', fontsize=11.5, weight='bold')
        for j, cl in enumerate(TARGETS):
            ax = axes[i, j]
            for series in ('multifoil_n9', 'mses_n9', 'multifoil_n6p25', 'experiment'):
                points = sorted((r for r in rows if r['dataset'] == dataset and r['target_cl'] == cl and r['series'] == series), key=lambda r: r['flap_deg'])
                style = SERIES[series]
                handle, = ax.plot([p['flap_deg'] for p in points], [p[ykey] for p in points],
                                  **style, ms=5.5, markerfacecolor='white' if series == 'mses_n9' else style['color'])
                handles[series] = handle
                for p in points:
                    if p['numerical_retry']:
                        ax.plot(p['flap_deg'], p[ykey], 'o', ms=13, mfc='none', mec=style['color'], mew=1.2, zorder=8)
            if delta:
                ax.axhline(0, color='#86929B', lw=1, zorder=0)
            ax.set_title(rf'$C_L={cl:g}$', fontsize=12, pad=10)
            ax.set_xlim(-2.8, 10.8)
            ax.set_xticks((-2, 0, 5, 10))
            ax.tick_params(labelbottom=True)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            ax.grid(axis='y', color='#DCE1E5', lw=.7)
            ax.set_axisbelow(True)
            ax.set_xlabel('Flap angle δ (degrees)', labelpad=6)
            ax.set_ylabel(r'$\Delta C_D$ (drag counts)' if delta else r'$C_D$ (drag counts)')
    for j, cl in enumerate(TARGETS):
        values = [r[ykey] for r in rows if r['target_cl'] == cl]
        if delta:
            values.append(0)
        pad = max(2, .12 * np.ptp(values))
        axes[0, j].set_ylim(min(values) - pad, max(values) + pad)
    fig.legend([handles[s] for s in SERIES], [SERIES[s]['label'] for s in SERIES],
               loc='upper left', bbox_to_anchor=(.068, .886), ncol=4,
               frameon=False, fontsize=10, handlelength=2.8, columnspacing=1.6)
    fig.text(.076, .097, 'Open rings: two calibrated MultiFoil points recovered with alternative numerical initialization.', fontsize=9.3)
    fig.text(.076, .071, 'MSES and experiment: digitized NASA/CR–20220014964 polars, interpolated in lift. MSES uses published Ncrit = 9.', fontsize=9.3)
    fig.text(.076, .045, 'Lines join sampled flap settings only; y-scales match within each lift column. No reserved data at −2°. MultiFoil geometry is reconstructed.', fontsize=9.1)
    name = 's702_flap_drag_change' if delta else 's702_flap_drag'
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=190, facecolor='white')
    plt.close(fig)


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, default=here / 'results/s702_ncrit_analysis')
    parser.add_argument('--data', type=Path, help='Use the portable 84-row plot-source CSV')
    parser.add_argument('--out', type=Path, default=None)
    args = parser.parse_args()
    data = args.data or (here / 's702_flap_drag_data.csv' if (here / 's702_flap_drag_data.csv').exists() else None)
    out = args.out or (here / 'figures' if data else here / 'results/s702_flap_effect')
    out.mkdir(parents=True, exist_ok=True)
    rows = validate(read_csv(data)) if data else from_analysis(args.analysis)
    write_csv(out / 's702_flap_drag_data.csv', rows)
    lookup = {(r['dataset'], r['target_cl'], r['flap_deg'], r['series']): r for r in rows}
    comparisons = []
    for dataset, re, mach, _, _ in CASES:
        for cl in TARGETS:
            row = dict(dataset=dataset, reynolds=re, mach=mach, target_cl=cl)
            for series in SERIES:
                row[series + '_delta_10_minus_5_counts'] = (lookup[(dataset, cl, 10, series)]['CD'] - lookup[(dataset, cl, 5, series)]['CD']) * 1e4
            comparisons.append(row)
    write_csv(out / 's702_flap_5_to_10.csv', comparisons)
    for delta in (False, True):
        make_figure(rows, out, delta)
    inputs = [data] if data else [args.analysis / c[-1] for c in CASES] + [args.analysis / 'frozen_selection.json']
    provenance = dict(
        inputs={str(p.name): sha(p) for p in inputs},
        script_sha256=sha(__file__), plot_rows=len(rows), numerical_rows=42,
        calibration_conditions=dict(reynolds=1000000, mach=0.1),
        reserved_conditions=dict(reynolds=700000, mach=0.07),
        recovered_numerical_points=[dict(dataset=r['dataset'], CL=r['target_cl'], flap=r['flap_deg'], result=r['result_path'], retry_of=r['retry_of']) for r in rows if r['numerical_retry']],
        reference=SOURCE, reference_report='NASA/CR–20220014964',
        note='Published MSES Ncrit=9; no new MSES execution. MultiFoil Ncrit=6.25 is effective calibration; geometry reconstructed. No CFD is rerun by this script.',
    )
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(dict(output=str(out.resolve()), rows=len(rows), comparisons=comparisons), indent=2))


if __name__ == '__main__':
    main()
