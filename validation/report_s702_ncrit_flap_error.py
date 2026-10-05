#!/usr/bin/env python3
"""Reanalyse archived S702 Ncrit sensitivity of the 10-minus-5 flap drag error.

No flow solves, calibration changes, or substitutions of finer-grid results.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = 'validation/results/s702_ncrit_analysis/'
MECH = 'validation/results/s702_drag_mechanism/'
MANIFESTS = {
    'handover/archive_manifest.json': '9102d74534effcc51c30872e0e1c5a2be1e5bc6e0d7fc5cd6d6c97c715fb4fc8',
    'handover/drag_mechanism_verification.json': '2c0f0041bed89c5261d27757da405b710772b68e80f571eedc05abc763fbbe66',
}
CLS = (.5, .8, 1.2)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Inputs:
    def __init__(self):
        self.expected, self.consumed = {}, {}
        for name, digest in MANIFESTS.items():
            require(sha(ROOT / name) == digest, 'Historical manifest changed: ' + name)
            self.consumed[name] = digest
            for row in json.loads((ROOT / name).read_text())['files']:
                self.expected[row['path']] = row['sha256']

    def read(self, name, kind='csv'):
        digest = sha(ROOT / name)
        require(name in self.expected and digest == self.expected[name], 'Archived input changed: ' + name)
        self.consumed[name] = digest
        with (ROOT / name).open() as stream:
            return list(csv.DictReader(stream)) if kind == 'csv' else json.load(stream)


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def close(a, b):
    return abs(float(a) - float(b)) < 1e-10


def read_points(inputs, name, dataset):
    rows = inputs.read(BASE + name)
    for row in rows:
        raw = inputs.read(row['result_path'], 'json')
        row.update(dataset=dataset, feasible=row['feasible'] == 'True')
        for key in ('ncrit', 'flap', 'target_cl', 'CD_experiment', 'CD_mses_n9'):
            row[key] = float(row[key])
        require(row['reference_conditions_matched'] == 'True', 'Unmatched reference conditions')
        require(row['feasible'] == raw['feasible'], 'Coverage differs from saved target')
        if row['feasible']:
            row['CD'], row['CL'] = float(row['CD']), float(row['CL'])
            require(np.isfinite(row['CD']) and row['CD'] > 0, 'Invalid converged drag')
            require(close(row['CD'], raw['CD']) and close(row['CL'], raw['CL']), 'Saved force differs')
            require(abs(row['CL'] - row['target_cl']) <= 2e-5, 'Lift mismatch')
        else:
            row['CD'], row['CL'] = None, None
    return rows


def pair_rows(points):
    rows = []
    for dataset in ('calibration', 'reserved'):
        group = [r for r in points if r['dataset'] == dataset]
        for ncrit in sorted({r['ncrit'] for r in group}):
            candidate = [r for r in group if r['ncrit'] == ncrit]
            expected = 12 if dataset == 'calibration' else 9
            require(len(candidate) == expected, 'Missing candidate entries')
            coverage = sum(r['feasible'] for r in candidate)
            for cl in CLS:
                p = {r['flap']: r for r in candidate if r['target_cl'] == cl and r['flap'] in (5, 10)}
                require(set(p) == {5, 10}, 'Missing flap pair')
                feasible = p[5]['feasible'] and p[10]['feasible']
                measured = (p[10]['CD_experiment'] - p[5]['CD_experiment']) * 1e4
                mses = (p[10]['CD_mses_n9'] - p[5]['CD_mses_n9']) * 1e4
                delta = (p[10]['CD'] - p[5]['CD']) * 1e4 if feasible else None
                rows.append(dict(dataset=dataset, reynolds=1e6 if dataset == 'calibration' else 7e5,
                    mach=.1 if dataset == 'calibration' else .07, target_cl=cl, ncrit=ncrit,
                    candidate_converged=coverage, candidate_required=expected,
                    pair_feasible=feasible, delta_counts=delta, experiment_delta_counts=measured,
                    error_counts=delta-measured if feasible else None, mses_n9_delta_counts=mses,
                    mses_n9_error_counts=mses-measured,
                    result_5=p[5]['result_path'], result_10=p[10]['result_path'],
                    retry_5=p[5].get('local_retry_of', ''), retry_10=p[10].get('local_retry_of', '')))
    return rows


def select(rows, **criteria):
    found = [r for r in rows if all(r[k] == v for k, v in criteria.items())]
    require(len(found) == 1, 'Selection is not unique: ' + str(criteria))
    return found[0]


def plot(out, points, pairs):
    colors = {5: '#3278a5', 10: '#dc862b'}
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), layout='constrained')
    for j, dataset in enumerate(('calibration', 'reserved')):
        rows = [r for r in points if r['dataset'] == dataset and r['target_cl'] == .8]
        ns = sorted({r['ncrit'] for r in rows})
        for flap in (5, 10):
            p = sorted((r for r in rows if r['flap'] == flap), key=lambda r: r['ncrit'])
            axes[0, j].plot(ns, [r['CD']*1e4 if r['feasible'] else np.nan for r in p],
                            'o-' if j == 0 else 'o--', color=colors[flap], ms=4, label=f'MultiFoil {flap}°')
            axes[0, j].axhline(p[0]['CD_experiment']*1e4, color=colors[flap], ls=':', label=f'Experiment {flap}°')
            axes[0, j].plot(9, p[0]['CD_mses_n9']*1e4, 'D', color=colors[flap], mec='black', ms=7)
        axes[0, j].plot([], [], 'D', color='white', mec='black', label='Published MSES, N=9 only')
        q = sorted((r for r in pairs if r['dataset'] == dataset and r['target_cl'] == .8), key=lambda r: r['ncrit'])
        axes[1, j].plot(ns, [r['delta_counts'] if r['pair_feasible'] else np.nan for r in q],
                        'o-' if j == 0 else 'o--', color='#3278a5', ms=4, label='MultiFoil')
        axes[1, j].axhline(q[0]['experiment_delta_counts'], color='#9f353b', ls='--', label='Experiment')
        axes[1, j].plot(9, q[0]['mses_n9_delta_counts'], 'D', color='#dc862b', mec='black', ms=7, label='Published MSES, N=9')
        for n in (6.25, 9.):
            p = select(q, ncrit=n)
            axes[1, j].annotate(f"{p['delta_counts']:+.2f}", (n, p['delta_counts']), xytext=(0, -17 if n == 6.25 else 9),
                                textcoords='offset points', ha='center', fontsize=9)
        axes[1, j].axhline(0, color='#999', lw=.8)
        axes[1, j].text(.025, .045, f"Experiment: {q[0]['experiment_delta_counts']:+.2f} counts",
                        transform=axes[1, j].transAxes, color='#9f353b', fontsize=9)
        axes[0, j].set_title('Re 1,000,000 · M 0.10 · calibration' if j == 0 else 'Re 700,000 · M 0.07 · reserved condition')
        axes[0, j].set_ylabel('Absolute CD, drag counts')
        axes[1, j].set_ylabel('CD(10°) − CD(5°), drag counts')
        for ax in axes[:, j]:
            ax.set_xlabel('Ncrit, common to both elements')
            ax.axvline(6.25, color='#777', ls=':', lw=.8)
            ax.grid(alpha=.18)
            ax.legend(fontsize=8, loc='best')
        axes[1, j].set_ylim(-37, 22)
        if j == 1:
            for ax in axes[:, j]:
                ax.set_xticks([6.25, 9])
                ax.set_xlim(5.9, 9.45)
            axes[0, j].set_title('Re 700,000 · M 0.07 · reserved condition\nOnly two Ncrit values; dashed joins are guides.', fontsize=11)
            axes[0, j].legend(fontsize=8, loc='lower left', bbox_to_anchor=(.015, .12))
        else:
            for ax in axes[:, j]:
                ax.set_xlim(2.7, 12.3)
            axes[1, j].text(.98, .96, 'N=12 pair failed; no value assigned.', ha='right', va='top',
                            transform=axes[1, j].transAxes, fontsize=8)
    fig.suptitle('S702 · CL = 0.8 · Ncrit changes the flap-drag error, but leaves a residual\n'
                 'Original 260/180 panel protocol held fixed; MSES sensitivity is not available.', fontsize=12)
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(out / ('ncrit_flap_error.' + ext), dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1), layout='constrained')
    for ax, cl in zip(axes, CLS):
        q = [r for r in pairs if r['dataset'] == 'calibration' and r['target_cl'] == cl]
        ax.plot([r['ncrit'] for r in q], [r['error_counts'] if r['pair_feasible'] else np.nan for r in q], 'o-', ms=4)
        partial = [r for r in q if r['pair_feasible'] and r['candidate_converged'] < r['candidate_required']]
        ax.plot([r['ncrit'] for r in partial], [r['error_counts'] for r in partial], 'o', mfc='white', mec='#3278a5', ms=7)
        for r in q:
            if not r['pair_feasible']:
                ax.plot(r['ncrit'], .03, 'x', color='#9f353b', transform=ax.get_xaxis_transform())
        ax.plot(9, q[0]['mses_n9_error_counts'], 'D', color='#dc862b', mec='black', ms=7)
        ax.axhline(0, color='#9f353b', ls='--')
        ax.axvline(6.25, color='#777', ls=':')
        ax.set(title=f'CL = {cl:g}', xlabel='Ncrit', ylabel='Error in flap drag change, counts')
        ax.grid(alpha=.18)
    fig.suptitle('Re 1,000,000 · M 0.10 · error = predicted flap change − measured flap change\n'
                 'Blue: MultiFoil; diamond: published MSES N=9; open circles: incomplete global candidate; ×: failed pair.', fontsize=11)
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(out / ('ncrit_other_lifts.' + ext), dpi=180)
    plt.close(fig)


def report(out):
    inputs = Inputs()
    points = read_points(inputs, 'points.csv', 'calibration') + read_points(inputs, 'holdout_after_retries.csv', 'reserved')
    pairs = pair_rows(points)
    archived_pairs = inputs.read(BASE + 'pairwise.csv')
    calibration_scores = inputs.read(BASE + 'scores.csv')
    fitted_score = select(calibration_scores, ncrit='6.25')
    for ncrit in ('4.0', '5.0'):
        require(float(select(calibration_scores, ncrit=ncrit)['mape_percent']) > float(fitted_score['mape_percent']),
                'Low-Ncrit global calibration comparison differs')
    for row in pairs:
        if row['dataset'] != 'calibration' or not row['pair_feasible']:
            continue
        original = select(archived_pairs, ncrit=str(row['ncrit']), first_flap='5', second_flap='10', target_cl=str(row['target_cl']))
        require(close(row['delta_counts'], original['predicted_delta_counts']), 'Archived pair mismatch')
        require(close(row['error_counts'], original['delta_error_counts']), 'Archived pair error mismatch')
    comparison = inputs.read(MECH + 'comparison_points.csv')
    states = inputs.read(MECH + 'selected_states.csv')
    surfaces = inputs.read(MECH + 'surface_summary.csv')
    transition_sweep = inputs.read(BASE + 'transition_locations.csv')
    element_deltas = inputs.read(MECH + 'element_drag_deltas.csv')
    element_focus = select(element_deltas, dataset='reserved', target_cl='0.8', ncrit='6.25')
    require(round(float(element_focus['fore_delta_wake_counts']), 2) == -21.72, 'Fore contribution differs')
    require(round(float(element_focus['aft_delta_wake_counts']), 2) == 20.38, 'Aft contribution differs')
    for ncrit, expected_x in (('4.0', .4085316685310256), ('6.0', .9045826720013942)):
        transition = select(transition_sweep, ncrit=ncrit, flap='10', target_cl='0.8', element='fore', side='upper')
        require(close(transition['x_over_local_chord'], expected_x), 'Low-Ncrit transition association differs')
    changes, transition_rows = [], []
    for dataset in ('calibration', 'reserved'):
        for cl in CLS:
            a = select(pairs, dataset=dataset, target_cl=cl, ncrit=9.)
            b = select(pairs, dataset=dataset, target_cl=cl, ncrit=6.25)
            shifts = {}
            for flap in (5., 10.):
                p = {n: select(points, dataset=dataset, target_cl=cl, flap=flap, ncrit=n) for n in (9., 6.25)}
                shifts[flap] = (p[6.25]['CD'] - p[9.]['CD'])*1e4
                for n, series in ((9., 'multifoil_n9'), (6.25, 'multifoil_n6p25')):
                    original = select(comparison, dataset=dataset, target_cl=str(cl), flap_deg=str(flap), series=series)
                    require(close(p[n]['CD'], original['CD']), 'Selected-state comparison mismatch')
                    if cl == .8:
                        state = select(states, dataset=dataset, target_cl=str(cl), flap=str(int(flap)), ncrit=str(n))
                        surface = select(surfaces, state_id=state['state_id'], element='fore', side='upper')
                        transition_rows.append(dict(dataset=dataset, ncrit=n, flap=flap,
                            fore_upper_transition_x_local=surface['transition_x_local'],
                            laminar_to_te=surface['laminar_to_te'], source_state_id=state['state_id']))
            change = b['delta_counts']-a['delta_counts']
            require(close(change, shifts[10.]-shifts[5.]), 'Ncrit-change decomposition fails')
            changes.append(dict(dataset=dataset, target_cl=cl, ncrit_from=9., ncrit_to=6.25,
                flap5_drag_change_counts=shifts[5.], flap10_drag_change_counts=shifts[10.],
                paired_change_counts=change, error_at_n9=a['error_counts'], error_at_n6p25=b['error_counts'],
                signed_error_reduction_percent=(a['error_counts']-b['error_counts'])/a['error_counts']*100))
    closest = {}
    for cl in CLS:
        q = [r for r in pairs if r['dataset'] == 'calibration' and r['target_cl'] == cl
             and r['pair_feasible'] and r['candidate_converged'] == r['candidate_required']]
        closest[str(cl)] = min(q, key=lambda r: abs(r['error_counts']))
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), scope='Archived-data Ncrit sensitivity; no re-fit or new flow solves',
        calibration=dict(targets=len([r for r in points if r['dataset'] == 'calibration']),
                         converged=sum(r['feasible'] for r in points if r['dataset'] == 'calibration'), ncrit_values=14),
        reserved=dict(targets=18, converged=18, ncrit_values=[6.25, 9.],
                      initialization_qualification='Ncrit 6.25 originally completed 7/9; two separately documented recoveries included.'),
        changes=changes, closest_sampled_calibration_pairs=closest, archived_calibration_scores=calibration_scores,
        closest_sampled_qualification='Post-hoc diagnostic among complete candidates; not a new chosen Ncrit, fitted root, or validation pass.',
        mses_sensitivity_available=False, numerical_protocol='Original 260 fore / 180 aft panels; no finer-grid mixing',
        new_flow_solves=0, source_files_checked=len(inputs.consumed))
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / 'flap_error_vs_ncrit.csv', pairs)
    write_csv(out / 'ncrit_change_attribution.csv', changes)
    write_csv(out / 'fore_transition_cl08.csv', transition_rows)
    write_csv(out / 'calibration_transition_sweep.csv', transition_sweep)
    write_csv(out / 'absolute_drag_cl08.csv', [r for r in points if r['target_cl'] == .8 and r['flap'] in (5., 10.)])
    json_write(out / 'summary.json', summary)
    plot(out, points, pairs)
    notes(out, pairs, changes, closest)
    json_write(out / 'input_manifest.json', dict(files=inputs.consumed,
        report_script_sha256=sha(Path(__file__)), guards='Original delivery and archived mechanism-audit receipt; every consumed result checked.'))
    return summary


def notes(out, pairs, changes, closest):
    table = []
    for dataset in ('calibration', 'reserved'):
        for cl in CLS:
            a = select(pairs, dataset=dataset, target_cl=cl, ncrit=9.)
            b = select(pairs, dataset=dataset, target_cl=cl, ncrit=6.25)
            table.append(f"| {a['reynolds']:,.0f} | {cl:g} | {a['experiment_delta_counts']:+.2f} | {a['mses_n9_delta_counts']:+.2f} | {a['delta_counts']:+.2f} | {b['delta_counts']:+.2f} | {a['error_counts']:+.2f} → {b['error_counts']:+.2f} |")
    focus = select(changes, dataset='reserved', target_cl=.8)
    (out / 'FINDINGS.md').write_text(
        '# Ncrit influence on the S702 flap-drag discrepancy\n\n'
        'Ncrit strongly affects MultiFoil’s fixed-lift flap response, but the saved data do not remove the experimental discrepancy. '
        'Published MSES is available only at Ncrit 9: its own Ncrit response cannot be inferred quantitatively from MultiFoil. '
        'No MSES or MultiFoil flow solve, model change, or parameter refit was performed.\n\n'
        '## Consistent comparison\n\n'
        'Delta means CD(10°) − CD(5°); negative means a flap-drag benefit. One count is 0.0001 CD. '
        'All MultiFoil sensitivities use the original 260/180 panel protocol with common Ncrit on both elements. '
        'Re=1 million uses M=0.10; Re=700,000 uses M=0.07. These are separate matched-condition comparisons, not a pure Reynolds-number experiment. '
        'The latest refined-grid Ncrit-6.25 values are deliberately excluded from Ncrit differences because no corresponding refined-grid Ncrit-9 baseline was run.\n\n'
        '| Re | CL | Experiment | Published MSES N=9 | MultiFoil N=9 | MultiFoil N=6.25 | MultiFoil error N=9 → 6.25 |\n'
        '|---:|---:|---:|---:|---:|---:|---:|\n' + '\n'.join(table) + '\n\n'
        '![Ncrit, individual flap drag and paired response](ncrit_flap_error.png)\n\n'
        '## What changes at the specific Re=700,000, CL=0.8 condition\n\n'
        f"Lowering Ncrit from 9 to 6.25 raises the 5° drag by {focus['flap5_drag_change_counts']:.2f} counts and changes the 10° drag by {focus['flap10_drag_change_counts']:.2f} counts. "
        f"The paired prediction shifts by {focus['paired_change_counts']:.2f} counts. Its signed error decreases from {focus['error_at_n9']:.2f} to {focus['error_at_n6p25']:.2f} counts ({focus['signed_error_reduction_percent']:.1f}% reduction). "
        'About 91% of that change comes from raising the 5° baseline drag. At Ncrit 6.25, the 5° absolute error is only +0.05 count, '
        'whereas the 10° error remains +19.43 counts. The ranking changes sign, but the measured benefit is still severely understated. '
        'These percentages describe MultiFoil, not the fraction of the MSES error explained by Ncrit.\n\n'
        'The audited fore upper transition moves from x/c_local≈0.953 to 0.571 at 5°, while the 10° state remains laminar almost to the trailing edge '
        '(transition at approximately 0.962 for Ncrit 6.25). Earlier transition raises the modeled 5° baseline drag. '
        'At Ncrit 6.25, the flap change gives a modeled fore wake benefit of −21.72 counts offset by an aft penalty of +20.38. '
        'These are modeled element terms, not experimentally measured allocations of the residual error.\n\n'
        '## Why further lowering Ncrit is not an established fix\n\n'
        'The 14-value sweep exists at Re=1 million, not at the reserved Re=700,000 condition. At CL=0.8 its response is non-monotonic: '
        'Ncrit 9 gives +9.32 counts, 6.25 gives −21.17, 6 gives −23.09, and 4 gives −7.71, versus experiment −29.86. '
        f"The closest complete sampled candidate there is Ncrit {closest['0.8']['ncrit']:g}, still {closest['0.8']['error_counts']:.2f} counts above measurement. "
        'The 10° fore upper transition shifts from x/c_local=0.905 at Ncrit 6 to 0.409 at Ncrit 4; both flap settings then transition relatively early. '
        'This is consistent with the loss of differential benefit at low thresholds. '
        'This does not exclude an untested value or prove the response at Re=700,000.\n\n'
        'At CL=1.2, Ncrit 6.25 retains the wrong sign at both Reynolds numbers. On the Re=1-million grid, Ncrit 4 or 5 can reverse the sign, '
        'but the closest sampled complete result is −1.58 versus measured −6.45 counts, and these settings degrade the broader drag calibration. '
        'A near-zero sign at Ncrit 5.5 is not reliable evidence of an accurate gain. No sampled complete Ncrit reproduces all three measured flap differences.\n\n'
        '![Sensitivity at other lift levels](ncrit_other_lifts.png)\n\n'
        '## Meaning and limits\n\n'
        'Ncrit is a logarithmic instability-amplification threshold: lower values generally permit earlier natural transition. '
        'Changing it can increase skin-friction losses or reduce separation-bubble losses; drag therefore need not vary monotonically. '
        'A single common threshold does not explicitly represent spatially varying fore-wake disturbances, hinge-film roughness or bypass transition. '
        'See [Drela’s transition documentation](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt). '
        'The [MSES manual](https://web.mit.edu/drela/Public/web/mses/mses.pdf) lists ACRIT and separate envelope/frequency-tracking options. '
        'The precise S702 input deck and an executable are unavailable; equal numeric Ncrit alone does not establish identical transition predictions.\n\n'
        'Ncrit 6.25 remains the frozen calibration choice; 6 was already nearly equivalent for overall error and better for absolute/pairwise RMSE. '
        'This analysis does not select a new value or establish a measured tunnel threshold. Quantifying MSES’s own sensitivity would require '
        'its matching geometry/input deck and actual fixed-lift runs at declared thresholds. A useful experimental discriminator is transition/bubble '
        'location at both flap settings; tuning total drag alone cannot identify the missing mechanism.\n\n'
        '## Coverage and verification\n\n'
        'All 168 calibration entries were retained: 165 converged, with failed targets at (N=3, flap10, CL0.5), '
        '(N=11, flap10, CL0.5), and (N=12, flap10, CL0.8). The reserved data contain 18 successful targets across two Ncrit values. '
        'The N=6.25 reserved set required two separately recorded initialization recoveries after the original 7/9 completion. '
        'Failed targets have no force value; incomplete global candidates cannot win by omission. '
        'Curve joins are visual guides, not additional solved points. These are raw differences from digitized reference curves; '
        'no experimental-uncertainty or digitization-based pass threshold is assigned here. Numerical convergence is still unestablished.\n\n'
        'Every consumed table/result was checked against the original delivery or archived mechanism-audit receipt. '
        'All successful table forces match their result JSON; paired values reproduce the existing independent pairwise table, '
        'and N=9/6.25 values reproduce the audited comparison. The Ncrit-change attribution closes algebraically. '
        'Existing full-state audits are reused; no new full-state audit or numerical solver test is claimed. '
        'See input_manifest.json and verification.json for receipts.\n\n'
        'Reproduce: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 validation/report_s702_ncrit_flap_error.py`. '
        'References and digitization: [S702 NASA report](https://ntrs.nasa.gov/citations/20220014964), '
        'Figures 25(a,b) and 26(b,c), retained in the original benchmark source archive.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'validation/results/s702_ncrit_flap_error_20261001')
    args = parser.parse_args()
    print(json.dumps(report(args.out), indent=2))
