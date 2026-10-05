#!/usr/bin/env python3
"""Compare the S702 trial with published MSES curves and wind-tunnel symbols."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from s702_geometry import DATA, geometry

DEGREES = [-2, 0, 5, 10]
COLORS = {'multifoil': '#2166ac', 'mses': '#c66b20', 'experiment': '#20262c'}
URL = 'https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reference(angle, kind, curve):
    return np.genfromtxt(DATA/f'd{angle:+d}'/f'{kind}_{curve}.csv', delimiter=',', names=True)


def interpolate_reference(angle, kind, target):
    rows = reference(angle, kind, 'drag')
    x, y = rows['CL'], rows['CD']
    if np.any(np.diff(x) <= 0):
        raise ValueError('Reference drag curve is not strictly increasing in lift')
    if not x[0] <= target <= x[-1]:
        return dict(CD=None, lower_CL=None, upper_CL=None)
    j = int(np.clip(np.searchsorted(x, target), 1, len(x)-1))
    return dict(CD=float(np.interp(target, x, y)), lower_CL=float(x[j-1]), upper_CL=float(x[j]))


def load_study(folder):
    spec = json.loads((folder/'study.json').read_text())
    for name, sha in spec['sha256'].items():
        if digest(ROOT/name) != sha:
            raise ValueError('Current source differs from archived study: '+name)
    if (spec['reynolds'], spec['mach'], spec['ncrit']) != (1e6, .1, 9.):
        raise ValueError('The extracted references require Re=1e6, Mach=.1, Ncrit=9')
    jobs = []
    for j in spec['jobs']:
        path = folder/f"d{j['deflection']:+d}__{j['setting']}"/'result.json'
        jobs.append(json.loads(path.read_text()))
    return spec, jobs


def make_plots(out, polar, fixed_rows):
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 12,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.5), sharex=True, sharey=True)
    for ax, angle in zip(axes.flat, DEGREES):
        for kind in ['experiment', 'mses']:
            ref = reference(angle, kind, 'drag')
            ax.plot(ref['CL'], ref['CD']*1e4, color=COLORS[kind],
                marker='o' if kind == 'experiment' else None, ms=4, mfc='white',
                lw=1 if kind == 'experiment' else 1.8,
                label='Wind tunnel (published)' if kind == 'experiment' else 'MSES (published)')
        for j in sorted(polar, key=lambda r: r['setting']):
            if j['deflection'] != angle or j['setting'] != 'base':
                continue
            rows = sorted(j['records'], key=lambda r: r['alpha_deg'])
            xy = np.array([[r['CL'], r['CD']*1e4] if r['admissible'] else [np.nan, np.nan] for r in rows])
            ax.plot(*xy.T, color=COLORS['multifoil'], ls='-', marker='s', ms=3.5, lw=2,
                label='MultiFoil 260/180')
        for setting in ['base', 'coarse']:
            chosen = [r for r in fixed_rows if r['deflection'] == angle and
                      r['setting'] == setting and r['feasible']]
            if chosen:
                ax.scatter([r['CL'] for r in chosen], [r['CD']*1e4 for r in chosen],
                    marker='D' if setting == 'base' else 's', s=32 if setting == 'base' else 60,
                    facecolors=COLORS['multifoil'] if setting == 'base' else 'none',
                    edgecolors=COLORS['multifoil'], zorder=5,
                    label=None if setting == 'base' else 'MultiFoil 180/120, fixed lift')
        ax.set_title(f'Flap deflection {angle:+d} degrees')
        ax.grid(alpha=.18)
        ax.set_xlim(0, 1.55)
        ax.set_ylim(40, 240)
        ax.set_xlabel('$C_L$')
        ax.set_ylabel('$C_D$ (drag counts)')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(.5, .935))
    fig.suptitle('S702: MultiFoil versus published MSES and experiment', y=.99, fontsize=16)
    fig.text(.5, .945, '$Re=10^6$, $M=0.10$, free transition; computational $N_{crit}=9$', ha='center')
    fig.text(.5, .028, 'Geometry: Figure 1 vector reconstruction. 1 drag count = $10^{-4}$. No force calibration.',
             ha='center', fontsize=9)
    fig.text(.5, .009, 'Blue line joins polar samples; diamonds are fixed-lift solves. Panel sensitivity is compared at common lift.',
             ha='center', fontsize=8.5)
    fig.tight_layout(rect=[0, .055, 1, .88])
    fig.savefig(out/'s702_drag_comparison.png', dpi=180)
    fig.savefig(out/'s702_drag_comparison.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True, sharey=True)
    for ax, angle in zip(axes.flat, DEGREES):
        for kind in ['experiment', 'mses']:
            ref = reference(angle, kind, 'lift')
            ax.plot(ref['alpha_deg'], ref['CL'], color=COLORS[kind],
                marker='o' if kind == 'experiment' else None, ms=3.5, mfc='white', lw=1.3,
                label='Wind tunnel' if kind == 'experiment' else 'MSES (published)')
        job = next(j for j in polar if j['deflection'] == angle and j['setting'] == 'base')
        rows = sorted(job['records'], key=lambda r: r['alpha_deg'])
        ax.plot([r['alpha_deg'] for r in rows], [r['CL'] if r['admissible'] else np.nan for r in rows],
                color=COLORS['multifoil'], marker='s', ms=3.5, label='MultiFoil')
        ax.set_title(f'Flap {angle:+d} degrees')
        ax.set_xlabel('Incidence (degrees)')
        ax.set_ylabel('$C_L$')
        ax.grid(alpha=.18)
        ax.set_xlim(-5, 18)
        ax.set_ylim(-.1, 2.3)
    fig.legend(*axes.flat[0].get_legend_handles_labels(), loc='upper center', ncol=3,
               frameon=False, bbox_to_anchor=(.5, .95))
    max_alpha = max(r['alpha_deg'] for j in polar if j['setting'] == 'base' for r in j['records'])
    fig.suptitle(f'S702 lift curves: MultiFoil sampled through {max_alpha:g} degrees', y=.99, fontsize=15)
    fig.text(.5, .015, 'No stall calibration or MultiFoil maximum-lift validation is included in this comparison.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .045, 1, .91])
    fig.savefig(out/'s702_lift_comparison.png', dpi=170)
    plt.close(fig)

    if fixed_rows:
        targets = sorted(set(r['target_cl'] for r in fixed_rows))
        fig, axes = plt.subplots(1, len(targets), figsize=(11.2, 3.9), squeeze=False)
        for ax, target in zip(axes[0], targets):
            rows = sorted([r for r in fixed_rows if r['target_cl'] == target and r['setting'] == 'base'],
                          key=lambda r: r['deflection'])
            for kind, field in [('experiment', 'experiment_CD'), ('mses', 'mses_CD'), ('multifoil', 'CD')]:
                ax.plot([r['deflection'] for r in rows],
                        [r[field]*1e4 if r[field] is not None else np.nan for r in rows],
                        marker='o', color=COLORS[kind], ms=5, mfc='white' if kind == 'experiment' else COLORS[kind],
                        label={'experiment': 'Wind tunnel', 'mses': 'MSES (published)', 'multifoil': 'MultiFoil'}[kind])
            ax.set_title(f'$C_L={target:g}$')
            ax.set_xticks(DEGREES)
            ax.set_xlabel('Flap deflection (degrees)')
            ax.set_ylabel('$C_D$ (drag counts)')
            ax.grid(alpha=.2)
        fig.suptitle('S702 drag at common lift', fontsize=15, y=.99)
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='upper center', ncol=3,
                   frameon=False, bbox_to_anchor=(.5, .925))
        fig.text(.5, .015, 'MultiFoil: converged fixed-lift solves. References: bracketed linear interpolation of published plots.',
                 ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .04, 1, .87])
        fig.savefig(out/'s702_fixed_lift.png', dpi=180)
        plt.close(fig)


def write_narrative(out, summary, rows, grid, rankings, polar):
    lines = ['# S702: MultiFoil, published MSES and experiment', '',
        'Executed 8 September 2026. The numerical flow solver is unchanged from the preceding package.', '',
        'This trial compares drag polars and the movement of the low-drag range with flap setting. '
        'Experimental configuration rankings remain a separate requirement from agreement with '
        'published MSES. These results do not establish reliable contour-shape '
        'gradients or an independently validated optimum.', '',
        '[Primary NASA report]('+URL+'). Geometry is extracted from vector Figure 1; '
        'MSES paths and experimental marker centers come from Figures 23(c), 24(a), 25(b) and 26(c). '
        'Original author geometry decks and raw tunnel data were not obtained. See '
        '[the benchmark protocol](../../S702_BENCHMARK.md) for source locations, assumptions and commands.', '',
        'Conditions: Re=1,000,000; Mach=0.10; free transition; computational Ncrit=9; '
        'flap deflections -2, 0, +5 and +10 degrees. Only the small flap on the aft element moves. '
        'There is no force calibration or fitted drag correction.', '',
        '![Three-way drag comparison](s702_drag_comparison.png)', '',
        '## Execution coverage', '',
        '| Study | Requested flow evaluations | Numerically admissible |',
        '| --- | ---: | ---: |']
    for m in summary['studies']:
        lines.append(f"| {Path(m['path']).name} | {m['requested']} | {m['admissible']} |")
    lines += ['', f"The comparison includes {summary['requested']} requested flow evaluations, "
        f"{summary['admissible']} admissible, and {summary['successful_fixed_targets']}/"
        f"{summary['fixed_targets']} successful fixed-lift targets. Internal continuation evaluations "
        'are additional and are archived in the flow diagnostics. A one-point initialization probe '
        'is retained separately in `s702_probe` and is not counted in this table.', '',
        'Base resolution is 260/180 panels; coarse resolution is 180/120. Fixed-lift solves use '
        'CL tolerance 2e-5 and residual tolerance 1e-7. The two polar grids use different incidence '
        'samples, so panel sensitivity is evaluated from the actual fixed-lift solutions. '
        'Lines on the polar plot join samples; diamonds show computed fixed-lift values.', '',
        '## Drag at common lift', '',
        'One drag count is 0.0001 in CD. References are linearly interpolated between recorded '
        'bracketing markers or MSES vertices; no extrapolation is used. Values below are rounded '
        'for reading; the [CSV](fixed_lift_comparison.csv) retains full precision and reference brackets.', '',
        '| CL | Flap (deg) | MultiFoil CD | Published MSES CD | Experimental CD | MF vs MSES (%) | MF vs experiment (%) |',
        '| ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    def f(value, fmt='.6f'):
        return 'unavailable' if value is None else format(value, fmt)
    for r in sorted([r for r in rows if r['setting'] == 'base'], key=lambda r: (r['target_cl'], r['deflection'])):
        lines.append(f"| {r['target_cl']:g} | {r['deflection']:+d} | {f(r['CD'])} | "
            f"{f(r['mses_CD'])} | {f(r['experiment_CD'])} | "
            f"{f(r['error_vs_mses_percent'], '+.1f')} | {f(r['error_vs_experiment_percent'], '+.1f')} |")
    lines += ['', '![Drag at common lift](s702_fixed_lift.png)', '']
    for kind in ['mses', 'experiment']:
        s = summary.get('base_error_vs_'+kind)
        if s:
            lines.append(f"Across these {s['points']} base-resolution targets, the mean absolute "
                f"percentage difference from {kind.upper() if kind == 'mses' else kind} is "
                f"{s['mean_absolute_percent']:.2f}%; signed differences range from "
                f"{s['min_percent']:+.2f}% to {s['max_percent']:+.2f}%. These are descriptive "
                'statistics for this selected matrix, not an uncertainty bound for other designs.')
            lines.append('')
    lines += ['## Configuration rankings and optimization', '',
        '| CL | Source / resolution | Ascending nominal drag: flap angles | Best-to-second margin (counts) |',
        '| ---: | --- | --- | ---: |']
    for r in rankings:
        if r['source'] != 'MultiFoil' and r['setting'] != 'base':
            continue
        label = r['source']+(' '+r['setting'] if r['source'] == 'MultiFoil' else '')
        order = ', '.join(f'{v:+d}' for v in r['ascending_drag_deflections'])
        lines.append(f"| {r['target_cl']:g} | {label} | {order} | {r['best_to_second_margin_counts']:.2f} |")
    focus = {r['deflection']: r for r in rows if r['setting'] == 'base' and r['target_cl'] == .8}
    if 5 in focus and 10 in focus and focus[5]['CD'] is not None and focus[10]['CD'] is not None:
        lines += ['', 'The CL=0.8 comparison exposes the most consequential disagreement in the '
            'selected matrix. Define the drag change as CD(+10 degrees) minus CD(+5 degrees):', '',
            '| Source | Drag change (counts) | Nominal preference between these two settings |',
            '| --- | ---: | --- |']
        for name, field in [('MultiFoil base', 'CD'), ('Published MSES', 'mses_CD'), ('Experiment', 'experiment_CD')]:
            delta = (focus[10][field]-focus[5][field])*1e4
            lines.append(f"| {name} | {delta:+.2f} | {'+5 degrees' if delta > 0 else '+10 degrees'} |")
        predicted_delta = focus[10]['CD']-focus[5]['CD']
        measured_delta = focus[10]['experiment_CD']-focus[5]['experiment_CD']
        if predicted_delta*measured_delta < 0:
            lines += ['', 'The computed and measured drag changes have opposite signs in this comparison. '
                'An optimizer choosing a flap setting solely from this prediction can therefore select a '
                'different configuration from the nominal experimental optimum. Agreement with MSES does '
                'not resolve that discrepancy.', '']
        else:
            lines += ['', 'The computed and measured drag changes have the same nominal sign in this comparison.', '']
        lines += ['The reference differences are based on the published plots; a formal statistical '
            'ranking claim would also require experimental uncertainty and exact geometry/transition matching.', '']
    if grid:
        lines += [f"The largest absolute coarse-minus-base fixed-lift drag difference is "
            f"{summary['max_panel_difference_counts']:.3f} counts across {len(grid)} matched pairs. "
            'See [panel sensitivity](panel_sensitivity.csv). Two resolutions check numerical sensitivity; '
            'they do not establish asymptotic grid convergence or bound model discrepancy.', '']
    lines += ['Tiny nominal ranking margins should remain unresolved when geometry, transition or '
        'measurement uncertainty is not bounded. This four-setting study tests configuration changes; '
        'it does not validate local contour perturbations. MultiFoil can support exploratory screening, '
        'with independent verification of shortlisted designs and particular attention to drag-bucket edges.', '',
        '## Lift, failures and physical limits', '',
        '![Lift comparison](s702_lift_comparison.png)', '',
        f"The main polar samples through {max(r['alpha_deg'] for j in polar if j['setting'] == 'base' for r in j['records']):g} degrees. "
        'No dedicated MultiFoil maximum-lift score is computed. The published zero-flap curves '
        'give MSES CLmax approximately 2.0088 and experimental CLmax '
        'approximately 1.8005, an 11.6% overprediction at Re=1e6. This checks the source extraction; '
        'it is not a new MultiFoil maximum-lift result.', '']
    failures = [(j, r) for j in polar for r in j['records'] if not r['admissible']]
    for j, r in failures:
        lines.append(f"- Retained polar failure: flap {j['deflection']:+d} degrees, {j['setting']}, "
            f"alpha {r['alpha_deg']:g} degrees: {r['failure_reason']}. Failed forces remain null/NaN.")
    for r in rows:
        if not r['feasible']:
            lines.append(f"- Unavailable fixed-lift target: CL={r['target_cl']:g}, flap {r['deflection']:+d}, "
                         f"{r['setting']}: {r['reason']}")
    lines += ['', 'Neither numerical convergence nor failure establishes a physical stall boundary. '
        'The existing model does not resolve massive separation, wake merging, unsteady stall, '
        'forced trips or the tunnel model\'s plastic-film hinge cover. Plot-vector precision is not '
        'a total experimental uncertainty bound. The plotted fore-element thickness is about '
        '12.943% in the retained vertical coordinate, compared with the report\'s stated airfoil '
        'thickness of 12.91%; the geometry has not been rescaled to force agreement.', '',
        '## Reproduction and saved evidence', '',
        'The archive includes source PDF, extracted geometry and reference CSVs, drawing-index '
        'audit, source overlays, requested-flow records, surface/wake data, continuation diagnostics, '
        'fixed-lift traces, manifests and the comparison scripts. Numerical source hashes are checked '
        'by the report. [The benchmark protocol](../../S702_BENCHMARK.md) gives the commands. '
        'Additional software, download and extraction checks are recorded in `execution.json`.', '']
    (out/'S702_RESULTS.md').write_text('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--polar', type=Path, action='append', required=True)
    p.add_argument('--fixed', type=Path, action='append', default=[])
    p.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_comparison')
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    polar, fixed, manifests = [], [], []
    for kind, folders in [('polar', args.polar), ('fixed', args.fixed)]:
        for folder in folders:
            folder = folder.resolve()
            spec, jobs = load_study(folder)
            if spec['mode'] != kind:
                raise ValueError('Study mode does not match report argument')
            (polar if kind == 'polar' else fixed).extend(jobs)
            manifests.append(dict(path=str(folder.relative_to(ROOT)), sha256=digest(folder/'study.json'),
                source_hashes_verified=True, requested=sum(j['flow_evaluations'] for j in jobs),
                admissible=sum(j['admissible_evaluations'] for j in jobs)))
    seen = [(j['deflection'], j['setting']) for j in polar]
    if len(seen) != len(set(seen)):
        raise ValueError('Duplicate polar setting/deflection')
    if set((d, 'base') for d in DEGREES)-set(seen):
        raise ValueError('A base polar is required for each published flap setting')
    rows = []
    for job in fixed:
        for t in job['targets']:
            row = dict(deflection=job['deflection'], setting=job['setting'], **t)
            row.pop('bounds', None)
            if t['feasible'] and abs(t['CL']-t['target_cl']) > 2e-5:
                raise ValueError('Fixed-lift tolerance violated')
            for kind in ['mses', 'experiment']:
                ref = interpolate_reference(job['deflection'], kind, t['target_cl'])
                row.update({kind+'_'+key: val for key, val in ref.items()})
            row['error_vs_experiment_counts'] = ((row['CD']-row['experiment_CD'])*1e4
                if row['CD'] is not None and row['experiment_CD'] is not None else None)
            row['error_vs_experiment_percent'] = (100*(row['CD']/row['experiment_CD']-1)
                if row['CD'] is not None and row['experiment_CD'] is not None else None)
            row['error_vs_mses_counts'] = ((row['CD']-row['mses_CD'])*1e4
                if row['CD'] is not None and row['mses_CD'] is not None else None)
            row['error_vs_mses_percent'] = (100*(row['CD']/row['mses_CD']-1)
                if row['CD'] is not None and row['mses_CD'] is not None else None)
            rows.append(row)
    keys = [(r['setting'], r['deflection'], r['target_cl']) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate fixed-lift result')
    write_csv(args.out/'fixed_lift_comparison.csv', rows)
    grid = []
    for r in rows:
        if r['setting'] != 'base' or not r['feasible']:
            continue
        other = next((s for s in rows if s['setting'] == 'coarse' and s['deflection'] == r['deflection']
                      and s['target_cl'] == r['target_cl'] and s['feasible']), None)
        if other:
            grid.append(dict(deflection=r['deflection'], target_cl=r['target_cl'],
                base_CD=r['CD'], coarse_CD=other['CD'],
                coarse_minus_base_counts=(other['CD']-r['CD'])*1e4))
    write_csv(args.out/'panel_sensitivity.csv', grid)
    rankings = []
    for target in sorted(set(r['target_cl'] for r in rows)):
        for setting in ['base', 'coarse']:
            chosen = [r for r in rows if r['target_cl'] == target and r['setting'] == setting]
            if len(chosen) != 4:
                continue
            for name, field in [('MultiFoil', 'CD'), ('MSES', 'mses_CD'), ('experiment', 'experiment_CD')]:
                if all(r[field] is not None for r in chosen):
                    ranks = sorted(chosen, key=lambda r: r[field])
                    rankings.append(dict(target_cl=target, setting=setting, source=name,
                        ascending_drag_deflections=[r['deflection'] for r in ranks],
                        best_to_second_margin_counts=(ranks[1][field]-ranks[0][field])*1e4))
    (args.out/'rankings.json').write_text(json.dumps(rankings, indent=2)+'\n')
    make_plots(args.out, polar, rows)
    summary = dict(validation_pass_claim=False, studies=manifests,
        requested=sum(m['requested'] for m in manifests), admissible=sum(m['admissible'] for m in manifests),
        fixed_targets=len(rows), successful_fixed_targets=sum(r['feasible'] for r in rows),
        max_panel_difference_counts=max((abs(r['coarse_minus_base_counts']) for r in grid), default=None),
        geometry_status='vector-figure reconstruction; author original deck not obtained',
        source_url=URL, report_code_sha256=digest(Path(__file__)),
        reference_files_sha256={str(f.relative_to(ROOT)): digest(f) for f in sorted(DATA.rglob('*')) if f.is_file()})
    for kind in ['mses', 'experiment']:
        errors = [r['error_vs_'+kind+'_percent'] for r in rows if r['setting'] == 'base'
                  and r['error_vs_'+kind+'_percent'] is not None]
        if errors:
            summary['base_error_vs_'+kind] = dict(points=len(errors),
                min_percent=min(errors), max_percent=max(errors),
                mean_absolute_percent=float(np.mean(np.abs(errors))))
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    write_narrative(args.out, summary, rows, grid, rankings, polar)
    print(json.dumps({k: v for k, v in summary.items() if k not in ['reference_files_sha256', 'studies']}, indent=2))


if __name__ == '__main__':
    main()
