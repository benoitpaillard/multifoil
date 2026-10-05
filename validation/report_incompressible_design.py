#!/usr/bin/env python3
"""Report the completed incompressible design trial without empirical fitting."""
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from shape_optimization import camber_mode

RESULTS = ROOT/'validation/results'
OUT = RESULTS/'s702_incompressible_design'


def read(path):
    return json.loads(Path(path).read_text())


def mapping(record):
    return {(r['scenario'], r['point']): r for r in record['flows']}


def paired_changes(record, baseline, cfg):
    current, base = mapping(record), mapping(baseline)
    points, means = [], {}
    for scenario in cfg['scenarios']:
        total = 0.
        for p in cfg['points']:
            key = scenario['name'], p['name']
            a, b = current.get(key), base.get(key)
            if not a or not b or not a['feasible'] or not b['feasible']:
                return [], {}
            d = (a['CD']-b['CD'])*1e4
            total += p['weight']*d
            points.append(dict(scenario=key[0], point=key[1], CL=p['cl'],
                               baseline_CD=b['CD'], candidate_CD=a['CD'], delta_cd_counts=d,
                               delta_cm=a['CM']-b['CM']))
        means[scenario['name']] = total/sum(p['weight'] for p in cfg['points'])
    return points, means


def main():
    spec, summary, search = (read(OUT/f) for f in ('study.json', 'summary.json', 'search.json'))
    for name, expected in spec['code_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Design source no longer matches: '+name)
    cfg = spec['config']
    history = search['history']
    base = history[0]
    candidate = search['best']
    changed = any(abs(v) > 1e-14 for v in summary['candidate_parameters'])
    if not changed:
        alternatives = [r for r in history[1:] if paired_changes(r, base, cfg)[1]]
        if alternatives:
            candidate = min(alternatives, key=lambda r: max(paired_changes(r, base, cfg)[1].values()))
    chosen_points, _ = paired_changes(candidate, base, cfg)
    fine_points = []
    if changed:
        fine_points, _ = paired_changes(read(OUT/'verification_candidate/result.json'),
                                        read(OUT/'verification_baseline/result.json'), cfg)
    mach_records = read(RESULTS/'s702_mach_check/summary.json')
    cd = {(r['flap'], r['mach']): r['CD'] for r in mach_records}
    legacy_path = RESULTS/'s702_comparison/fixed_lift_comparison.csv'
    with legacy_path.open() as f:
        reference = {int(r['deflection']): r for r in csv.DictReader(f)
                     if r['setting'] == 'base' and float(r['target_cl']) == .8}
    gaps = [(cd[10, 0]-cd[5, 0])*1e4, (cd[10, .1]-cd[5, .1])*1e4,
            (float(reference[10]['mses_CD'])-float(reference[5]['mses_CD']))*1e4,
            (float(reference[10]['experiment_CD'])-float(reference[5]['experiment_CD']))*1e4]
    native = read(RESULTS/'incompressible_shape_check/summary.json')
    entries = []
    for r in history:
        points, means = paired_changes(r, base, cfg)
        entries.append(dict(folder=r['folder'], parameters=r['parameters'],
                            feasible=r['feasible'], reason=r['reason'], points=points, means=means))
    write = lambda p, v: p.write_text(json.dumps(v, indent=2, allow_nan=False)+'\n')
    write(OUT/'report_data.json', dict(mach_ranking_gap_counts=gaps, designs=entries,
          selected_for_display=candidate['folder'], base_point_changes=chosen_points,
          fine_point_changes=fine_points, native_check=native,
          reference_csv_sha256=hashlib.sha256(legacy_path.read_bytes()).hexdigest()))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.labelcolor': '#25364a', 'text.color': '#25364a',
                         'axes.titleweight': 'bold'})
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.8), layout='constrained')
    ax = axes[0, 0]
    labels = ['MultiFoil\nM=0', 'MultiFoil\nM=0.1', 'MSES\nM=0.1', 'Experiment\nM=0.1']
    bars = ax.bar(np.arange(4), gaps, color=['#116b87', '#63a7b4', '#8894a2', '#cf713d'], width=.6)
    ax.axhline(0, color='#566372', lw=.8)
    for bar, v in zip(bars, gaps):
        ax.annotate(f'{v:+.2f}', (bar.get_x()+bar.get_width()/2, v),
                    xytext=(0, 4 if v > 0 else -4), textcoords='offset points',
                    ha='center', va='bottom' if v > 0 else 'top')
    ax.set(xticks=np.arange(4), xticklabels=labels, ylim=(-38, 17),
           ylabel=r'$C_{D,10°}-C_{D,5°}$  (drag counts)',
           title='A  Mach zero leaves the flap ranking unchanged')
    ax.text(.02, .05, 'Positive: 5° preferred; negative: 10° preferred\nRe=1 million; CL=0.8; published references at M=0.1',
            transform=ax.transAxes, fontsize=8)
    ax = axes[0, 1]
    ordered = sorted(entries, key=lambda r: r['parameters'][0])
    for scenario, color in [('ncrit9', '#116b87'), ('ncrit6', '#bb6c34')]:
        values = [(r['parameters'][0]*1e3, r['means'][scenario]) for r in ordered if scenario in r['means']]
        ax.plot(*np.array(values).T, 'o-', color=color, label=scenario.replace('ncrit', 'Ncrit='))
    for r in ordered:
        if not r['feasible'] and r['means']:
            ax.scatter([r['parameters'][0]*1e3]*len(r['means']), list(r['means'].values()),
                       marker='x', s=95, color='#8e3440', zorder=5, label='Point restriction violated')
    ax.axhline(0, color='#566372', lw=.8)
    ax.axhline(-cfg['constraints']['acceptance_counts'], color='#566372', lw=.9, ls='--',
               label='Required mean change < −0.5 counts')
    ax.set(xlabel='Camber amplitude / local chord × 1000', ylabel='Mean ΔCD at CL=0.5 and 0.8 (counts)',
           title='B  Paired incompressible shape changes')
    ax.set_ylim(bottom=-.62)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=.15)
    ax = axes[1, 0]
    x = np.linspace(0, 1, 501)
    amp = candidate['parameters'][0]
    ax.plot(x, np.zeros_like(x), color='#a3aab3', label='Baseline')
    ax.plot(x, amp*1e3*camber_mode(x, *cfg['modes'][0]['support']), color='#116b87', lw=2,
            label=f'Displayed trial: a={amp:+.4f}')
    ax.set(xlabel='x / fore-element chord', ylabel='Normal displacement / local chord × 1000',
           title='C  Smooth camber change; ends and thickness fixed')
    ax.grid(alpha=.15)
    ax.legend(frameon=False)
    ax = axes[1, 1]
    positions = np.arange(len(chosen_points))
    ax.axhline(0, color='#566372', lw=.8)
    ax.axhline(cfg['constraints']['max_point_drag_increase_counts'], color='#8e3440', ls=':',
               label='Allowed pointwise increase')
    if fine_points:
        ax.bar(positions-.18, [p['delta_cd_counts'] for p in chosen_points], .34, color='#116b87', label='260/180 panels')
        ax.bar(positions+.18, [p['delta_cd_counts'] for p in fine_points], .34, color='#63a7b4', label='340/240 panels')
    else:
        ax.bar(positions, [p['delta_cd_counts'] for p in chosen_points], .6, color='#116b87', label='260/180 panels')
    labels = [f"Ncrit={p['scenario'].replace('ncrit', '')}\nCL={p['CL']:g}" for p in chosen_points]
    ax.set(xticks=positions, xticklabels=labels, ylabel='Candidate − own-grid baseline CD (counts)',
           title='D  Best small trial: changes at each lift')
    ax.legend(frameon=False, fontsize=8)
    passed = summary['numerical_screen_passed']
    fig.suptitle('S702 incompressible contour design: '+('numerical screen passed' if passed else 'baseline retained'),
                 fontsize=16, fontweight='bold')
    fig.supxlabel('Natural transition; figure-reconstructed geometry. Numerical screening does not establish experimental shape improvement.', fontsize=10)
    fig.savefig(OUT/'incompressible_design.png', dpi=180)
    fig.savefig(OUT/'incompressible_design.svg')
    plt.close(fig)
    flows = [read(p) for p in OUT.glob('design_*/*/result.json')]
    flows += [read(p) for p in OUT.glob('verification_*/*/result.json')]
    flow_count = sum(r['flow_evaluations'] for r in flows)
    converged = sum(r['admissible_evaluations'] for r in flows)
    interruption_path = OUT/'interruptions/INTERRUPTION.json'
    interrupted = read(interruption_path) if interruption_path.exists() else {'records': []}
    interrupted_count = sum(r['requested_completed'] for r in interrupted['records'])
    baseline_by = mapping(base)
    transition_delta = (baseline_by['ncrit6', 'cl08']['CD']-baseline_by['ncrit9', 'cl08']['CD'])*1e4
    md = ['# Executed incompressible contour-design trial', '',
          'Reproduced 8 September 2026. See [the design protocol](../../INCOMPRESSIBLE_DESIGN.md) for all settings and limits.', '',
          '**Outcome: '+('the changed contour passed the declared numerical screen.' if passed else 'the baseline was retained by the numerical acceptance rules.')+'**', '',
          f"The bounded search evaluated {len(history)} designs and stopped with `{search['status']}`. "
           'This is a local trial with a finite evaluation budget, not a proven optimum. No aerodynamic closure or empirical drag correction was changed.', '',
          '## Compressibility diagnosis', '',
          '| Fixed CL=0.8 comparison | CD(10°) − CD(5°), counts |', '|---|---:|']
    gap_names = ['MultiFoil, Mach 0', 'MultiFoil, Mach 0.1', 'Published MSES, Mach 0.1', 'Experiment, Mach 0.1']
    md += [f'| {label} | {v:+.6f} |' for label, v in zip(gap_names, gaps)]
    md += ['', f'The Mach change alters this ranking margin by {gaps[0]-gaps[1]:+.6f} counts. '
           'It does not explain the disagreement with experiment. Mach-zero values are sensitivity calculations; '
           'the published comparisons remain at Mach 0.1. Positive means 5° is preferred.', '',
           f'At fixed CL=0.8 on the 5-degree geometry, changing Ncrit from 9 to 6 changes '
           f'the incompressible baseline drag by {transition_delta:+.4f} counts. This model is far more '
           'sensitive to the transition assumption here than to Mach 0.1 versus zero. '
           'Ncrit=6 was not selected to fit the experiment and is not a universal drag correction.', '',
           '## Actual contour search', '',
           'Mach=0, Re=1 million, S702 5-degree flap. Only the fore-element camber changes. '
           'Each transition scenario gives equal weight to CL=0.5 and 0.8.', '',
           '| Fore camber amplitude / local chord | Mean ΔCD, Ncrit 9 (counts) | Mean ΔCD, Ncrit 6 (counts) | Design admissibility |',
           '|---:|---:|---:|---|']
    for r in ordered:
        fmt = lambda key: f"{r['means'][key]:+.6f}" if key in r['means'] else 'unavailable'
        md.append(f"| {r['parameters'][0]:+.5f} | {fmt('ncrit9')} | {fmt('ncrit6')} | {'admissible' if r['feasible'] else r['reason']} |")
    md += ['', 'The search objective is the worse of the two scenario means, each paired with '
           'its own baseline. Every operating point also has a one-count deterioration limit '
           'and a 0.015 moment-change limit. Geometry and separation restrictions apply before acceptance.', '',
           f"Search candidate: `{summary['candidate_parameters']}`. Delivered amplitudes: `{summary['delivered_parameters']}`. "
           f"Verification status: `{summary['verification']['status']}`.", '',
           (f"The table and plot below show trial `{candidate['folder']}`, amplitude `{candidate['parameters']}`. "
            'Its mean gains did not exceed the predeclared 0.5-count threshold in both scenarios, '
            'so it was not accepted and no aerodynamic finer-grid run was triggered.' if not changed else
            f"The table below compares the selected trial `{candidate['folder']}` on both grids."), '',
           '| Scenario | CL | Base-grid ΔCD, counts | Fine-grid ΔCD, counts |', '|---|---:|---:|---:|']
    for i, p in enumerate(chosen_points):
        fine_value = f"{fine_points[i]['delta_cd_counts']:+.6f}" if fine_points else 'not run'
        md.append(f"| {p['scenario']} | {p['CL']:g} | {p['delta_cd_counts']:+.6f} | {fine_value} |")
    md += ['', 'If a trial is displayed after baseline retention, this table describes that trial, '
           'not an accepted improvement. All coefficients and geometry files are retained in the trial folders.', '',
           f'The design run contains {len(flows)} fixed-lift point evaluations and {flow_count} requested flow solves; '
           f'{converged} of those requested solves were admissible. Internal angle-continuation steps are additional. '
           'A design can be rejected for performance restrictions even when all its flow solutions converge.', '',
           (f'An execution-runtime interruption required replaying {len(interrupted["records"])} unfinished targets. '
            f'An additional {interrupted_count} completed requested solves from the interrupted attempts '
            'are preserved under `interruptions/`, separately from the completed study count above. '
            'This interruption was not classified as aerodynamic nonconvergence.' if interrupted_count else
            'No interrupted-attempt records were present in this run.'), '',
           '## Single-element implementation check', '',
           'NACA0012, Mach=0, Re=1 million, Ncrit=9, CL=0.5, 200 panels. '
           'Native mfoil uses its own coupled Newton/lift-constraint system on identical actual surface coordinates. '
           'The Python reference needed its existing sparse reallocation adapter plus two array-layout repairs '
           '(the Nx2 edge-velocity basis and scalar lift-residual concatenation). '
           'The velocity-basis derivative was checked against native edge velocities. '
           'The unmodified failures and runner snapshots are preserved in sibling result folders. '
           'The vendored solver source is unchanged.', '',
           '| Camber amplitude | MultiFoil CD | Native mfoil CD | Difference, counts |', '|---:|---:|---:|---:|']
    for r in native['records']:
        md.append(f"| {r['amplitude']:+.4f} | {r['multifoil']['CD']:.9f} | {r['native_mfoil']['CD']:.9f} | {r['delta_cd_counts']:+.6f} |")
    d = native['central_dcd_damplitude']
    relative = abs(d['multifoil']/d['native_mfoil']-1)*100
    md += ['', f"Central drag sensitivity dCD/da: MultiFoil {d['multifoil']:+.9f}, native mfoil {d['native_mfoil']:+.9f}. "
           f"Same sign: {native['same_derivative_sign']}; magnitude difference {relative:.2f}%. This is implementation verification with shared closures, "
           'not independent aerodynamic validation or an MSES run.', '',
           '## What remains uncertain', '',
           'The S702 source is a vector-figure reconstruction, not the original coordinate deck. '
           'No measured polar exists for these new contours. Two Ncrit settings are sensitivity choices; '
           'they do not model a specified trip or all disturbances. Bubble drag, hinge-film effects and '
           'wake interaction remain limitations. A two-grid screening spread measures numerical sensitivity, not a physical error bound.', '',
           'This workflow makes low-Mach shape trials concrete and reviewable. It does not correct '
           'the demonstrated experimental flap-ranking error or establish reliable multi-element shape gradients.', '',
           '## Reproduce', '', '```bash',
           'python shape_optimization.py examples/s702_incompressible_design.json --out validation/results/s702_incompressible_design --workers 4',
           'python validation/s702_mach_check.py --workers 2',
           'python validation/incompressible_shape_check.py',
           'python validation/report_incompressible_design.py', '```', '',
           'Scientific context: [Drela, LINDOP section 14](https://web.mit.edu/drela/Public/web/mses/lindop.pdf) '
           'supports multipoint, low-dimensional design with transition precautions. '
           '[Fidkowski, mfoil](https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf) documents the underlying method. '
           'The [NASA S702 report](https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf) '
           'provides the geometry drawing and the original MSES/experimental reference.', '']
    (OUT/'INCOMPRESSIBLE_RESULTS.md').write_text('\n'.join(md))


if __name__ == '__main__':
    main()
