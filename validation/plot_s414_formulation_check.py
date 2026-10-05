#!/usr/bin/env python3
"""Figures for the S414 inviscid cross-formulation check (V3); reads its CSV outputs only."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_formulation_check_20261004'
TAPS = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv'
# Fixed entity colors (validated light-surface categorical slots); text stays in ink.
PRODUCTION, HESS, BETWEEN = '#2a78d6', '#eb6834', '#1baf7a'
SURFACE, INK, SECONDARY, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def style(ax):
    ax.set_facecolor(SURFACE); ax.grid(color=GRID, linewidth=.6); ax.set_axisbelow(True)
    for side in ('top', 'right'): ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'): ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8.5); ax.xaxis.label.set_color(SECONDARY); ax.yaxis.label.set_color(SECONDARY)


def save(fig, out, name):
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(out/f'{name}.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def nose(out):
    data = rows(out/'fore_nose_profiles.csv')
    taps = sorted({float(r['x_over_c']) for r in rows(TAPS) if r['case'] == 'cl0056' and r['method'] == 'multifoil'
                   and r['coordinate'] == 'table' and r['included'] == 'True' and r['surface'] == 'fore_lower'
                   and float(r['x_over_c']) < .06})
    panels = [('measured_cl0056', 'Measured incidence'), ('matched_cl0056', 'Matched-lift incidence, CL 0.056')]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]), layout='constrained')
    fig.patch.set_facecolor(SURFACE)
    for column, (label, title) in enumerate(panels):
        r = [x for x in data if x['scenario'] == 'nominal' and x['incidence'] == label and x['panels_fore'] == '340']
        r.sort(key=lambda x: float(x['s'])); x = np.array([float(v['x']) for v in r])
        p = np.array([float(v['production_Cp']) for v in r]); h = np.array([float(v['hess_Cp']) for v in r])
        keep = x < .06; top, low = axes[0, column], axes[1, column]
        top.plot(x[keep], p[keep], color=PRODUCTION, linewidth=1.6, solid_capstyle='round', label='Production streamfunction')
        top.plot(x[keep], h[keep], color=HESS, linewidth=1.6, solid_capstyle='round', label='Hess–Smith')
        top.plot(taps, np.full(len(taps), top.get_ylim()[1]), linestyle='none', marker='|', markersize=9,
                 markeredgewidth=1.2, color=MUTED, clip_on=False, label='Supported fore-lower taps')
        top.invert_yaxis(); style(top)
        top.set_title(f"{title}, α = {float(r[0]['alpha_deg']):.3f}°", fontsize=10, color=INK, loc='left')
        k = int(np.argmin(np.where(keep, p, np.inf)))
        top.annotate(f'minimum {p[k]:.3f}\nbetween taps', (x[k], p[k]), xytext=(.018, p[k]+.12), fontsize=8.5, color=SECONDARY,
                     arrowprops=dict(arrowstyle='-', color=MUTED, linewidth=.7))
        low.plot(x[keep], (h-p)[keep], color=SECONDARY, linewidth=1.4, solid_capstyle='round')
        low.axhline(0, color=AXIS, linewidth=.8); style(low); low.set_xlabel('x / total chord, fore-element nose')
        low.set_xlim(-.001, .06)
        # Clockwise node order: the lower surface precedes the leading edge.
        le = int(np.argmin(x)); sides = (('lower surface', np.arange(len(x)) < le), ('upper surface', np.arange(len(x)) > le))
        for name, side in sides:
            j = np.flatnonzero(side & (x > .044) & (x < .05))[0]
            top.annotate(name, (x[j], p[j]), xytext=(0, 9 if p[j] < 0 else -13), textcoords='offset points',
                         fontsize=8.5, color=SECONDARY, ha='center')
            low.annotate(name, (x[j], (h-p)[j]), xytext=(0, 7 if name.startswith('upper') else -13), textcoords='offset points',
                         fontsize=8.5, color=SECONDARY, ha='center')
    axes[0, 0].set_ylabel('Inviscid Cp at Mach 0.10 (suction upward)'); axes[1, 0].set_ylabel('Hess–Smith − production')
    handles, names = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, names, frameon=False, fontsize=8.5, loc='outside lower center', ncols=3, labelcolor=SECONDARY)
    fig.suptitle('Two inviscid formulations on the identical 340/240 polygon', fontsize=11.5, color=INK, x=.01, ha='left')
    save(fig, out, 'nose_pressure_formulations')


def refinement(out):
    data = rows(out/'refinement.csv'); grids = ('340', '480', '680'); labels = ('340 / 240', '480 / 340', '680 / 480')
    series = [('between_methods_nose_Cp_', 'Between the two methods', BETWEEN),
              ('production_nose_Cp_change_to_finest_', 'Production, change to finest grid', PRODUCTION),
              ('hess_nose_Cp_change_to_finest_', 'Hess–Smith, change to finest grid', HESS)]
    fig, ax = plt.subplots(figsize=(7.4, 4.4), layout='constrained'); fig.patch.set_facecolor(SURFACE)
    for prefix, label, color in series:
        worst = [max(float(r[prefix+g]) for r in data) for g in grids]
        ax.plot(range(3), worst, color=color, linewidth=1.6, marker='o', markersize=7, markeredgecolor=SURFACE,
                markeredgewidth=1.5, solid_capstyle='round', label=label)
        ax.annotate(f'{worst[0]:.4f}', (0, worst[0]), xytext=(-8, 0), textcoords='offset points', ha='right', va='center',
                    fontsize=8.5, color=SECONDARY)
        if prefix.startswith('between'):
            ax.annotate(f'{worst[-1]:.4f}', (2, worst[-1]), xytext=(8, 0), textcoords='offset points', va='center',
                        fontsize=8.5, color=SECONDARY)
    style(ax); ax.set_xticks(range(3), labels); ax.set_xlim(-.45, 2.45); ax.set_ylim(bottom=0)
    ax.set_xlabel('Panels, fore / aft (same continuous contour)')
    ax.set_ylabel('Largest fore-nose |ΔCp| over all ten cases')
    ax.legend(frameon=False, fontsize=8.5, labelcolor=SECONDARY, loc='upper right')
    ax.set_title('Formulation difference against each method’s own grid sensitivity', fontsize=11, color=INK, loc='left')
    save(fig, out, 'refinement_trend')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    nose(out); refinement(out)


if __name__ == '__main__':
    main()
