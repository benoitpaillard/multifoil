#!/usr/bin/env python3
"""Figures for V5 tranche 1 from saved outputs only; no flow or boundary-layer calculation."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_laminar_replay_20261004'
BRANCH = ROOT/'validation/results/s414_direct_branch_diagnostic_20261004'
TAPS = ROOT/'validation/results/s414_fore_nose_20261003/pressure_at_taps.csv'
CASES = (('cl0056', 'CL 0.056'), ('cl0175', 'CL 0.175'), ('cl028', 'CL 0.284'))
NU, X_END = 1e-6, .03424
SURFACE, INK, SECONDARY, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e6e5e0'
SERIES = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100')  # validated categorical order, light surface


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID)
    ax.grid(color=GRID, linewidth=.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=SECONDARY, labelsize=9, length=0)


def after_leading_edge(x):
    return slice(int(np.argmin(x)), None)


def save(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(OUT/f'{name}.{ext}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def branches(z):
    """One prescribed pressure, two boundary layers."""
    rows = list(csv.DictReader((BRANCH/'stations.csv').open()))
    roots = list(csv.DictReader((BRANCH/'station_roots.csv').open()))
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 9.6), sharex=True)
    fig.patch.set_facecolor(SURFACE)
    for column, (case, title) in enumerate(CASES):
        r = [q for q in rows if q['case'] == case]
        x = np.array([float(q['x']) for q in r]); keep = after_leading_edge(x)
        laminar = np.array([q['archive_laminar'] == 'True' for q in r])
        value = lambda key: np.array([float(q[key]) if q[key] else np.nan for q in r])
        double = [float(q['x']) for q in roots if q['case'] == case and q['two_distinct_roots'] == 'True']
        ax = axes[0, column]; style(ax)
        mx, mcp = z[f'native_{case}/x'], z[f'native_{case}/cp']; k = after_leading_edge(mx)
        ax.plot(mx[k], mcp[k], color=SECONDARY, lw=2)
        ax.invert_yaxis(); ax.set_title(title+(' (control: one root)' if not double else ''), color=INK, fontsize=12, loc='left')
        ax = axes[1, column]; style(ax)
        ax.plot(x[keep][laminar[keep]], value('archive_Hk')[keep][laminar[keep]], color=SERIES[0], lw=2, marker='o', ms=6.5,
                mec=SURFACE, mew=1.5, label='Archived coupled state')
        ax.plot(x[keep], value('march_Hk')[keep], color=SERIES[1], lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5,
                label='Same equations and pressure,\nstarted from the upstream station')
        fx, fh = z[f'native_{case}/fd_x'], z[f'native_{case}/fd_H']; k = after_leading_edge(fx)
        ax.plot(fx[k], fh[k], color=SERIES[2], lw=2, label='Independent finite-difference solution')
        ax = axes[2, column]; style(ax)
        ax.plot(x[keep][laminar[keep]], value('archive_N')[keep][laminar[keep]], color=SERIES[0], lw=2, marker='o', ms=6.5,
                mec=SURFACE, mew=1.5)
        ax.plot(x[keep], value('march_N')[keep], color=SERIES[1], lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5)
        ax.axhline(9, color=SECONDARY, lw=1); ax.set_ylim(-.4, 10)
        ax.text(.0335, 9.15, 'Ncrit = 9', color=SECONDARY, fontsize=9, ha='right', va='bottom')
        ax.set_xlabel('x / total chord', color=SECONDARY, fontsize=10)
        for row in (1, 2):
            for position in double[:1]:
                axes[row, column].axvline(position, color=SECONDARY, lw=1)
        if double:
            axes[2, column].annotate('two roots at\nthis station', (double[0], 8.6), xytext=(6, 0),
                                     textcoords='offset points', color=SECONDARY, fontsize=9, va='top')
        ends = np.where(laminar[keep])[0]
        last = x[keep][ends[-1]]
        if double:
            axes[2, column].annotate(f'{value("archive_N")[keep][ends[-1]]:.1f}', (last, value('archive_N')[keep][ends[-1]]),
                                     xytext=(7, 0), textcoords='offset points', color=INK, fontsize=9, va='center')
            axes[2, column].annotate(f'{value("march_N")[keep][ends[-1]]:.1f}', (last, value('march_N')[keep][ends[-1]]),
                                     xytext=(7, 0), textcoords='offset points', color=INK, fontsize=9, va='center')
    axes[0, 0].set_ylabel('Prescribed Cp\n(the model\'s own)', color=SECONDARY, fontsize=10)
    axes[1, 0].set_ylabel('Shape factor H', color=SECONDARY, fontsize=10)
    axes[2, 0].set_ylabel('Amplification N', color=SECONDARY, fontsize=10)
    axes[0, 0].set_xlim(0, X_END)
    handles, labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .925), ncol=3, frameon=False, fontsize=10, labelcolor=INK)
    fig.suptitle('One prescribed pressure, two boundary layers', color=INK, fontsize=15, x=.065, ha='left', y=.985)
    fig.text(.065, .94, 'Fore-lower nose. Both the separated and the attached state satisfy the inherited equations to 1e−9 '
             'with the same edge velocity and the same upstream state.', color=SECONDARY, fontsize=10, ha='left')
    fig.tight_layout(rect=(0, 0, 1, .895))
    save(fig, 'two_branches')


def taps(z):
    """Measured taps with four between-tap reconstructions and the finite-difference response."""
    measured = {}
    for row in csv.DictReader(TAPS.open()):
        if row['coordinate'] == 'table':
            measured.setdefault(row['case'], []).append((float(row['x']), float(row['experiment_Cp'])))
    experiment = {(r['case'], r['family'], r['upstream']): r for r in csv.DictReader((OUT/'experiment.csv').open())
                  if r['coordinate'] == 'table'}
    families = (('linear', 'Taps only, linear', '-'), ('monotone', 'Taps only, monotone cubic', (0, (4, 2))),
                ('additive', 'Model shape + tap differences', '-'), ('ratio', 'Model shape × tap ratios', (0, (4, 2))))
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.6), sharex=True)
    fig.patch.set_facecolor(SURFACE)
    for column, (case, title) in enumerate(CASES):
        top, bottom = axes[0, column], axes[1, column]; style(top); style(bottom)
        mx, mcp = z[f'native_{case}/x'], z[f'native_{case}/cp']; k = after_leading_edge(mx)
        top.plot(mx[k], mcp[k], color='#b5b3aa', lw=1.5, label='Model pressure (reference)')
        for (family, label, dash), color in zip(families, SERIES):
            key = f'exp_{case}_{family}_hold/'
            x, cp = z[key+'x'], z[key+'cp']; k = after_leading_edge(x)
            top.plot(x[k], cp[k], color=color, lw=2, ls=dash, label=label)
            fx = z[key+'fd_x']; k = after_leading_edge(fx)
            shear = z[key+'fd_cf']*z[key+'fd_speed']*z[key+'fd_theta']/NU
            bottom.plot(fx[k], shear[k], color=color, lw=2, ls=dash)
            record = experiment[(case, family, 'hold')]
            if record['fd_status'] == 'separated':
                position = float(record['fd_x_separation'])
                bottom.plot([position], [0], marker='v', ms=9, color=color, mec=SURFACE, mew=1.5, clip_on=False, zorder=5)
                top.plot([position], [np.interp(position, x[after_leading_edge(x)], cp[after_leading_edge(x)])], marker='v', ms=9,
                         color=color, mec=SURFACE, mew=1.5, zorder=5)
        points = np.array(measured[case][:4])
        top.plot(points[:, 0], points[:, 1], ls='none', marker='o', ms=8, mfc=SURFACE, mec=INK, mew=1.6, zorder=6,
                 label='Measured taps')
        top.invert_yaxis(); top.set_title(title, color=INK, fontsize=12, loc='left')
        bottom.axhline(0, color=SECONDARY, lw=1); bottom.set_ylim(-.03, .8)
        bottom.set_xlabel('x / total chord', color=SECONDARY, fontsize=10)
        bottom.text(.0338, .71, 'taps only:\nattached', color=INK, fontsize=9, ha='right', va='center')
        bottom.text(.0075, .055, 'model shape added:\nseparates (▼)', color=INK, fontsize=9, ha='left', va='center')
    axes[0, 0].set_ylabel('Reconstructed Cp', color=SECONDARY, fontsize=10)
    axes[1, 0].set_ylabel('Wall-shear parameter cf·Re_θ\n(0 = laminar separation)', color=SECONDARY, fontsize=10)
    axes[0, 0].set_xlim(0, .035)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .915), ncol=3, frameon=False, fontsize=10, labelcolor=INK)
    fig.suptitle('The measured taps do not decide whether the laminar layer separates', color=INK, fontsize=15, x=.065,
                 ha='left', y=.985)
    fig.text(.065, .935, 'Four ways to fill the gaps between taps (upstream treatment: hold), each with the finite-difference '
             'laminar solution.', color=SECONDARY, fontsize=10, ha='left')
    fig.tight_layout(rect=(0, 0, 1, .865))
    save(fig, 'measured_taps')


def main():
    z = np.load(OUT/'profiles.npz')
    branches(z); taps(z)
    print('Wrote two_branches and measured_taps (png, pdf, svg)')


if __name__ == '__main__':
    main()
