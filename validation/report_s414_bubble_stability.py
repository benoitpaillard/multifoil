#!/usr/bin/env python3
"""Figures for V5 tranche 2 from saved outputs only; no boundary-layer or stability calculation."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_s414_laminar_replay import SURFACE, INK, SECONDARY, SERIES, style, after_leading_edge

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_bubble_stability_20261005'
CASES = (('cl0056', 'CL 0.056'), ('cl0175', 'CL 0.175'), ('cl028', 'CL 0.284 (no bubble, control)'))
X_END = .02082853867960402
ARCHIVE, INHERITED, INDEPENDENT = SERIES[0], SERIES[1], SERIES[2]


def table(name):
    with (OUT/name).open(newline='') as stream:
        return list(csv.DictReader(stream))


def save(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(OUT/f'{name}.{ext}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def closure(z, nodes):
    """Archived closure against the independent boundary layer on the same displacement thickness."""
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 9.4), sharex=True)
    fig.patch.set_facecolor(SURFACE)
    for column, (case, title) in enumerate(CASES):
        rows = [r for r in nodes if r['case'] == case and r['variant'] == 'primary']
        x = np.array([float(r['x']) for r in rows])
        switch = max(float(r['x']) for r in rows if r['mode'] == 'direct')
        fx = z[f'{case}/primary/x']; k = after_leading_edge(fx)
        keep = fx[k] >= x[0]
        for row, (archive, independent, key) in enumerate((('Cp_archive', 'Cp_independent', 'cp'), ('H_archive', 'H_independent', 'H'),
                ('cf_Re_theta_archive', 'cf_Re_theta_independent', 'cf_re_theta'))):
            ax = axes[row, column]; style(ax)
            ax.plot(x, [float(r[archive]) for r in rows], color=ARCHIVE, lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5,
                    label='Archived state (inherited closure)')
            ax.plot(fx[k][keep], z[f'{case}/primary/{key}'][k][keep], color=INDEPENDENT, lw=2,
                    label='Independent finite-difference layer')
            ax.axvline(switch, color=SECONDARY, lw=1)
        axes[0, column].invert_yaxis(); axes[0, column].set_title(title, color=INK, fontsize=12, loc='left')
        axes[2, column].axhline(0, color=SECONDARY, lw=1)
        axes[2, column].set_xlabel('x / total chord', color=SECONDARY, fontsize=10)
        bottom, top = axes[0, column].get_ylim()  # inverted axis: the first value is at the bottom
        axes[0, column].annotate('same\npressure', (switch, bottom), xytext=(-6, 6), textcoords='offset points',
                                 color=SECONDARY, fontsize=9, va='bottom', ha='right')
        axes[0, column].annotate('same displacement\nthickness downstream', (switch, top), xytext=(6, -4),
                                 textcoords='offset points', color=SECONDARY, fontsize=9, va='top', ha='left')
    axes[0, 0].set_ylabel('Cp', color=SECONDARY, fontsize=10)
    axes[1, 0].set_ylabel('Shape factor H', color=SECONDARY, fontsize=10)
    axes[2, 0].set_ylabel('Wall-shear parameter cf·Re_θ', color=SECONDARY, fontsize=10)
    axes[0, 0].set_xlim(.004, .0215)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .925), ncol=2, frameon=False, fontsize=10, labelcolor=INK)
    fig.suptitle('The nose layer recomputed with independent laminar equations', color=INK, fontsize=15, x=.065, ha='left', y=.985)
    fig.text(.065, .94, 'Fore-lower nose. Downstream of the vertical line the archived displacement thickness is prescribed and '
             'the pressure is solved.', color=SECONDARY, fontsize=10, ha='left')
    fig.tight_layout(rect=(0, 0, 1, .895))
    save(fig, 'bubble_closure')


def amplification(z):
    """Archived N, the inherited law on the independent state, and Orr-Sommerfeld N."""
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.9), sharex=True, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (case, title) in zip(axes, CASES):
        style(ax)
        x, state, turbulent = z[f'{case}/archive_x'], z[f'{case}/archive_U'], z[f'{case}/archive_turbulent']
        k = np.arange(len(x)) >= int(np.argmin(x)); k &= (x >= .004) & (x <= X_END+1e-9) & (turbulent == 0)
        ax.plot(x[k], state[2][k], color=ARCHIVE, lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5, label='Archived state')
        ix, inherited = z[f'{case}/primary/inherited_x'], z[f'{case}/primary/inherited_N']
        m = np.arange(len(ix)) >= int(np.argmin(ix)); m &= ix >= .004
        ax.plot(ix[m], inherited[m], color=INHERITED, lw=2, label='Inherited growth law on the independent layer')
        sx, stability = z[f'{case}/primary/stability_x'], z[f'{case}/primary/stability_N']
        ax.plot(sx, stability, color=INDEPENDENT, lw=2, label='Linear stability (Orr–Sommerfeld) on the independent layer')
        ax.axhline(9, color=SECONDARY, lw=1)
        ax.text(.0046, 9.2, 'Ncrit = 9', color=SECONDARY, fontsize=9, va='bottom')
        ends = sorted(((state[2][k][-1], 'archive'), (inherited[m][-1], 'inherited'), (stability[-1], 'stability')))
        placed = []
        for value, _ in ends:
            position = value if not placed or value-placed[-1] > .55 else placed[-1]+.55
            placed.append(position)
            ax.annotate(f'{value:.1f}', (X_END, position), xytext=(7, 0), textcoords='offset points', color=INK, fontsize=9,
                        va='center')
        ax.set_title(title, color=INK, fontsize=12, loc='left')
        ax.set_xlabel('x / total chord', color=SECONDARY, fontsize=10)
    axes[0].set_ylabel('Amplification N', color=SECONDARY, fontsize=10)
    axes[0].set_xlim(.004, .0228); axes[0].set_ylim(-.4, 11); axes[0].set_xticks([.005, .010, .015, .020])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .875), ncol=3, frameon=False, fontsize=9.5, labelcolor=INK)
    return fig


def similar(rows):
    """e^N on similar profiles: inherited correlation against Orr-Sommerfeld."""
    names = list(dict.fromkeys(r['profile'] for r in rows))
    fig, axes = plt.subplots(1, len(names), figsize=(12.5, 3.9))
    fig.patch.set_facecolor(SURFACE)
    for ax, name in zip(axes, names):
        style(ax)
        r = [q for q in rows if q['profile'] == name]
        re = [float(q['re_theta']) for q in r]
        ax.plot(re, [float(q['N_inherited']) for q in r], color=INHERITED, lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5,
                label='Inherited correlation')
        ax.plot(re, [float(q['N_orr_sommerfeld']) for q in r], color=INDEPENDENT, lw=2, marker='o', ms=6.5, mec=SURFACE, mew=1.5,
                label='Linear stability (Orr–Sommerfeld)')
        ax.set_title(f'H = {float(r[0]["H"]):.1f}'+(' (Blasius)' if name == 'Blasius' else ''), color=INK, fontsize=11, loc='left')
        ax.set_xlabel('Re_θ', color=SECONDARY, fontsize=10); ax.set_ylim(bottom=0)
    axes[0].set_ylabel('Amplification N', color=SECONDARY, fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .86), ncol=2, frameon=False, fontsize=10, labelcolor=INK)
    return fig


def main():
    titles = dict(amplification='Linear stability amplifies slightly more than the inherited growth law in the bubble',
                  similar='On similar profiles the inherited correlation grows faster than linear stability as H rises')
    z = np.load(OUT/'profiles.npz')
    closure(z, table('nodes.csv'))
    fig = amplification(z)
    fig.suptitle(titles['amplification'], color=INK, fontsize=15, x=.065, ha='left', y=.985)
    fig.text(.065, .905, 'Fore-lower nose, same displacement thickness in all three curves downstream of the switch station.',
             color=SECONDARY, fontsize=10, ha='left')
    fig.tight_layout(rect=(0, 0, 1, .82)); save(fig, 'bubble_amplification')
    fig = similar(table('similar_profiles.csv'))
    fig.suptitle(titles['similar'], color=INK, fontsize=15, x=.065, ha='left', y=.985)
    fig.tight_layout(rect=(0, 0, 1, .82)); save(fig, 'similar_profiles')
    print('Wrote bubble_closure, bubble_amplification and similar_profiles (png, pdf, svg)')


if __name__ == '__main__':
    main()
