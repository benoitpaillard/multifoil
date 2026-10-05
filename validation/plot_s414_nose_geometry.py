#!/usr/bin/env python3
"""Figures for the S414 fore-nose geometry trace and inviscid comparison (V4); reads its CSV outputs only."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_nose_geometry_20261004'
NOMINAL, CONVEX = '#2a78d6', '#eb6834'  # validated categorical slots 1 and 2; text stays in ink
SURFACE, INK, SECONDARY, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
NAMES = {'nominal': 'Nominal (production PCHIP)', 'convex_nose': 'Convex-nose alternative'}
CHORD_MM, THRESHOLD = 457.2, .05


def rows(name, out):
    with (out/name).open(newline='') as stream: return list(csv.DictReader(stream))


def style(ax):
    ax.set_facecolor(SURFACE); ax.grid(color=GRID, linewidth=.6); ax.set_axisbelow(True)
    for side in ('top', 'right'): ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'): ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8.5); ax.xaxis.label.set_color(SECONDARY); ax.yaxis.label.set_color(SECONDARY)


def save(fig, out, name):
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/f'{name}.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def distance(points, polyline):
    a, b = polyline[:-1], polyline[1:]; d = b-a; offset = points[:, None, :]-a[None]
    t = np.clip(np.sum(offset*d[None], axis=2)/np.sum(d*d, axis=1), 0, 1)
    return np.min(np.linalg.norm(offset-t[:, :, None]*d[None], axis=2), axis=1)


def geometry(out):
    curvature = rows('lower_nose_curvature.csv', out); nodes = rows('final_nodes_nose.csv', out)
    reported = [float(r['x']) for r in rows('reported_nose_vertices.csv', out) if r['surface'] == 'lower' and float(r['x']) <= .06]
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 8.6), sharex=True, layout='constrained'); fig.patch.set_facecolor(SURFACE)
    lower = {v: sorted(([float(r['x']), float(r['y']), float(r['turn_degrees'])] for r in nodes if r['variant'] == v and r['surface'] == 'lower'),
                       key=lambda p: p[0]) for v in NAMES}
    for variant, color in (('nominal', NOMINAL), ('convex_nose', CONVEX)):
        data = sorted((float(r['x']), float(r['curvature'])) for r in curvature if r['curve'] == variant)
        axes[0].plot(*zip(*data), color=color, linewidth=1.6, solid_capstyle='round', label=NAMES[variant])
        p = np.array(lower[variant])
        axes[1].plot(p[:, 0], p[:, 2], color=color, linewidth=1.2, marker='o', markersize=5.5, markeredgecolor=SURFACE, markeredgewidth=1.)
    axes[0].plot(reported, np.full(len(reported), .04), transform=axes[0].get_xaxis_transform(), linestyle='none', marker='|',
                 markersize=11, markeredgewidth=1.3, color=MUTED, clip_on=False, label='Reported vertices')
    worst = min(lower['nominal'], key=lambda p: p[2] if p[0] > .005 else 9)
    axes[1].annotate(f'{worst[2]:.2f}° at x/c {worst[0]:.4f}', (worst[0], worst[2]), xytext=(18, -4), textcoords='offset points',
                     fontsize=8.5, color=SECONDARY, va='center', arrowprops=dict(arrowstyle='-', color=MUTED, linewidth=.7))
    axes[0].annotate('concave', (.0158, -19), xytext=(.021, -17), fontsize=8.5, color=SECONDARY, va='center',
                     arrowprops=dict(arrowstyle='-', color=MUTED, linewidth=.7))
    for ax in axes[:2]: ax.axhline(0, color=AXIS, linewidth=.9)
    # Distances are taken to the whole nominal nose polygon (both surfaces), in node order.
    nominal = np.array([[float(r['x']), float(r['y'])] for r in sorted((r for r in nodes if r['variant'] == 'nominal'), key=lambda r: int(r['node']))])
    convex = np.array(lower['convex_nose'])[:, :2]
    axes[2].plot(convex[:, 0], distance(convex, nominal)*CHORD_MM, color=SECONDARY, linewidth=1.5, solid_capstyle='round')
    axes[2].axhline(.13, color=AXIS, linewidth=.9)
    axes[2].annotate('reported contour tolerance 0.13 mm', (.059, .13), xytext=(0, 3), textcoords='offset points', ha='right',
                     fontsize=8.5, color=SECONDARY)
    for ax, label in zip(axes, ('Curvature of the reconstructed\ncurve [1 / chord]', 'Turn at each final\nsolver node [deg]',
                                'Convex-nose nodes: distance\nfrom nominal polygon [mm]')):
        style(ax); ax.set_ylabel(label)
    axes[0].set_xlim(-.001, .06); axes[0].set_ylim(-30, 100); axes[2].set_ylim(0, .15)
    axes[2].set_xlabel('x / total chord, fore-element lower surface')
    # The production curve has a tangent break at the branch-split vertex; its curvature there is off scale.
    axes[0].annotate('tangent break at the\nbranch-split vertex (off scale)', (.0005, 98), xytext=(.006, 88), fontsize=8.5,
                     color=SECONDARY, va='center', arrowprops=dict(arrowstyle='-', color=MUTED, linewidth=.7))
    axes[0].set_title('Where the reversal comes from: the interpolated curve between sparse reported vertices', fontsize=10.5, color=INK, loc='left')
    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, frameon=False, fontsize=8.5, loc='outside lower center', ncols=3, labelcolor=SECONDARY)
    save(fig, out, 'nose_geometry_trace')


def pressure(out):
    data = rows('nose_pressure_profiles.csv', out)
    taps = [float(r['x']) for r in rows('reported_nose_vertices.csv', out)][:0]  # taps are drawn from the pressure table below
    cases = [('matched_cl0056', 'Matched CL 0.056 incidence'), ('matched_cl0175', 'Matched CL 0.175 incidence')]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[2.1, 1]), layout='constrained')
    fig.patch.set_facecolor(SURFACE); query = np.linspace(.0035, .06, 400)
    for column, (label, title) in enumerate(cases):
        curves = {}
        for variant, color in (('nominal', NOMINAL), ('convex_nose', CONVEX)):
            r = sorted(((float(v['x']), float(v['Cp'])) for v in data if v['variant'] == variant and v['incidence'] == label
                        and v['surface'] == 'lower' and float(v['x']) <= .062), key=lambda p: p[0])
            x, cp = map(np.array, zip(*r)); curves[variant] = (x, cp)
            axes[0, column].plot(x, cp, color=color, linewidth=1.6, solid_capstyle='round', label=NAMES[variant])
        alpha = next(float(v['alpha_deg']) for v in data if v['incidence'] == label)
        axes[0, column].invert_yaxis(); style(axes[0, column])
        axes[0, column].set_title(f'{title}, α = {alpha:.3f}°', fontsize=10, color=INK, loc='left')
        delta = np.interp(query, *curves['convex_nose'])-np.interp(query, *curves['nominal'])
        low = axes[1, column]; low.plot(query, delta, color=SECONDARY, linewidth=1.5, solid_capstyle='round'); style(low)
        low.axhline(0, color=AXIS, linewidth=.9)
        for level in (THRESHOLD, -THRESHOLD): low.axhline(level, color=MUTED, linewidth=.8)
        low.axvspan(.008, .06, color=GRID, alpha=.35, linewidth=0)
        low.set_xlabel('x / total chord, fore-element lower surface'); low.set_xlim(-.001, .06); low.set_ylim(-.13, .13)
        if column == 0:
            low.annotate('declared ±0.05 threshold', (.059, THRESHOLD), xytext=(0, 3), textcoords='offset points', ha='right', fontsize=8.5, color=SECONDARY)
            low.annotate('recovery window', (.034, -.115), ha='center', fontsize=8.5, color=SECONDARY)
    axes[0, 0].set_ylabel('Inviscid Cp at Mach 0.10 (suction upward)'); axes[1, 0].set_ylabel('Convex nose − nominal')
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, frameon=False, fontsize=8.5, loc='outside lower center', ncols=2, labelcolor=SECONDARY)
    fig.suptitle('Inviscid lower-nose pressure with and without the reversal (production operator, 340/240 panels)',
                 fontsize=11, color=INK, x=.01, ha='left')
    save(fig, out, 'nose_pressure_effect')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve(); geometry(out); pressure(out)


if __name__ == '__main__':
    main()
