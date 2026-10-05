#!/usr/bin/env python3
"""Render source overlays and numerical geometry checks for S702 extraction.

Optional dependency: PyMuPDF. Does not modify the archived source PDF.
"""
import json
import contextlib
import io
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitz
from s702_geometry import geometry, DATA
from optimization import minimum_clearance
from extract_s702 import extract


def main():
    with tempfile.TemporaryDirectory() as name:
        regenerated = Path(name)
        with contextlib.redirect_stdout(io.StringIO()):
            extract(DATA/'NASA-CR-20220014964.pdf', regenerated)
        count = 0
        for path in regenerated.rglob('*'):
            if path.is_file():
                if path.read_bytes() != (DATA/path.relative_to(regenerated)).read_bytes():
                    raise ValueError('Distributed extraction differs: '+str(path.relative_to(regenerated)))
                count += 1
    doc = fitz.open(DATA/'NASA-CR-20220014964.pdf')
    source = json.loads((DATA/'source.json').read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.2))
    for ax, angle in zip(axes.flat, [-2, 0, 5, 10]):
        case = next(c for c in source['reference_cases'] if c['deflection_deg'] == angle)
        page = doc[case['page']-1]
        pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
        ax.imshow(np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3))
        for kind, color in [('experiment', '#00866b'), ('mses', '#dc7625')]:
            data = np.genfromtxt(DATA/f'd{angle:+d}'/f'{kind}_drag.csv', delimiter=',', names=True)
            if kind == 'experiment':
                ax.scatter((792-data['pdf_y'])*2, data['pdf_x']*2,
                           s=45, marker='+', color=color, linewidths=1.3)
            else:
                ax.plot((792-data['pdf_y'])*2, data['pdf_x']*2, color=color, lw=1)
        c = case['calibration']
        ax.set_xlim((792-c['cd_zero_pdf_y']+.002*c['points_per_cd'])*2,
                    (792-c['cd_zero_pdf_y']+.019*c['points_per_cd'])*2)
        ax.set_ylim((c['cl_zero_pdf_x']+.04*c['points_per_cl'])*2,
                    (c['cl_zero_pdf_x']-1.55*c['points_per_cl'])*2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"Figure {case['figure']}, flap {angle:+d} degrees")
    fig.suptitle('Source overlay: extracted MSES path (orange), experimental centers (green)', fontsize=13)
    fig.text(.5, .01, 'Original PDF raster behind the extracted vectors. Drawing IDs and transforms are archived in source.json.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .025, 1, .965])
    fig.savefig(DATA/'extraction_audit.png', dpi=170)
    plt.close(fig)
    fig = plt.figure(figsize=(11, 5.8))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.4])
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    audits = []
    baseline = geometry(0)
    hinge = np.array([.9001, -.0114])
    for angle, color in zip([-2, 0, 5, 10], ['#666666', '#2166ac', '#00866b', '#c66b20']):
        elements = geometry(angle)
        for ax in axes:
            if angle == 0:
                ax.plot(*elements[0][1].T, color='#20262c', lw=1.5)
            ax.plot(*elements[1][1].T, color=color, lw=1.1, label=f'{angle:+d} deg')
        a = np.deg2rad(-angle)
        r = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        te = (baseline[1][1][[0, -1]]-hinge) @ r.T+hinge
        drag = np.genfromtxt(DATA/f'd{angle:+d}'/'experiment_drag.csv', delimiter=',', names=True)
        lift = np.genfromtxt(DATA/f'd{angle:+d}'/'experiment_lift.csv', delimiter=',', names=True)
        audits.append(dict(deflection_deg=angle,
            minimum_clearance_c=minimum_clearance(elements[0][1], elements[1][1]),
            te_rotation_residual_c=float(np.max(np.linalg.norm(te-elements[1][1][[0, -1]], axis=1))),
            max_drag_vs_lift_marker_CL_difference=float(np.max(abs(drag['CL']-lift['CL'][:len(drag)])))))
    for ax in axes:
        ax.set_aspect('equal')
        ax.set_xlabel('$x/c$')
        ax.set_ylabel('$z/c$')
        ax.grid(alpha=.15)
    axes[0].set_xlim(-.02, 1.02)
    axes[0].set_ylim(-.065, .10)
    axes[0].legend(ncol=4, loc='upper right', fontsize=9)
    axes[1].set_xlim(.67, .84)
    axes[1].set_ylim(-.045, .065)
    axes[1].set_title('Fixed slot and aft-element nose')
    axes[2].set_xlim(.855, 1.015)
    axes[2].set_ylim(-.04, .025)
    axes[2].scatter(*hinge, marker='+', c='red', s=80)
    axes[2].set_title('Small aft flap about reported hinge')
    fig.suptitle('S702 geometry recovered from Figure 1 vector paths', fontsize=15)
    fig.tight_layout(rect=[0, .025, 1, .94])
    fig.text(.5, .01, 'Original common reference frame, finite trailing edges and all four drawn configurations retained.',
             ha='center', fontsize=9)
    fig.savefig(DATA/'geometry_audit.png', dpi=180)
    plt.close(fig)
    f = baseline[0][1]
    le = np.argmin(f[:, 0])
    x = np.linspace(0, .80, 20001)
    thickness = np.interp(x, f[:le+1, 0][::-1], f[:le+1, 1][::-1])-np.interp(x, f[le:, 0], f[le:, 1])
    audit = dict(cases=audits, reported_thickness_ratio=.1291,
        max_vertical_linear_interpolated_fore_thickness=float(thickness.max()),
        thickness_qualification='This plotted vertical thickness is not used to rescale the geometry; original author coordinates are still needed for exact matching.',
        reference_roundtrip=dict(files_compared=count, all_byte_identical=True),
        source_pdf_sha256=source['pdf_sha256'])
    (DATA/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
