#!/usr/bin/env python3
"""Qualify S414 source geometry; no flow solve and no aerodynamic fitting.

Use /tmp/multifoil-pressure-extraction-env/bin/python in the present workspace.
Run acquire_s414_sources.py first. Sources and selectors are SHA-256 pinned.
"""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from multifoil import validate_geometry
from extract_s702 import vertices

OUT = ROOT / 'validation/results/s414_geometry_qualification_20261002'
PINNED = {
    's414_report_2010.pdf': '52756797afa2f7089be9e7d2ac807fff024b77769751b84b81ad7d6312fc0c5c',
    'archive_s414_report_2010.pdf': 'd9c37324cf4492ba6eac9b8f26c1b107f11b9308878489e067c932add5ea6799',
    's414.pts': '6b4a1f1a5fbdfc3efeaeba2fa7129fed251a52a19c4095b555f4577c168afce3',
    's414-1.pts': '1e03af896c83ee65260397ac6a0b6741f4bc8baf58f836fb4191a4cc53a23146',
}
REPORT_URL = 'https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_pts(path):
    """Preserve all typed boundary blocks; a .pts file is not just an airfoil."""
    lines = Path(path).read_text().splitlines()
    blocks = []
    for i, line in enumerate(lines):
        if line.strip() != 'NEWBND':
            continue
        meta = {}
        j = i + 1
        while lines[j].strip() != 'BNDEXY':
            meta[lines[j].strip()] = int(lines[j + 1].strip())
            j += 2
        count = meta['NRBNDE']
        points = np.array([[float(x) for x in s.split()] for s in lines[j+1:j+1+count]])
        require(points.shape == (count, 2) and np.isfinite(points).all(), 'Invalid .pts block')
        blocks.append(dict(metadata=meta, first_coordinate_line=j+2, points=points))
    require(blocks, 'No .pts blocks found')
    return blocks


def comparison_geometry(blocks, name):
    if name == 's414.pts':
        fore = np.vstack((blocks[0]['points'][::-1], blocks[1]['points'][::-1][1:]))
        aft = np.vstack((blocks[2]['points'][::-1], blocks[3]['points'][::-1][1:]))
    elif name == 's414-1.pts':
        fore, aft = blocks[1]['points'], blocks[2]['points']
    else:
        raise ValueError('Unknown public file selector')
    # Only undo apparent fore-unit scaling. Never move the elements separately.
    scale = float(aft[:, 0].max())
    return {'fore': fore / scale, 'aft': aft / scale}, scale


def point_distances(points, contour):
    """Distance to the closed piecewise-linear contour, ignoring duplicate closure."""
    a = np.asarray(contour)
    delta = np.roll(a, -1, axis=0) - a
    norm2 = np.sum(delta**2, axis=1)
    keep = norm2 > 1e-28
    a, delta, norm2 = a[keep], delta[keep], norm2[keep]
    require(len(a) >= 3, 'Degenerate contour')
    diff = np.asarray(points)[:, None, :] - a[None, :, :]
    t = np.clip(np.sum(diff*delta, axis=2)/norm2, 0, 1)
    return np.sqrt(np.min(np.sum((diff-t[:, :, None]*delta)**2, axis=2), axis=1))


def symmetric_vertex_distance(a, b):
    # This is a sampled diagnostic, not the exact continuous Hausdorff distance.
    return float(max(point_distances(a, b).max(), point_distances(b, a).max()))


def lower_branch(contour):
    return contour[np.argmin(contour[:, 0]):]


def geometry_metrics(elements):
    fore, aft = elements['fore'], elements['aft']
    validate_geometry(list(elements.items()))
    tip = (fore[0] + fore[-1]) / 2
    return dict(nonintersecting=True, fore_vertices=len(fore), aft_vertices=len(aft),
                fore_min_x=float(fore[:, 0].min()), aft_min_x=float(aft[:, 0].min()),
                fore_TE_midpoint=tip.tolist(), aft_TE_midpoint=((aft[0]+aft[-1])/2).tolist(),
                fore_TE_thickness=float(np.linalg.norm(fore[0]-fore[-1])),
                aft_TE_thickness=float(np.linalg.norm(aft[0]-aft[-1])),
                fore_TE_to_aft_surface_gap=float(point_distances([tip], aft)[0]),
                projected_overlap=float(tip[0]-aft[:, 0].min()))


def extract_taps(doc):
    rows = []
    for page, element, counts in [(19, 'fore', (25, 35)), (20, 'aft', (20, 10))]:
        columns = (142.5, 358.5) if element == 'fore' else (124.5, 340.5)
        minimum_y = 200 if element == 'fore' else 175
        for side, column, count in zip(('upper', 'lower'), columns, counts):
            words = [w for w in doc[page-1].get_text('words')
                     if column-1 <= w[0] <= column+10 and w[1] > minimum_y
                     and re.fullmatch(r'(?:[01])?\.\d{5}', w[4])]
            words.sort(key=lambda w: w[1])
            require(len(words) == count, f'Table II count changed: {element}/{side}')
            for index, w in enumerate(words):
                rows.append(dict(element=element, surface=side, index=index, x_over_c=float(w[4]),
                                 printed_value=w[4], pdf_page=page, pdf_x=w[0], pdf_y=w[1]))
    return rows


def save_geometry(prefix, elements):
    for name, xy in elements.items():
        np.savetxt(OUT/f'{prefix}_{name}.dat', xy, fmt='%.16g', comments='',
                   header=f'S414 {prefix} {name}; figure/public-file reconstruction; not author/as-built coordinates')


def main():
    import fitz
    for name, digest in PINNED.items():
        require(sha(OUT/'sources'/name) == digest, f'Changed source: {name}')
    manifest = json.loads((OUT/'source_manifest.json').read_text())
    for record in manifest['sources']:
        if record.get('path'):
            require(sha(ROOT/record['path']) == record['sha256'], 'Acquisition receipt mismatch')
    doc = fitz.open(OUT/'sources/s414_report_2010.pdf')
    archive = fitz.open(OUT/'sources/archive_s414_report_2010.pdf')
    require(len(doc) == len(archive) == 97, 'Report page count changed')
    metadata = json.loads((OUT/'sources/archive_metadata.json').read_text())
    supplied = next(f for f in metadata['files'] if f['name'] == 'DTIC_ADA545821.pdf')
    body = (OUT/'sources/archive_s414_report_2010.pdf').read_bytes()
    require(hashlib.sha1(body).hexdigest() == supplied['sha1'], 'Archive supplied SHA-1 differs')
    require(hashlib.md5(body).hexdigest() == supplied['md5'], 'Archive supplied MD5 differs')
    ds = doc[21].get_drawings()
    require(ds == archive[21].get_drawings(), 'Independent report copies differ in Figure 1')
    a, b = ds[2]['items'][0][1:]
    scale = b.y - a.y
    require(scale > 600 and a.x == b.x, 'Figure 1 chord calibration changed')
    provenance = dict(figure='1', pdf_page=22, printed_page=20, chord_path=2,
                      reference_chord_mm=457.2, origin_pdf_y=b.y, zero_z_pdf_x=a.x,
                      points_per_chord=scale,
                      transform='x/c=(origin_pdf_y-PDF_y)/scale; z/c=(zero_z_pdf_x-PDF_x)/scale',
                      as_built_coordinates=False, original_author_coordinates=False, drawings=[])
    shapes = {}
    for label, i in [('aft', 0), ('design_fore', 3), ('test_fore', 9)]:
        raw = np.array(vertices(ds[i]))
        xy = np.column_stack(((b.y-raw[:, 1])/scale, (a.x-raw[:, 0])/scale))
        shapes[label] = xy
        provenance['drawings'].append(dict(label=label, drawing_index=i, vertices=len(xy),
                                           raw_pdf_vertices=raw.tolist()))
    require(np.array_equal(np.array(vertices(ds[0])), np.array(vertices(ds[6]))),
            'Design and tested aft elements differ')
    tested = dict(fore=shapes['test_fore'], aft=shapes['aft'])
    design = dict(fore=shapes['design_fore'], aft=shapes['aft'])
    metrics = {'report_test': geometry_metrics(tested), 'report_design': geometry_metrics(design)}
    save_geometry('report_test', tested)
    save_geometry('report_design', design)
    lo_design, lo_test = lower_branch(design['fore']), lower_branch(tested['fore'])
    x = np.linspace(.4, .68, 2801)
    difference = np.interp(x, lo_test[:, 0], lo_test[:, 1])-np.interp(x, lo_design[:, 0], lo_design[:, 1])
    k = int(np.argmax(difference))
    modification = dict(max_lower_surface_thinning_chord=float(difference[k]), x_over_c=float(x[k]),
                        max_lower_surface_thinning_mm=float(difference[k]*457.2),
                        interpretation='Vertical separation between published piecewise-linear outlines; not an as-built measurement')
    # Independent x-axis calibration; one common z translation from fore vertices.
    repeated = []
    for page, x0path, x1path in [(34, 33, 47), (38, 33, 47), (66, 26, 40), (67, 40, 54)]:
        paths = doc[page-1].get_drawings()
        x0, x1 = paths[x0path]['rect'].x1, paths[x1path]['rect'].x1
        sc = x1-x0
        raw_fore = np.array(vertices(paths[1]))
        require(raw_fore.shape == tested['fore'].shape, 'Repeated fore contour count changed')
        offset = float(np.mean(tested['fore'][:, 1]+raw_fore[:, 1]/sc))
        contours = {}
        for name, index in [('fore', 1), ('aft', 0)]:
            raw = np.array(vertices(paths[index]))
            contours[name] = np.column_stack(((raw[:, 0]-x0)/sc, offset-raw[:, 1]/sc))
        err = symmetric_vertex_distance(contours['aft'], tested['aft'])
        fore_err = float(abs(contours['fore']-tested['fore']).max())
        require(fore_err < 2e-5 and err < 2e-5, 'Repeated test outline mismatch: check axes/selectors')
        repeated.append(dict(pdf_page=page, x0_path=x0path, x1_path=x1path, x0=x0, scale=sc,
                             common_z_offset=offset, fore_max_coordinate_difference=fore_err,
                             aft_max_vertex_to_segment_distance=err,
                             fore_vertices=len(contours['fore']), aft_vertices=len(contours['aft'])))
    provenance['pressure_figure_crosschecks'] = repeated
    candidates = {}
    public_comparisons = []
    pts_info = []
    for name in ['s414.pts', 's414-1.pts']:
        blocks = parse_pts(OUT/'sources'/name)
        elems, sc = comparison_geometry(blocks, name)
        candidates[name] = elems
        metrics[name] = geometry_metrics(elems)
        save_geometry(name.replace('.pts', '').replace('-', '_')+'_public', elems)
        pts_info.append(dict(file=name, scale_divisor=sc,
                             transform='Both coordinates divided by maximum aft x; no translation, rotation or separate placement adjustment',
                             author_origin_and_provenance='Not established by the public teaching catalog',
                             blocks=[dict(b['metadata'], first_coordinate_line=b['first_coordinate_line']) for b in blocks]))
        for element in ('fore', 'aft'):
            public_comparisons.append(dict(file=name, element=element,
                symmetric_max_vertex_to_test_segments=symmetric_vertex_distance(elems[element], tested[element]),
                symmetric_max_vertex_to_design_segments=symmetric_vertex_distance(elems[element], design[element]),
                maximum_vertex_distance_mm=symmetric_vertex_distance(elems[element], tested[element])*457.2,
                TE_thickness=metrics[name][element+'_TE_thickness']))
    taps = extract_taps(doc)
    for row in taps:
        c = tested[row['element']]
        i = int(np.argmin(c[:, 0]))
        branch = c[:i+1] if row['surface'] == 'upper' else c[i:]
        low, high = float(branch[:, 0].min()), float(branch[:, 0].max())
        excess = max(low-row['x_over_c'], row['x_over_c']-high, 0.)
        row.update(source_branch_xmin=low, source_branch_xmax=high,
                   unsupported_center=int(excess > 1e-12), outside_extent=excess,
                   outside_extent_mm=excess*457.2)
    write_csv(OUT/'pressure_tap_geometry_support.csv', taps)
    write_csv(OUT/'public_geometry_comparison.csv', public_comparisons)
    write_json(OUT/'geometry_extraction.json', provenance)
    write_json(OUT/'public_pts_provenance.json', pts_info)
    write_json(OUT/'geometry_metrics.json', metrics)
    # Save scoped primary evidence and readable renderings, not full-article quotations.
    for page in (19, 20, 22, 66, 67):
        pg = doc[page-1]
        pg.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(OUT/f'source_page_{page}.png')
    summary = dict(source_geometry_reconstruction_pass=True, original_author_coordinates_recovered=False,
                   tested_geometry_benchmark_admitted=False, new_flow_solves=0,
                   status='Nominal tested figure reconstruction qualified; exact/as-built benchmark admission remains open',
                   source_copies_geometry_identical=True, archive_supplied_checksums_pass=True,
                   vertices=dict(fore=len(tested['fore']), aft=len(tested['aft'])),
                   geometry_nonintersecting=True, modification=modification,
                   repeated_outline_checks=repeated, public_file_comparisons=public_comparisons,
                   tap_count=len(taps), unsupported_tap_centers=[t for t in taps if t['unsupported_center']],
                   reported_model_contour_tolerance_mm=.13, reported_tap_diameter_mm=.51,
                   qualification='Matching figure paths do not establish exact manufactured contour or original MSES input. Public .pts provenance/registration is undocumented and both close the trailing edges to points.',
                   next_step='Figure 12 pressure extraction is provided separately; Figure 13 drag extraction remains. Obtain original modified tested contour/placement before strict matched validation. Any reconstruction sensitivity study must be explicitly labeled exploratory.',
                   pymupdf=fitz.VersionBind)
    write_json(OUT/'summary.json', summary)
    plot(tested, design, candidates, modification)
    inputs = {str((OUT/'sources'/n).relative_to(ROOT)):sha(OUT/'sources'/n) for n in PINNED}
    for path in (Path(__file__), ROOT/'validation/extract_s702.py', ROOT/'multifoil.py', OUT/'source_manifest.json'):
        inputs[str(path.relative_to(ROOT))] = sha(path)
    write_json(OUT/'input_manifest.json', dict(input_sha256=inputs, outputs_sha256={
        p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.suffix in ('.dat','.csv','.png','.pdf','.svg','.json')
        and p.name not in ('input_manifest.json','verification.json')}))
    print(json.dumps(summary, indent=2))


def plot(tested, design, candidates, modification):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    a = axes[0, 0]
    a.plot(*design['fore'].T, color='#c56d27', label='Rotorcraft design (Figure 1)')
    for i, xy in enumerate(tested.values()):
        a.plot(*xy.T, color='#185e85', label='Modified tunnel shape' if i == 0 else None)
    a.axhline(0, color='gray', lw=.5)
    a.set(title='Published design and tested geometry', xlabel='x / complete chord', ylabel='z / complete chord')
    a.set_aspect('equal', adjustable='box');a.legend(fontsize=8, loc='lower center', bbox_to_anchor=(.5,1.6))
    a = axes[0, 1]
    for xy, color, label in [(design['fore'], '#c56d27', 'Design'), (tested['fore'], '#185e85', 'Test')]:
        branch = lower_branch(xy)
        a.plot(*branch.T, color=color, label=label)
    a.set(xlim=(.4,.69), ylim=(-.061,.019), title=f"Fore lower-surface change: up to {modification['max_lower_surface_thinning_mm']:.2f} mm", xlabel='x / complete chord', ylabel='z / complete chord')
    a.legend(fontsize=8)
    for a, element in zip(axes[1], ('fore', 'aft')):
        a.plot(*tested[element].T, color='black', lw=1.7, label='Report tested outline')
        for (name, elems), color in zip(candidates.items(), ('#9462a8', '#259982')):
            a.plot(*elems[element].T, color=color, ls='--', lw=1.2, label=name+' (global rescale only)')
        a.set(title=f'Public coordinate candidates: {element}', xlabel='x / complete chord', ylabel='z / complete chord')
        a.set_aspect('equal', adjustable='box');a.legend(fontsize=8)
    for a in axes.flat:a.grid(alpha=.2)
    fig.suptitle('S414 source-geometry qualification\nFigure reconstruction and public files; no aerodynamic solve or fitted placement')
    for ext in ('png','svg','pdf'):fig.savefig(OUT/f'geometry_qualification.{ext}', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
