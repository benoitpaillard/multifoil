#!/usr/bin/env python3
"""Extract S702 plotted geometry and references from NASA's vector PDF.

Requires PyMuPDF only for re-extraction. The distributed .dat/.csv files are
sufficient to run the benchmark. Drawing IDs refer to the pinned PDF below;
the coordinates are plotted vertices, not the author's original geometry deck.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'validation/data/s702'
PDF_SHA256 = '01b8bb6fc9b4fb2bc3246611106d64f980c03c41b18dcd1740a59a769e66cf03'
URL = 'https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf'
# PDF page numbers are one-based; drawing indices are zero-based, unrotated.
CASES = {
    0: dict(page=85, figure='23(c)', theory=[240, 241, 242],
            markers=[[267, 315], [335, 399], [421, 485]],
            cl_zero_path=49, cl_two_path=37, cd_twenty_path=74, alpha_zero_path=87),
    -2: dict(page=87, figure='24(a)', theory=[224, 225, 226],
             markers=[[251, 297], [315, 375], [393, 453]],
             cl_zero_path=41, cl_two_path=37, cd_twenty_path=66, alpha_zero_path=79),
    5: dict(page=90, figure='25(b)', theory=[224, 225, 226],
            markers=[[251, 293], [310, 366], [386, 442]],
            cl_zero_path=40, cl_two_path=28, cd_twenty_path=66, alpha_zero_path=79),
    10: dict(page=94, figure='26(c)', theory=[224, 225, 226],
             markers=[[251, 299], [320, 388], [408, 476]],
             cl_zero_path=40, cl_two_path=28, cd_twenty_path=66, alpha_zero_path=79),
}


def vertices(drawing):
    points = []
    for item in drawing['items']:
        if item[0] != 'l':
            raise ValueError('Expected a line-only published polyline')
        a, b = (list(p) for p in item[1:])
        if not points:
            points.append(a)
        elif max(abs(x-y) for x, y in zip(a, points[-1])) > .002:
            raise ValueError('Disconnected line segments in selected path')
        points.append(b)
    return points


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def extract(pdf, out):
    import fitz
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError('PDF hash changed: audit the drawing selectors before extraction')
    out.mkdir(parents=True, exist_ok=True)
    dest = out / 'NASA-CR-20220014964.pdf'
    if pdf.resolve() != dest.resolve():
        shutil.copyfile(pdf, dest)
    doc = fitz.open(pdf)
    paths = doc[26].get_drawings()
    chord = paths[2]['items'][0]
    z_zero = chord[1].x
    x_origin = max(chord[1].y, chord[2].y)
    scale = abs(chord[1].y-chord[2].y)
    geo = dict(pdf_page=27, figure='1', chord_path=2,
               x_origin_pdf_y=x_origin, z_zero_pdf_x=z_zero,
               points_per_reference_chord=scale,
               transformation='x/c=(x_origin_pdf_y-PDF_y)/scale; z/c=(z_zero_pdf_x-PDF_x)/scale',
               status='vector-figure reconstruction; original author coordinates not obtained',
               hinge_reported=[.9001, -.0114], elements=[])
    for name, index, angle in [('fore', 24, None), ('aft_d0', 0, 0),
                               ('aft_dm2', 3, -2), ('aft_d5', 10, 5), ('aft_d10', 17, 10)]:
        raw = vertices(paths[index])
        xy = [[(x_origin-y)/scale, (z_zero-x)/scale] for x, y in raw]
        with (out / (name+'.dat')).open('w') as f:
            f.write('S702 '+name+' - NASA Figure 1 vector reconstruction\n')
            f.writelines(f'{x:.16g} {y:.16g}\n' for x, y in xy)
        geo['elements'].append(dict(name=name, drawing_index=index,
            deflection_deg=angle, vertex_count=len(xy), pdf_vertices=raw,
            trailing_edge_upper=xy[0], trailing_edge_lower=xy[-1]))
    # The hinge cross is independently checkable against the printed coordinates.
    hinge = paths[5]['items'][0][1]
    geo['hinge_from_cross'] = [(x_origin-hinge.y)/scale, (z_zero-hinge.x)/scale]
    (out / 'geometry_extraction.json').write_text(json.dumps(geo, indent=2)+'\n')
    refs = []
    for angle, cfg in CASES.items():
        ds = doc[cfg['page']-1].get_drawings()
        cl0 = ds[cfg['cl_zero_path']]['rect'].x0
        cl_scale = (cl0-ds[cfg['cl_two_path']]['rect'].x0)/2
        cd0 = ds[0]['rect'].y0
        cd_scale = (cd0-ds[cfg['cd_twenty_path']]['rect'].y0)/.02
        alpha0 = ds[cfg['alpha_zero_path']]['rect'].y0
        # The alpha axis spans 20 degrees over the same distance as one CL unit.
        alpha_scale = cl_scale/20
        calibration = dict(cl_zero_pdf_x=cl0, points_per_cl=cl_scale,
            cd_zero_pdf_y=cd0, points_per_cd=cd_scale,
            alpha_zero_pdf_y=alpha0, points_per_degree=alpha_scale,
            points_per_cm=4*cl_scale)
        folder = out / f'd{angle:+d}'
        folder.mkdir(exist_ok=True)
        for curve, tid, bounds in zip(['drag', 'lift', 'moment'], cfg['theory'], cfg['markers']):
            for kind in ['mses', 'experiment']:
                if kind == 'mses':
                    raw = [(tid, *p) for p in vertices(ds[tid])]
                else:
                    raw = []
                    for idx in range(bounds[0], bounds[1]+1, 2):
                        d = ds[idx]
                        r = d['rect']
                        if len(d['items']) != 24 or abs(r.width-r.height) > .003:
                            raise ValueError('Selected experimental symbol is not a circular marker')
                        raw.append((idx, .5*(r.x0+r.x1), .5*(r.y0+r.y1)))
                rows = []
                for idx, px, py in raw:
                    row = dict(drawing_index=idx, pdf_x=px, pdf_y=py)
                    if curve == 'drag':
                        row.update(CL=(cl0-px)/cl_scale, CD=(cd0-py)/cd_scale)
                    elif curve == 'lift':
                        row.update(alpha_deg=(alpha0-py)/alpha_scale, CL=(cl0-px)/cl_scale)
                    else:
                        row.update(alpha_deg=(alpha0-py)/alpha_scale, CM=(px-cl0)/(4*cl_scale))
                    rows.append(row)
                write_csv(folder / f'{kind}_{curve}.csv', rows)
        refs.append(dict(deflection_deg=angle, **cfg, calibration=calibration))
    source = dict(schema_version=1,
        title='Design and Experimental Results for a Low-Speed, Flapped, Slotted, Natural-Laminar-Flow Airfoil',
        authors=['Dan M. Somers', 'Mark D. Maughmer'], year=2022,
        report='NASA/CR-20220014964', url=URL, pdf_sha256=PDF_SHA256,
        checked_utc='2026-09-08', geometry=geo['status'], reference_cases=refs,
        conditions=dict(reynolds=1e6, mach=.10, ncrit=9., transition='free',
            reference_chord_m=.45748, moment_reference=[.25, 0.],
            incidence_reference='original report x axis; aft trailing edge lies below chord line'),
        qualifications=[
            'MSES curves are published results extracted from the PDF; MSES was not run here.',
            'Original author geometry decks and tunnel raw data have not been obtained.',
            'Vector centerlines reduce reading error; they do not bound source or model-form uncertainty.',
            'The plastic-film hinge cover and measured model deviations are not represented.',
            'Only natural-transition cases are used; forced trips are not approximated by low Ncrit.',
            'MSES uses a subcritical Euler formulation; MultiFoil uses potential flow with a compressibility correction.',
            'These are fixed-configuration polars, not measured local contour-design derivatives.'],
        extraction_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'source.json').write_text(json.dumps(source, indent=2)+'\n')
    print(json.dumps(dict(geometry_vertices=[e['vertex_count'] for e in geo['elements']],
        hinge=geo['hinge_from_cross'], cases=len(refs)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', type=Path, default=DATA/'NASA-CR-20220014964.pdf')
    parser.add_argument('--out', type=Path, default=DATA)
    args = parser.parse_args()
    extract(args.pdf, args.out)
