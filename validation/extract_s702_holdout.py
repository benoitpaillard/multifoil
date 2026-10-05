#!/usr/bin/env python3
"""Extract previously unused Re=700000, M=0.07 S702 references.

These three figures were reserved before choosing the sweep winner. They must
not be used to refit Ncrit after the transfer check. Drawing selectors and PDF
hash are pinned, and rendered source pages are retained for visual inspection.
"""
import argparse
import hashlib
import json
from pathlib import Path

from extract_s702 import DATA, PDF_SHA256, URL, vertices, write_csv

CASES = {
    0: dict(page=84, figure='23(b)', theory=[224, 225, 226],
            markers=[[251, 297], [314, 380], [399, 465]], cl_zero_path=41, cl_two_path=37),
    5: dict(page=89, figure='25(a)', theory=[224, 225, 226],
            markers=[[251, 295], [313, 375], [392, 454]], cl_zero_path=40, cl_two_path=28),
    10: dict(page=93, figure='26(b)', theory=[224, 225, 226],
             markers=[[251, 297], [315, 381], [401, 467]], cl_zero_path=40, cl_two_path=28),
}


def extract(out):
    import fitz
    pdf = DATA/'NASA-CR-20220014964.pdf'
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError('Source PDF changed; drawing selectors must be audited')
    out.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(pdf)
    cases = []
    for flap, cfg in CASES.items():
        page = doc[cfg['page']-1]
        if '0.70' not in page.get_text() or '0.07' not in page.get_text():
            raise ValueError('Caption conditions differ from the reserved case')
        drawings = page.get_drawings()
        cl0 = drawings[cfg['cl_zero_path']]['rect'].x0
        cl_scale = (cl0-drawings[cfg['cl_two_path']]['rect'].x0)/2.
        cd0 = drawings[0]['rect'].y0
        cd_scale = (cd0-drawings[66]['rect'].y0)/.02
        alpha0 = drawings[79]['rect'].y0
        alpha_scale = cl_scale/20.
        folder = out/f'd{flap:+d}'
        folder.mkdir(exist_ok=True)
        for curve, tid, bounds in zip(('drag', 'lift', 'moment'), cfg['theory'], cfg['markers']):
            for source in ('mses', 'experiment'):
                if source == 'mses':
                    raw = [(tid, *p) for p in vertices(drawings[tid])]
                else:
                    raw = []
                    for index in range(bounds[0], bounds[1]+1, 2):
                        drawing = drawings[index]
                        rect = drawing['rect']
                        if len(drawing['items']) != 24 or abs(rect.width-rect.height) > .003:
                            raise ValueError('Experimental selector is not a circular marker')
                        raw.append((index, .5*(rect.x0+rect.x1), .5*(rect.y0+rect.y1)))
                rows = []
                for index, px, py in raw:
                    row = dict(drawing_index=index, pdf_x=px, pdf_y=py)
                    if curve == 'drag':
                        row.update(CL=(cl0-px)/cl_scale, CD=(cd0-py)/cd_scale)
                    elif curve == 'lift':
                        row.update(alpha_deg=(alpha0-py)/alpha_scale, CL=(cl0-px)/cl_scale)
                    else:
                        row.update(alpha_deg=(alpha0-py)/alpha_scale, CM=(px-cl0)/(4*cl_scale))
                    rows.append(row)
                write_csv(folder/f'{source}_{curve}.csv', rows)
        page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(str(folder/'source_page.png'))
        cases.append(dict(flap=flap, **cfg, reynolds=700000., mach=.07,
            calibration=dict(cl_zero_pdf_x=cl0, points_per_cl=cl_scale,
                cd_zero_pdf_y=cd0, points_per_cd=cd_scale,
                alpha_zero_pdf_y=alpha0, points_per_degree=alpha_scale),
            caption=page.get_text()))
    source = dict(schema=1, reference='NASA/CR-20220014964', url=URL,
        pdf_sha256=PDF_SHA256, reynolds=700000., mach=.07, mses_ncrit=9.,
        transition='free', flaps=[0, 5, 10], targets=[.5, .8, 1.2], cases=cases,
        role='Reserved Reynolds-number transfer check; not a parameter-fitting dataset',
        limitations=['No minus-two-degree flap case is published at this Reynolds number.',
                    'Vector-figure values; experimental and reconstruction uncertainties are not supplied.',
                    'Published MSES N=9 results; no new MSES runs.'],
        sha256={str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(out.glob('d*/*'))})
    (out/'source.json').write_text(json.dumps(source, indent=2)+'\n')
    print(json.dumps(dict(cases=len(cases), out=str(out)), indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=DATA/'re700k_holdout')
    extract(ap.parse_args().out)
