#!/usr/bin/env python3
"""Pinned vector extraction of the independent Figure 22(b) Cp check.

The report labels the common lift as 1.04 (two decimals). No exact underlying
experimental lift or original pressure file is claimed. Drawing selectors and
axis calibration are retained for visual and numerical audits.
"""
import argparse
import hashlib
import json
from pathlib import Path

from extract_s702 import DATA, PDF_SHA256, URL, vertices, write_csv

SELECTORS = {
    'fore_upper': (69, 97, 167), 'fore_lower': (70, 195, 251),
    'aft_upper': (71, 268, 302), 'aft_lower': (72, 314, 336),
}


def extract(pdf, out):
    import fitz
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError('Source PDF differs from audited NASA report')
    out.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(pdf)
    page = doc[80]
    ds = page.get_drawings()
    x0 = ds[38]['items'][0][1].x
    x1 = ds[52]['items'][0][2].x
    cp0 = ds[13]['items'][0][2].y
    cp1 = ds[7]['items'][0][2].y
    assert abs(x0-117.084) < .002 and abs(x1-537.064) < .002
    assert abs(cp0-409.991) < .002 and abs(cp1-514.986) < .002
    manifest = dict(pdf_sha256=PDF_SHA256, source_url=URL, pdf_page=81,
        figure='22(b)', flap=0, reynolds=1e6, mach=.1, target_cl=1.04,
        target_lift_precision='Two decimal places in the published caption',
        mses_ncrit=9., transition='free',
        calibration=dict(x0=x0, x1=x1, cp0=cp0, cp1=cp1), selectors=SELECTORS,
        extraction_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for label, (theory, lo, hi) in SELECTORS.items():
        for kind in ('mses', 'experiment'):
            if kind == 'mses':
                points = [(theory, *p) for p in vertices(ds[theory])]
            else:
                points = []
                for i in range(lo, hi+1, 2):
                    d = ds[i]
                    rect = d['rect']
                    if len(d['items']) != 24 or abs(rect.width-rect.height) > .004:
                        raise ValueError('Expected a circular pressure marker')
                    points.append((i, .5*(rect.x0+rect.x1), .5*(rect.y0+rect.y1)))
            rows = [dict(drawing_index=i, pdf_x=px, pdf_y=py,
                         x_over_cref=(px-x0)/(x1-x0), Cp=(py-cp0)/(cp1-cp0))
                    for i, px, py in points]
            write_csv(out/(kind+'_'+label+'.csv'), rows)
    manifest['files_sha256'] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.glob('*.csv'))}
    (out/'source.json').write_text(json.dumps(manifest, indent=2)+'\n')
    page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(out/'source_page.png')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pdf', type=Path, default=DATA/'NASA-CR-20220014964.pdf')
    ap.add_argument('--out', type=Path, default=DATA/'pressure_cl1.04')
    a = ap.parse_args()
    extract(a.pdf, a.out)
