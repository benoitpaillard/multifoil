#!/usr/bin/env python3
"""Reproduce COBEM2005 Fig.3 MSES drag samples from archived PDF vectors.

Default checks the reference CSV. --pdf SOURCE additionally verifies the
archived vectors against the original, hash-checked PDF (requires PyMuPDF).
--write updates the CSV. Coordinates are graphical readings, not native MSES data.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

DIRECTORY = Path(__file__).resolve().parent/'data/naca0012_mses_drag'


def pdf_segments(path, cfg):
    import fitz
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != cfg['source_pdf_sha256']:
        raise ValueError('Original PDF hash mismatch')
    with fitz.open(path) as document:
        page = document[cfg['pdf_page_index']]
        segments = []
        for drawing in page.get_drawings():
            color = drawing['color']
            # Red MSES line within the data region; excludes legend and Fig.5.
            if not color or not (color[0] > .8 and color[1] < .01 and color[2] < .01):
                continue
            if not fitz.Rect(*cfg['curve_selection_rect_pt']).contains(drawing['rect']):
                continue
            for item in drawing['items']:
                if item[0] == 'l':
                    segments.append([list(item[1]), list(item[2])])
    return segments


def extract(directory=DIRECTORY, pdf=None):
    cfg=json.loads((directory/'digitization.json').read_text())
    vectors=directory/'source_vectors.json'
    if hashlib.sha256(vectors.read_bytes()).hexdigest() != cfg['vectors_sha256']:
        raise ValueError('Archived PDF vector hash mismatch')
    segments=json.loads(vectors.read_text())['mses_segments_pt']
    if pdf is not None and pdf_segments(pdf,cfg) != segments:
        raise ValueError('Archived vectors differ from original PDF')
    points=[segments[0][0]]
    for start,end in segments:
        if start != points[-1] or end[0] <= start[0]:
            raise ValueError('MSES curve must be continuous and increasing in incidence')
        points.append(end)
    (x0,a0),(x1,a1)=cfg['alpha_calibration_pt']
    (y0,c0),(y1,c1)=cfg['cd_calibration_pt']
    rows=[]
    for x,y in points:
        alpha=round(a0+(x-x0)*(a1-a0)/(x1-x0),6)
        cd=round(c0+(y-y0)*(c1-c0)/(y1-y0),6)
        allowance=cfg['cd_reading_allowance']
        rows.append(dict(alpha_deg=alpha,CD=cd,CD_low=round(cd-allowance,6),
            CD_high=round(cd+allowance,6),
            alpha_reading_allowance_deg=cfg['alpha_reading_allowance_deg'],
            reference_kind='published_solver',reference_solver='MSES'))
    if [r['alpha_deg'] for r in rows] != cfg['sampling_alpha_deg']:
        raise ValueError('Published point locations differ from expected incidences')
    return rows


def csv_text(rows):
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return stream.getvalue()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf',type=Path)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    rows=extract(pdf=args.pdf);target=DIRECTORY/'polar.csv'
    if args.write:
        target.write_text(csv_text(rows))
    elif target.read_text() != csv_text(rows):
        raise SystemExit('Distributed drag reference differs from extraction')
    print(f'{len(rows)} published MSES CD samples '+('written' if args.write else 'reproduced exactly'))


if __name__=='__main__': main()
