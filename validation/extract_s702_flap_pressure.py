#!/usr/bin/env python3
"""Qualify four S702 pressure snapshots without interpolating pressure curves.

Re-extraction uses PyMuPDF 1.26.7; the numerical runner uses only saved CSV/JSON.
Selectors refer to the SHA-pinned NASA PDF, with one-based PDF page numbers.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from extract_s702 import DATA, PDF_SHA256, URL, write_csv

ROOT = DATA.parents[2]
OUT = ROOT/'validation/results/s702_pressure_qualification_20261002'
CASES = [
    dict(name='d5_a2p03', flap=5, alpha=2.03, page=43, figure='7(b)', marker='diamond',
         x1_path=45, starts=[842, 940, 1014, 1062], counts=[36, 29, 19, 12]),
    dict(name='d5_a6p11', flap=5, alpha=6.11, page=44, figure='7(c)', marker='diamond',
         x1_path=60, starts=[870, 969, 1045, 1095], counts=[36, 29, 19, 12]),
    dict(name='d10_am0p01', flap=10, alpha=-.01, page=48, figure='8(b)', marker='triangle',
         x1_path=45, starts=[584, 679, 755, 804], counts=[36, 29, 20, 12]),
    dict(name='d10_a4p07', flap=10, alpha=4.07, page=49, figure='8(c)', marker='triangle',
         x1_path=52, starts=[603, 701, 776, 826], counts=[36, 29, 20, 12]),
]
SURFACES = ('fore_upper', 'fore_lower', 'aft_upper', 'aft_lower')
# Independent Table II(a), not derived from pressure-figure coordinates.
FORE_UPPER = np.array([.00026,.00198,.00682,.01462,.02564,.03961,.05661,.07680,
    .09926,.12465,.15099,.18048,.21118,.24364,.27736,.31247,.34791,.38394,
    .41894,.45569,.49061,.52582,.55998,.59183,.62247,.65136,.67858,.70403,
    .72663,.74583,.76395,.77766,.78828,.79672,.80154,.80270])
AFT_UPPER_PREHINGE = np.array([.69941,.70194,.70930,.72177,.73844,.75893,
    .78231,.80843,.83542,.84940,.86303,.87212,.88138,.88928])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def match_lift(rows, alpha):
    """Associate with a measured polar marker, within declared 0.02-degree limit."""
    row = min(rows, key=lambda r: abs(float(r['alpha_deg'])-alpha))
    delta = float(row['alpha_deg'])-alpha
    if abs(delta) > .02:
        raise ValueError('No experimental lift marker within 0.02 degrees')
    return dict(target_cl=float(row['CL']), polar_alpha_deg=float(row['alpha_deg']),
                alpha_difference_deg=delta, polar_drawing_index=int(row['drawing_index']))


def extract(out):
    import fitz
    pdf = DATA/'NASA-CR-20220014964.pdf'
    if sha(pdf) != PDF_SHA256:
        raise ValueError('Source PDF differs from audited NASA report')
    out.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(pdf)
    archive=ROOT/'handover/archive_manifest.json'
    if sha(archive)!='9102d74534effcc51c30872e0e1c5a2be1e5bc6e0d7fc5cd6d6c97c715fb4fc8':
        raise ValueError('Original data receipt changed')
    original={r['path']:r['sha256'] for r in json.loads(archive.read_text())['files']}
    source = dict(source_url=URL, pdf_sha256=PDF_SHA256, product='NASA/CR-20220014964',
        access_mode='public anonymous; existing local source snapshot', download_status='verified',
        reynolds=1e6, mach=.1, transition='free', reference_chord_m=.45748,
        extraction_runtime=dict(pymupdf=fitz.VersionBind), cases=[],
        extraction_method='Vector marker central dots, not triangle bounding-box centers',
        lift_matching='Nearest measured lift marker within 0.02 deg of rounded pressure caption; no Cp interpolation',
        pressure_uncertainty='No case-specific experimental Cp uncertainty budget recovered',
        transition_observation='Qualitative pressure features and author discussion; no exact per-snapshot transition coordinate',
        mses_pressure_available_for_these_cases=False, files_sha256={},
        input_sha256={str(pdf.relative_to(ROOT)):sha(pdf),str(archive.relative_to(ROOT)):sha(archive),
                      str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)})
    features=[]
    for case in CASES:
        page = doc[case['page']-1]; ds = page.get_drawings()
        x0 = ds[0]['items'][0][1].x
        x1 = ds[case['x1_path']]['items'][0][2].x
        cp0 = ds[13]['items'][0][2].y; cp1 = ds[7]['items'][0][2].y
        scale = x1-x0
        if abs((cp1-cp0)/scale-.25) > 1e-5:
            raise ValueError('Pressure and chord axis scales inconsistent')
        folder = out/case['name']; folder.mkdir(exist_ok=True)
        curves = {}
        for label, start, count in zip(SURFACES, case['starts'], case['counts']):
            rows = []
            for j in range(count):
                i = start+2*j; outer=ds[i]['rect']; dot=ds[i+1]['rect']
                if not (.7 < dot.width < 1.1 and .7 < dot.height < 1.1 and
                        4 < outer.width < 7.3 and 4 < outer.height < 7.3):
                    raise ValueError(f'Unexpected marker structure: {case["page"]}:{i}')
                px = (dot.x0+dot.x1)/2; py=(dot.y0+dot.y1)/2
                if not (outer.x0 < px < outer.x1 and outer.y0 < py < outer.y1):
                    raise ValueError('Central dot outside marker')
                # Lower polylines repeat the upper leading/trailing endpoints.
                # Preserve those plotted values but exclude them from lower-tap scores.
                duplicate = label.endswith('lower') and j in (0, count-1)
                rows.append(dict(drawing_index=i, central_dot_index=i+1, pdf_x=px, pdf_y=py,
                    x_over_cref=(px-x0)/scale, Cp=(py-cp0)/(cp1-cp0),
                    endpoint_duplicate=int(duplicate)))
            curves[label] = rows
            write_csv(folder/(label+'.csv'), rows)
        upper_x = np.array([r['x_over_cref'] for r in curves['fore_upper']])
        fore_error = float(np.max(abs(upper_x-FORE_UPPER)))
        aft_x = np.array([r['x_over_cref'] for r in curves['aft_upper'] if r['x_over_cref']<.89])
        aft_error = float(abs(aft_x[:,None]-AFT_UPPER_PREHINGE[None,:]).min(axis=1).max())
        if fore_error > 8e-5 or aft_error > 8e-5:
            raise ValueError('Pressure taps disagree with independent Table II')
        endpoint_error = 0.
        for element in ('fore','aft'):
            for j in (0,-1):
                a,b=curves[element+'_upper'][j],curves[element+'_lower'][j]
                endpoint_error=max(endpoint_error,abs(a['Cp']-b['Cp']))
        if endpoint_error > .001:
            raise ValueError('Expected duplicated pressure endpoints differ')
        polar = DATA/f'd+{case["flap"]}'/'experiment_lift.csv'
        if sha(polar)!=original[str(polar.relative_to(ROOT))]:
            raise ValueError('Original experimental lift source changed')
        source['input_sha256'][str(polar.relative_to(ROOT))] = sha(polar)
        matched = match_lift(read_csv(polar),case['alpha'])
        source['cases'].append(dict(case, **matched,
            calibration=dict(x0=x0,x1=x1,cp0=cp0,cp1=cp1),
            audit=dict(fore_tap_max_error=fore_error,aft_prehinge_tap_max_error=aft_error,
                       duplicated_endpoint_max_cp_difference=endpoint_error),
            caption=page.get_text().strip()))
        values={}
        for kind in ('experiment','mses'):
            polar=DATA/f'd+{case["flap"]}'/(kind+'_drag.csv')
            name=str(polar.relative_to(ROOT))
            if sha(polar)!=original[name]:raise ValueError('Original drag context source changed')
            source['input_sha256'][name]=sha(polar)
            drag=sorted(read_csv(polar),key=lambda r:float(r['CL']))
            cl=np.array([float(r['CL']) for r in drag]);cd=np.array([float(r['CD']) for r in drag])
            if not cl[0]<=matched['target_cl']<=cl[-1]:raise ValueError('Drag context outside support')
            values[kind+'_CD']=float(np.interp(matched['target_cl'],cl,cd))
        up=curves['aft_upper']
        left=min(up,key=lambda r:abs(r['x_over_cref']-.87212))
        right=min(up,key=lambda r:abs(r['x_over_cref']-.88928))
        features.append(dict(case=case['name'],flap=case['flap'],alpha_deg=case['alpha'],
            target_cl=matched['target_cl'],prehinge_left_x=left['x_over_cref'],
            prehinge_right_x=right['x_over_cref'],prehinge_delta_Cp=right['Cp']-left['Cp'],
            **values,mses_minus_experiment_drag_counts=(values['mses_CD']-values['experiment_CD'])*1e4))
        page.get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(folder/'source_page.png')
    write_csv(out/'snapshot_features_and_drag_context.csv',features)
    plot(out,source)
    for path in sorted(out.rglob('*')):
        if path.is_file() and path.suffix in ('.csv','.png','.pdf','.svg'):
            source['files_sha256'][str(path.relative_to(out))]=sha(path)
    (out/'source.json').write_text(json.dumps(source,indent=2)+'\n')
    print(json.dumps(dict(cases=len(CASES),markers=sum(sum(c['counts']) for c in CASES),
                         matched_lifts=[c['target_cl'] for c in source['cases']]),indent=2))


def plot(out, source):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for element in ('aft','fore'):
        fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        for ax,case in zip(axes.flat,source['cases']):
            for side,color in (('upper','#1f77b4'),('lower','#c66b19')):
                rows=read_csv(out/case['name']/(element+'_'+side+'.csv'))
                rows=[r for r in rows if r['endpoint_duplicate']=='0']
                ax.plot([float(r['x_over_cref']) for r in rows],[float(r['Cp']) for r in rows],
                    'o-',color=color,ms=4,lw=.9,label=side+' pressure taps')
            if element=='aft':
                ax.axvline(.9001,color='gray',ls=':',lw=1,label='hinge x/c')
            ax.invert_yaxis();ax.grid(alpha=.2);ax.set_xlabel('x / reference chord');ax.set_ylabel('$C_p$')
            ax.set_title(f"Flap {case['flap']}° · α {case['alpha']:.2f}° · measured CL ≈ {case['target_cl']:.3f}\nFigure {case['figure']}",fontsize=11)
            ax.legend(fontsize=8)
        fig.suptitle(f'S702 {element} element — qualified experimental snapshots\nRe = 1,000,000 · M = 0.10 · free transition',fontsize=14)
        for ext in ('png','pdf','svg'):fig.savefig(out/f'{element}_pressure_snapshots.{ext}',dpi=170)
        plt.close(fig)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=OUT)
    extract(ap.parse_args().out)
