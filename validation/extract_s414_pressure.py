#!/usr/bin/env python3
"""Extract two published S414 MSES/experimental pressure comparisons.

No CFD. The 2010 approximate-envelope MSES curves remain distinct from 2014
full-e^N curves. Geometry and tap source qualification must be run first.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qualify_s414_geometry import OUT, ROOT, PINNED, sha, require, write_json, write_csv, extract_taps
from extract_s702 import vertices

CASES = [dict(name='cl028', page=66, figure='12(a)', nominal_cl=.28, experimental_cl=.284,
              experimental_alpha_deg=-1.04, snapshot_page=36, snapshot_cd=.00655,
              x0path=26, x1path=40, theory=[57,58,59,60], starts=[85,155,409,461]),
         dict(name='cl104', page=67, figure='12(b)', nominal_cl=1.04, experimental_cl=1.041,
              experimental_alpha_deg=6.09, snapshot_page=43, snapshot_cd=.01449,
              x0path=40, x1path=54, theory=[71,72,73,74], starts=[99,169,433,485])]
SURFACES = ('fore_upper','fore_lower','aft_upper','aft_lower')
COUNTS = (25,37,20,12)


def pressure_coordinates(px, py, calibration):
    return ((px-calibration['x0'])/calibration['scale'],
            (py-calibration['cp0'])/calibration['cp_scale'])


def main():
    import fitz
    pdf = OUT/'sources/s414_report_2010.pdf'
    require(sha(pdf) == PINNED[pdf.name], 'S414 report changed')
    article=ROOT/'validation/results/reference_review_20260924/sources/s414_article.pdf'
    require(sha(article)=='82b7e6292278e459dc11bfc9d9d64259c1129bc4d99e0f4f4f9870af0d91763b',
            '2014 transition-method attribution source changed')
    source_geometry = json.loads((OUT/'geometry_extraction.json').read_text())
    require(source_geometry['as_built_coordinates'] is False, 'Unexpected geometry qualification')
    doc = fitz.open(pdf)
    taps = extract_taps(doc)
    output = OUT/'pressure';output.mkdir(exist_ok=True)
    exp_rows, theory_rows, audits, tap_checks = [], [], [], []
    for case in CASES:
        ds = doc[case['page']-1].get_drawings()
        x0 = ds[case['x0path']]['rect'].x1
        scale = ds[case['x1path']]['rect'].x1-x0
        cp0 = ds[15]['rect'].y0
        cp_scale = ds[9]['rect'].y0-cp0
        require(abs(cp_scale/scale-.25) < 1e-5, 'Pressure calibration changed')
        cal = dict(x0=x0, scale=scale, cp0=cp0, cp_scale=cp_scale)
        max_tap_error = 0.
        curves = {}
        for surface, start, count, theory in zip(SURFACES, case['starts'], COUNTS, case['theory']):
            element, side = surface.split('_')
            rows = []
            table = [r['x_over_c'] for r in taps if r['element']==element and r['surface']==side]
            for j in range(count):
                i = start+2*j
                outer, dot = ds[i]['rect'], ds[i+1]['rect']
                require(len(ds[i]['items'])==24 and 4<outer.width<6 and .7<dot.width<1.1,
                        'Pressure marker selector changed')
                require(.7<dot.height<1.1 and 4<outer.height<6, 'Unexpected marker height')
                px,py = (dot.x0+dot.x1)/2, (dot.y0+dot.y1)/2
                require(outer.contains(dot), 'Central marker outside circle')
                x,cp = pressure_coordinates(px,py,cal)
                duplicate = side=='lower' and j in (0,count-1)
                rows.append(dict(case=case['name'], surface=surface, index=j, pdf_page=case['page'],
                                 drawing_index=i, central_dot_index=i+1, pdf_x=px,pdf_y=py,
                                 x_over_c=x,Cp=cp,endpoint_duplicate=int(duplicate)))
            actual = [r['x_over_c'] for r in rows if not r['endpoint_duplicate']]
            require(len(actual)==len(table), 'Table II marker count mismatch')
            error = float(np.max(abs(np.array(actual)-table)))
            # The source fails a 2e-5c exact-coordinate check. Preserve the mismatch;
            # require unambiguous one-to-one nearest-table association, not a fitted axis.
            nearest = np.argmin(abs(np.array(actual)[:,None]-np.array(table)[None,:]),axis=1)
            require(np.array_equal(nearest,np.arange(len(table))), 'Ambiguous tap association')
            for j,(measured,printed) in enumerate(zip(actual,table)):
                tap_checks.append(dict(case=case['name'],surface=surface,index=j,
                    plotted_x=measured,table_x=printed,plotted_minus_table_x=measured-printed,
                    passes_2e_minus_5_check=int(abs(measured-printed)<2e-5)))
            max_tap_error = max(error,max_tap_error)
            curves[surface] = rows
            exp_rows.extend(rows)
            raw = np.array(vertices(ds[theory]))
            for j,(px,py) in enumerate(raw):
                x,cp=pressure_coordinates(px,py,cal)
                theory_rows.append(dict(case=case['name'],surface=surface,index=j,
                    pdf_page=case['page'],drawing_index=theory,pdf_x=px,pdf_y=py,x_over_c=x,Cp=cp))
        duplicate_error=0.
        for element in ('fore','aft'):
            for j in (0,-1):
                duplicate_error=max(duplicate_error,abs(curves[element+'_upper'][j]['Cp']-curves[element+'_lower'][j]['Cp']))
        # Lower curves repeat upper endpoint taps, with small source plotting offsets.
        # Keep both values and flag the repeats rather than forcing equality.
        for element in ('fore','aft'):
            for j in (0,-1):
                require(abs(curves[element+'_upper'][j]['x_over_c']-curves[element+'_lower'][j]['x_over_c'])<1e-4,
                        'Repeated endpoint identity is ambiguous')
        # Verify the experimental snapshot association against the caption values.
        text=doc[case['snapshot_page']-1].get_text()
        require(f"{case['experimental_cl']:.3f}" in text and f"{case['experimental_alpha_deg']:.2f}" in text.replace('−','-')
                and f"{case['snapshot_cd']:.5f}" in text,'Experimental snapshot caption differs')
        audits.append(dict(case, calibration=cal,maximum_tap_x_error=max_tap_error,
                          duplicate_endpoint_max_cp_error=duplicate_error,
                          exact_table_coordinate_check_pass=max_tap_error<2e-5,
                          ordered_nearest_table_association_pass=True,
                          matching='Caption nominal CL is rounded; experimental Figure 7 CL retained separately'))
    write_csv(output/'experiment.csv',exp_rows)
    write_csv(output/'mses_2010_envelope.csv',theory_rows)
    write_csv(output/'tap_coordinate_checks.csv',tap_checks)
    plot(exp_rows,theory_rows)
    receipt=dict(source_url='https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf',
        product='RDECOM TR 10-D-112 (2010), Figure 12(a,b)',reynolds=1e6,mach=.1,
        reference_chord_m=.4572,transition='free',mses_ncrit=9,
        mses_transition='Approximate envelope e^N; identified by Coder et al. 2014, section III.A',
        mses_executable_run=False, measured_transition_endpoints_available=False,
        pressure_uncertainty='No case-specific Cp uncertainty budget recovered',
        actual_geometry='Modified tested shape; original solver geometry file remains unavailable',
        status='Traceable plotted curves with unambiguous tap association; exact Table II coordinate check fails, offsets retained',
        cases=audits,experiment_markers=len(exp_rows),repeated_endpoints=sum(r['endpoint_duplicate'] for r in exp_rows),
        independent_tap_locations=len(taps),theory_vertices=len(theory_rows),
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in (pdf,Path(__file__),
            ROOT/'validation/qualify_s414_geometry.py',ROOT/'validation/extract_s702.py',OUT/'geometry_extraction.json',
            ROOT/'validation/results/reference_review_20260924/sources/s414_article.pdf')},
        output_sha256={str(p.relative_to(output)):sha(p) for p in sorted(output.iterdir()) if p.suffix in ('.csv','.png','.svg','.pdf')})
    write_json(output/'source.json',receipt)
    print(json.dumps({k:receipt[k] for k in ('experiment_markers','repeated_endpoints','independent_tap_locations','theory_vertices')},indent=2))


def plot(exp,theory):
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for row,case in enumerate(CASES):
        for col,element in enumerate(('fore','aft')):
            ax=axes[row,col]
            for side,color in [('upper','#1b6485'),('lower','#b66b24')]:
                surface=element+'_'+side
                e=[r for r in exp if r['case']==case['name'] and r['surface']==surface and not r['endpoint_duplicate']]
                t=[r for r in theory if r['case']==case['name'] and r['surface']==surface]
                ax.plot([r['x_over_c'] for r in t],[r['Cp'] for r in t],color=color,label=f'MSES {side}')
                ax.plot([r['x_over_c'] for r in e],[r['Cp'] for r in e],'o',mfc='white',ms=4,color=color,label=f'Experiment {side}')
            ax.invert_yaxis();ax.grid(alpha=.2)
            ax.set(xlabel='x / complete chord',ylabel='$C_p$',title=f"{element.capitalize()} · nominal CL {case['nominal_cl']:.2f} (measured {case['experimental_cl']:.3f})")
            if element=='aft':ax.set_xlim(.665,1.005)
            ax.legend(fontsize=8)
    fig.suptitle('S414 published pressure comparison · modified tested shape\nRe = 1,000,000 · M = 0.10 · 2010 MSES envelope method, Ncrit 9 · no MultiFoil run')
    for ext in ('png','svg','pdf'):fig.savefig(OUT/'pressure'/f'published_pressure.{ext}',dpi=170)
    plt.close(fig)


if __name__=='__main__':main()
