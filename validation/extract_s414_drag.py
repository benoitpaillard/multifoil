#!/usr/bin/env python3
"""Extract S414 Figure 13(c), retaining measured lift/drag plotting offsets."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qualify_s414_geometry import ROOT, OUT as GEOMETRY, PINNED, sha, require, write_json, write_csv
from extract_s702 import vertices

OUT=ROOT/'validation/results/s414_drag_reference_20261002'
TARGETS=[dict(case='cl028',CL=.284,alpha_deg=-1.04,caption_CD=.00655,nominal_mses_pressure_CL=.28),
         dict(case='cl104',CL=1.041,alpha_deg=6.09,caption_CD=.01449,nominal_mses_pressure_CL=1.04)]


def interpolate(rows,target,x='CL',y='CD'):
    a=np.array([[r[x],r[y]] for r in rows],float)
    if not np.isfinite(a).all() or np.any(np.diff(a[:,0])<=0):
        raise ValueError('Require finite, strictly increasing independent coordinates')
    if not a[0,0]<=target<=a[-1,0]:raise ValueError('Reference extrapolation is forbidden')
    return float(np.interp(target,a[:,0],a[:,1]))


def ascending_branch(rows):
    """Keep the declared increasing-alpha branch up to first maximum CL."""
    i=int(np.argmax([r['CL'] for r in rows]))
    branch=rows[:i+1]
    if len(branch)<2 or any(b['CL']<=a['CL'] or b['alpha_deg']<=a['alpha_deg'] for a,b in zip(branch,branch[1:])):
        raise ValueError('Ambiguous published pre-stall lift branch')
    return branch


def extract():
    import fitz
    pdf=GEOMETRY/'sources/s414_report_2010.pdf'
    require(sha(pdf)==PINNED[pdf.name],'Primary report changed')
    receipt=json.loads((GEOMETRY/'verification.json').read_text())
    for name in ('pressure/source.json','pressure/experiment.csv','pressure/mses_2010_envelope.csv'):
        require(sha(GEOMETRY/name)==receipt['output_sha256'][name], 'Prior pressure qualification changed')
    doc=fitz.open(pdf);page=doc[70];ds=page.get_drawings()
    cl0=ds[49]['rect'].x0;cl_scale=(cl0-ds[37]['rect'].x0)/2
    cd0=ds[0]['rect'].y0;cd_scale=(cd0-ds[74]['rect'].y0)/.02
    alpha0=ds[87]['rect'].y0;alpha_scale=cl_scale/20
    require(abs(cd_scale/cl_scale-100)<.002,'Axis scale ratio changed')
    cal=dict(cl_zero_pdf_x=cl0,points_per_CL=cl_scale,cd_zero_pdf_y=cd0,points_per_CD=cd_scale,
             alpha_zero_pdf_y=alpha0,points_per_degree=alpha_scale)
    curves={}
    for kind,theory,limits in [('drag',232,(259,302)),('lift',233,(319,384))]:
        for method in ('experiment','mses'):
            raw=[]
            if method=='mses':raw=[(theory,None,*p) for p in vertices(ds[theory])]
            else:
                for i in range(*limits,2):
                    dot=ds[i+1]['rect'];circle=ds[i]['rect']
                    require(len(ds[i]['items'])==24 and .6<dot.width<.7 and .6<dot.height<.7
                            and 3.4<circle.width<3.6 and circle.contains(dot),'Marker selector changed')
                    raw.append((i,i+1,(dot.x0+dot.x1)/2,(dot.y0+dot.y1)/2))
            key=method+'_'+kind;rows=[]
            for j,(i,dot,px,py) in enumerate(raw):
                row=dict(index=j,pdf_page=71,drawing_index=i,central_dot_index=dot,pdf_x=px,pdf_y=py,CL=(cl0-px)/cl_scale)
                row['CD' if kind=='drag' else 'alpha_deg']=(cd0-py)/cd_scale if kind=='drag' else (alpha0-py)/alpha_scale
                rows.append(row)
            curves[key]=rows
    OUT.mkdir(exist_ok=True)
    for key,rows in curves.items():write_csv(OUT/(key+'.csv'),rows)
    # The original plots have a small constant lift-axis offset; never align them by fitting.
    pairs=[]
    for drag,lift in zip(curves['experiment_drag'],curves['experiment_lift']):
        delta=drag['CL']-lift['CL']
        require(abs(delta)<1e-4,'Drag/lift markers no longer associate by acquisition order')
        pairs.append(dict(index=drag['index'],drag_CL=drag['CL'],lift_CL=lift['CL'],CL_difference=delta,alpha_deg=lift['alpha_deg']))
    refs=[]
    for t in TARGETS:
        k=min(range(len(curves['experiment_drag'])),key=lambda i:abs(curves['experiment_lift'][i]['alpha_deg']-t['alpha_deg']))
        drag,lift=curves['experiment_drag'][k],curves['experiment_lift'][k]
        require(abs(lift['alpha_deg']-t['alpha_deg'])<.01 and abs(drag['CL']-t['CL'])<.0006
                and abs(drag['CD']-t['caption_CD'])<6e-6,'Caption association failed')
        exp=interpolate(curves['experiment_drag'],t['CL']);mses=interpolate(curves['mses_drag'],t['CL'])
        refs.append(dict(t,experiment_CD_interpolated=exp,mses_CD_interpolated=mses,
                         mses_minus_experiment_counts=(mses-exp)*1e4,
                         experimental_marker_index=k,experimental_marker_CL=drag['CL'],experimental_marker_CD=drag['CD'],
                         experimental_marker_alpha=lift['alpha_deg'],
                         mses_seed_alpha=interpolate(ascending_branch(curves['mses_lift']),t['CL'],y='alpha_deg')))
    write_csv(OUT/'target_references.csv',refs);write_csv(OUT/'lift_drag_marker_associations.csv',pairs)
    page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(OUT/'source_page_71.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for method,color in [('experiment','#222222'),('mses','#327da8')]:
        rows=curves[method+'_drag'];axes[0].plot([r['CD']*1e4 for r in rows],[r['CL'] for r in rows],
            'o-' if method=='experiment' else '-',color=color,ms=4,label='Experiment' if method=='experiment' else 'Published MSES N9 (2010 envelope)')
    for row in refs:axes[0].axhline(row['CL'],color='gray',ls=':',lw=.8)
    axes[0].set(xlabel='CD (drag counts)',ylabel='CL',title='Re 1,000,000 · M 0.10 · free transition');axes[0].legend(fontsize=8)
    bars=axes[1].bar([f"CL {r['CL']:.3f}" for r in refs],[r['mses_minus_experiment_counts'] for r in refs],color=['#bb7733','#327da8'])
    axes[1].bar_label(bars,fmt='%+.2f',padding=4);axes[1].axhline(0,color='black',lw=.7)
    axes[1].set(ylabel='MSES minus experiment (drag counts)',title='Same-lift interpolation at pressure targets',ylim=(-17,8))
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('S414 matching drag reference · Figure 13(c)\nPublished curves; no new MSES run or MultiFoil result')
    for ext in ('png','svg','pdf'):fig.savefig(OUT/f'drag_reference.{ext}',dpi=170)
    plt.close(fig)
    inputs=[pdf,GEOMETRY/'verification.json',GEOMETRY/'pressure/source.json',GEOMETRY/'pressure/experiment.csv',
            GEOMETRY/'pressure/mses_2010_envelope.csv',Path(__file__),ROOT/'validation/qualify_s414_geometry.py',ROOT/'validation/extract_s702.py']
    source=dict(product='RDECOM TR 10-D-112 (2010), Figure 13(c), PDF p.71 / printed p.69',
        url='https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf',access_mode='Previously acquired public PDF, hash verified',
        reynolds=1e6,mach=.1,ncrit_mses=9,mses_transition='Approximate envelope e^N, not 2014 full e^N',reference_chord_m=.4572,
        calibration=cal,experiment_drag_markers=len(curves['experiment_drag']),experiment_lift_markers=len(curves['experiment_lift']),
        mses_drag_vertices=len(curves['mses_drag']),targets=refs,
        maximum_drag_lift_CL_offset=max(abs(p['CL_difference']) for p in pairs),
        qualification='Pressure target CL uses rounded experimental captions. Interpolation is piecewise linear without extrapolation. CL plotting offsets and direct caption CD are retained separately. Figure 13 drag coverage ends before experimental stall.',
        extraction_precision='PDF coordinate precision is not physical measurement uncertainty; no exact Cp/bubble-location uncertainty recovered',
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        files_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.suffix in ('.csv','.png','.svg','.pdf')},pymupdf=fitz.VersionBind)
    write_json(OUT/'source.json',source);print(json.dumps(refs,indent=2))


if __name__=='__main__':extract()
