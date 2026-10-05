#!/usr/bin/env python3
"""Qualify S414 trip comparisons and free-transition boundary pressures. No CFD."""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qualify_s414_geometry import ROOT,OUT as GEOMETRY,PINNED,sha,require,write_json,write_csv,extract_taps
from extract_s702 import vertices
from extract_s414_drag import interpolate,ascending_branch

OUT=ROOT/'validation/results/s414_trip_qualification_20261003'
PDF=GEOMETRY/'sources/s414_report_2010.pdf'
STATES=('free','fore_fixed','both_fixed')
LABELS={'free':'Free transition','fore_fixed':'Fore trips; aft free','both_fixed':'Fore + aft trips'}
COLORS={'free':'#236890','fore_fixed':'#c16a2b','both_fixed':'#6b5191'}
CONDITIONS=[dict(reynolds=700000,mach=.07,tables=[88,89,90],figures=[70,74,78],
                table_counts=[31,25,25],marker_counts=[22,15,15]),
            dict(reynolds=1000000,mach=.10,tables=[91,92,93],figures=[71,75,79],
                table_counts=[33,29,26],marker_counts=[22,18,16])]
PRESSURES=[dict(page=33,figure='7(a)',CL=-.072,alpha=-4.09,CD=.01069,starts=[53,123,387,439]),
           dict(page=34,figure='7(b)',CL=.056,alpha=-3.07,CD=.00618,starts=[53,123,381,433]),
           dict(page=35,figure='7(c)',CL=.175,alpha=-2.06,CD=.00669,starts=[53,122,374,426]),
           dict(page=36,figure='7(d)',CL=.284,alpha=-1.04,CD=.00655,starts=[53,123,377,429])]
SURFACES=('fore_upper','fore_lower','aft_upper','aft_lower')
TARGETS=(0.,.284,.5,1.041)


def table_tokens(text):
    rows=[]
    for line in text.splitlines():
        tokens=line.replace('−','-').split()
        if len(tokens)==4 and all(re.fullmatch(r'-?(?:\d+(?:\.\d+)?|\.\d+)|—',x) for x in tokens):
            rows.append(tokens)
    return rows


def parse_table(text,condition,state,page):
    expected='transition free' if state=='free' else 'transition fixed on fore'+(' and aft elements' if state=='both_fixed' else ' element')
    require(expected in text,'Table condition caption mismatch')
    require(f"{condition['reynolds']/1e6:.2f}" in text and f"{condition['mach']:.2f}" in text,'Re/M caption mismatch')
    rows=[]
    for i,tokens in enumerate(table_tokens(text)):
        values=[None if t=='—' else float(t) for t in tokens]
        require(all(values[k] is not None for k in (0,1,3)),'Missing non-drag table value')
        rows.append(dict(reynolds=condition['reynolds'],mach=condition['mach'],configuration=state,
            pdf_page=page,printed_page=page-2,index=i,raw_tokens=' '.join(tokens),
            alpha_deg=values[0],CL=values[1],CD=values[2],CM=values[3],drag_available=values[2] is not None))
    branch=ascending_branch(rows)
    for i,r in enumerate(rows):r['prestall_branch']=i<len(branch)
    return rows


def extract_polar(page,condition,state,pdf_page,expected_markers):
    ds=page.get_drawings();cl0=ds[49]['rect'].x0;cl_scale=(cl0-ds[37]['rect'].x0)/2
    cd0=ds[0]['rect'].y0;cd_scale=(cd0-ds[74]['rect'].y0)/.02
    require(abs(cd_scale/cl_scale-100)<.002,'Polar calibration ratio changed')
    require(ds[232]['type']=='s' and ds[232]['dashes']=='[] 0','MSES curve style changed')
    base=dict(reynolds=condition['reynolds'],mach=condition['mach'],configuration=state,pdf_page=pdf_page)
    theory=[dict(base,index=i,drawing_index=232,pdf_x=x,pdf_y=y,CL=(cl0-x)/cl_scale,CD=(cd0-y)/cd_scale)
            for i,(x,y) in enumerate(vertices(ds[232]))]
    require(all(b['CL']>a['CL'] for a,b in zip(theory,theory[1:])),'MSES drag curve folds')
    markers=[];i=259
    while len(ds[i]['items'])==24 and 3.4<ds[i]['rect'].width<3.6:
        outer,dot=ds[i]['rect'],ds[i+1]['rect']
        require(.6<dot.width<.7 and .6<dot.height<.7 and outer.contains(dot),'Drag marker selector changed')
        x,y=(dot.x0+dot.x1)/2,(dot.y0+dot.y1)/2
        markers.append(dict(base,index=len(markers),drawing_index=i,central_dot_index=i+1,
            pdf_x=x,pdf_y=y,CL=(cl0-x)/cl_scale,CD=(cd0-y)/cd_scale));i+=2
    require(len(markers)==expected_markers,'Polar marker coverage changed')
    return theory,markers,dict(base,cl_zero_pdf_x=cl0,points_per_CL=cl_scale,
        cd_zero_pdf_y=cd0,points_per_CD=cd_scale,theory_path=232,theory_vertices=len(theory),experimental_markers=len(markers))


def delta_at(curves,reynolds,target):
    row=dict(reynolds=reynolds,mach=.1 if reynolds==1000000 else .07,CL=target)
    for method in ('experiment','mses'):
        for state in STATES:
            row[method+'_'+state+'_CD']=interpolate(curves[(reynolds,method,state)],target)
        for change,a,b in [('fore_minus_free','fore_fixed','free'),('both_minus_fore','both_fixed','fore_fixed')]:
            row[method+'_'+change+'_counts']=1e4*(row[method+'_'+a+'_CD']-row[method+'_'+b+'_CD'])
    for change in ('fore_minus_free','both_minus_fore'):
        row['mses_minus_experiment_'+change+'_counts']=row['mses_'+change+'_counts']-row['experiment_'+change+'_counts']
    for state in STATES:
        row['mses_minus_experiment_'+state+'_counts']=1e4*(row['mses_'+state+'_CD']-row['experiment_'+state+'_CD'])
        row['experiment_'+state+'_alpha_deg']=interpolate(curves[(reynolds,'experiment',state)],target,y='alpha_deg')
        row['experiment_'+state+'_CM']=interpolate(curves[(reynolds,'experiment',state)],target,y='CM')
    return row


def extract_pressure(doc):
    taps=extract_taps(doc);rows=[];checks=[];recovery=[]
    for case in PRESSURES:
        page=doc[case['page']-1];text=page.get_text().replace('−','-')
        require(f"{case['CL']:.3f}" in text and f"{case['alpha']:.2f}" in text and f"{case['CD']:.5f}" in text,'Pressure caption changed')
        ds=page.get_drawings();x0=ds[33]['rect'].x0;scale=ds[47]['rect'].x1-x0
        cp0=ds[15]['rect'].y0;cp_scale=ds[9]['rect'].y0-cp0
        require(abs(cp_scale/scale-.25)<1e-5,'Pressure scale changed')
        for surface,start,count in zip(SURFACES,case['starts'],(25,37,20,12)):
            element,side=surface.split('_');table=[t['x_over_c'] for t in taps if t['element']==element and t['surface']==side]
            selected=[];k=0
            for j in range(count):
                i=start+2*j;outer,dot=ds[i]['rect'],ds[i+1]['rect']
                require(len(ds[i]['items'])==24 and 4<outer.width<6 and 4<outer.height<6 and
                        .7<dot.width<1.1 and .7<dot.height<1.1 and outer.contains(dot),'Pressure marker selector changed')
                x,y=(dot.x0+dot.x1)/2,(dot.y0+dot.y1)/2;duplicate=side=='lower' and j in (0,count-1)
                tx=None if duplicate else table[k]
                if not duplicate:k+=1
                r=dict(figure=case['figure'],pdf_page=case['page'],CL=case['CL'],alpha_deg=case['alpha'],
                    surface=surface,index=j,drawing_index=i,central_dot_index=i+1,pdf_x=x,pdf_y=y,
                    plotted_x=(x-x0)/scale,table_x=tx,Cp=(y-cp0)/cp_scale,endpoint_duplicate=duplicate)
                rows.append(r)
                if not duplicate:selected.append(r)
            actual=np.array([r['plotted_x'] for r in selected]);table=np.array(table)
            require(k==len(table) and np.array_equal(np.argmin(abs(actual[:,None]-table),axis=1),np.arange(k)),
                    'Ambiguous pressure tap association')
            checks.append(dict(figure=case['figure'],surface=surface,count=k,
                maximum_plotted_table_x_difference=float(abs(actual-table).max()),
                exact_2e_minus_5_coordinate_check=bool(np.all(abs(actual-table)<2e-5))))
            if surface=='aft_upper':
                for convention in ('table_x','plotted_x'):
                    rr=[r for r in selected if .8<=r[convention]<=.98]
                    for left,right in zip(rr,rr[1:]):
                        dc=right['Cp']-left['Cp'];dx=right[convention]-left[convention]
                        recovery.append(dict(figure=case['figure'],CL=case['CL'],coordinate=convention,
                            left_x=left[convention],right_x=right[convention],delta_Cp=dc,slope=dc/dx))
    old=list(csv.DictReader((GEOMETRY/'pressure/experiment.csv').open()))
    duplicate_errors=[]
    for row in [r for r in rows if r['figure']=='7(d)']:
        ref=next(r for r in old if r['case']=='cl028' and r['surface']==row['surface'] and int(r['index'])==row['index'])
        duplicate_errors.append(dict(surface=row['surface'],index=row['index'],
            delta_Cp=row['Cp']-float(ref['Cp']),delta_plotted_x=row['plotted_x']-float(ref['x_over_c'])))
    # Two separately plotted copies of one measurement: retain drawing offsets.
    for r in duplicate_errors:
        r['strict_copy_check_pass']=abs(r['delta_Cp'])<3e-5 and abs(r['delta_plotted_x'])<1e-5
    # The strict copy check fails: two plotted copies have small coordinate/Cp offsets.
    # Preserve failures, never align by fitting. This separate guard only rejects
    # a grossly different curve; it does not replace the failed precision check.
    require(max(abs(r['delta_Cp']) for r in duplicate_errors)<1e-3 and
            max(abs(r['delta_plotted_x']) for r in duplicate_errors)<5e-4,'Pressure copies fail gross identity screen')
    return rows,checks,recovery,duplicate_errors


def save(fig,name):
    for ext in ('png','svg','pdf'):fig.savefig(OUT/(name+'.'+ext),dpi=170)
    plt.close(fig)


def plots(curves,continuous,pressure):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,condition in enumerate(CONDITIONS):
        re=condition['reynolds'];ax=axes[0,col]
        for state in STATES:
            for method in ('experiment','mses'):
                rows=curves[(re,method,state)];rows=[r for r in rows if -.10<=r['CL']<=1.35]
                ax.plot([r['CL'] for r in rows],[r['CD']*1e4 for r in rows],
                    'o-' if method=='experiment' else '--',ms=3,color=COLORS[state],
                    label=LABELS[state]+(' · experiment' if method=='experiment' else ' · MSES'))
        ax.set(xlim=(-.05,1.2),xlabel='CL',ylabel='CD (counts)',title=f'Re {re/1e6:g} million · M {condition["mach"]:.2f}')
        if col==0:ax.legend(fontsize=7)
        ax.grid(alpha=.2);ax=axes[1,col];rows=[r for r in continuous if r['reynolds']==re]
        for method,style in [('experiment','-'),('mses','--')]:
            for change,color,label in [('fore_minus_free',COLORS['fore_fixed'],'Fore trips − free'),
                                       ('both_minus_fore',COLORS['both_fixed'],'Add aft trips')]:
                ax.plot([r['CL'] for r in rows],[r[method+'_'+change+'_counts'] for r in rows],style,
                    color=color,label=label+' · '+('experiment' if method=='experiment' else 'MSES'))
        ax.axhline(0,color='.3',lw=.7);ax.set(xlim=(-.05,1.2),xlabel='Matched CL',ylabel='Change in CD (counts)');ax.grid(alpha=.2)
        if col==0:ax.legend(fontsize=8)
    fig.suptitle('S414 trip comparison · experiment tables and published MSES (2010)\nPhysical trips alter displacement; low-lift fore-upper transition may remain free')
    save(fig,'trip_drag_comparison')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    colors=('#6b5191','#236890','#c16a2b','#3a845b')
    for case,color in zip(PRESSURES,colors):
        for ax,surface in zip(axes,('fore_lower','aft_upper')):
            rows=[r for r in pressure if r['figure']==case['figure'] and r['surface']==surface and not r['endpoint_duplicate']]
            ax.plot([r['table_x'] for r in rows],[r['Cp'] for r in rows],'o-',ms=4,color=color,
                label=f"CL {case['CL']:.3f} · α {case['alpha']:.2f}°")
    axes[0].set(xlim=(0,.81),title='Fore lower surface');axes[1].set(xlim=(.80,.98),title='Aft upper recovery region')
    visible=[r['Cp'] for r in pressure if r['surface']=='aft_upper' and .8<=r['table_x']<=.98]
    axes[1].set_ylim(min(visible)-.06,max(visible)+.06)
    for ax in axes:ax.invert_yaxis();ax.grid(alpha=.2);ax.set(xlabel='Table-II x / total chord',ylabel='Cp');ax.legend(fontsize=8)
    fig.suptitle('S414 free-transition approach to the lower bucket boundary\nRe 1 million · M 0.10 · measured taps; no exact bubble endpoints inferred')
    save(fig,'free_boundary_pressure')


def main():
    OUT.mkdir(exist_ok=True)
    require(sha(PDF)==PINNED[PDF.name],'Primary report changed')
    prior=json.loads((GEOMETRY/'verification.json').read_text())
    for name in ('pressure/experiment.csv','pressure/source.json'):
        require(sha(GEOMETRY/name)==prior['output_sha256'][name],'Prior pressure reference changed')
    textfile=GEOMETRY/'sources/s414_report_2010.txt';textpages=textfile.read_text().split('\f')
    doc=fitz.open(PDF);require(len(doc)==97,'Unexpected report length')
    curves={};tables=[];theory=[];markers=[];calibrations=[];crosschecks=[]
    for condition in CONDITIONS:
        for state,table_page,figure_page,count,mcount in zip(STATES,condition['tables'],condition['figures'],condition['table_counts'],condition['marker_counts']):
            text=doc[table_page-1].get_text(sort=True)
            rows=parse_table(text,condition,state,table_page);require(len(rows)==count,'Appendix row count changed')
            require(table_tokens(text)==table_tokens(textpages[table_page-1]),'Independent PDF text extractions disagree')
            th,mk,cal=extract_polar(doc[figure_page-1],condition,state,figure_page,mcount)
            for i,marker in enumerate(mk):
                reference=rows[i];dl=marker['CL']-reference['CL'];dd=marker['CD']-reference['CD']
                require(abs(dl)<6e-4 and abs(dd)<5e-6,'Figure/table association failed')
                crosschecks.append(dict(reynolds=condition['reynolds'],configuration=state,index=i,
                    table_page=table_page,figure_page=figure_page,delta_CL=dl,delta_CD=dd))
            tables+=rows;theory+=th;markers+=mk;calibrations.append(cal)
            branch=ascending_branch(rows);require(all(r['CD'] is not None for r in branch),'Missing pre-stall drag forbids interpolation')
            curves[(condition['reynolds'],'experiment',state)]=branch;curves[(condition['reynolds'],'mses',state)]=th
    pressure,pchecks,recovery,duplicates=extract_pressure(doc)
    targets=[delta_at(curves,c['reynolds'],cl) for c in CONDITIONS for cl in TARGETS]
    continuous=[];coverage=[]
    for condition in CONDITIONS:
        re=condition['reynolds'];cc=[curves[(re,m,s)] for m in ('experiment','mses') for s in STATES]
        low=max(-.05,max(r[0]['CL'] for r in cc));high=min(1.2,min(r[-1]['CL'] for r in cc))
        # Union of original knots preserves every piecewise-linear difference without smoothing.
        knots=sorted({low,high,*[r['CL'] for c in cc for r in c if low<=r['CL']<=high]})
        continuous.extend(delta_at(curves,re,cl) for cl in knots)
        coverage.append(dict(reynolds=re,common_plot_CL_min=low,common_plot_CL_max=high,
            configurations=3,methods=2,targets=len(TARGETS)))
    output=[('experimental_tables.csv',tables),('published_mses_drag.csv',theory),('experimental_plot_markers.csv',markers),
            ('table_plot_crosschecks.csv',crosschecks),('matched_lift_targets.csv',targets),('matched_lift_curves.csv',continuous),
            ('free_boundary_pressure.csv',pressure),('pressure_tap_checks.csv',pchecks),('pressure_recovery_intervals.csv',recovery),
            ('repeated_pressure_crosscheck.csv',duplicates)]
    for name,rows in output:write_csv(OUT/name,rows)
    plots(curves,continuous,pressure)
    # Inventory the complete report so an absent fixed-transition pressure comparison is scoped.
    inventory=[]
    for i,page in enumerate(doc):
        text=' '.join(page.get_text().split())
        inventory.append(dict(pdf_page=i+1,printed_page=i-1,caption_text=text if i>=18 else '',
            has_pressure_caption='pressure distribution' in text.lower(),has_fixed_label='transition fixed' in text.lower()))
    write_json(OUT/'report_page_inventory.json',inventory)
    evidence=dict(primary_source_url='https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf',
        report='Somers and Maughmer, RDECOM TR 10-D-112, August 2010',
        primary_condition=dict(reynolds=1000000,mach=.1,reference_chord_m=.4572),secondary_condition=dict(reynolds=700000,mach=.07),
        trip_table_pdf_page=21,methods_pdf_page=11,trip_confounds_pdf_page=13,comparison_qualification_pdf_page=14,
        trips=[dict(element='fore',side='upper',x_over_total_chord=.02,type='grit',
                    size_mm_by_Re={'700000':.178,'1000000':.124},grit_number_by_Re={'700000':90,'1000000':120}),
               dict(element='fore',side='lower',x_over_total_chord=.10,type='serrated tape',size_mm_by_Re={'700000':.457,'1000000':.457}),
               dict(element='aft',side='upper',x_over_aft_chord=.05,type='serrated tape',size_mm_by_Re={'700000':.191,'1000000':.114}),
               dict(element='aft',side='lower',x_over_aft_chord=.10,type='serrated tape',size_mm_by_Re={'700000':.686,'1000000':.610})],
        confounds=['Author reports excessive displacement caused by large trips.',
                   'At low lift the fore-upper grit can be too small to force transition; threshold CL not tabulated.',
                   'At fixed lift the incidences and mean loading differ.',
                   'Adding aft trips changes its friction, pressure/displacement and bubble; delta is not a measured bubble-drag component.'],
        pressure_coverage=dict(experimental_free_figures='7(a–x)',published_MSES_free_figures='12(a–c)',
            fore_only_fixed_Cp_found_in_report=False,both_fixed_Cp_found_in_report=False,
            measured_case_specific_bubble_endpoints_found=False,wake_turbulence_profiles_found=False,
            oil_flow='Qualitative aft-bubble shortening near the lower free-transition bucket boundary; no paired tripped endpoints.'),
        excluded_conditions=[dict(reynolds=500000,reason='Outside bounded primary/secondary qualification'),
            dict(reynolds=1500000,reason='Free M0.17 versus fixed M0.16; not an exactly matched Mach intervention')],
        mses=dict(figures='13–15',transition_method='2010 envelope e^N',ncrit=9,
            original_run_deck_recovered=False,new_executable_run=False,actual_trip_displacement_equivalence_established=False),
        admission=dict(matched_global_force_context=True,clean_wake_disturbance_intervention=False,
            matched_tripped_aft_pressure_benchmark=False,exact_geometry=False,closure_fit_justified=False))
    write_json(OUT/'evidence_qualification.json',evidence)
    peaks=[]
    for case in PRESSURES:
        for coordinate in ('table_x','plotted_x'):
            rr=[r for r in recovery if r['figure']==case['figure'] and r['coordinate']==coordinate]
            peaks.append(max(rr,key=lambda r:r['delta_Cp']))
    summary=dict(complete=True,experiment_table_rows=len(tables),experimental_curves=6,published_mses_curves=6,
        published_mses_vertices=len(theory),independent_plot_table_checks=len(crosschecks),
        maximum_plot_table_CL_difference=max(abs(r['delta_CL']) for r in crosschecks),
        maximum_plot_table_CD_difference=max(abs(r['delta_CD']) for r in crosschecks),
        free_pressure_markers=len(pressure),repeated_pressure_checks=len(duplicates),
        maximum_repeated_pressure_Cp_difference=max(abs(r['delta_Cp']) for r in duplicates),
        maximum_repeated_pressure_x_difference=max(abs(r['delta_plotted_x']) for r in duplicates),
        strict_repeated_pressure_copy_check_pass=all(r['strict_copy_check_pass'] for r in duplicates),
        strongest_sampled_pressure_recovery=peaks,coverage=coverage,targets=targets,
        new_flow_solves=0,causal_wake_transfer_benchmark_admitted=False,
        qualification='Observed total-force contrasts with documented trip confounds; no tripped aft pressure/endpoints in the report.')
    write_json(OUT/'summary.json',summary)
    paths=[PDF,textfile,Path(__file__),ROOT/'validation/S414_TRIP_QUALIFICATION_PROTOCOL.md',
        ROOT/'validation/qualify_s414_geometry.py',ROOT/'validation/extract_s702.py',ROOT/'validation/extract_s414_drag.py',
        GEOMETRY/'pressure/experiment.csv',GEOMETRY/'pressure/source.json',GEOMETRY/'verification.json']
    write_json(OUT/'source_manifest.json',dict(product=evidence['report'],version='2010 primary report',
        authoritative_url=evidence['primary_source_url'],access_mode='Reused public PDF; pinned SHA-256 verified',
        download_status='Previously completed; no new acquisition',expected_sha256=PINNED[PDF.name],source_supplied_sha256=None,
        scope=dict(table_pages=[88,89,90,91,92,93],polar_pages=[70,71,74,75,78,79],pressure_pages=[33,34,35,36]),
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},runtime=dict(python=sys.version,numpy=np.__version__,pymupdf=fitz.VersionBind),
        calibration=calibrations,verified_blockers=evidence['admission'],
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.suffix in ('.csv','.png','.svg','.pdf')}))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('targets','strongest_sampled_pressure_recovery')},indent=2))
    for r in targets:
        print(r['reynolds'],r['CL'],'fore-free experiment/MSES counts',r['experiment_fore_minus_free_counts'],r['mses_fore_minus_free_counts'])


if __name__=='__main__':main()
