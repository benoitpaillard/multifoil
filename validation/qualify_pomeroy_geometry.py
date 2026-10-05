#!/usr/bin/env python3
"""Extract published Pomeroy coordinates and check placement consistency; no CFD."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from multifoil import repanel, transform, validate_geometry

SOURCE=ROOT/'validation/results/reference_review_20260924/sources/pomeroy2014.pdf'
SOURCE_HASH='c6235093e0812173d50f04987dc5ade1c17b8e2b61dafdb41767b194cdc4bff3'
SOURCE_URL='https://m-selig.ae.illinois.edu/pubs/PomeroyDieboldAnselSelig-2014-JofAC-BurstWakes-FINAL.pdf'
COMPANION=SOURCE.with_name('pomeroy2015.pdf')
COMPANION_HASH='96ef7d349d8be463ebf61936acc0e8ef395eb2d6f12c17c13038d871e124604f'
NAMES=('main','flap1','flap2')
PLACEMENT=dict(main=dict(chord=.7028,origin=[0.,0.],deflection_deg=0.),
               flap1=dict(chord=.2139,origin=[.6860,-.0344],deflection_deg=26.43),
               flap2=dict(chord=.1834,origin=[.8644,-.1510],deflection_deg=42.32))
# Poppler page coordinates, inspected against the rendered primary tables.
BLOCKS=dict(main=[(0,'left',80,750,74),(0,'right',60,140,7)],
            flap1=[(0,'right',200,750,60),(1,'left',60,255,21)],
            flap2=[(1,'left',335,750,46),(1,'right',60,380,35)])
NS={'h':'http://www.w3.org/1999/xhtml'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def require(ok,message):
    if not ok:raise ValueError(message)


def extract(bbox):
    pages=ET.parse(bbox).findall('.//h:page',NS);require(len(pages)==2,'Expected PDF pages 9 and 10')
    local,provenance={},[]
    for name in NAMES:
        points=[]
        for page,side,low,high,count in BLOCKS[name]:
            xr,yr=((100,140),(175,220)) if side=='left' else ((365,402),(440,477))
            words=[]
            for w in pages[page].findall('.//h:word',NS):
                if re.fullmatch(r'[-−]?\d\.\d{4}',w.text or ''):
                    words.append((float(w.attrib['xMin']),float(w.attrib['yMin']),w.text))
            x=sorted((y,t,px) for px,y,t in words if xr[0]<px<xr[1] and low<y<high)
            y=sorted((py,t,px) for px,py,t in words if yr[0]<px<yr[1] and low<py<high)
            require(len(x)==len(y)==count,f'Coordinate block length differs: {name}/{page}/{side}')
            for a,b in zip(x,y):
                # A minus glyph expands the word bounding box upward by ~2.13 pt.
                require(abs(a[0]-b[0])<2.2,'Coordinate row pairing differs')
                xx,yy=float(a[1].replace('−','-')),float(b[1].replace('−','-'))
                provenance.append(dict(element=name,index=len(points),pdf_page=page+9,printed_page=page+829,
                    column=side,x_string=a[1],y_string=b[1],x=xx,y=yy,x_word_left=a[2],x_word_top=a[0],
                    y_word_left=b[2],y_word_top=b[0]))
                points.append([xx,yy])
        xy=np.asarray(points);require(xy.shape==(81,2),'Expected 81 coordinates per element')
        require(xy[0,0]==xy[-1,0]==1. and xy[0,1]>0 and xy[-1,1]<0,'Table ordering/TE differs')
        require(np.all(np.isfinite(xy)) and len(np.unique(xy,axis=0))==81,'Invalid/duplicated coordinate')
        local[name]=xy
    return local,provenance


def assemble(local):
    return {n:transform(local[n],scale=p['chord'],angle_deg=-p['deflection_deg'],translation=p['origin'])
            for n,p in PLACEMENT.items()}


def point_to_contour(point,xy):
    delta=np.roll(xy,-1,axis=0)-xy;size=np.sum(delta*delta,axis=1)
    require(np.all(size>0),'Degenerate segment in gap check')
    t=np.clip(np.sum((point-xy)*delta,axis=1)/size,0.,1.)
    feet=xy+t[:,None]*delta;dist=np.linalg.norm(point-feet,axis=1);i=int(np.argmin(dist))
    return float(dist[i]),feet[i]


def rigging_checks(local):
    assembled=assemble(local);rows=[]
    for i in (0,1):
        first,second=NAMES[i:i+2];body=assembled[second]
        tips=dict(upper=assembled[first][0],midpoint=.5*(assembled[first][0]+assembled[first][-1]),lower=assembled[first][-1])
        angle=np.deg2rad(-PLACEMENT[first]['deflection_deg']);axis=np.array([np.cos(angle),np.sin(angle)])
        overhang=float((np.asarray(PLACEMENT[second]['origin'])-tips['midpoint'])@axis)
        for samples in (80,800,1600):
            target=body if samples==80 else assemble({n:repanel(xy,samples) for n,xy in local.items()})[second]
            for endpoint,tip in tips.items():
                gap,foot=point_to_contour(tip,target)
                rows.append(dict(first=first,second=second,target_panels=samples,TE_definition=endpoint,
                    gap=gap,table2_gap=.0100,gap_minus_table2=gap-.0100,
                    nominal_overhang=overhang,table2_overhang=-.0150,overhang_minus_table2=overhang+.0150,
                    closest_x=float(foot[0]),closest_y=float(foot[1])))
    return rows


def qualify(out,render=False):
    require(sha(SOURCE)==SOURCE_HASH,'Primary PDF changed')
    require(sha(COMPANION)==COMPANION_HASH,'Companion PDF changed')
    out.mkdir(parents=True,exist_ok=True)
    bbox=out/'coordinates_bbox.html'
    subprocess.run(['pdftotext','-f','9','-l','10','-bbox-layout',str(SOURCE),str(bbox)],check=True)
    local,provenance=extract(bbox);assembled=assemble(local)
    # Local contours share an origin before assembly; validate them individually.
    for n,xy in local.items():validate_geometry([(n,xy)])
    validate_geometry(list(assembled.items()))
    for name in NAMES:
        np.savetxt(out/(name+'_local.dat'),local[name],header=name+' local chord; Pomeroy 2014 Table A'+str(NAMES.index(name)+1),comments='',fmt='%.4f')
        np.savetxt(out/(name+'_assembled.dat'),assembled[name],header=name+' / csys; Pomeroy 2014 Table 1 placement; NOT QUALIFIED',comments='',fmt='%.12g')
    checks=rigging_checks(local)
    for filename,rows in (('coordinate_provenance.csv',provenance),('rigging_checks.csv',checks)):
        with (out/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    selected=[r for r in checks if r['target_panels']==80 and r['TE_definition']=='midpoint']
    chord=float(.5*(assembled['flap2'][0,0]+assembled['flap2'][-1,0]))
    summary=dict(coordinates_per_element={n:len(local[n]) for n in NAMES},coordinate_extraction_pass=True,
        assembled_geometry_nonintersecting=True,assembly_transform='xy/csys = nominal_origin + chord_ratio * R(-deflection) * xy/local_chord',
        table1_placement=PLACEMENT,projected_TE_midpoint_chord=chord,primary_rigging_checks=selected,
        minimum_gap_across_TE_and_sampling_checks=min(r['gap'] for r in checks),
        maximum_gap_interpolation_change=max(abs(next(v['gap'] for v in checks if v['first']==r['first'] and v['TE_definition']==r['TE_definition'] and v['target_panels']==80)-r['gap']) for r in checks),
        geometry_benchmark_admitted=False,flow_runs=0,
        reason='Published absolute placement plus tabulated contours does not reproduce published relative gaps/overhangs; no geometry was fitted to fix it.',
        companion_2015_config10=dict(flap1_origin=[.6907,-.0178],flap1_absolute_angle=26.4,
            flap2_origin=[.8619,-.1352],flap2_absolute_angle=42.3,
            qualification='Different absolute coordinates; its positive overhang is upstream, opposite to the 2014 sign. Not substituted into the 2014 test.'),
        experimental_conditions=dict(reynolds=1e6,flowfield_alpha_deg=0.,reference_chord_body_inches=16.36,
            reference_chord_body_m=16.36*.0254,mach='not specified in reviewed methods',
            transition='No Ncrit or trip specification located in reviewed methods; requires qualification'),
        observable_compatibility=dict(global_CL_CD='Potentially comparable after geometry/conditions and figure digitization are qualified.',
            surface_attachment='Qualitative oil visualization only; use its limitations.',
            wake_core_velocity_profiles='No directly matching reconstructed viscous wake cross-profile in present solver.',
            reynolds_stress_TKE_spectra='Not predicted by the steady integral model.',
            companion_static_total_pressure='2015 off-body maps; do not relabel them as measured surface Cp or compare to wake-edge speed.'),
        next_action='Reconcile paper-specific placement from original conference/thesis data or an author geometry file; independently qualify S103/S414 tested geometry for MSES/experiment comparison.')
    write(out/'summary.json',summary)
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for name,color in zip(NAMES,('#3c79a4','#c57925','#7655a0')):
        xy=assembled[name];axes[0].plot(xy[:,0],xy[:,1],'.-',ms=2,label=name,color=color)
    for row in selected:
        body=assembled[row['first']];tip=.5*(body[0]+body[-1]);foot=np.array([row['closest_x'],row['closest_y']])
        axes[0].plot([tip[0],foot[0]],[tip[1],foot[1]],color='red',lw=2)
    axes[0].set_aspect('equal',adjustable='box');axes[0].set_xlabel('x / published csys');axes[0].set_ylabel('y / published csys')
    axes[0].set_title('Table A coordinates + Table 1 placement');axes[0].legend();axes[0].grid(alpha=.2)
    x=np.arange(2);bars=axes[1].bar(x,[r['gap'] for r in selected],color='#3c79a4',width=.55)
    axes[1].bar_label(bars,labels=[f"{r['gap']:.5f}" for r in selected],padding=4)
    axes[1].axhline(.01,color='#b54536',ls='--',label='Table 2: 0.0100');axes[1].legend()
    axes[1].set_xticks(x,['main → flap 1','flap 1 → flap 2']);axes[1].set_ylabel('TE-midpoint to next contour / csys');axes[1].set_ylim(0,.028)
    axes[1].set_title('Unresolved gap discrepancy');axes[1].grid(axis='y',alpha=.2)
    fig.suptitle('Pomeroy 2014 geometry qualification · 81 published points per element\n'
                 'No aerodynamic solve; original placement retained without fitting.',fontsize=11)
    for ext in ('png','svg','pdf'):fig.savefig(out/('geometry_qualification.'+ext),dpi=180)
    plt.close(fig)
    if render:
        for page,label in ((2,822),(3,823),(9,829),(10,830)):
            subprocess.run(['pdftoppm','-f',str(page),'-l',str(page),'-scale-to','1600','-png','-singlefile',str(SOURCE),str(out/f'source_page_{label}')],check=True)
        subprocess.run(['pdftoppm','-f','4','-l','4','-scale-to','1600','-png','-singlefile',str(COMPANION),str(out/'companion_page_4')],check=True)
    notes(out,summary)
    original=json.loads((SOURCE.parent.parent/'source_manifest.json').read_text())
    sources=[v for v in original['sources'] if v.get('path') in (str(SOURCE.relative_to(ROOT)),str(COMPANION.relative_to(ROOT)))]
    write(out/'source_manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),product='Pomeroy 2014 Table A1–A3 extraction and Table 1/Table 2 consistency check',
        sources=sources,access_mode='Previously downloaded anonymous public author PDFs; local hash verification and extraction',
        scope='2014: pages 822–823 and 829–830; 2015: geometry context only, Tables 1–2 page 4',
        new_downloads=0,supplied_checksums=None,extraction_tool=subprocess.run(['pdftotext','-v'],text=True,capture_output=True).stderr.splitlines()[0],
        extraction_script_sha256=sha(Path(__file__)),geometry_helper_sha256=sha(ROOT/'multifoil.py'),
        files={str(p.relative_to(ROOT)):sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name not in ('source_manifest.json','verification.json','checks.log')}))
    return summary


def notes(out,s):
    rows='\n'.join(f"| {r['first']} → {r['second']} | {r['gap']:.8f} | {r['table2_gap']:.4f} | {r['nominal_overhang']:.8f} | {r['table2_overhang']:.4f} |" for r in s['primary_rigging_checks'])
    (out/'FINDINGS.md').write_text(
        '# Pomeroy geometry and observable qualification\n\n'
        'Recovered all 243 tabulated coordinates (81 per element) from Tables A1–A3 of the '
        '[author-hosted 2014 paper]('+SOURCE_URL+'). PDF pages 9–10 correspond to printed pages 829–830. '
        'Numeric word positions, signs, row order and split-table continuations are retained in coordinate_provenance.csv. '
        'The local and assembled contours pass self-intersection and inter-element intersection checks. '
        'The assembled geometry is a traceable candidate, not an admitted experimental benchmark.\n\n'
        '## Placement inconsistency\n\n'
        'Apply the tabulated chord ratios, downward-positive angles (clockwise rotations) and leading-edge origins '
        'directly; do not shift each local minimum-x point to zero or normalize the assembled chord afterward. The '
        f"projected final TE-midpoint chord is {s['projected_TE_midpoint_chord']:.9f}, consistent with unit normalization to tabular precision. "
        'Gap is the distance from the preceding TE midpoint to the closest segment on the next element; signed overhang '
        'uses the nominal next leading-edge origin projected along the preceding nominal chord, as in the 2014 definition.\n\n'
        '| Pair | Reconstructed gap/csys | Table 2 gap/csys | Reconstructed overhang/csys | Table 2 overhang/csys |\n'
        '|---|---:|---:|---:|---:|\n'+rows+'\n\n'
        'These are calculated inconsistencies between the printed absolute and relative descriptions under the stated '
        'transform, not evidence that a particular table is the correct tunnel rigging. Both alternative TE endpoints '
        'and 800/1600-panel PCHIP resampling were checked. The smallest gap across those checks is '
        f"{s['minimum_gap_across_TE_and_sampling_checks']:.8f}c and the largest resampling change is {s['maximum_gap_interpolation_change']:.3g}c. "
        'Neither TE thickness nor this contour-faceting check reconciles the gaps. No geometry was adjusted to force '
        'agreement and no aerodynamic solve was launched. Nominal local origins and rounded coordinate tables remain explicit assumptions.\n\n'
        'The [2015 companion](https://m-selig.ae.illinois.edu/pubs/PomeroySelig-2015-AIAA-2015-2569-ThreeElementAirfoilWakes.pdf) '
        'uses different absolute coordinates for baseline Configuration 10: flap 1 (0.6907,-0.0178), 26.4°; flap 2 '
        '(0.8619,-0.1352), 42.3°. It also defines positive overhang as upstream, whereas the 2014 table/diagram uses '
        'negative upstream. The 2015 placement must not silently replace the 2014 placement. Reconcile each paper\'s '
        'case with its original conference/thesis geometry or author coordinate deck before quantitative scoring.\n\n'
        '## What can be compared\n\n'
        'The 2014 paper reports Re=1,000,000 and alpha=0° for the detailed wake measurements. Its body defines csys '
        'as 16.36 in (0.415544 m); the abstract uses a different rounded chord (1.35 ft / 0.411 m). Preserve the body\'s '
        'explicit geometry convention rather than deriving normalization from the abstract. Mach and a trip/Ncrit '
        'specification were not established from the reviewed methods. Do not infer a measured Ncrit from low tunnel turbulence.\n\n'
        '* Global lift and profile drag in Figure 8 are potential comparison quantities after rigging, conditions and '
        'digitization are qualified; they have not been extracted or scored here.\n'
        '* Oil-flow attachment is qualitative evidence, not a measured force or transition curve.\n'
        '* The 2014 split-film velocity profiles/core widths, Reynolds stresses and spectra do not directly match the '
        'present steady integral solver\'s wake-edge speed or Squire–Young drag output. No viscous wake cross-profile '
        'reconstruction is provided by the current model.\n'
        '* The 2015 static/total-pressure maps are off-body probe measurements. They are not surface Cp; matching '
        'spatial locations, viscous deficits and pressure definitions must be established before comparison.\n\n'
        'Geometry extraction is complete; experimental benchmark admission remains unresolved. The new evidence lowers '
        'Pomeroy\'s readiness despite its valuable coordinate tables. Continue exact tested-geometry qualification for '
        'S103/S414 to retain a direct experiment/MultiFoil/published-MSES path while investigating Pomeroy\'s rigging. '
        'No physical closure or transition constant has been tuned.\n\n'
        '![Reconstructed geometry and gap check](geometry_qualification.png)\n\n'
        'Reproduce: `OPENBLAS_NUM_THREADS=1 python3 validation/qualify_pomeroy_geometry.py --render`. '
        'Requires existing Poppler pdftotext/pdftoppm and the project Python dependencies. Original PDFs remain unchanged; '
        'source_manifest.json records original URLs/access records, hashes, extraction scope and derived files.\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/pomeroy_geometry_qualification_20261001')
    parser.add_argument('--render',action='store_true')
    args=parser.parse_args();print(json.dumps(qualify(args.out,args.render),indent=2))
