#!/usr/bin/env python3
"""Independent scalar quadrature and provenance checks; no flow solves."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_pressure_budget_20261003'
PRIOR = ROOT/'validation/results/s414_incidence_analysis_20261003'
STUDY = ROOT/'validation/results/s414_incidence_20261003'


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    if abs(float(a)-float(b)) > 1e-12:
        raise ValueError(f'Budget mismatch: {a} != {b}')


def trapezoid(points, field):
    return math.fsum(.5*(float(a[field])+float(b[field]))*(float(b['x_over_c'])-float(a['x_over_c']))
                     for a,b in zip(points,points[1:]))


def main():
    inputs=read(OUT/'input_manifest.json')['files']
    for name,digest in inputs.items():
        if sha(ROOT/name)!=digest: raise ValueError('Input changed: '+name)
    previous=read(PRIOR/'verification.json')
    for name,digest in previous['output_sha256'].items():
        if sha(PRIOR/name)!=digest: raise ValueError('Prior output changed: '+name)
    surface=rows(OUT/'surface_budgets.csv'); intervals=rows(OUT/'interval_budgets.csv')
    budgets=rows(OUT/'budgets.csv'); missing=rows(OUT/'uncovered_regions.csv')
    taps=rows(PRIOR/'pressure_at_taps.csv'); dense=rows(OUT/'dense_force_checks.csv')
    assert len(surface)==32 and len(budgets)==8 and len(missing)==80 and len(dense)==8
    for r in surface:
        selected=[s for s in taps if all(s[k]==r[k] for k in ('job','coordinate','surface')) and s['included']=='True']
        selected.sort(key=lambda s:float(s['x_over_c']))
        assert len(selected)==int(r['taps'])
        sign=-1 if r['surface'].endswith('upper') else 1
        for prefix,field in [('experiment','experiment_Cp'),('model','model_Cp')]:
            close(sign*trapezoid(selected,field),r[prefix+'_CN'])
        close(float(r['model_CN'])-float(r['experiment_CN']),r['delta_CN'])
        parts=[s for s in intervals if all(s[k]==r[k] for k in ('job','coordinate','surface'))]
        for key in ('experiment_CN','model_CN','delta_CN','experiment_CC_pressure','model_CC_pressure'):
            close(math.fsum(float(s[key]) for s in parts),r[key])
    for r in dense:
        folder=STUDY/r['job']; result=read(folder/'result.json'); i=result['selected_flow_index']
        s=rows(folder/f'flow/surface_{i:03d}_{r["element"]}.csv')
        # Explicit independent loop including the blunt TE face.
        cn=cc=0.
        for a,b in zip(s,s[1:]+s[:1]):
            cp=.5*(float(a['Cp'])+float(b['Cp']))
            cn+=cp*(float(b['x_over_cref'])-float(a['x_over_cref']))
            cc-=cp*(float(b['y_over_cref'])-float(a['y_over_cref']))
        close(cn,r['CN']); close(cc,r['CC_pressure'])
        alpha=math.radians(result['alpha_deg'])
        close(cn*math.cos(alpha)-cc*math.sin(alpha),result['element_forces'][r['element']]['CL'])
    targets={t['job']:t for t in read(PRIOR/'summary.json')['targets']}
    for r in budgets:
        t=targets[r['job']]; group=[s for s in surface if s['job']==r['job'] and s['coordinate']==r['coordinate']]
        gaps=[s for s in missing if s['job']==r['job'] and s['coordinate']==r['coordinate']]
        assert sum(int(s['taps']) for s in group)==int(r['taps'])==87
        close(sum(float(s['delta_CN']) for s in group),r['delta_partial_CN'])
        close(sum(float(s['experiment_CN']) for s in group),r['experiment_partial_CN'])
        close(sum(float(s['model_CN']) for s in gaps),r['model_uncovered_CN'])
        close(sum(float(s['model_sampling_correction_CN']) for s in group),r['model_sampling_correction_CN'])
        close(sum(float(s['CN']) for s in dense if s['job']==r['job']),r['model_full_CN'])
        close(float(r['model_partial_CN'])+float(r['model_sampling_correction_CN'])+float(r['model_uncovered_CN']),r['model_full_CN'])
        a=math.radians(t['reference_appendix_alpha'])
        expected=t['reference_appendix_CL']*math.cos(a)+t['reference_appendix_CD']*math.sin(a)
        close(expected,r['reported_appendix_CN'])
        close(float(r['experiment_partial_CN'])+float(r['experimental_unresolved_CN']),expected)
        close(float(r['model_full_CN'])-expected,r['full_model_minus_reported_CN'])
        for endpoint in ('low','high'):
            bound=sum(float(s['conditional_experiment_CN_'+endpoint]) for s in gaps)
            close(bound,r['conditional_missing_CN_'+endpoint])
        assert float(r['conditional_missing_CN_low'])<=float(r['experimental_unresolved_CN'])<=float(r['conditional_missing_CN_high'])
        if r['comparison_basis']=='fixed_alpha':
            assert float(r['model_full_CN'])>float(r['conditional_completed_exp_CN_high'])
    mses=rows(OUT/'published_mses_partial.csv')
    old=rows(ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv')
    assert len(mses)==8 and {r['case'] for r in mses}=={'cl028'}
    for r in mses:
        s=[s for s in old if s['method']=='mses' and s['included']=='True'
           and all(s[k]==r[k] for k in ('case','coordinate','surface'))]
        s.sort(key=lambda s:float(s['x_over_c']))
        sign=-1 if r['surface'].endswith('upper') else 1
        close(sign*trapezoid(s,'model_Cp'),r['published_mses_CN'])
    summary=read(OUT/'summary.json')
    assert summary['complete'] and summary['new_nonlinear_flow_solves']==0 and not summary['total_force_validation']
    visual=read(OUT/'visual_review.json')
    assert visual['passed'] and len(visual['images'])==4
    for name,digest in visual['images'].items(): assert sha(OUT/name)==digest
    tests=read(OUT/'tests.json')
    assert tests['returncode']==0 and tests['tests_run']==12
    assert sha(OUT/'tests.log')==tests['log_sha256']
    log=(OUT/'tests.log').read_text()
    assert 'Ran 12 tests in ' in log and log.rstrip().endswith('OK')
    links=0
    for name in re.findall(r'\]\(([^)]+)\)',(OUT/'FINDINGS.md').read_text()):
        if name.startswith(('http:','https:')) or name=='verification.json':continue
        assert (OUT/name.split('#')[0]).resolve().exists(),name
        links+=1
    receipt=dict(complete=True,verified_utc=datetime.now(timezone.utc).isoformat(),
        inputs_verified=len(inputs),prior_outputs_unchanged=len(previous['output_sha256']),
        independent_surface_quadratures=32,independent_mses_quadratures=8,
        dense_element_force_checks=8,budget_partitions_checked=8,tests_passed=12,
        visual_review_pass=True,local_report_links_checked=links,new_nonlinear_flow_solves=0,
        verifier_sha256=sha(Path(__file__)),
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'})
    write=OUT/'verification.json'; write.write_text(json.dumps(receipt,indent=2)+'\n')
    print({k:v for k,v in receipt.items() if k!='output_sha256'})


if __name__=='__main__': main()
