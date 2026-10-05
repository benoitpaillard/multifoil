#!/usr/bin/env python3
"""Verify saved linear equations and independently recompute pressure/load splits."""
import ast
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'validation/results/s414_inviscid_budget_20261003'
PRIOR=ROOT/'validation/results/s414_pressure_budget_20261003'
STUDY=ROOT/'validation/results/s414_incidence_20261003'


def read(path): return json.loads(path.read_text())
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):
    with path.open() as stream: return list(csv.DictReader(stream))
def close(a,b,tol=1e-12):
    if not (math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a)-float(b))<=tol):
        raise ValueError(f'Comparison differs: {a} vs {b}')
def trap(data,field):
    return math.fsum(.5*(float(a[field])+float(b[field]))*(float(b['x'])-float(a['x']))
                     for a,b in zip(data,data[1:]))
def cp(u,mach=.1):
    c=1-u*u;b=math.sqrt(1-mach*mach)
    return c/(b+mach*mach*c/(2*(1+b)))


def main():
    manifest=read(OUT/'input_manifest.json')
    for name,digest in manifest['files'].items():
        if sha(ROOT/name)!=digest: raise ValueError('Frozen input changed: '+name)
    prior=read(PRIOR/'verification.json')
    for name,digest in prior['output_sha256'].items():
        if sha(PRIOR/name)!=digest: raise ValueError('Previous output changed: '+name)
    summary=read(OUT/'summary.json')
    assert summary['complete'] and summary['cases']==2 and summary['operator_rebuilds']==1
    assert summary['inviscid_angle_evaluations']==2 and summary['new_nonlinear_flow_solves']==0
    assert not summary['numerical_convergence_established'] and not summary['pure_inviscid_experimental_reference']
    assert not summary['isolated_transition_or_wake_intervention'] and summary['new_mses_runs']==0
    nodes=rows(OUT/'nodes.csv'); taps=rows(OUT/'pressure_at_taps.csv')
    surfaces=rows(OUT/'surface_budgets.csv'); totals=rows(OUT/'totals.csv')
    regions=rows(OUT/'region_budgets.csv'); dense=rows(OUT/'dense_forces.csv')
    assert len(nodes)==1164 and len(taps)==348 and len(surfaces)==16 and len(totals)==4
    assert len(regions)==56 and len(dense)==4
    with np.load(OUT/'inviscid_operator.npz',allow_pickle=False) as operator:
        assert operator['matrix'].shape==(584,584) and operator['rhs'].shape==(584,2)
        residual=float(np.max(abs(operator['matrix']@operator['reference']-operator['rhs'])))
        close(residual,summary['matrix_residual']); assert residual<1e-12
        assert np.max(abs(operator['reference'][:582]-operator['production_gref']))<2e-11
        for job in sorted({r['job'] for r in nodes}):
            alpha=math.radians(next(float(r['alpha_deg']) for r in nodes if r['job']==job))
            expected=operator['production_gref']@np.array([math.cos(alpha),math.sin(alpha)])
            for e,lo,hi in [('fore',0,341),('aft',341,582)]:
                selected=[r for r in nodes if r['job']==job and r['element']==e]
                assert [int(r['node']) for r in selected]==list(range(hi-lo))
                xy=operator[e+'_xy'][:,::-1]
                with np.load(STUDY/job/f'{e}_state.npz',allow_pickle=False) as saved:
                    assert np.array_equal(operator[e+'_xy'],saved['xy'])
                    for i,r in enumerate(selected):
                        close(r['x'],xy[0,i]);close(r['y'],xy[1,i])
                        close(r['inviscid_signed_speed'],expected[lo:hi][::-1][i])
                        close(r['viscous_raw_speed'],saved['U'][3,:hi-lo][::-1][i])
                        close(r['inviscid_Cp'],cp(float(r['inviscid_signed_speed'])))
                        close(r['viscous_Cp'],cp(float(r['viscous_raw_speed'])))
                        close(float(r['viscous_Cp'])-float(r['inviscid_Cp']),r['viscous_change_Cp'])
                result=read(STUDY/job/'result.json'); d=next(r for r in dense if r['job']==job and r['element']==e)
                for mode in ('inviscid','viscous'):
                    cn=cc=0.
                    for a,b in zip(selected,selected[1:]+selected[:1]):
                        mean=.5*(float(a[mode+'_Cp'])+float(b[mode+'_Cp']))
                        cn+=mean*(float(b['x'])-float(a['x']))
                        cc-=mean*(float(b['y'])-float(a['y']))
                    close(cn,d[mode+'_CN']);close(cc,d[mode+'_CC_pressure'])
                    close(cn*math.cos(alpha)-cc*math.sin(alpha),d[mode+'_CL_pressure'])
                close(d['viscous_CL_pressure'],result['element_forces'][e]['CL'])
    for r in taps:
        for field,a,b in [('inviscid_error_Cp','inviscid_Cp','experiment_Cp'),
                          ('viscous_change_Cp','viscous_Cp','inviscid_Cp'),('viscous_error_Cp','viscous_Cp','experiment_Cp')]:
            close(float(r[a])-float(r[b]),r[field])
        close(float(r['inviscid_error_Cp'])+float(r['viscous_change_Cp']),r['viscous_error_Cp'])
    for r in surfaces:
        points=[t for t in taps if all(t[k]==r[k] for k in ('job','coordinate','surface'))]
        points.sort(key=lambda t:float(t['x'])); sign=-1 if r['surface'].endswith('upper') else 1
        assert len(points)==int(r['taps'])
        for mode in ('experiment','inviscid','viscous'):close(sign*trap(points,mode+'_Cp'),r[mode+'_CN'])
        for field,a,b in [('inviscid_error_CN','inviscid_CN','experiment_CN'),
                          ('viscous_change_CN','viscous_CN','inviscid_CN'),('viscous_error_CN','viscous_CN','experiment_CN')]:
            close(float(r[a])-float(r[b]),r[field])
            parts=[x for x in regions if all(x[k]==r[k] for k in ('job','coordinate','surface'))]
            close(math.fsum(float(x[field]) for x in parts),r[field])
    old=rows(PRIOR/'budgets.csv')
    for r in totals:
        s=[s for s in surfaces if all(s[k]==r[k] for k in ('job','coordinate'))]
        assert sum(int(x['taps']) for x in s)==87
        previous=next(x for x in old if all(x[k]==r[k] for k in ('job','coordinate')))
        for mode in ('inviscid','viscous'):
            close(math.fsum(float(x[mode+'_CN']) for x in s),r[mode+'_partial_CN'])
            close(math.fsum(float(x[mode+'_CN']) for x in dense if x['job']==r['job']),r[mode+'_full_CN'])
            close(float(r[mode+'_partial_CN'])+float(r[mode+'_sampling_correction_CN'])+float(r[mode+'_uncovered_CN']),r[mode+'_full_CN'])
            close(float(r[mode+'_full_CN'])-float(r['reported_appendix_CN']),r[mode+'_error_full_CN'])
        close(float(r['inviscid_error_partial_CN'])+float(r['viscous_change_partial_CN']),r['viscous_error_partial_CN'])
        close(float(r['inviscid_error_full_CN'])+float(r['viscous_change_full_CN']),r['viscous_error_full_CN'])
        close(r['viscous_error_partial_CN'],previous['delta_partial_CN'])
        close(r['viscous_error_full_CN'],previous['full_model_minus_reported_CN'])
        close(r['experimental_unresolved_CN'],previous['experimental_unresolved_CN'])
    for r in rows(OUT/'pressure_metrics.csv'):
        selected=[t for t in taps if all(t[k]==r[k] for k in ('job','coordinate')) and (r['surface']=='all_surfaces' or t['surface']==r['surface'])]
        errors=[float(t[r['method']+'_Cp'])-float(t['experiment_Cp']) for t in selected]
        assert len(selected)==int(r['taps'])
        close(math.sqrt(math.fsum(e*e for e in errors)/len(errors)),r['rmse_Cp'])
    test=read(OUT/'tests.json'); assert test['returncode']==0 and test['tests_run']==18
    assert sha(OUT/'tests.log')==test['log_sha256']
    assert 'Ran 18 tests in ' in (OUT/'tests.log').read_text() and (OUT/'tests.log').read_text().rstrip().endswith('OK')
    visual=read(OUT/'visual_review.json'); assert visual['passed'] and len(visual['images'])==3
    for name,digest in visual['images'].items(): assert sha(OUT/name)==digest
    tree=ast.parse((ROOT/'validation/report_s414_inviscid_budget.py').read_text())
    calls=[node.func for node in ast.walk(tree) if isinstance(node,ast.Call)]
    assert not any(isinstance(f,ast.Attribute) and f.attr in ('_newton','_solve_point','_initialize','_build_wakes') for f in calls)
    assert not any(isinstance(f,ast.Attribute) and f.attr=='solve' and isinstance(f.value,ast.Name) and f.value.id=='solver' for f in calls)
    links=0
    for name in re.findall(r'\]\(([^)]+)\)',(OUT/'FINDINGS.md').read_text()):
        if name.startswith(('http:','https:')) or name=='verification.json':continue
        assert (OUT/name.split('#')[0]).resolve().exists(),name
        links+=1
    receipt=dict(complete=True,verified_utc=datetime.now(timezone.utc).isoformat(),
        frozen_inputs_verified=len(manifest['files']),previous_budget_outputs_unchanged=len(prior['output_sha256']),
        linear_matrix_residual=residual,independent_nodal_pressure_checks=len(nodes),
        independent_surface_quadratures=16,independent_dense_element_checks=4,complete_budgets=4,
        tests_passed=18,visual_review_pass=True,local_report_links_checked=links,new_nonlinear_flow_solves=0,
        verifier_sha256=sha(Path(__file__)),
        output_sha256={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='verification.json'})
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print({k:v for k,v in receipt.items() if k!='output_sha256'})


if __name__=='__main__':main()
