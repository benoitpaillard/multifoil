#!/usr/bin/env python3
"""Independent check of the V3 cross-formulation outputs through the public solver entry points."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import csv
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from report_design_initialization import sha, write
from report_s702_drag_mechanism import require
from multifoil import PanelSolver, HessSmithSolver
from viscous import CoupledSolver
from s414_geometry import elements

OUT = ROOT/'validation/results/s414_formulation_check_20261004'


def table(name):
    with (OUT/name).open(newline='') as stream: return list(csv.DictReader(stream))


def main():
    require(not (OUT/'verification.json').exists(), 'Preserve completed verification')
    summary = json.loads((OUT/'summary.json').read_text()); require(summary['complete'], 'Incomplete check')
    manifest = json.loads((OUT/'input_manifest.json').read_text())['files']
    changed = [name for name, digest in manifest.items() if sha(ROOT/name) != digest]
    require(not changed, 'Inputs changed since the check: '+str(changed))
    # Public entry points, not the helper used by the report: full production constructor and PanelSolver.
    solver = CoupledSolver(elements('nominal', [340, 240]), reynolds=1e6, mach=.1, ncrit=9.)
    panel = PanelSolver(solver.elements); require(type(panel) is HessSmithSolver, 'Finite trailing edges must select Hess–Smith')
    profiles = [r for r in table('fore_nose_profiles.csv') if r['scenario'] == 'nominal' and r['panels_fore'] == '340']
    checks = []
    for label in ('measured_cl0056', 'matched_cl0056', 'bracket_lower'):
        rows = [r for r in profiles if r['incidence'] == label]; alpha = float(rows[0]['alpha_deg'])
        u = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
        gamma = solver.models[0].isol.gamref@u; mid = .5*(gamma[:-1]+gamma[1:])
        vt = panel.solve(alpha).surfaces[0]['ue_signed']; n = solver.models[0].foil.N; hess = -vt[:n-1][::-1]
        index = np.array([int(r['panel']) for r in rows])
        a = float(np.max(abs(mid[index]-np.array([float(r['production_speed']) for r in rows]))))
        b = float(np.max(abs(hess[index]-np.array([float(r['hess_speed']) for r in rows]))))
        checks.append(dict(incidence=label, alpha_deg=alpha, panels_compared=len(rows),
                           production_speed_maximum_difference=a, hess_speed_maximum_difference=b))
        require(max(a, b) < 1e-11, 'Independent entry points differ from reported speeds')
    differences = table('s414_differences.csv'); refinement = table('refinement.csv'); controls = table('controls.csv')
    base = [d for d in differences if d['panels_fore'] == '340']
    require(abs(max(float(d['fore_nose_max_Cp_difference']) for d in base)-summary['baseline']['max_fore_nose_Cp_difference']) < 1e-15,
            'Summary differs from table')
    require(sum(r['agree_within_resolved_spread'] == 'True' for r in refinement) == summary['cases_agreeing_within_resolved_spread']
            and sum(r['between_methods_decreasing'] == 'True' for r in refinement) == summary['cases_with_decreasing_difference'],
            'Refinement counts differ')
    require(len(differences) == 30 and len(refinement) == 10 and len(table('s414_cases.csv')) == 60, 'Coverage differs')
    kt = {(r['method'], r['panels']): r for r in controls if r['control'] == 'karman_trefftz'}
    order = {m: float(np.log(float(kt[(m, '240')]['rms_speed_error'])/float(kt[(m, '480')]['rms_speed_error']))/np.log(2))
             for m in ('production', 'hess')}
    te = table('te_treatment_diagnostic.csv')
    sharp = [r for r in te if float(r['te_factor']) == 0]; blunt = [r for r in te if float(r['te_factor']) == 1]
    require(abs(float(sharp[1]['relative_CL_difference'])) < abs(float(sharp[0]['relative_CL_difference']))
            and abs(float(blunt[1]['relative_CL_difference'])) > abs(float(blunt[0]['relative_CL_difference'])),
            'Trailing-edge diagnostic trend differs from the finding')
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', 'test_s414_formulation_check'], cwd=ROOT, text=True, capture_output=True)
    (OUT/'tests.log').write_text(tests.stdout+tests.stderr); require(tests.returncode == 0, 'Helper tests failed')
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, controls_qualified=summary['controls_qualified'],
        predeclared_agreement_criterion_met=summary['cases_agreeing_within_resolved_spread'] == summary['cases'],
        gate='open: the predeclared within-spread criterion is not met; disagreement is localized and bounded, not eliminated',
        independent_entry_point_checks=checks, observed_control_order=order, inputs_verified=len(manifest),
        tests=dict(module='test_s414_formulation_check', passed=True, count=5),
        figures_reviewed=['nose_pressure_formulations.png', 'refinement_trend.png'],
        figure_review='Rendered and inspected: no label collision; taps and unmeasured peak visible; legend outside the plot area',
        initial_run='initial_run_stopped/: the first run stopped on an over-strict guard before any S414 solve',
        new_nonlinear_flow_solves=0, output_sha256={p.name: sha(p) for p in outputs})
    write(OUT/'verification.json', receipt); print(json.dumps({k: receipt[k] for k in ('controls_qualified',
        'predeclared_agreement_criterion_met', 'independent_entry_point_checks', 'observed_control_order')}, indent=1))


if __name__ == '__main__':
    main()
