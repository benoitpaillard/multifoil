#!/usr/bin/env python3
"""Bounded inviscid incidence diagnostic using qualified report angles, without fitting."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from viscous import CoupledSolver, mf
from s414_geometry import elements
from report_s702_wake_extent import ROOT, Inputs, sha, require, csv_write, write

OUT = ROOT/'validation/results/s414_incidence_preflight_20261003'
PRIOR = ROOT/'validation/results/s414_lower_bucket_analysis_20261003'
CURRENT = ROOT/'validation/results/s414_jacobian_analysis_20261003'


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs()
    summary = inputs.json(CURRENT/'summary.json'); require(summary['full_state_audit_pass'], 'Current audit incomplete')
    reference = {r['case']: r for r in inputs.csv(PRIOR/'targets.csv')}
    taps = inputs.csv(PRIOR/'pressure_at_taps.csv')
    solver = CoupledSolver(elements('nominal'), reynolds=1e6, mach=.1, ncrit=9)
    comparisons, metrics = [], []
    for case in ('cl0056', 'cl028'):
        row = next(t for t in summary['targets'] if t['case'] == case and t['mode'] == 'corrected_af')
        ref = reference[case]
        angles = [('matched_lift', row['alpha_deg']), ('measured_caption', float(ref['experiment_caption_alpha'])),
                  ('measured_appendix', float(ref['appendix_alpha']))]
        for label, angle in angles:
            direction = np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))])
            for e, model in zip(('fore', 'aft'), solver.models):
                x = model.foil.x[0, ::-1]; le = int(np.argmin(x)); cp = mf.get_cp(model.isol.gamref@direction, model.param)[0][::-1]
                for side in ('upper', 'lower'):
                    bx, bc = (x[:le+1][::-1], cp[:le+1][::-1]) if side == 'upper' else (x[le:], cp[le:])
                    require(np.all(np.diff(bx) > 0), 'Pressure branch folded')
                    chosen = [t for t in taps if t['case'] == case and t['surface'] == e+'_'+side
                        and t['coordinate'] == 'table' and t['included'] == 'True' and t['method'] == 'multifoil']
                    values = []
                    for tap in chosen:
                        tx = float(tap['x_over_c']); require(bx[0] <= tx <= bx[-1], 'Unsupported tap')
                        value = float(np.interp(tx, bx, bc)); error = value-float(tap['experiment_Cp'])
                        entry = dict(case=case, incidence=label, alpha_deg=angle, surface=e+'_'+side,
                            table_index=int(tap['table_index']), x=tx, inviscid_Cp=value,
                            experiment_Cp=float(tap['experiment_Cp']), error_Cp=error)
                        comparisons.append(entry); values.append(entry)
                    for region, selected in [('all_supported', values), ('nose_x_le_006', [v for v in values if v['x'] <= .06])]:
                        if not selected: continue
                        metrics.append(dict(case=case, incidence=label, alpha_deg=angle, surface=e+'_'+side,
                            region=region, taps=len(selected), rmse_Cp=float(np.sqrt(np.mean([v['error_Cp']**2 for v in selected])))))
    csv_write(OUT/'pressure_at_taps.csv', comparisons); csv_write(OUT/'pressure_metrics.csv', metrics)
    for name in ('viscous.py', 'multifoil.py', 'paneling.py', 'vendor/mfoil.py', 'validation/s414_geometry.py',
                 'validation/check_s414_incidence.py', 'validation/results/s414_geometry_qualification_20261002/report_test_fore.dat',
                 'validation/results/s414_geometry_qualification_20261002/report_test_aft.dat',
                 'validation/results/s414_geometry_qualification_20261002/sources/s414_report_2010.txt'):
        inputs.add(ROOT/name)
    write(OUT/'input_manifest.json', dict(files=inputs.files))
    write(OUT/'summary.json', dict(complete=True, new_nonlinear_flow_solves=0, inviscid_operator_rebuilds=1,
        incidence_evaluations=6, comparison_rows=len(comparisons), metrics=metrics,
        interpretation='Changing to the measured incidence is a diagnostic, not a correction or a matched-lift validation; all supported taps retained.'))
    print('Incidence preflight complete; six linear evaluations, no nonlinear solves.')


if __name__ == '__main__': main()
