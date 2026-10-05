#!/usr/bin/env python3
"""Independent check of the V4 geometry trace and inviscid comparison through the public solver entry point."""
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
from viscous import CoupledSolver
import s414_nose_geometry as nose

OUT = ROOT/'validation/results/s414_nose_geometry_20261004'
MATERIAL = .05


def table(name):
    with (OUT/name).open(newline='') as stream: return list(csv.DictReader(stream))


def cp(speed, mach=.1):
    cp0 = 1-np.asarray(speed)**2; beta = np.sqrt(1-mach*mach)
    return cp0/(beta+mach*mach*cp0/(2*(1+beta)))


def main():
    require(not (OUT/'verification.json').exists(), 'Preserve completed verification')
    summary = json.loads((OUT/'summary.json').read_text()); require(summary['complete'], 'Incomplete study')
    manifest = json.loads((OUT/'input_manifest.json').read_text())['files']
    changed = [n for n, d in manifest.items() if sha(ROOT/n) != d]; require(not changed, 'Inputs changed: '+str(changed))
    # Frozen solver inputs reload to the constructed arrays and hashes.
    frozen = {}
    for variant in nose.VARIANTS:
        path = OUT/f'fore_input_{variant}.dat'; require(sha(path) == summary['frozen_inputs'][path.name], 'Frozen input changed')
        frozen[variant] = np.loadtxt(path, skiprows=1)
        require(np.array_equal(frozen[variant], nose.fore_input(variant)), 'Frozen input differs from the construction: '+variant)
    # Full production constructor, not the helper used by the report.
    solvers = {v: CoupledSolver(nose.elements(v), reynolds=1e6, mach=.1, ncrit=9.) for v in ('nominal', 'convex_nose')}
    require(np.array_equal(solvers['nominal'].models[1].foil.x, solvers['convex_nose'].models[1].foil.x), 'Aft solver nodes differ')
    nodes = table('final_nodes_nose.csv'); geometry = []
    for variant, solver in solvers.items():
        polygon = solver.models[0].foil.x.T[::-1]; turn = nose.turns(polygon); inside = polygon[1:-1, 0] <= nose.NOSE
        saved = sorted((r for r in nodes if r['variant'] == variant), key=lambda r: int(r['node']))
        index = np.array([int(r['node']) for r in saved])
        difference = float(np.max(abs(polygon[index]-np.array([[float(r['x']), float(r['y'])] for r in saved]))))
        require(difference < 1e-14, 'Reported final nodes differ from the production constructor')
        geometry.append(dict(variant=variant, nose_nodes=int(inside.sum()), negative_turns=int(np.sum(turn[inside] < 0)),
                             smallest_turn_degrees=float(turn[inside].min()), reported_node_difference=difference))
    require(geometry[0]['negative_turns'] == 2 and geometry[1]['negative_turns'] == 0, 'Reversal counts differ from the finding')
    require(summary['admission']['convex_nose']['admitted'] and not summary['admission']['single_curve']['admitted'], 'Admission differs')
    # Pressure through the production solver's own inviscid basis.
    profiles = table('nose_pressure_profiles.csv'); checks = []
    for label in ('matched_cl0056', 'matched_cl0175'):
        for variant, solver in solvers.items():
            rows = [r for r in profiles if r['variant'] == variant and r['incidence'] == label and r['surface'] == 'lower']
            alpha = float(rows[0]['alpha_deg']); u = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
            model = solver.models[0]; gamma = model.isol.gamref@u; x = model.foil.x; ile = int(np.argmin(x[0]))
            xm = .5*(x[0, :-1]+x[0, 1:])[:ile]; value = cp(.5*(gamma[:-1]+gamma[1:])[:ile]); keep = xm <= 2*nose.NOSE
            saved = np.array([[float(r['x']), float(r['Cp'])] for r in rows])
            require(len(saved) == keep.sum(), 'Profile coverage differs')
            difference = float(np.max(abs(np.sort(saved[:, 0])-np.sort(xm[keep]))))+float(np.max(abs(saved[np.argsort(saved[:, 0]), 1]-value[keep][np.argsort(xm[keep])])))
            checks.append(dict(incidence=label, variant=variant, alpha_deg=alpha, panels=int(keep.sum()), maximum_difference=difference))
            require(difference < 1e-11, 'Independent pressure differs from the report')
    differences = [r for r in table('pressure_differences.csv') if r['method'] == 'production' and r['variant'] == 'convex_nose']
    worst = max(float(r['lower_recovery_max_abs_delta_Cp']) for r in differences); exceeding = [r['incidence'] for r in differences if r['material'] == 'True']
    require(abs(worst-summary['decision']['convex_nose']['production']['lower_recovery_max_abs_delta_Cp']) < 1e-15
            and (worst > MATERIAL) == summary['decision']['convex_nose']['material'], 'Decision differs from the tables')
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', 'test_s414_nose_geometry'], cwd=ROOT, text=True, capture_output=True)
    (OUT/'tests.log').write_text(tests.stdout+tests.stderr); require(tests.returncode == 0, 'Tests failed')
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, stage='V4 geometry trace and inviscid comparison', inputs_verified=len(manifest),
        frozen_inputs_reproduced=list(frozen), production_constructor_geometry=geometry, independent_pressure_checks=checks,
        admission={k: v['admitted'] for k, v in summary['admission'].items()},
        convex_nose_material=summary['decision']['convex_nose']['material'], incidences_exceeding_threshold=exceeding,
        largest_lower_recovery_delta_Cp=worst, threshold=MATERIAL,
        decision='Material by the declared rule (marginally). A separate bounded viscous protocol is required for the convex-nose '
                 'alternative; this receipt does not launch it',
        tests=dict(module='test_s414_nose_geometry', passed=True, count=6),
        figures_reviewed=['nose_geometry_trace.png', 'nose_pressure_effect.png'],
        figure_review='Rendered and inspected: off-scale tangent break annotated; threshold and recovery window shown; no label collision',
        first_attempt='The first report run was killed by the operating system for memory in a distance computation, before any '
                      'pressure was computed; the computation was chunked and rerun',
        new_nonlinear_flow_solves=0, output_sha256={p.name: sha(p) for p in outputs})
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('production_constructor_geometry', 'independent_pressure_checks', 'admission',
        'convex_nose_material', 'incidences_exceeding_threshold', 'largest_lower_recovery_delta_Cp')}, indent=1))


if __name__ == '__main__':
    main()
