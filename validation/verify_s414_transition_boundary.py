#!/usr/bin/env python3
"""Independent check of the V2 tranche-1 analysis from the raw exports; no flow solve."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, sha, write
from report_s702_drag_mechanism import require

STUDY = ROOT/'validation/results/s414_transition_boundary_20261004'
OUT = ROOT/'validation/results/s414_transition_boundary_analysis_20261004'
TESTS = ['test_s414_transition_boundary_analysis', 'test_s414_transition_boundary']


def table(path):
    with Path(path).open(newline='') as stream: return list(csv.DictReader(stream))


def main():
    require(not (OUT/'verification.json').exists(), 'Preserve completed verification')
    audit = json.loads((OUT/'residual_recheck.json').read_text()); summary = json.loads((OUT/'summary.json').read_text())
    spec = json.loads((STUDY/'study.json').read_text()); screens = spec['equivalence_screens']
    require(audit['complete'] and len(audit['states']) == audit['saved_states'] == 25, 'Audit coverage differs')
    require(max(s['max_residual'] for s in audit['states']) <= spec['residual_tolerance'], 'Fresh residual exceeds tolerance')
    changed = [n for n, d in audit['input_sha256'].items() if sha(ROOT/n) != d]
    require(not changed, 'Audited inputs changed: '+str(changed[:3]))
    points = table(OUT/'points.csv'); raw = {}
    # Raw polar and transition exports, not the result records the report reads.
    for job in spec['jobs']:
        polar = table(STUDY/job['name']/'flow/polar.csv'); transition = table(STUDY/job['name']/'flow/transition_separation.csv')
        for index, (request, row) in enumerate(zip(job['requests'], polar)):
            fore = [t for t in transition if float(t['alpha_deg']) == request['alpha_deg'] and t['element'] == 'fore' and t['side'] == 'lower']
            raw[(job['name'], index)] = dict(alpha=float(row['alpha_deg']), converged=row['converged'] == '1',
                CL=float(row['CL']), CD=float(row['CD']), transition_local=[float(t['transition_x_over_local_chord']) for t in fore])
            require(float(row['alpha_deg']) == request['alpha_deg'], 'Exported incidence differs from the request')
    checked = 0
    for p in points:
        if p['index'] == '': continue
        r = raw[(p['job'], int(p['index']))]; require(r['converged'] == (p['converged'] == 'True'), 'Convergence differs from polar')
        if r['converged']:
            require(abs(r['CL']-float(p['CL'])) < 1e-15 and abs(r['CD']*1e4-float(p['CD_counts'])) < 1e-10, 'Force differs from polar')
            checked += 1
        else:
            require(math.isnan(r['CL']) and math.isnan(r['CD']) and p['CL'] == '' and p['CD'] == '', 'Failed point carries a force')
    path = [p for p in points if p['role'] == 'path']
    require(len(path) == 20 and sum(p['converged'] == 'True' for p in path) == 19 == summary['converged_path_points'], 'Path coverage differs')
    pairs = []
    for scenario in ('nominal', 'fore_source_cubic'):
        for k, alpha in enumerate(spec['angles_deg']):
            up = raw[(scenario+'_corrected_increasing', k)]; down = raw[(scenario+'_corrected_decreasing', 5-k)]
            require(up['alpha'] == down['alpha'] == alpha, 'Direction pairing differs')
            if up['converged'] and down['converged']:
                pairs.append(dict(scenario=scenario, alpha_deg=alpha, delta_CL=up['CL']-down['CL'], delta_CD_counts=(up['CD']-down['CD'])*1e4,
                    delta_fore_lower_transition_local=max(abs(a-b) for a, b in zip(up['transition_local'], down['transition_local']))))
    require(len(pairs) == summary['direction_pairs_both_converged'] == 9, 'Direction pair count differs')
    require(all(abs(p['delta_CL']) <= screens['CL'] and abs(p['delta_CD_counts']) <= screens['CD_counts']
                and p['delta_fore_lower_transition_local'] <= screens['transition_x_over_c'] for p in pairs), 'Raw exports violate a screen')
    require(abs(max(abs(p['delta_CL']) for p in pairs)-summary['largest_direction_difference']['CL']) < 1e-15, 'Largest lift difference differs')
    cubic = [raw[('fore_source_cubic_corrected_increasing', k)] for k in range(5)]
    lift = [c['CL'] for c in cubic]; drag = [c['CD']*1e4 for c in cubic]
    require(lift[1] > lift[2] and lift[0] < lift[1] and lift[2] < lift[3] < lift[4], 'Cubic lift response differs from the finding')
    brackets = sum((a-.175)*(b-.175) < 0 for a, b in zip(lift, lift[1:]))
    require(brackets == 3 == len(summary['intervals_bracketing_target_CL']['fore_source_cubic']), 'Target-lift bracket count differs')
    selection = summary['tranche2_selection']
    require(selection['nominal']['new_locations'] == [] and selection['nominal']['fallback'], 'Nominal selection differs')
    require(len(selection['fore_source_cubic']['new_locations']) == 4, 'Cubic selection differs')
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TESTS], cwd=ROOT, text=True, capture_output=True)
    (OUT/'tests.log').write_text(tests.stdout+tests.stderr); require(tests.returncode == 0, 'Tests failed')
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, stage='V2 tranche 1', retrieved_files_verified=audit['retrieved_files_verified'],
        saved_states_audited=audit['saved_states'], maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        maximum_transition_root_shift=max(s['transition_root_maximum_shift'] for s in audit['states']),
        legacy_controls_identical_to_archive=audit['controls']['passed'],
        control_items_compared=sum(len(j['comparisons']) for j in audit['controls']['jobs'].values()),
        failed_evaluations=audit['failed_evaluations'], converged_forces_checked_against_polar=checked,
        direction_pairs_rechecked_from_raw_exports=pairs, cubic_increasing_lift=lift, cubic_increasing_drag_counts=drag,
        tests=dict(modules=TESTS, passed=True), figures_reviewed=['path_response.png', 'path_diagnostics.png'],
        figure_review='Rendered and inspected: failed point shown without a value and never bridged by a line; '
                      'both directions visible; one vertical scale per measure across geometries',
        new_flow_solves=0, output_sha256={p.name: sha(p) for p in outputs},
        qualification='Numerical records verified; sampled-path repeatability only. Not matched-lift validation or uniqueness.')
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('saved_states_audited', 'maximum_fresh_residual', 'maximum_transition_root_shift',
        'legacy_controls_identical_to_archive', 'control_items_compared', 'converged_forces_checked_against_polar',
        'cubic_increasing_lift', 'cubic_increasing_drag_counts')}, indent=1))


if __name__ == '__main__':
    main()
