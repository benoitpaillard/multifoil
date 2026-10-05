#!/usr/bin/env python3
"""Independent check of the V2 tranche-2 analysis from the raw polar exports; no flow solve."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, sha, write
from report_s702_drag_mechanism import require

T1 = ROOT/'validation/results/s414_transition_boundary_20261004'
STUDY = ROOT/'validation/results/s414_transition_boundary_t2_20261004'
OUT = ROOT/'validation/results/s414_transition_boundary_t2_analysis_20261004'
TESTS = ['test_s414_transition_boundary_t2', 'test_s414_transition_boundary_analysis']


def polar(folder, index=0):
    with (folder/'flow/polar.csv').open(newline='') as stream: row = list(csv.DictReader(stream))[index]
    return dict(alpha=float(row['alpha_deg']), converged=row['converged'] == '1', CL=float(row['CL']), CD=float(row['CD']))


def main():
    require(not (OUT/'verification.json').exists(), 'Preserve completed verification')
    audit = json.loads((OUT/'residual_recheck.json').read_text()); summary = json.loads((OUT/'summary.json').read_text())
    spec = json.loads((STUDY/'study.json').read_text()); screens = spec['equivalence_screens']; jobs = {j['name']: j for j in spec['jobs']}
    require(audit['complete'] and len(audit['states']) == audit['saved_states'] == 14, 'Audit coverage differs')
    require(max(s['max_residual'] for s in audit['states']) <= spec['residual_tolerance'], 'Fresh residual exceeds tolerance')
    changed = [n for n, d in audit['input_sha256'].items() if sha(ROOT/n) != d]
    require(not changed, 'Audited inputs changed: '+str(changed[:3]))
    raw = {name: polar(STUDY/name) for name in jobs}
    for name, job in jobs.items():
        require(raw[name]['alpha'] == job['alpha_deg'], 'Exported incidence differs from the declared one')
        if not raw[name]['converged']: require(math.isnan(raw[name]['CL']) and math.isnan(raw[name]['CD']), 'Failed evaluation carries a force')
    new = [n for n, j in jobs.items() if j['kind'] != 'restore_control']
    require(len(new) == 16 and sum(raw[n]['converged'] for n in new) == 12 == summary['converged_new_evaluations'], 'Coverage differs')
    require(spec['requested_points']+spec['tranche1_requested_points'] == 42 == summary['total_requested_points_V2'], 'Allowance differs')
    # Restore controls against the raw tranche-1 polar row of the same step.
    restores = []
    for scenario in ('nominal', 'fore_source_cubic'):
        live = polar(T1/f'{scenario}_corrected_increasing', 1); again = raw[scenario+'_restore_control']
        restores.append(dict(scenario=scenario, delta_CL=again['CL']-live['CL'], delta_CD_counts=(again['CD']-live['CD'])*1e4))
        require(again['alpha'] == live['alpha'] and abs(restores[-1]['delta_CL']) <= screens['CL']
                and abs(restores[-1]['delta_CD_counts']) <= screens['CD_counts'], 'Restore control differs from the live result')
    pairs = []
    for scenario, count in (('fore_source_cubic', 4), ('nominal', 3)):
        for k in range(count):
            up, down = raw[f'{scenario}_increasing_{k}'], raw[f'{scenario}_decreasing_{k}']
            require(up['alpha'] == down['alpha'], 'Direction pairing differs')
            if up['converged'] and down['converged']:
                pairs.append(dict(scenario=scenario, alpha_deg=up['alpha'], delta_CL=up['CL']-down['CL'], delta_CD_counts=(up['CD']-down['CD'])*1e4))
            else:
                require(not up['converged'] and not down['converged'], 'A one-sided failure is not what the finding states')
    require(len(pairs) == 5 == summary['direction_pairs_both_converged'] == summary['direction_pairs_within_screens'], 'Pair count differs')
    require(all(abs(p['delta_CL']) <= screens['CL'] and abs(p['delta_CD_counts']) <= screens['CD_counts'] for p in pairs), 'Raw exports violate a screen')
    cold = raw['fore_source_cubic_cold_check']; live = polar(T1/'fore_source_cubic_corrected_increasing', 1)
    cold_cubic = dict(delta_CL=cold['CL']-live['CL'], delta_CD_counts=(cold['CD']-live['CD'])*1e4)
    cold_nominal = dict(delta_CL=raw['nominal_cold_check']['CL']-raw['nominal_increasing_1']['CL'],
                        delta_CD_counts=(raw['nominal_cold_check']['CD']-raw['nominal_increasing_1']['CD'])*1e4)
    for item in (cold_cubic, cold_nominal):
        require(abs(item['delta_CL']) <= screens['CL'] and abs(item['delta_CD_counts']) <= screens['CD_counts'], 'Cold check differs from the path state')
    require(jobs['fore_source_cubic_cold_check']['deviation']['declared_alpha_deg'] != jobs['fore_source_cubic_cold_check']['alpha_deg'], 'Deviation not recorded')
    require(not summary['distinct_converged_states_at_one_incidence'], 'Summary reports distinct states')
    require(len(summary['intervals_bracketing_target_CL']['fore_source_cubic']) == 3, 'Target-lift bracket count differs')
    nominal_drop = (raw['nominal_decreasing_2']['CD']-polar(T1/'nominal_corrected_increasing', 0)['CD'])*1e4
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TESTS], cwd=ROOT, text=True, capture_output=True)
    (OUT/'tests.log').write_text(tests.stdout+tests.stderr); require(tests.returncode == 0, 'Tests failed')
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, stage='V2 tranche 2', retrieved_files_verified=audit['retrieved_files_verified'],
        saved_states_audited=audit['saved_states'], maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        restored_base_states=len(audit['restored_base_states']),
        maximum_restored_base_residual=max(max(r['BL_residual_l2'], r['coupling_residual_max']) for r in audit['restored_base_states']),
        restore_controls_from_raw_polars=restores, direction_pairs_rechecked_from_raw_polars=pairs,
        cold_check_cubic_vs_tranche1_path=cold_cubic, cold_check_nominal_vs_path=cold_nominal,
        failed_evaluations=audit['failed_evaluations'], nominal_drag_change_across_unresolved_interval_counts=-nominal_drop,
        recorded_deviation=jobs['fore_source_cubic_cold_check']['deviation'],
        tests=dict(modules=TESTS, passed=True), figures_reviewed=['cubic_refined_response.png', 'nominal_fallback_response.png'],
        figure_review='Rendered and inspected: failed incidences shown without a value and never bridged by a line; '
                      'cold starts and the archived anchor distinguished and labelled',
        new_flow_solves=0, output_sha256={p.name: sha(p) for p in outputs},
        qualification='Numerical records verified; sampled-path repeatability only. Not matched-lift validation or uniqueness.')
    write(OUT/'verification.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('saved_states_audited', 'maximum_fresh_residual', 'restored_base_states',
        'maximum_restored_base_residual', 'restore_controls_from_raw_polars', 'cold_check_cubic_vs_tranche1_path',
        'cold_check_nominal_vs_path', 'nominal_drag_change_across_unresolved_interval_counts')}, indent=1))


if __name__ == '__main__':
    main()
