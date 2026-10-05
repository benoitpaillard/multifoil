#!/usr/bin/env python3
"""V1 gate: implementation receipt and source snapshot for the integrated default.

Combines the archived-state replay, the targeted tests, a classification of the
complete suite against the same suite on pre-integration sources, and the known
far-field limitation evaluated in both derivative modes. No S414 flow solve.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, sha, write, csv_write
from report_s702_drag_mechanism import require

OUT = ROOT/'validation/results/s414_default_integration_20261004'
FROZEN = ROOT/'validation/results/s414_fore_reconstruction_20261004/source_snapshot'
INTEGRATED = ('amplification_jacobian.py', 'multifoil.py', 'shape_optimization.py', 'validate.py', 'viscous.py')
TARGETED = ['test_default_integration', 'test_amplification_jacobian', 'test_viscous', 'test_multifoil',
            'test_optimization', 'test_shape_optimization', 'test_fixed_lift_policy', 'test_paneling', 'test_validation']
FAR_FIELD = 'test_viscous.ViscousPhysics.test_multiple_elements_order_and_far_limit'
GUARD = re.compile(r'ValueError: [^\n]*(changed|differs)[^\n]*: ('+'|'.join(re.escape(n) for n in INTEGRATED)+r')\s*$', re.M)
BLOCK = re.compile(r'^={70}\n(ERROR|FAIL): (\S+) \(([\w.]+)\)\n-{70}\n(.*?)(?=^={70}\n|^-{70}\nRan )', re.M | re.S)


def parse(path):
    text = Path(path).read_text()
    ran = re.search(r'^Ran (\d+) tests', text, re.M)
    require(ran is not None, 'Incomplete test log: '+str(path))
    failures = {m.group(3) if m.group(2) != 'setUpClass' else m.group(3)+'.setUpClass': (m.group(1), m.group(4))
                for m in BLOCK.finditer(text)}
    counts = re.search(r'FAILED \((?:failures=(\d+))?(?:, )?(?:errors=(\d+))?', text)
    declared = sum(int(v) for v in counts.groups() if v) if counts else 0
    require(declared == len(failures), 'Could not parse every failure block in '+str(path))
    return int(ran.group(1)), failures


def classify(test, body, hashes):
    last = [line for line in body.splitlines() if re.match(r'^[A-Za-z_.]+(Error|Exception): ', line)]
    message = last[-1] if last else ''
    if test == FAR_FIELD and 'not less than 5e-07' in message:
        return 'known_far_field_limitation', message
    if GUARD.search(body):
        return 'historical_source_guard', GUARD.search(body).group(0).strip()
    pair = re.search(r"AssertionError: '([0-9a-f]{64})' != '([0-9a-f]{64})'", body)
    if pair and any({pair.group(1), pair.group(2)} == set(v) for v in hashes.values()):
        return 'historical_source_guard', message
    return 'unexplained', message


def far_field(mode):
    """The preserved limitation, evaluated exactly as test_viscous does, in one explicit mode."""
    from multifoil import naca4
    from viscous import CoupledSolver
    xy = naca4('0012', 81); e = [('a', xy), ('b', xy+[0, 100])]
    forward = CoupledSolver(e, tolerance=1e-9, amplification_jacobian=mode).solve(2)
    reverse = CoupledSolver(e[::-1], tolerance=1e-9, amplification_jacobian=mode).solve(2)
    require(forward.converged and reverse.converged, 'Far-field check did not converge')
    return dict(mode=mode, forward_CL=forward.cl, reverse_CL=reverse.cl, CL_difference=abs(forward.cl-reverse.cl),
                CD_difference=abs(forward.cd-reverse.cd), exceeds_5e_7=bool(abs(forward.cl-reverse.cl) >= 5e-7),
                iterations=[forward.diagnostics['iterations'], reverse.diagnostics['iterations']])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-log', type=Path, required=True,
                        help='Complete-suite log from a copy with the five pre-integration files restored')
    args = parser.parse_args(); out = OUT
    require(not (out/'verification.json').exists(), 'Preserve completed receipt; choose a new output directory')
    replay = json.loads((out/'replay_summary.json').read_text())
    require(replay['complete'] and replay['passed'], 'Archived-state replay has not passed')
    hashes = {name: (sha(FROZEN/name), sha(ROOT/name)) for name in INTEGRATED}
    require(all(a != b for a, b in hashes.values()), 'Integrated file equals its frozen version')
    baseline = out/'full_suite_pre_integration_source.log'
    if not baseline.exists(): shutil.copyfile(args.baseline_log, baseline)
    require(sha(baseline) == sha(args.baseline_log), 'Baseline log differs from the preserved copy')
    ran_before, failed_before = parse(baseline); ran_now, failed_now = parse(out/'full_suite_current_source.log')
    rows = []
    for test, (kind, body) in sorted(failed_now.items()):
        category, message = classify(test, body, hashes)
        rows.append(dict(test=test, outcome=kind, failed_before_integration=test in failed_before,
                         classification=category, message=message))
    csv_write(out/'test_classification.csv', rows)
    before = {test: classify(test, body, hashes) for test, (_, body) in failed_before.items()}
    require(set(before) == {FAR_FIELD} and before[FAR_FIELD][0] == 'known_far_field_limitation',
            'Pre-integration baseline has failures other than the known far-field limitation')
    unexplained = [r['test'] for r in rows if r['classification'] == 'unexplained']
    tests = subprocess.run([sys.executable, '-m', 'unittest', '-v', *TARGETED], cwd=ROOT, text=True, capture_output=True)
    (out/'targeted_tests.log').write_text(tests.stdout+tests.stderr)
    ran_targeted, failed_targeted = parse(out/'targeted_tests.log')
    targeted_ok = set(failed_targeted) == {FAR_FIELD}
    value = lambda body: float(re.search(r'AssertionError: (\S+) not less than 5e-07', body).group(1))
    limitation = dict(test=FAR_FIELD, threshold=5e-7,
        pre_integration_value=value(failed_before[FAR_FIELD][1]), current_default_value=value(failed_now[FAR_FIELD][1]),
        explicit_modes=[far_field('legacy'), far_field('corrected_af')])
    legacy, corrected = limitation['explicit_modes']
    limitation.update(
        legacy_reproduces_pre_integration_value_exactly=legacy['CL_difference'] == limitation['pre_integration_value'],
        corrected_reproduces_current_default_value_exactly=corrected['CL_difference'] == limitation['current_default_value'],
        change_between_modes=abs(corrected['CL_difference']-legacy['CL_difference']),
        interpretation='The limitation is preserved. The default mode changes the Newton path, moving this converged '
                       'lift difference by far less than the solver tolerance; explicit legacy restores the archived value.')
    limitation_ok = bool(legacy['exceeds_5e_7'] and corrected['exceeds_5e_7']
        and limitation['legacy_reproduces_pre_integration_value_exactly']
        and limitation['corrected_reproduces_current_default_value_exactly'] and limitation['change_between_modes'] < 1e-12)
    snapshot = out/'source_snapshot'
    names = sorted([p.name for p in ROOT.glob('*.py')]+['vendor/mfoil.py', 'vendor/__init__.py', 'requirements.txt',
        'requirements-handover.txt', 'validation/verify_s414_default_integration.py',
        'validation/s414_default_integration_frozen_child.py', 'validation/receipt_s414_default_integration.py'])
    manifest = {}
    for name in names:
        saved = snapshot/name; saved.parent.mkdir(parents=True, exist_ok=True)
        if not saved.exists(): shutil.copyfile(ROOT/name, saved)
        require(sha(saved) == sha(ROOT/name), 'Source changed while writing the snapshot: '+name)
        manifest[name] = sha(saved)
    write(out/'source_manifest.json', dict(files=manifest, integrated={n: dict(frozen_pre_integration=a, current=b)
          for n, (a, b) in hashes.items()}))
    gate = bool(replay['passed'] and not unexplained and targeted_ok and limitation_ok)
    outputs = sorted(p for p in out.iterdir() if p.is_file() and p.name not in ('verification.json', 'FINDINGS.md'))
    receipt = dict(complete=True, gate_passed=gate, stage='V1: integrated default amplification Jacobian',
        replay=dict(passed=replay['passed'], state_replay_comparisons=replay['state_replay_comparisons'],
            state_replay_exact=replay['state_replay_exact'], derivative_comparisons=replay['derivative_comparisons'],
            archived_rule_comparisons=replay['archived_rule_comparisons'],
            archived_rule_maximum_corrected_relative_error=replay['archived_rule_maximum_corrected_relative_error'],
            declared_step_comparisons_at_or_above_criterion=replay['declared_step_comparisons_at_or_above_criterion'],
            exceedances_classified_as_truncation=sum(i['classified_as_finite_difference_truncation'] for i in replay['exceedance_investigations']),
            maximum_corrected_relative_error_at_1e_5=replay['maximum_corrected_relative_error_at_1e_5'],
            maximum_legacy_relative_error=replay['maximum_legacy_relative_error'],
            first_attempt='initial_replay_failed/: both steps of all 80 probes were held to 2e-5; two coarse-step rows exceeded it'),
        targeted_tests=dict(modules=TARGETED, ran=ran_targeted, failed=sorted(failed_targeted), only_known_limitation_fails=targeted_ok),
        complete_suite=dict(pre_integration=dict(ran=ran_before, failed=sorted(failed_before)),
            current=dict(ran=ran_now, failed=len(failed_now),
                by_classification={c: sum(r['classification'] == c for r in rows) for c in sorted({r['classification'] for r in rows})}),
            unexplained=unexplained,
            note='Guard failures are historical study manifests refusing the changed source, as designed. '
                 'They are not relaxed; archived studies replay from their frozen snapshots. A class-level '
                 'guard error replaces several tests by one entry, hence the smaller number of tests run.'),
        far_field_limitation=limitation, source_snapshot_files=len(manifest),
        new_nonlinear_s414_flow_solves=0,
        runtime=dict(python=platform.python_version(), versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')},
            note='Local analysis runtime; exact archived-array replay under the pinned 3.12.13 runtime is the V2 control gate.'),
        output_sha256={p.name: sha(p) for p in outputs},
        qualification='Implementation verification only. No aerodynamic benefit, physical validation or solution uniqueness is claimed.')
    write(out/'verification.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('gate_passed', 'targeted_tests', 'complete_suite')}, indent=2))
    require(gate, 'V1 gate failed; stop new nonlinear interpretation')


if __name__ == '__main__':
    main()
