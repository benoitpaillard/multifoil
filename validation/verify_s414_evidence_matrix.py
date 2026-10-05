#!/usr/bin/env python3
"""Independent check of the S414 evidence matrix (stage V0); no flow solves.

Recomputes every hash of input_manifest.json, checks that each cited path exists
and still has the hash cited beside it, re-reads every `data` value and quotation
from its archived file, and checks that no row lacks an admitted use, a stated
uncertainty entry or an explicit unknown entry. Writes verification.json.
"""
import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_evidence_matrix_20261004'
CASES = ('CL -0.072', 'CL 0.056', 'CL 0.175', 'CL 0.284', 'CL 1.041')
KINDS = ('experiment_pressure', 'experiment_drag', 'mses2010_drag', 'mses2010_pressure', 'mses_transition',
         'experiment_transition_reattachment_location', 'multifoil_archived_audited_state')
ENVELOPE_2010 = ('mses_2010_envelope.csv', 'published_mses_drag.csv', 'mses_drag.csv', 'mses_lift.csv')
ARTICLE_2014 = 'reference_review_20260924/sources/s414_article'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(text):
    return re.sub(r'\s+', ' ', text.replace('*', '')).strip()


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def reread(record, cache):
    """Value of a `data` record read again from its archived file, by its locator."""
    path = ROOT/record['source']
    if record['source'].endswith('.json'):
        node = cache.setdefault(path, json.loads(path.read_text()))
        for key in record['locator'].split('/'):
            node = node[int(key)] if isinstance(node, list) else node[key]
        return node
    where, field = re.fullmatch(r'row\[(.*)\]\.(\w+)', record['locator']).groups()
    filters = dict(item.split('=', 1) for item in where.split(','))
    table = cache.setdefault(path, rows(path))
    hits = [r for r in table if all(r[k] == v for k, v in filters.items())]
    if len(hits) != 1:
        raise ValueError(f'{record["source"]}: {len(hits)} rows for {filters}')
    return None if hits[0][field] == '' else hits[0][field]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    failures, cache = [], {}

    def check(ok, message):
        if not ok:
            failures.append(message)
        return ok

    manifest = json.loads((out/'input_manifest.json').read_text())
    summary = json.loads((out/'summary.json').read_text())
    matrix = json.loads((out/'evidence_matrix.json').read_text())
    flat = rows(out/'evidence_matrix.csv')
    coverage = rows(out/'case_coverage.csv')
    classes = rows(out/'uncertainty_classes.csv')
    unknowns = rows(out/'unknown_quantities.csv')
    current = {}

    def digest(rel):
        if rel not in current:
            path = ROOT/rel
            current[rel] = sha(path) if path.is_file() else None
        return current[rel]

    # 1. every hash of the input manifest
    changed = {}
    for category in ('files', 'living_documents', 'reconciliation_only_files'):
        changed[category] = sorted(rel for rel, recorded in manifest[category].items() if digest(rel) != recorded)
    check(not changed['files'], 'Evidence files changed or missing: '+', '.join(changed['files']))
    check(summary['complete'] and summary['new_flow_solves'] == 0, 'Summary incomplete or reports a flow solve')

    # 2. every cited path exists, is in the manifest and still has the cited hash
    known = dict(manifest['files'], **manifest['living_documents'])
    cited = 0
    tables = [(r['row_id'], r['source_paths'], r['source_sha256']) for r in matrix['rows']]
    tables += [(c['case']+' '+c['evidence'], c['source_paths'].split('; '), c['source_sha256'].split('; ')) for c in coverage]
    tables += [(c['item'], [c['source_path']], [c['source_sha256']]) for c in classes]
    tables += [(u['id'], u['source_paths'].split('; '), u['source_sha256'].split('; ')) for u in unknowns]
    for label, paths, hashes in tables:
        check(paths and len(paths) == len(hashes) and all(paths), f'{label}: missing source path or hash')
        for rel, cited_hash in zip(paths, hashes):
            cited += 1
            check((ROOT/rel).is_file(), f'{label}: cited path missing: {rel}')
            check(known.get(rel) == cited_hash, f'{label}: cited hash not in input manifest: {rel}')
            if rel in manifest['files']:
                check(digest(rel) == cited_hash, f'{label}: cited hash differs from the current file: {rel}')

    # 3. no row without admitted use, uncertainty entry or explicit unknown entry
    statuses = set(matrix['statuses'])
    check([r['row_id'] for r in matrix['rows']] == [r['row_id'] for r in flat], 'CSV and JSON matrices list different rows')
    check(len({r['row_id'] for r in flat}) == len(flat) == summary['rows'], 'Row identifiers not unique or row count differs')
    for r in flat:
        for field in ('admitted_use', 'unknown_quantities', 'reported_uncertainty', 'conditions', 'measurement_resolution',
                      'source_page_figure', 'strict_matched_geometry_validation', 'source_paths', 'source_sha256'):
            check(r[field].strip(), f'{r["row_id"]}: empty {field}')
        check(r['admission_status'] in statuses, f'{r["row_id"]}: unknown status')
        if r['admission_status'] in ('missing', 'not_extracted'):
            check('unknown U' in r['unknown_quantities'], f'{r["row_id"]}: missing evidence without a registered unknown')
    registry = matrix['unknown_registry']
    for r in matrix['rows']:
        for item in r['unknown_quantities']:
            check(item['id'] in registry and r['row_id'] in registry[item['id']]['rows'], f'{r["row_id"]}: unregistered unknown {item["id"]}')
    for status in ('admitted', 'conditional', 'not_admitted', 'not_extracted', 'missing'):
        check(summary[status] == sum(r['admission_status'] == status for r in flat), 'Summary count differs: '+status)
    check(summary['unknown_items'] == len(registry) == len(unknowns), 'Unknown registry count differs')

    # 4. re-read every data value and quotation from its archived file
    reread_values = quotations = 0
    records = [(r['row_id'], e) for r in matrix['rows'] for e in r['evidence']]
    records += [(key, e) for key, item in registry.items() for e in item['basis']]
    for label, e in records:
        try:
            if e['kind'] == 'data':
                reread_values += 1
                check(reread(e, cache) == e['value'], f'{label}: value differs from {e["source"]}#{e["locator"]}')
            elif e['kind'] == 'quote':
                quotations += 1
                text = cache.setdefault(('text', e['source']), norm((ROOT/e['source']).read_text(encoding='utf-8', errors='replace')))
                check(e['value'] in text, f'{label}: quotation no longer in {e["source"]}')
            elif e['kind'] == 'source_text':
                quotations += 1
                lines = cache.setdefault(('lines', e['source']), (ROOT/e['source']).read_text(encoding='utf-8', errors='replace').split('\n'))
                check(e['value'] in lines[int(e['locator'].split()[1])-1], f'{label}: source text not on {e["locator"]} of {e["source"]}')
        except (KeyError, IndexError, ValueError, AttributeError, OSError) as error:
            check(False, f'{label}: cannot re-read {e["source"]}#{e["locator"]}: {error}')

    # 5. series separation and coverage completeness
    for r in matrix['rows']:
        sources = ' '.join(r['source_paths'])
        if r['reference'].startswith(('mses_2014', 'experiment_2014', 'overflow_2014')):
            check(not any(name in sources for name in ENVELOPE_2010), f'{r["row_id"]}: 2014 row cites 2010 envelope-MSES data')
        if r['reference'].startswith('mses_2010'):
            check(ARTICLE_2014 not in sources, f'{r["row_id"]}: 2010 envelope-MSES row cites the 2014 article')
        if r['reference'].startswith('experiment_2010') and ' and ' not in r['reference']:
            check(not any(e['kind'] == 'data' and e['source'].endswith(ENVELOPE_2010) for e in r['evidence']),
                  f'{r["row_id"]}: experimental row carries a published-MSES value')
        if r['reference'].startswith('multifoil'):
            check('not a reference measurement' in r['reference'], f'{r["row_id"]}: model state not labelled as model')
    entries = {(c['case'], c['evidence']): c for c in coverage}
    check(len(entries) == len(coverage) == len(CASES)*len(KINDS), 'Coverage table is not 5 cases x 7 evidence kinds')
    for case in CASES:
        for kind in KINDS:
            entry = entries.get((case, kind))
            if check(entry is not None, f'No coverage entry: {case} {kind}'):
                check(entry['status'].startswith(('exists', 'missing')) and entry['matrix_row_ids'], f'Coverage status unclear: {case} {kind}')
    for c in classes:
        check(c['value'].strip() and c['status'].strip(), 'Uncertainty item without value/status: '+c['item'])
    for name in ('extraction repeatability', 'reported manufacturing tolerance', 'reconstruction spread', 'actual as-built uncertainty'):
        check(any(c['uncertainty_class'] == name for c in classes), 'Missing uncertainty class: '+name)

    outputs = sorted(p.name for p in out.iterdir() if p.is_file() and p.name != 'verification.json')
    receipt = dict(
        verified_utc=datetime.now(timezone.utc).isoformat(), passed=not failures, failures=failures,
        evidence_files_checked=len(manifest['files']), evidence_files_changed=changed['files'],
        living_documents_changed=changed['living_documents'],
        reconciliation_only_files_checked=len(manifest['reconciliation_only_files']),
        reconciliation_only_files_changed=changed['reconciliation_only_files'],
        cited_paths_checked=cited, matrix_rows=len(flat), coverage_entries=len(coverage),
        uncertainty_class_items=len(classes), unknown_items=len(unknowns),
        data_values_reread=reread_values, quotations_rechecked=quotations, new_flow_solves=0,
        note='A change in a living document or a reconciliation-only file is reported, not failed: those files are '
             'expected to be edited after this matrix; an evidence-file change fails the check.',
        code_sha256={name: sha(ROOT/name) for name in ('validation/report_s414_evidence_matrix.py',
                                                       'validation/verify_s414_evidence_matrix.py')},
        report_script_unchanged=sha(ROOT/'validation/report_s414_evidence_matrix.py') == manifest['code_sha256']['validation/report_s414_evidence_matrix.py'],
        output_sha256={name: sha(out/name) for name in outputs})
    (out/'verification.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False)+'\n')
    print(f'{"PASS" if not failures else "FAIL"}: {len(manifest["files"])} evidence files, {cited} cited paths, '
          f'{reread_values} values re-read, {quotations} quotations re-checked, {len(flat)} rows')
    for failure in failures[:40]:
        print('  '+failure)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
