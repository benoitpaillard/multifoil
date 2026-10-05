#!/usr/bin/env python3
"""Stage V0: consolidate the existing S414 reference manifests into one evidence matrix.

Consolidation only: no flow solve, download, literature search or PDF extraction.
Every value is (a) read from an archived CSV/JSON file, (b) a `derived` count,
extremum or difference of archived columns, (c) a literal quotation whose presence
in the cited FINDINGS/protocol/register file is checked here, or (d) a line of an
archived source-text rendition located by a literal fragment. A quantity recorded
nowhere is written as `unknown`; nothing is estimated.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/results/s414_evidence_matrix_20261004'
RES = 'validation/results/'
GEO = RES+'s414_geometry_qualification_20261002/'
DRAG = RES+'s414_drag_reference_20261002/'
TRIP = RES+'s414_trip_qualification_20261003/'
LBA = RES+'s414_lower_bucket_analysis_20261003/'
RCA = RES+'s414_reconstruction_analysis_20261003/'
JCA = RES+'s414_jacobian_analysis_20261003/'
INA = RES+'s414_incidence_analysis_20261003/'
FRA = RES+'s414_fore_reconstruction_analysis_20261004/'
NOSE = RES+'s414_fore_nose_20261003/'
PB = RES+'s414_pressure_budget_20261003/'
WAKE = RES+'s414_wake_transfer_20261003/'
REVIEW = RES+'reference_review_20260924/'
PLAN = 'DEVELOPMENT_PLAN.md'
REGISTER = 'validation/REFERENCES.md'
LB_PROTOCOL = 'validation/S414_LOWER_BUCKET_PROTOCOL.md'
RC_PROTOCOL = 'validation/S414_RECONSTRUCTION_PROTOCOL.md'
REPORT_TXT = GEO+'sources/s414_report_2010.txt'
REPORT_PDF = GEO+'sources/s414_report_2010.pdf'
ARTICLE_TXT = REVIEW+'sources/s414_article.txt'
ARTICLE_PDF = REVIEW+'sources/s414_article.pdf'
# Documents that are edited as the project advances; hashes recorded, changes expected.
LIVING = {PLAN, REGISTER}
# Archived solver studies: (analysis directory, state directory).
STUDIES = [(RCA, RES+'s414_reconstruction_20261002/'), (LBA, RES+'s414_lower_bucket_20261003/'),
           (JCA, RES+'s414_jacobian_20261003/'), (INA, RES+'s414_incidence_20261003/'),
           (FRA, RES+'s414_fore_reconstruction_20261004/')]
RECEIPTS = ('verification.json', 'input_manifest.json', 'source.json', 'pressure/source.json',
            'source_manifest.json')
RECONCILED = [GEO, DRAG, TRIP, LBA, RCA, JCA, INA, FRA, NOSE, PB, WAKE,
              RES+'s414_inviscid_budget_20261003/', RES+'s414_prescribed_pressure_20261004/',
              RES+'s414_jacobian_check_20261003/']
GROUPS = {
    'A': '2010 free-transition experiment, CL -0.072/0.056/0.175/0.284',
    'B': 'Published 2010 envelope-MSES',
    'C': 'CL 1.041 pressure/drag and other completed S414 states (regression controls)',
    'D': '2014 full-e^N MSES and related measurements',
    'E': '2010 tripped and secondary-Re series (supporting context; not a row of the V0 plan table)',
}
EXP = 'experiment_2010 (Somers and Maughmer, RDECOM TR 10-D-112)'
MSES10 = 'mses_2010_envelope (published curves; no MSES executable run)'
MSES14 = 'mses_2014_full_eN (Coder, Maughmer and Somers, J. Aircraft 51(6); not digitized)'
EXP14 = 'experiment_2014_article (overlaps the 2010 campaign; not independent)'
MODEL = 'multifoil_archived_state (model result; not a reference measurement)'
STRICT = 'conditional (geometry not qualified)'
NA = 'not applicable'
STATUSES = ('admitted', 'conditional', 'not_admitted', 'not_extracted', 'missing')
HEX = re.compile(r'^[0-9a-f]{64}$')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(text):
    return re.sub(r'\s+', ' ', text.replace('*', '')).strip()


def show(value):
    return value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)


class Evidence:
    """Registers every file read and builds traceable value/quotation records."""

    def __init__(self):
        self.files, self.living, self.reconciled = {}, {}, {}
        self._json, self._rows, self._text = {}, {}, {}

    def path(self, rel):
        target = self.living if rel in LIVING else self.files
        if rel not in target:
            path = ROOT/rel
            require(path.is_file(), 'Missing evidence file: '+rel)
            target[rel] = sha(path)
        return ROOT/rel

    def digest(self, rel):
        self.path(rel)
        return self.living[rel] if rel in LIVING else self.files[rel]

    def json(self, rel):
        if rel not in self._json:
            self._json[rel] = json.loads(self.path(rel).read_text())
        return self._json[rel]

    def rows(self, rel):
        if rel not in self._rows:
            with self.path(rel).open(newline='') as stream:
                self._rows[rel] = list(csv.DictReader(stream))
        return self._rows[rel]

    def text(self, rel):
        if rel not in self._text:
            self._text[rel] = self.path(rel).read_text(encoding='utf-8', errors='replace')
        return self._text[rel]

    # --- traceable records -------------------------------------------------
    def val(self, name, rel, *keys, unit=''):
        value = self.json(rel)
        for key in keys:
            value = value[key]
        return dict(kind='data', name=name, value=value, unit=unit, source=rel,
                    locator='/'.join(str(k) for k in keys))

    def find(self, rel, keys, **where):
        """Index of the single list item matching `where` below JSON `keys`."""
        items = self.json(rel)
        for key in keys:
            items = items[key]
        hits = [i for i, item in enumerate(items) if all(item[k] == v for k, v in where.items())]
        require(len(hits) == 1, f'Expected one JSON item in {rel}: {where}')
        return hits[0]

    def row(self, rel, **where):
        hits = [r for r in self.rows(rel) if all(r[k] == str(v) for k, v in where.items())]
        require(len(hits) == 1, f'Expected one CSV row in {rel}: {where}, found {len(hits)}')
        return hits[0]

    def cell(self, name, rel, field, unit='', **where):
        value = self.row(rel, **where)[field]
        value = None if value == '' else value   # an empty archived field is an explicit null
        locator = 'row['+','.join(f'{k}={v}' for k, v in where.items())+'].'+field
        return dict(kind='data', name=name, value=value, unit=unit, source=rel, locator=locator)

    def derived(self, name, value, rel, how, unit=''):
        self.path(rel)
        return dict(kind='derived', name=name, value=value, unit=unit, source=rel, locator=how)

    def quote(self, rel, text):
        require(norm(text) in norm(self.text(rel)), f'Quotation not found in {rel}: {text}')
        return dict(kind='quote', name='quotation', value=norm(text), unit='', source=rel,
                    locator='literal text (whitespace and emphasis markers normalised)')

    def srcline(self, rel, fragment):
        """Single line of an archived source-text rendition containing `fragment`."""
        hits = [(i+1, line) for i, line in enumerate(self.text(rel).split('\n')) if fragment in line]
        require(len(hits) == 1, f'Expected one line with {fragment!r} in {rel}, found {len(hits)}')
        return dict(kind='source_text', name='archived source text', value=fragment, unit='',
                    source=rel, locator=f'line {hits[0][0]}')


def fmt(record):
    unit = (' '+record['unit']) if record['unit'] else ''
    if record['kind'] in ('quote', 'source_text'):
        return f'"{record["value"]}" [{record["source"]}#{record["locator"]}]'
    return f'{record["name"]}={show(record["value"])}{unit} {{{record["kind"]}}} [{record["source"]}#{record["locator"]}]'


class Matrix:
    def __init__(self, ev, unknowns):
        self.ev, self.unknowns, self.rows, self.count = ev, unknowns, [], {}

    def add(self, group, reference, case, quantity, status, admitted_use, not_for, page,
            conditions, resolution, uncertainty, unknown, evidence=(), cite=(), notes='',
            strict=STRICT):
        require(status in STATUSES, 'Unknown status '+status)
        require(admitted_use and uncertainty and resolution and conditions, 'Incomplete row: '+quantity)
        for key in unknown:
            require(key in self.unknowns, 'Unregistered unknown '+key)
        self.count[group] = self.count.get(group, 0)+1
        row_id = f'{group}{self.count[group]:02d}'
        evidence = [e for e in evidence if e]
        sources = []
        primary = [REPORT_PDF] if reference.startswith(('experiment_2010', 'mses_2010', 'report_2010')) else []
        for rel in list(cite)+[e['source'] for e in evidence]+primary:
            if rel not in sources:
                sources.append(rel)
        require(sources, 'Row without a source: '+quantity)
        for key in unknown:
            self.unknowns[key]['rows'].append(row_id)
        self.rows.append(dict(
            row_id=row_id, reference_group=GROUPS[group], reference=reference, case=case,
            quantity=quantity, admission_status=status, admitted_use=admitted_use,
            not_admitted_for=not_for, strict_matched_geometry_validation=strict,
            conditions=conditions, source_page_figure=page, measurement_resolution=resolution,
            reported_uncertainty=uncertainty,
            unknown_quantities=[dict(id=k, text=self.unknowns[k]['text']) for k in unknown],
            evidence=evidence, source_paths=sources,
            source_sha256=[self.ev.digest(s) for s in sources], notes=notes))
        return row_id

    def ordered(self):
        return sorted(self.rows, key=lambda r: r['row_id'])

    def flat(self):
        out = []
        for r in self.ordered():
            flat = dict(r)
            flat['unknown_quantities'] = ('; '.join(f'unknown {u["id"]}: {u["text"]}' for u in r['unknown_quantities'])
                                          or 'none identified in the cited files')
            flat['values'] = ' | '.join(fmt(e) for e in r['evidence'] if e['kind'] in ('data', 'derived'))
            flat['evidence_quotes'] = ' | '.join(fmt(e) for e in r['evidence'] if e['kind'] in ('quote', 'source_text'))
            flat['source_paths'] = '; '.join(r['source_paths'])
            flat['source_sha256'] = '; '.join(r['source_sha256'])
            del flat['evidence']
            out.append(flat)
        return out


def write_csv(path, rows, fields=None):
    fields = fields or list(rows[0].keys())
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')


def neighbours(values, target):
    """Nearest archived abscissae below/above a target (no interpolation performed here)."""
    below = [v for v in values if v <= target]
    above = [v for v in values if v >= target]
    return (max(below) if below else None, min(above) if above else None)


def build_unknowns(ev):
    """Registry of explicitly unknown quantities, each tied to an archived statement or flag."""
    U = {}

    def add(key, text, *basis):
        U[key] = dict(text=text, basis=[b for b in basis if b], rows=[])

    add('U01', 'exact as-built tested contour and original modified author coordinates',
        ev.val('original_author_coordinates_recovered', GEO+'summary.json', 'original_author_coordinates_recovered'),
        ev.val('tested_geometry_benchmark_admitted', GEO+'summary.json', 'tested_geometry_benchmark_admitted'))
    add('U02', 'measured element placement/rigging (gap, overlap, deflection) of the tested model',
        ev.quote(GEO+'FINDINGS.md', 'The gap and overlap are our explicit geometric diagnostics, not a recovered rigging table.'))
    add('U03', 'published-MSES calculation geometry, input deck and detailed solver settings',
        ev.val('actual_geometry', GEO+'pressure/source.json', 'actual_geometry'),
        ev.val('original_run_deck_recovered', TRIP+'evidence_qualification.json', 'mses', 'original_run_deck_recovered'))
    add('U04', 'case-specific measured transition, separation, reattachment or bubble-endpoint locations',
        ev.val('measured_transition_endpoints_available', GEO+'pressure/source.json', 'measured_transition_endpoints_available'),
        ev.val('measured_case_specific_bubble_endpoints_found', TRIP+'evidence_qualification.json', 'pressure_coverage', 'measured_case_specific_bubble_endpoints_found'))
    add('U05', 'experimental Cp measurement uncertainty',
        ev.val('pressure_uncertainty', GEO+'pressure/source.json', 'pressure_uncertainty'),
        ev.quote(GEO+'FINDINGS.md', 'Exact experimental Cp uncertainty and bubble endpoints were not recovered.'))
    add('U06', 'magnitude and position of the fore-nose pressure peak between the first supported taps',
        ev.quote(LBA+'FINDINGS.md', 'The modeled nose peak is displaced relative to the sparse taps, so comparing the continuous computed minimum directly with the smallest measured marker would exaggerate what is known.'),
        ev.quote(PLAN, "The model's unresolved between-tap nose peak must remain visibly unmeasured."))
    add('U07', 'published-MSES pressure at CL -0.072, 0.056 and 0.175',
        ev.val('mses_pressure_cases', LBA+'summary.json', 'mses_pressure_cases'),
        ev.quote(LB_PROTOCOL, 'Only the 0.284 case has a matching published MSES pressure reference (nominal CL 0.28); mark the other three pressure references missing.'))
    add('U08', 'published-MSES drag at CL -0.072 (outside the published curve)',
        ev.quote(LBA+'FINDINGS.md', 'MSES drag begins at CL −0.01220, so the first reference remains missing.'))
    add('U09', 'published-MSES transition/separation/reattachment locations at any case',
        ev.quote(LBA+'FINDINGS.md', 'No MSES transition locations or new pressure curves are invented for the other conditions.'))
    add('U10', 'measured wake turbulence/velocity profiles',
        ev.val('wake_turbulence_profiles_found', TRIP+'evidence_qualification.json', 'pressure_coverage', 'wake_turbulence_profiles_found'))
    add('U11', 'aft or fore Cp with fore-only or fore+aft trips',
        ev.val('fore_only_fixed_Cp_found_in_report', TRIP+'evidence_qualification.json', 'pressure_coverage', 'fore_only_fixed_Cp_found_in_report'),
        ev.val('both_fixed_Cp_found_in_report', TRIP+'evidence_qualification.json', 'pressure_coverage', 'both_fixed_Cp_found_in_report'))
    add('U12', 'trip effectiveness: transition onset per case and threshold CL of the fore-upper grit',
        ev.val('confound', TRIP+'evidence_qualification.json', 'confounds', 1))
    add('U13', 'which tap-coordinate convention (printed Table II or plotted) is correct, and the true position of three tap centres outside the reconstructed branches',
        ev.val('retained_exact_table_check_failure', GEO+'verification.json', 'checks', 'retained_exact_table_check_failure'),
        ev.quote(GEO+'FINDINGS.md', 'Finite tap size, source faceting and coordinate conventions are possible contributors; these comparisons do not identify one.'))
    offsets = [abs(float(r['caption_minus_appendix_CL'])) for r in ev.rows(LBA+'targets.csv')]
    add('U14', f'unrounded experimental CL/alpha of the pressure snapshots (captions are rounded; appendix rows differ by up to {max(offsets):.4f} CL)',
        ev.derived('maximum |caption CL - appendix CL| over the four cases', round(max(offsets), 6), LBA+'targets.csv', 'max(|caption_minus_appendix_CL|)'),
        ev.quote(LBA+'FINDINGS.md', 'Its full ±0.0005 caption-rounding interval is not covered; the unsupported endpoint is null in the output.'))
    add('U15', 'numerical values of the 2014 full-e^N MSES curves (lift, drag, pressure, transition)',
        ev.quote(GEO+'FINDINGS.md', 'The 2014 Figures 11–12 offer another pressure reference at Re=0.7/1.5 million and CL=0.493/0.518; their full-e^N curves have not been digitized here.'))
    add('U16', 'numerical values of the 2014 OVERFLOW results',
        ev.quote(REGISTER, 'the paper compares lift, drag and pressure with MSES and OVERFLOW.'))
    add('U17', 'grid convergence and solution uniqueness of the archived MultiFoil states',
        ev.val('numerical_convergence_established', FRA+'summary.json', 'numerical_convergence_established'),
        ev.val('solution_uniqueness_established', FRA+'summary.json', 'solution_uniqueness_established'))
    add('U18', 'fore-only cubic reconstruction result at CL 0.175 and its paired CL 0.056 -> 0.175 drag change (target failed)',
        ev.val('paired_drag_changes_counts.fore_source_cubic', FRA+'summary.json', 'paired_drag_changes_counts', 'fore_source_cubic'))
    add('U19', 'measurement uncertainty stated by the 2010 report itself for CL, CD, CM, alpha',
        ev.derived("occurrences of 'uncertaint' in the archived 2010 report text",
                   ev.text(REPORT_TXT).lower().count('uncertaint'), REPORT_TXT, "case-insensitive count of 'uncertaint'"))
    inventory = ev.json(TRIP+'report_page_inventory.json')
    last12 = max(i for i, page in enumerate(inventory) if 'Figure 12.-' in page['caption_text'])
    add('U20', 'pressure markers of Figure 7(e-x) (including the 7(k) copy at CL 1.041) and of the third Figure 12 panel',
        ev.val('pressure_pages extracted', TRIP+'source_manifest.json', 'scope', 'pressure_pages'),
        ev.val('caption of the third Figure 12 panel', TRIP+'report_page_inventory.json', last12, 'caption_text'))
    add('U21', 'actual as-built contour error and measured rigging uncertainty',
        ev.quote(RC_PROTOCOL, 'it is not measured rigging uncertainty, a probability distribution, or a confidence bound'),
        ev.quote(RCA+'FINDINGS.md', 'They do not exhaust geometry uncertainty: original modified coordinates, as-built shape, nose/edge details, horizontal placement and tunnel reference conditions remain unverified.'))
    add('U22', 'conditions, geometry and data of the Twiss (2019) thesis campaign',
        ev.val('local_download_status', TRIP+'research_search.json', 'additional_lead', 'local_download_status'))
    add('U23', 'measured per-element drag',
        ev.quote(LBA+'FINDINGS.md', 'measured element drag is unavailable'))
    add('U24', 'quantitative oil-flow bubble lengths',
        ev.quote(TRIP+'FINDINGS.md', 'The oil-flow statement cannot be converted into a quantitative bubble-shortening curve from these sparse taps alone.'))
    add('U25', 'provenance, registration and coordinate frames of the public .pts files',
        ev.val('author_origin_and_provenance', GEO+'public_pts_provenance.json', 0, 'author_origin_and_provenance'))
    add('U26', 'magnitude of the tunnel-boundary and probe corrections already applied by the authors',
        ev.quote(PB+'FINDINGS.md', 'The source already reports standard tunnel-boundary corrections; none was applied again here.'))
    add('U27', 'polars at Re 0.5 and 1.5 million (present in the report, outside the bounded qualification; never extracted)',
        ev.val('excluded_conditions', TRIP+'evidence_qualification.json', 'excluded_conditions'))
    return U


def tap_resolution(ev):
    """Resolution facts derived from the printed Table II tap list."""
    rel = GEO+'pressure_tap_geometry_support.csv'
    taps = ev.rows(rel)
    chord_mm = ev.json(GEO+'pressure/source.json')['reference_chord_m']*1000

    def xs(element, surface, supported=True):
        return [r['x_over_c'] for r in taps if r['element'] == element and r['surface'] == surface
                and (not supported or r['unsupported_center'] == '0')]

    lower, upper = xs('fore', 'lower'), xs('fore', 'upper')
    spacing = float(lower[1])-float(lower[0])
    counts = {f'{e}_{s}': (len(xs(e, s, False)), len(xs(e, s)))
              for e in ('fore', 'aft') for s in ('upper', 'lower')}
    unsupported = [r for r in taps if r['unsupported_center'] == '1']
    records = [
        ev.derived('printed Table II taps', len(taps), rel, 'row count'),
        ev.derived('taps (listed, supported) per surface', counts, rel, 'count by element/surface; supported = unsupported_center 0'),
        ev.derived('first two supported fore-lower taps x/c', lower[:2], rel, 'first rows with element=fore,surface=lower'),
        ev.derived('first two supported fore-upper taps x/c', upper[:2], rel, 'first rows with element=fore,surface=upper,unsupported_center=0'),
        ev.derived('spacing of the first two fore-lower taps', round(spacing, 5), rel, 'difference of x_over_c', 'c'),
        ev.derived('same spacing', round(spacing*chord_mm, 3), GEO+'pressure/source.json', 'spacing * reference_chord_m * 1000', 'mm'),
        ev.derived('supported fore taps with x/c <= 0.06 (upper, lower)',
                   [sum(float(x) <= .06 for x in upper), sum(float(x) <= .06 for x in lower)], rel, 'count'),
        ev.derived('unsupported tap centres (element surface x/c: mm outside branch)',
                   [f'{r["element"]} {r["surface"]} {r["printed_value"]}: {float(r["outside_extent_mm"]):.3f} mm'
                    for r in unsupported], rel, 'rows with unsupported_center=1'),
        ev.val('reported_tap_diameter_mm', GEO+'summary.json', 'reported_tap_diameter_mm', unit='mm'),
    ]
    text = (f'{sum(c[1] for c in counts.values())} supported taps of {len(taps)} printed '
            f'({counts["fore_upper"][1]} fore upper, {counts["fore_lower"][1]} fore lower, '
            f'{counts["aft_upper"][1]} aft upper, {counts["aft_lower"][1]} aft lower); first supported fore-lower taps '
            f'x/c {lower[0]} and {lower[1]} (spacing {spacing:.5f}c = {spacing*chord_mm:.2f} mm); first supported '
            f'fore-upper tap x/c {upper[0]} (printed tap 0.00000 lies outside the reconstructed branch); '
            f'orifice diameter {ev.json(GEO+"summary.json")["reported_tap_diameter_mm"]} mm')
    return records, text, lower, upper


def reconcile(ev):
    """Compare every sha256 map of the earlier S414 receipts with the files as they are now."""
    cache, rows = {}, []

    def digest(path):
        if path not in cache:
            cache[path] = sha(path)
        return cache[path]

    def maps(node, trail):
        if isinstance(node, dict):
            if node and all(isinstance(v, str) and HEX.match(v) for v in node.values()):
                yield trail, node
            else:
                for key, value in node.items():
                    yield from maps(value, trail+[str(key)])
        elif isinstance(node, list):
            for i, value in enumerate(node):
                yield from maps(value, trail+[str(i)])

    for directory in RECONCILED:
        for receipt in RECEIPTS:
            rel = directory+receipt
            if not (ROOT/rel).is_file():
                continue
            for trail, mapping in maps(ev.json(rel), []):
                matched, changed, unresolved = 0, [], []
                for name, recorded in mapping.items():
                    candidates = [ROOT/name, ROOT/directory/name, (ROOT/rel).parent/name]
                    path = next((c for c in candidates if c.is_file()), None)
                    if path is None:
                        unresolved.append(name)
                        continue
                    current = digest(path)
                    key = str(path.relative_to(ROOT))
                    if key not in ev.files and key not in ev.living:
                        ev.reconciled[key] = current
                    if current == recorded:
                        matched += 1
                    else:
                        changed.append(key)
                rows.append(dict(receipt=rel, sha256_map='/'.join(trail), entries=len(mapping),
                                 matched=matched, changed=len(changed), unresolved=len(unresolved),
                                 changed_files='; '.join(changed), unresolved_names='; '.join(unresolved)))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    require(not (out/'summary.json').exists(), 'Preserve completed matrix; choose a new output directory')
    ev = Evidence()
    unknowns = build_unknowns(ev)
    M = Matrix(ev, unknowns)
    coverage, classes, notes = [], [], []

    # ------------------------------------------------------------------ shared facts
    psrc = GEO+'pressure/source.json'
    quali = TRIP+'evidence_qualification.json'
    tsum = TRIP+'summary.json'
    primary = ev.json(quali)['primary_condition']
    secondary = ev.json(quali)['secondary_condition']
    chord_mm = ev.json(psrc)['reference_chord_m']*1000
    cond_free = (f'Re {primary["reynolds"]:,}; Mach {primary["mach"]}; transition free; modified tested S414 model, '
                 f'reference chord {primary["reference_chord_m"]} m')
    cond_mses = (f'Re {ev.json(psrc)["reynolds"]:,.0f}; Mach {ev.json(psrc)["mach"]}; transition '
                 f'{ev.json(psrc)["transition"]}; Ncrit {ev.json(psrc)["mses_ncrit"]}; '
                 f'{ev.json(psrc)["mses_transition"]}; mses_executable_run={ev.json(psrc)["mses_executable_run"]}')
    gate = ev.quote(PLAN, 'Strict matched-geometry validation stays conditional until geometry is qualified.')
    tap_records, tap_text, lower_taps, upper_taps = tap_resolution(ev)
    inventory = {p['pdf_page']: p for p in ev.json(TRIP+'report_page_inventory.json')}
    report_uncertainty_hits = ev.text(REPORT_TXT).lower().count('uncertaint')
    article_drag = [ev.quote(REGISTER, 'Reported drag uncertainty is 0.5 count in the low-drag range, rising toward stall.'),
                    ev.srcline(ARTICLE_TXT, '0.00005 (0.5'), ev.srcline(ARTICLE_TXT, '0.00015 (1.5 drag counts) as the'),
                    ev.srcline(ARTICLE_TXT, 'efficients measured in the low-drag range is'),
                    ev.srcline(ARTICLE_TXT, 'angle of attack approaches stall.')]
    article_lift = [ev.srcline(ARTICLE_TXT, 'increasing angle of attack [23]. At higher lift, for which the un-'),
                    ev.srcline(ARTICLE_TXT, '0.005, and the pitching-moment coefficients have an'),
                    ev.srcline(ARTICLE_TXT, '0.002. The uncertainty of the drag co-')]
    drag_uncertainty = (
        f"2010 report: not reported (archived report text has {report_uncertainty_hits} occurrences of 'uncertaint'). "
        'Overlapping 2014 article: 0.00005 (0.5 drag counts) in the low-drag range, increasing to 0.00015 '
        '(1.5 drag counts) approaching stall [archived article text; the relation/sign glyphs are unreadable in '
        'the text rendition]. Extraction: plotted markers and appendix rows agree within '
        f'{ev.json(tsum)["maximum_plot_table_CD_difference"]:.3e} CD; this is not measurement uncertainty.')
    cp_uncertainty = ('not reported (no case-specific Cp uncertainty budget recovered); extraction offsets are '
                      'recorded separately and are not measurement uncertainty')
    experiments = ev.rows(TRIP+'experimental_tables.csv')
    mses_curves = ev.rows(TRIP+'published_mses_drag.csv')
    free_pressure = ev.rows(TRIP+'free_boundary_pressure.csv')
    lb_targets = ev.rows(LBA+'targets.csv')
    mses_free = sorted(float(r['CL']) for r in mses_curves if r['reynolds'] == '1000000' and r['configuration'] == 'free')
    table_free = sorted(float(r['CL']) for r in experiments if r['reynolds'] == '1000000'
                        and r['configuration'] == 'free' and r['prestall_branch'] == 'True')
    study_of = {analysis: ev.json(study+'study.json') for analysis, study in STUDIES}
    fig13 = ev.val('product (figure and pages)', DRAG+'source.json', 'product')
    page13 = fig13['value'].split(', ', 1)[1]
    pages_report = ev.json(TRIP+'report_page_inventory.json')
    fig12_captions = [ev.val('caption', TRIP+'report_page_inventory.json', i, 'caption_text')
                      for i, p in enumerate(pages_report) if 'Figure 12.-' in p['caption_text']]
    fig12_pages = [p['pdf_page'] for p in pages_report if 'Figure 12.-' in p['caption_text']]
    discussion = ev.quote(WAKE+'FINDINGS.md', 'printed pp. 9–11 (PDF pp. 11–13), describes oil-flow and acoustic checks')
    aft_upper = [float(r['x_over_c']) for r in ev.rows(GEO+'pressure_tap_geometry_support.csv')
                 if r['element'] == 'aft' and r['surface'] == 'upper' and .86 <= float(r['x_over_c']) <= .94]
    aft_spacing = [b-a for a, b in zip(aft_upper, aft_upper[1:])]
    aft_record = ev.derived('aft-upper tap spacing for 0.86 <= x/c <= 0.94 [min, max]',
                            [round(min(aft_spacing), 5), round(max(aft_spacing), 5)],
                            GEO+'pressure_tap_geometry_support.csv', 'differences of consecutive x_over_c', 'c')
    bubble_quote = ev.quote(GEO+'FINDINGS.md', 'The 2010 report interprets the predicted upper-aft bubble as downstream of the measured one at nominal CL 0.28 and 1.04.')
    cover = {}

    def mark(case, kind, status, detail, sources, page, row_id):
        cover.setdefault(case, {})[kind] = dict(status=status, detail=detail, sources=list(sources),
                                                page=page, row_id=row_id)

    # ------------------------------------------------------------------ A/B: four lower-bucket cases
    for t in lb_targets:
        key, fig, cl = t['case'], t['figure'], t['target_cl']
        case = f'CL {cl} (Figure {fig})'
        markers = [r for r in free_pressure if r['figure'] == fig]
        page = int(markers[0]['pdf_page'])
        require(all(int(r['pdf_page']) == page for r in markers), 'Mixed pages for Figure '+fig)
        printed = inventory[page]['printed_page']
        where_fig = f'Figure {fig}, PDF p. {page} / printed p. {printed}, RDECOM TR 10-D-112'
        caption = ev.val('figure caption', TRIP+'report_page_inventory.json', page-1, 'caption_text')
        require(ev.json(TRIP+'report_page_inventory.json')[page-1]['pdf_page'] == page, 'Page inventory order changed')
        taps = [r for r in ev.rows(LBA+'pressure_at_taps.csv')
                if r['case'] == key and r['coordinate'] == 'table' and r['method'] == 'multifoil']
        supported = [r for r in taps if r['included'] == 'True']
        excluded = [r for r in ev.rows(LBA+'excluded_taps.csv') if r['case'] == key]
        checks = [r for r in ev.rows(TRIP+'pressure_tap_checks.csv') if r['figure'] == fig]
        max_dx = max(float(r['maximum_plotted_table_x_difference']) for r in checks)
        evidence = [
            caption,
            ev.derived('pressure markers extracted', len(markers), TRIP+'free_boundary_pressure.csv', f'count(figure={fig})'),
            ev.derived('repeated plotting endpoints among them', sum(r['endpoint_duplicate'] == 'True' for r in markers),
                       TRIP+'free_boundary_pressure.csv', 'count(endpoint_duplicate=True)'),
            ev.derived('supported independent taps scored', len(supported), LBA+'pressure_at_taps.csv',
                       f'count(case={key},coordinate=table,method=multifoil,included=True)'),
            ev.derived('excluded observations (reason: count)',
                       {reason: sum(r['reason'] == reason for r in excluded) for reason in sorted({r['reason'] for r in excluded})},
                       LBA+'excluded_taps.csv', f'count by reason (case={key})'),
            ev.derived('maximum |plotted - printed| tap x in this figure', max_dx, TRIP+'pressure_tap_checks.csv',
                       f'max over surfaces (figure={fig})', 'c'),
            ev.quote(LB_PROTOCOL, 'Score all commonly supported independent taps, with Table-II coordinates primary and original plotted coordinates as a sensitivity.'),
            ev.quote(LB_PROTOCOL, 'Exclude repeated endpoints and the existing unsupported tap centers; disclose masks and preserve every source marker.'),
            ev.quote(LBA+'FINDINGS.md', 'No measured case-specific pressure uncertainty, transition map or wake profile has been recovered to separate all remaining alternatives.'),
            ev.quote(PLAN, 'Qualified pressure markers and matching drag records, with original masks and conditions'),
        ]
        cite = [TRIP+'free_boundary_pressure.csv', GEO+'pressure_tap_geometry_support.csv']
        extra = ''
        if key == 'cl028':
            evidence += [
                ev.val('strict Figure 7(d)/12(a) copy check passed', tsum, 'strict_repeated_pressure_copy_check_pass'),
                ev.val('maximum Cp difference between the two copies', tsum, 'maximum_repeated_pressure_Cp_difference'),
                ev.val('maximum plotted-x difference between the two copies', tsum, 'maximum_repeated_pressure_x_difference', unit='c'),
                ev.derived('markers in the Figure 12(a) copy', sum(r['case'] == 'cl028' for r in ev.rows(GEO+'pressure/experiment.csv')),
                           GEO+'pressure/experiment.csv', 'count(case=cl028)'),
            ]
            cite += [GEO+'pressure/experiment.csv', TRIP+'repeated_pressure_crosscheck.csv',
                     TRIP+'pressure_duplicate_strict_check_failed.log']
            extra = (' A second copy of the same measurement exists as Figure 12(a) (PDF p. 66); the strict '
                     'copy check fails and both copies are retained without alignment.')
        row_id = M.add('A', EXP, case, 'pressure (surface Cp at taps)', 'admitted',
                       'Qualified free-transition pressure markers at the supported taps, with the original masks '
                       '(repeated endpoints and three unsupported tap centres excluded); printed Table II coordinates '
                       'primary, plotted coordinates as a sensitivity; pressure diagnostic on the nominal reconstructed outline',
                       'Exact matched-geometry benchmark; bubble/transition endpoints; magnitude of the nose peak between '
                       'taps; any fitted alignment or pressure-averaging correction',
                       where_fig, cond_free, tap_text, cp_uncertainty, ['U05', 'U06', 'U13', 'U01'], evidence, cite,
                       notes='Strongest-rise intervals are tap-sampling statistics, not bubble endpoints.'+extra)
        mark(key, 'experiment_pressure', 'exists', f'{len(markers)} markers; {len(supported)} supported taps',
             cite, where_fig, row_id)

        # experimental drag
        appendix = dict(reynolds='1000000', configuration='free', index=t['appendix_row_index'])
        arow = ev.row(TRIP+'experimental_tables.csv', **appendix)
        low, high = neighbours(table_free, float(cl))
        evidence = [
            caption,
            ev.cell('appendix row (alpha CL CD CM as printed)', TRIP+'experimental_tables.csv', 'raw_tokens', **appendix),
            ev.cell('experiment CD interpolated at caption CL', LBA+'targets.csv', 'experiment_CD', case=key),
            ev.cell('caption CL minus appendix CL', LBA+'targets.csv', 'caption_minus_appendix_CL', case=key),
            ev.cell('CD at caption CL - 0.0005', LBA+'targets.csv', 'experiment_rounding_CD_low', case=key),
            ev.cell('CD at caption CL + 0.0005', LBA+'targets.csv', 'experiment_rounding_CD_high', case=key),
            ev.cell('rounding interval fully inside the table', LBA+'targets.csv', 'experiment_rounding_fully_supported', case=key),
            ev.derived('nearest pre-stall appendix CL below/above the caption CL', [low, high],
                       TRIP+'experimental_tables.csv', 'reynolds=1000000,configuration=free,prestall_branch=True'),
            ev.val('maximum plotted-marker vs appendix CD difference', tsum, 'maximum_plot_table_CD_difference'),
            ev.quote(DRAG+'FINDINGS.md', 'Caption CL rounding and piecewise-linear reference interpolation do not define measurement uncertainty.'),
            fig13,
        ]+article_drag
        cite = [TRIP+'experimental_tables.csv', DRAG+'experiment_drag.csv']
        unknown = ['U19', 'U23']
        extra = ''
        if key == 'cl028':
            evidence += [
                ev.cell('experiment CD interpolated from Figure 13(c) markers', DRAG+'target_references.csv', 'experiment_CD_interpolated', case='cl028'),
                ev.quote(LBA+'FINDINGS.md', 'Using appendix values rather than the earlier Figure-13 marker extraction changes the 0.284 experimental reference by only about 0.00059 count.'),
            ]
            extra = ' Two archived reference conventions exist (appendix table and Figure 13(c) markers); both are kept.'
        if t['experiment_rounding_fully_supported'] != 'True':
            unknown.append('U14')
            evidence.append(ev.quote(LBA+'FINDINGS.md', 'Its full ±0.0005 caption-rounding interval is not covered; the unsupported endpoint is null in the output.'))
            extra += ' The lower caption-rounding endpoint lies outside the appendix table and is null.'
        where = (f'Appendix table PDF p. {arow["pdf_page"]} / printed p. {arow["printed_page"]}; caption of Figure {fig}, '
                 f'PDF p. {page}; plotted markers {page13}')
        row_id = M.add('A', EXP, case, 'drag (wake-survey section CD)', 'admitted',
                       'Matched-lift drag reference: printed appendix table interpolated piecewise-linearly at the rounded '
                       'caption CL on the increasing-alpha branch, without extrapolation',
                       'Per-element drag; a statistical confidence band; any fitted drag correction',
                       where, cond_free,
                       f'appendix row printed as "{arow["raw_tokens"]}" (alpha, CL, CD, CM); nearest pre-stall appendix CL '
                       f'{low} / {high}; one caption value per pressure snapshot',
                       drag_uncertainty, unknown, evidence, cite,
                       notes='Caption rounding and interpolation offsets are not measurement uncertainty.'+extra)
        mark(key, 'experiment_drag', 'exists', f'CD {t["experiment_CD"]} at CL {cl} (appendix interpolation)', cite, where, row_id)

        # lift / incidence association
        evidence = [
            caption,
            ev.cell('caption alpha', LBA+'targets.csv', 'experiment_caption_alpha', unit='deg', case=key),
            ev.cell('appendix alpha', LBA+'targets.csv', 'appendix_alpha', unit='deg', case=key),
            ev.cell('appendix CL', LBA+'targets.csv', 'appendix_CL', case=key),
            ev.cell('caption CL minus appendix CL', LBA+'targets.csv', 'caption_minus_appendix_CL', case=key),
            ev.quote(PLAN, 'Maintain three separate comparison bases: identical incidence for mechanism isolation; identical lift for drag validation; identical prescribed edge velocity for a conditional boundary-layer comparison.'),
            ev.quote(PLAN, 'do not apply another fitted angle correction'),
        ]+article_lift
        M.add('A', EXP, case, 'lift and incidence association (caption versus appendix)', 'admitted',
              'Case identification; matched-lift target (rounded caption CL) and measured incidence for same-angle '
              'mechanism comparisons, kept as two separate comparison bases',
              'A fitted angle or lift correction; treating fixed-angle drag as matched-lift validation',
              f'caption of Figure {fig}, PDF p. {page}; appendix table PDF p. {arow["pdf_page"]} / printed p. {arow["printed_page"]}', cond_free,
              f'caption CL to 3 decimals and alpha to 0.01 deg; appendix CL to 4 decimals and alpha to 0.001 deg '
              f'(row "{arow["raw_tokens"]}")',
              "2010 report: not reported. Overlapping 2014 article: lift-coefficient uncertainty 0.005 'at higher lift' "
              '[archived article text]; alpha uncertainty: not reported',
              ['U14', 'U19', 'U26'], evidence, [TRIP+'experimental_tables.csv', LBA+'targets.csv'])

        # experimental transition / reattachment
        index = ev.find(tsum, ['strongest_sampled_pressure_recovery'], figure=fig, coordinate='table_x')
        evidence = [
            ev.val('measured_case_specific_bubble_endpoints_found', quali, 'pressure_coverage', 'measured_case_specific_bubble_endpoints_found'),
            ev.val('largest sampled aft-upper pressure-rise interval, left tap x/c (sampling statistic, not an endpoint)', tsum, 'strongest_sampled_pressure_recovery', index, 'left_x'),
            ev.val('same interval, right tap x/c', tsum, 'strongest_sampled_pressure_recovery', index, 'right_x'),
            ev.val('Cp rise across that interval', tsum, 'strongest_sampled_pressure_recovery', index, 'delta_Cp'),
            ev.quote(TRIP+'FINDINGS.md', 'The location is not monotonic across these four snapshots; a tapwise pressure rise is not an exact bubble length or reattachment measurement.'),
            ev.quote(LB_PROTOCOL, 'no measured case-specific transition coordinates or wake turbulence are available'),
            ev.quote(PLAN, 'Exact case-specific transition/reattachment locations and sufficiently resolved nose pressure'),
            discussion, aft_record,
        ]
        qualitative = 'no case-specific qualitative statement recorded in the existing files'
        if key in ('clm0072', 'cl0056'):
            evidence.append(ev.quote(TRIP+'FINDINGS.md', 'The author describes turbulent fore-lower flow for 7(a), and laminar flow surviving its pressure peak near the lower boundary at 7(b).'))
            qualitative = 'qualitative author description only (see quotation)'
        elif key == 'cl028':
            evidence += [ev.quote(WAKE+'FINDINGS.md', 'confirms laminar flow on both fore surfaces within the low-drag range'),
                         ev.quote(WAKE+'FINDINGS.md', 'Its Figure 7(d), our CL 0.284 pressure case, lies inside that range.')]
            qualitative = 'qualitative statement only: fore surfaces laminar within the low-drag range (see quotations)'
        else:
            evidence.append(ev.quote(LBA+'FINDINGS.md', 'inside the bucket, oil-flow and acoustic observations support laminar fore surfaces'))
            qualitative = ('general in-bucket statement only; the existing files do not state explicitly whether '
                           'CL 0.175 is inside the bucket: unknown')
        row_id = M.add('A', EXP, case, 'transition / separation / reattachment location', 'missing',
                       'none as a quantitative reference; '+qualitative,
                       'Any numerical transition, bubble-length or reattachment target',
                       where_fig+'; report discussion printed pp. 9-11 (PDF pp. 11-13)', cond_free,
                       f'aft-upper tap spacing {min(aft_spacing):.5f}c to {max(aft_spacing):.5f}c for 0.86 <= x/c <= 0.94; a pressure rise between two taps only',
                       'not reported', ['U04', 'U24', 'U10'], evidence, [TRIP+'pressure_recovery_intervals.csv'],
                       notes='The sampled largest-rise interval is a tap-sampling statistic and is not an endpoint.')
        mark(key, 'experiment_transition_reattachment_location', 'missing',
             'no measured coordinates; '+qualitative, [quali, tsum], where_fig, row_id)

        # published MSES drag
        low, high = neighbours(mses_free, float(cl))
        base = [
            ev.derived('published MSES drag curve vertices (Re 1e6, free)', len(mses_free), TRIP+'published_mses_drag.csv',
                       'count(reynolds=1000000,configuration=free)'),
            ev.derived('lowest CL of the published curve', min(mses_free), TRIP+'published_mses_drag.csv', 'min(CL)'),
            ev.derived('nearest curve vertices CL below/above the target', [low, high], TRIP+'published_mses_drag.csv', 'neighbours of target CL'),
            ev.quote(LB_PROTOCOL, 'Published 2010 envelope-MSES drag covers 0.056, 0.175 and 0.284, not −0.072.'),
        ]
        mses_page = page13+', RDECOM TR 10-D-112'
        base += [fig13, ev.val('extraction_precision', DRAG+'source.json', 'extraction_precision')]
        if t['mses_CD']:
            evidence = [
                ev.cell('published MSES CD interpolated at the target CL', LBA+'targets.csv', 'mses_CD', case=key),
                ev.cell('MSES minus experiment', LBA+'targets.csv', 'mses_error_CD_counts', unit='counts', case=key),
                ev.cell('rounding interval fully inside the curve', LBA+'targets.csv', 'mses_rounding_fully_supported', case=key),
            ]+base+[ev.quote(PLAN, 'Drag at 0.056/0.175/0.284; pressure at nominal 0.28 as the 0.284 control')]
            extra = ''
            if key == 'cl028':
                evidence.append(ev.cell('MSES minus experiment using Figure 13(c) experimental markers', DRAG+'target_references.csv',
                                        'mses_minus_experiment_counts', unit='counts', case='cl028'))
                extra = ' Two archived bias values differ only through the experimental reference convention.'
            row_id = M.add('B', MSES10, case, 'drag (published curve)', 'admitted',
                           'Published-MSES drag at the same lift, as a second reference curve',
                           'Experiment; a run of an installed MSES; the 2014 full-e^N series',
                           mses_page, cond_mses,
                           f'piecewise-linear published curve; nearest vertices at CL {low:.4f} and {high:.4f}',
                           'not reported (computation; PDF vector coordinates carry no physical uncertainty)',
                           ['U03'], evidence, [TRIP+'published_mses_drag.csv', DRAG+'mses_drag.csv'],
                           notes='Interpolated between sparse published vertices; no extrapolation.'+extra)
            mark(key, 'mses2010_drag', 'exists', f'CD {t["mses_CD"]} at CL {cl}', [TRIP+'published_mses_drag.csv'], mses_page, row_id)
        else:
            evidence = [ev.cell('published MSES CD at the target CL (empty = null)', LBA+'targets.csv', 'mses_CD', case=key)]+base+[
                ev.quote(LBA+'FINDINGS.md', 'MSES drag begins at CL −0.01220, so the first reference remains missing.'),
                ev.quote(LB_PROTOCOL, 'No extrapolated MSES value is a validation reference.')]
            row_id = M.add('B', MSES10, case, 'drag (published curve)', 'missing',
                           'none: the published curve does not reach this lift',
                           'Extrapolation of the published curve', mses_page, cond_mses,
                           f'published curve starts at CL {min(mses_free):.5f}', 'not reported', ['U08', 'U03'],
                           evidence, [TRIP+'published_mses_drag.csv'])
            mark(key, 'mses2010_drag', 'missing', f'published curve starts at CL {min(mses_free):.5f}',
                 [TRIP+'published_mses_drag.csv'], mses_page, row_id)

        # published MSES pressure
        if t['mses_pressure_available'] == 'True':
            index = ev.find(psrc, ['cases'], name=key)
            item = ev.json(psrc)['cases'][index]
            vertices = sum(r['case'] == key for r in ev.rows(GEO+'pressure/mses_2010_envelope.csv'))
            where_p = f'Figure {item["figure"]}, PDF p. {item["page"]}, RDECOM TR 10-D-112'
            evidence = [
                ev.val('nominal theoretical CL of the published curve', psrc, 'cases', index, 'nominal_cl'),
                ev.val('associated experimental CL', psrc, 'cases', index, 'experimental_cl'),
                ev.val('matching', psrc, 'cases', index, 'matching'),
                ev.derived('theory curve vertices', vertices, GEO+'pressure/mses_2010_envelope.csv', f'count(case={key})'),
                ev.val('mses_pressure_cases', LBA+'summary.json', 'mses_pressure_cases'),
                ev.quote(PLAN, 'Drag at 0.056/0.175/0.284; pressure at nominal 0.28 as the 0.284 control'),
            ]
            row_id = M.add('B', MSES10, case, 'pressure (published curve)', 'conditional',
                           'Published-MSES pressure at nominal CL 0.28, used as the CL 0.284 control',
                           'Experiment; an exact-lift match (caption CL is rounded); the 2014 full-e^N series',
                           where_p, cond_mses, f'{vertices} published curve vertices over four surface branches',
                           'not reported', ['U03', 'U09'], evidence, [GEO+'pressure/mses_2010_envelope.csv'],
                           notes='Nominal theoretical CL 0.28 versus measured 0.284: a disclosed mismatch.')
            mark(key, 'mses2010_pressure', 'exists (nominal CL 0.28)', f'{vertices} vertices', [GEO+'pressure/mses_2010_envelope.csv'], where_p, row_id)
        else:
            evidence = [
                ev.cell('mses_pressure_available', LBA+'targets.csv', 'mses_pressure_available', case=key),
                ev.val('mses_pressure_cases', LBA+'summary.json', 'mses_pressure_cases'),
                ev.quote(LB_PROTOCOL, 'Only the 0.284 case has a matching published MSES pressure reference (nominal CL 0.28); mark the other three pressure references missing.'),
            ]
            row_id = M.add('B', MSES10, case, 'pressure (published curve)', 'missing',
                           'none: no published MSES pressure curve at this lift',
                           'Substituting the nominal CL 0.28 curve or the 2014 curves',
                           f'none at this lift (Figure 12 panels on PDF pp. {fig12_pages[0]}-{fig12_pages[-1]}: see archived captions)',
                           cond_mses, 'no curve', 'not reported', ['U07'], evidence+fig12_captions)
            mark(key, 'mses2010_pressure', 'missing', 'no published curve at this lift', [LBA+'targets.csv'], 'none', row_id)

        # published MSES transition
        evidence = [
            ev.derived('columns of the archived MSES pressure file', list(ev.rows(GEO+'pressure/mses_2010_envelope.csv')[0].keys()),
                       GEO+'pressure/mses_2010_envelope.csv', 'header (no transition field)'),
            ev.derived('columns of the archived MSES drag file', list(mses_curves[0].keys()), TRIP+'published_mses_drag.csv', 'header (no transition field)'),
            ev.quote(LBA+'FINDINGS.md', 'No MSES transition locations or new pressure curves are invented for the other conditions.'),
            ev.quote(PLAN, 'Pressure and transition for 0.056/0.175, exact calculation geometry/deck and detailed solver settings'),
        ]
        if key == 'cl028':
            evidence.append(bubble_quote)
        row_id = M.add('B', MSES10, case, 'transition / separation / reattachment location', 'missing',
                       'none: no MSES transition or bubble location is recorded in the existing files'
                       + ('; qualitative author interpretation only (see quotation)' if key == 'cl028' else ''),
                       'Inferring an MSES transition location from a Cp kink', 'none recorded', cond_mses, 'no data',
                       'not reported', ['U09', 'U03'], evidence)
        mark(key, 'mses_transition', 'missing', 'not recorded in any archived file', [GEO+'pressure/mses_2010_envelope.csv', TRIP+'published_mses_drag.csv'], 'none', row_id)

    # ------------------------------------------------------------------ A: series-level experiment rows
    summary = GEO+'summary.json'
    outlines = ev.json(summary)['repeated_outline_checks']
    geometry_evidence = [
        ev.val('source_geometry_reconstruction_pass', summary, 'source_geometry_reconstruction_pass'),
        ev.val('tested_geometry_benchmark_admitted', summary, 'tested_geometry_benchmark_admitted'),
        ev.val('original_author_coordinates_recovered', summary, 'original_author_coordinates_recovered'),
        ev.val('fore vertices', GEO+'geometry_metrics.json', 'report_test', 'fore_vertices'),
        ev.val('aft vertices', GEO+'geometry_metrics.json', 'report_test', 'aft_vertices'),
        ev.val('fore TE thickness', GEO+'geometry_metrics.json', 'report_test', 'fore_TE_thickness', unit='c'),
        ev.val('aft TE thickness', GEO+'geometry_metrics.json', 'report_test', 'aft_TE_thickness', unit='c'),
        ev.val('fore TE midpoint to nearest aft segment', GEO+'geometry_metrics.json', 'report_test', 'fore_TE_to_aft_surface_gap', unit='c'),
        ev.val('projected overlap', GEO+'geometry_metrics.json', 'report_test', 'projected_overlap', unit='c'),
        ev.val('maximum lower-surface thinning of the test shape versus design', summary, 'modification', 'max_lower_surface_thinning_chord', unit='c'),
        ev.derived('maximum fore coordinate difference over four repeated outlines', max(o['fore_max_coordinate_difference'] for o in outlines),
                   summary, 'max(repeated_outline_checks[].fore_max_coordinate_difference)', 'c'),
        ev.derived('maximum aft vertex-to-segment distance over four repeated outlines', max(o['aft_max_vertex_to_segment_distance'] for o in outlines),
                   summary, 'max(repeated_outline_checks[].aft_max_vertex_to_segment_distance)', 'c'),
        ev.val('reported_model_contour_tolerance_mm', summary, 'reported_model_contour_tolerance_mm', unit='mm'),
        ev.val('reference_chord_m', psrc, 'reference_chord_m', unit='m'),
        ev.quote(GEO+'FINDINGS.md', 'Admit the result as a traceable nominal test-shape reconstruction and a published-pressure diagnostic, not an exact matched-geometry benchmark.'),
        ev.quote(GEO+'FINDINGS.md', "This verifies extraction and repeated placement, not manufactured accuracy or the authors' original MSES input deck."),
        ev.quote(GEO+'FINDINGS.md', 'Figure 1, PDF page 22 / printed page 20, paths 9 and 0 provide the tested fore and aft outlines (83 and 68 points).'),
        ev.srcline(REPORT_TXT, 'to ensure an aerodynamically smooth finish. The measured model contour was within'),
        ev.srcline(REPORT_TXT, '0.13 mm (0.005 in.) of the prescribed shape.'),
        gate,
    ]
    geometry_row = M.add(
        'A', EXP, 'all cases (tested model)', 'geometry: nominal modified test outline (finite trailing edges)', 'conditional',
        'Traceable nominal test-shape reconstruction for exploratory comparisons that declare a reconstruction sensitivity',
        'Exact matched-geometry benchmark; certification of small design gains; as-built contour or rigging',
        'Figure 1, PDF p. 22 / printed p. 20 (paths 9 and 0); repeated outlines PDF pp. '
        + ', '.join(str(o['pdf_page']) for o in outlines), cond_free,
        f'{ev.json(GEO+"geometry_metrics.json")["report_test"]["fore_vertices"]} fore and '
        f'{ev.json(GEO+"geometry_metrics.json")["report_test"]["aft_vertices"]} aft published polyline vertices; '
        f'{ev.json(NOSE+"summary.json")["geometry"][0]["nose_vertices"]} source vertices at x/c <= 0.06',
        f'reported manufacturing tolerance: measured model contour within {ev.json(summary)["reported_model_contour_tolerance_mm"]} mm '
        'of the prescribed shape; extraction repeatability and reconstruction spread are separate classes '
        '(see uncertainty_classes.csv); as-built uncertainty: unknown',
        ['U01', 'U02', 'U21'], geometry_evidence
        + [ev.val('source nose vertices at x/c <= 0.06', NOSE+'summary.json', 'geometry', 0, 'nose_vertices')],
        [GEO+'report_test_fore.dat', GEO+'report_test_aft.dat', GEO+'geometry_extraction.json'],
        notes='Gap and overlap are geometric diagnostics of the drawing, not a recovered rigging table.')

    public = ev.rows(GEO+'public_geometry_comparison.csv')
    M.add('A', 'public_pts_files (unqualified geometry candidates)', 'all cases (tested model)',
          'geometry: public sharp-trailing-edge .pts coordinate files', 'not_admitted',
          'none as tested or MSES geometry; retained only as unqualified candidates',
          'Substitute for the modified finite-TE test model',
          'public teaching catalogue files s414.pts and s414-1.pts (no report page)', cond_free,
          'vertex-to-segment distances after a declared scale normalisation only',
          'not reported', ['U25', 'U01'],
          [ev.derived('file / element / max vertex-to-test-segment distance [c] / TE thickness [c]',
                      [f'{r["file"]} {r["element"]} {float(r["symmetric_max_vertex_to_test_segments"]):.6f} {r["TE_thickness"]}' for r in public],
                      GEO+'public_geometry_comparison.csv', 'all rows'),
           ev.val('tested fore TE thickness', GEO+'geometry_metrics.json', 'report_test', 'fore_TE_thickness', unit='c'),
           ev.quote(GEO+'FINDINGS.md', 'Neither file is adopted as the exact tested/MSES geometry.'),
           ev.quote(PLAN, 'Public sharp-TE coordinate files are not substitutes for the modified finite-TE test model.')],
          [GEO+'sources/s414.pts', GEO+'sources/s414-1.pts', GEO+'public_pts_provenance.json'], strict=NA)

    tap_pages = sorted({int(r['pdf_page']) for r in ev.rows(GEO+'pressure_tap_geometry_support.csv')})
    M.add('A', EXP, 'all cases (tested model)', 'tap coordinates: printed Table II', 'admitted',
          'Primary tap coordinates for scoring; three unsupported tap centres stay flagged and excluded',
          'Contour shift or pressure extrapolation to accommodate the unsupported centres',
          f'Table II, PDF pp. {tap_pages[0]}-{tap_pages[-1]} / printed pp. {inventory[tap_pages[0]]["printed_page"]}-{inventory[tap_pages[-1]]["printed_page"]}',
          cond_free, tap_text,
          'not reported (x/c printed to 5 decimals); printed-versus-plotted difference recorded as its own item',
          ['U13'], tap_records+[
              ev.quote(GEO+'FINDINGS.md', 'Keep the three taps explicitly flagged; no contour shift or pressure extrapolation is justified from them.'),
              ev.quote(GEO+'FINDINGS.md', 'Finite tap size, source faceting and coordinate conventions are possible contributors; these comparisons do not identify one.'),
              ev.srcline(REPORT_TXT, 'span at the staggered positions listed in table II. All the orifices were 0.51 mm (0.020 in.) in')],
          [GEO+'pressure_tap_geometry_support.csv'])

    checks12 = ev.rows(GEO+'pressure/tap_coordinate_checks.csv')
    checks7 = ev.rows(TRIP+'pressure_tap_checks.csv')
    fail12 = {c: sum(r['case'] == c and r['passes_2e_minus_5_check'] == '0' for r in checks12) for c in ('cl028', 'cl104')}
    max7 = max(float(r['maximum_plotted_table_x_difference']) for r in checks7)
    worst7 = max(checks7, key=lambda r: float(r['maximum_plotted_table_x_difference']))
    worst12 = max(checks12, key=lambda r: abs(float(r['plotted_minus_table_x'])))
    pages7 = sorted({int(r['pdf_page']) for r in free_pressure})
    pages12 = [c['page'] for c in ev.json(psrc)['cases']]
    fore_tap_page = sorted({r['pdf_page'] for r in ev.rows(GEO+'pressure_tap_geometry_support.csv') if r['element'] == 'fore'})
    max12 = max(ev.json(psrc)['cases'][i]['maximum_tap_x_error'] for i in (0, 1))
    plotted_evidence = [
        ev.val('maximum |plotted - printed| tap x, Figure 12(a)', psrc, 'cases', 0, 'maximum_tap_x_error', unit='c'),
        ev.val('maximum |plotted - printed| tap x, Figure 12(b)', psrc, 'cases', 1, 'maximum_tap_x_error', unit='c'),
        ev.val('exact 2e-5c table check passed, Figure 12(a)', psrc, 'cases', 0, 'exact_table_coordinate_check_pass'),
        ev.val('ordered nearest-tap association passed, Figure 12(a)', psrc, 'cases', 0, 'ordered_nearest_table_association_pass'),
        ev.derived('taps failing the 2e-5c check in Figures 12(a)/12(b), of 90 each', fail12, GEO+'pressure/tap_coordinate_checks.csv', 'count(passes_2e_minus_5_check=0) by case'),
        ev.derived('maximum |plotted - printed| tap x over Figure 7(a-d)', max7, TRIP+'pressure_tap_checks.csv', 'max(maximum_plotted_table_x_difference)', 'c'),
        ev.derived('surface groups of Figure 7(a-d) passing the 2e-5c check, of 16',
                   sum(r['exact_2e_minus_5_coordinate_check'] == 'True' for r in checks7), TRIP+'pressure_tap_checks.csv', 'count(exact_2e_minus_5_coordinate_check=True)'),
        ev.derived('largest plotted-printed difference in mm', round(max(max7, max12)*chord_mm, 4), psrc, 'max difference * reference_chord_m * 1000', 'mm'),
        ev.derived('surface and figure of the largest Figure 7 offset', [worst7['surface'], worst7['figure']], TRIP+'pressure_tap_checks.csv', 'argmax(maximum_plotted_table_x_difference)'),
        ev.derived('surface and case of the largest Figure 12 offset', [worst12['surface'], worst12['case']], GEO+'pressure/tap_coordinate_checks.csv', 'argmax(|plotted_minus_table_x|)'),
        ev.quote(GEO+'FINDINGS.md', 'maximum offset is 0.000260c (0.119 mm). Offsets show systematic trends, so they are not dismissed as PDF rounding.'),
        ev.quote(LBA+'FINDINGS.md', 'Fore-front coordinate sensitivity changes RMSE by at most 0.00223 and aft bubble recovery by at most 0.00118.'),
    ]
    M.add('A', EXP, 'all cases (tested model)', 'tap coordinates: plotted marker positions', 'conditional',
          'Sensitivity coordinate convention, reported beside the printed-table result; retained without fitted axes',
          'Replacing Table II; dismissing the offsets as PDF rounding',
          f'Figures 7(a-d), PDF pp. {pages7[0]}-{pages7[-1]}, and 12(a,b), PDF pp. {pages12[0]}-{pages12[-1]}', cond_free,
          f'offsets up to {max(max7, max12):.6f}c ({max(max7, max12)*chord_mm:.3f} mm); largest in Figure 7 on {worst7["surface"]} '
          f'({worst7["figure"]}), in Figure 12 on {worst12["surface"]} ({worst12["case"]})',
          'not reported; the printed-versus-plotted offset is an unresolved source inconsistency, not a measurement uncertainty',
          ['U13'], plotted_evidence, [GEO+'pressure_table_check_initial.log'])

    nose_rows = ev.rows(NOSE+'fore_lower_pressure.csv')
    nose_evidence = list(tap_records[2:7])
    for key in ('cl0056', 'cl0175', 'cl028'):
        nose = [r for r in nose_rows if r['case'] == key and float(r['x']) <= .06]
        peak = min(nose, key=lambda r: float(r['viscous_Cp']))
        nose_evidence += [
            ev.derived(f'MODEL (not measured) minimum viscous Cp on the fore-lower nose, {key}: [x/c, Cp]', [peak['x'], peak['viscous_Cp']],
                       NOSE+'fore_lower_pressure.csv', f'argmin viscous_Cp over x<=0.06 (case={key})'),
            ev.cell(f'measured Cp at first fore-lower tap x/c {lower_taps[0]}, {key}', NOSE+'pressure_at_taps.csv', 'experiment_Cp', case=key, coordinate='table', table_index=0),
            ev.cell(f'measured Cp at second fore-lower tap x/c {lower_taps[1]}, {key}', NOSE+'pressure_at_taps.csv', 'experiment_Cp', case=key, coordinate='table', table_index=1),
        ]
        require(float(lower_taps[0]) < float(peak['x']) < float(lower_taps[1]), 'Model nose peak no longer between the first two taps')
    nose_evidence += [
        ev.val('maximum nominal/cubic model Cp difference on the nose (model-only)', FRA+'summary.json', 'low_lift_comparison', 'maximum_nose_Cp_change'),
        ev.val('x/c of that difference', FRA+'summary.json', 'low_lift_comparison', 'x_at_maximum_nose_Cp_change'),
        ev.quote(FRA+'FINDINGS.md', 'The largest interpolated Cp difference over the full nose overlap is 0.298 at x/c = 0.000364, ahead of the first tap; it must not be presented as a measured pressure correction.'),
        ev.quote(LBA+'FINDINGS.md', 'The modeled nose peak is displaced relative to the sparse taps, so comparing the continuous computed minimum directly with the smallest measured marker would exaggerate what is known.'),
        ev.quote(PLAN, "The model's unresolved between-tap nose peak must remain visibly unmeasured."),
        ev.quote(PLAN, 'Finite tap size alone does not justify a fitted pressure-averaging correction.'),
    ]
    nose_row = M.add(
        'A', EXP, 'CL 0.056 / 0.175 / 0.284 (fore-lower nose)', 'fore-nose pressure peak between taps', 'missing',
        'none: compare model and experiment only at the supported taps',
        'A measured peak value; a fitted pressure-averaging correction; comparing the model minimum with the lowest marker',
        f'Table II fore-element taps, PDF p. {", ".join(fore_tap_page)}; Figures '+', '.join(t['figure'] for t in lb_targets if t['case'] != 'clm0072'), cond_free,
        f'no tap between x/c {lower_taps[0]} and {lower_taps[1]} on the fore lower surface; the archived model minimum lies inside that gap in all three states',
        'not reported', ['U06', 'U05'], nose_evidence,
        notes='Model values are shown only to locate the unmeasured interval; they are not measurements. No archived nose-pressure file exists for CL -0.072 in the fore-nose diagnostic.')

    M.add('A', EXP, 'all cases (data reduction)', 'lift / normal-force definition and applied tunnel corrections', 'admitted',
          "Use the report's stated reduction explicitly: normal force from surface pressure, lift using wake-survey drag; "
          'authors already applied standard boundary corrections',
          'A further fitted angle or boundary correction; explaining the measured-angle lift excess by the force definition',
          f'report printed p. {ev.json(PB+"source_manifest.json")["scope"]["definitions_printed_page"]} (definitions) and '
          f'p. {ev.json(PB+"source_manifest.json")["scope"]["methods_printed_page"]} (methods)', cond_free,
          'definition only; no additional measurement',
          "2010 report: not reported. Overlapping 2014 article: lift-coefficient uncertainty 0.005 'at higher lift', "
          'pitching-moment 0.002 [archived article text]',
          ['U26', 'U19'],
          [ev.val('experimental_reduction', PB+'summary.json', 'experimental_reduction'),
           ev.val('definitions_printed_page', PB+'source_manifest.json', 'scope', 'definitions_printed_page'),
           ev.val('methods_printed_page', PB+'source_manifest.json', 'scope', 'methods_printed_page'),
           ev.quote(PB+'FINDINGS.md', 'defines normal force from surface pressure and lift using wake-survey drag'),
           ev.quote(PB+'FINDINGS.md', 'CL_report = CN / cos(alpha) - CD_wake tan(alpha)'),
           ev.quote(PB+'FINDINGS.md', "Applying the report's lift formula to the model changes the measured-angle CL by only +0.000539 / +0.000123."),
           ev.quote(PB+'FINDINGS.md', 'The source already reports standard tunnel-boundary corrections; none was applied again here.'),
           ev.srcline(REPORT_TXT, 'Standard, low-speed, wind-tunnel boundary corrections (ref. 18) have been applied to'),
           ev.srcline(REPORT_TXT, 'the data. The wake-survey-probe total-pressure-tube displacement correction (ref. 13) has'),
           ev.quote(PLAN, 'do not apply another fitted angle correction')]+article_lift,
          [REPORT_PDF], strict=NA)

    figure7 = [p for p in pages_report if 'Figure 7.-' in p['caption_text']]
    extracted_pages = ev.json(TRIP+'source_manifest.json')['scope']['pressure_pages']
    other7 = [p['pdf_page'] for p in figure7 if p['pdf_page'] not in extracted_pages]
    index_k = next(i for i, p in enumerate(pages_report) if p['pdf_page'] == ev.json(psrc)['cases'][1]['snapshot_page'])
    M.add('A', EXP, f'other Figure 7 snapshots ({len(other7)} panels, PDF pp. {min(other7)}-{max(other7)})',
          'pressure (surface Cp at taps)', 'not_extracted',
          'none yet: only the caption values (alpha, cl, cd, cm) are archived; markers were not extracted',
          'Any pressure comparison at these lifts', f'Figure 7(e-x), PDF pp. {min(other7)}-{max(other7)}', cond_free,
          'captions only', 'not reported', ['U20'],
          [ev.derived('Figure 7 panels listed in the page inventory', len(figure7), TRIP+'report_page_inventory.json', "count(caption_text contains 'Figure 7.-')"),
           ev.val('pressure pages extracted', TRIP+'source_manifest.json', 'scope', 'pressure_pages'),
           ev.val('caption of Figure 7(k)', TRIP+'report_page_inventory.json', index_k, 'caption_text')],
          notes='The CL 1.041 measurement was extracted from its Figure 12(b) copy, not from Figure 7(k).')

    # ------------------------------------------------------------------ B: series-level MSES rows
    M.add('B', MSES10, 'all cases (calculation setup)', 'transition model, Ncrit and calculation inputs', 'conditional',
          'Identify the published curves as approximate-envelope e^N, Ncrit 9, free transition',
          'Assuming identical geometry, closures or constants to MultiFoil; treating the curves as an installed-MSES run',
          'report Figures 12-13; transition method identified in the 2014 article, section III.A', cond_mses,
          'settings as stated only', 'not reported', ['U03'],
          [ev.val('mses_ncrit', psrc, 'mses_ncrit'), ev.val('mses_transition', psrc, 'mses_transition'),
           ev.val('mses_executable_run', psrc, 'mses_executable_run'), ev.val('actual_geometry', psrc, 'actual_geometry'),
           ev.val('original_run_deck_recovered', quali, 'mses', 'original_run_deck_recovered'),
           ev.val('new_executable_run', quali, 'mses', 'new_executable_run'),
           ev.quote(REGISTER, 'MSES uses full e^N at Ncrit 9; its earlier report uses an approximate envelope method.')])
    mses_lift = ev.rows(DRAG+'mses_lift.csv')
    M.add('B', MSES10, 'Re 1e6 free-transition polar', 'lift curve (alpha versus CL)', 'conditional',
          'Initialisation seed for solver incidence only (increasing-lift branch)',
          'Validation reference for incidence; extrapolation below the first vertex',
          page13, cond_mses,
          f'{len(mses_lift)} published vertices; lowest CL {min(float(r["CL"]) for r in mses_lift):.5f}', 'not reported', ['U03'],
          [ev.derived('vertices', len(mses_lift), DRAG+'mses_lift.csv', 'row count'), fig13,
           ev.quote(LB_PROTOCOL, 'New CL 0.056 and 0.175 centers come from the previously qualified increasing-lift branch of the published MSES lift curve.'),
           ev.quote(LB_PROTOCOL, 'No extrapolated MSES value is a validation reference.')])
    last12 = pages_report[fig12_pages[-1]-1]
    require(last12['pdf_page'] == fig12_pages[-1], 'Page inventory order changed')
    M.add('B', MSES10, 'third Figure 12 panel (see archived caption)', 'pressure (published curve and experiment)', 'not_extracted',
          'none: the panel exists in the report but was not extracted',
          'Any comparison at this lift', f'Figure 12(c), PDF p. {last12["pdf_page"]} / printed p. {last12["printed_page"]}', cond_mses, 'caption only',
          'not reported', ['U20'],
          [fig12_captions[-1],
           ev.val('published_MSES_free_figures', quali, 'pressure_coverage', 'published_MSES_free_figures'),
           ev.derived('cases extracted from Figure 12', [c['figure'] for c in ev.json(psrc)['cases']], psrc, 'cases[].figure')])

    # ------------------------------------------------------------------ C: CL 1.041 control
    index = ev.find(psrc, ['cases'], name='cl104')
    item = ev.json(psrc)['cases'][index]
    control_case = f'CL {item["experimental_cl"]} (Figure {item["figure"]} copy)'
    markers = [r for r in ev.rows(GEO+'pressure/experiment.csv') if r['case'] == 'cl104']
    supported = [r for r in ev.rows(RCA+'pressure_at_taps.csv') if r['case'] == 'cl104' and r['coordinate'] == 'table'
                 and r['method'] == 'nominal' and r['included'] == 'True']
    where = f'Figure {item["figure"]}, PDF p. {item["page"]}; measurement identified by Figure 7(k), PDF p. {item["snapshot_page"]}'
    control_quotes = [ev.quote(PLAN, 'Regression controls after a proposed remedy'),
                      ev.quote(PLAN, 'They are already inspected, not blind validation')]
    row_id = M.add('C', EXP, control_case, 'pressure (surface Cp at taps)', 'admitted',
                   'Regression control after a proposed remedy (already inspected)',
                   'Blind validation; exact matched-geometry benchmark; bubble endpoints',
                   where, cond_free, tap_text, cp_uncertainty, ['U05', 'U13', 'U01', 'U20'],
                   [ev.derived('pressure markers extracted', len(markers), GEO+'pressure/experiment.csv', 'count(case=cl104)'),
                    ev.derived('repeated plotting endpoints among them', sum(r['endpoint_duplicate'] == '1' for r in markers), GEO+'pressure/experiment.csv', 'count(endpoint_duplicate=1)'),
                    ev.derived('supported independent taps scored', len(supported), RCA+'pressure_at_taps.csv', 'count(case=cl104,coordinate=table,method=nominal,included=True)'),
                    ev.val('experimental alpha', psrc, 'cases', index, 'experimental_alpha_deg', unit='deg'),
                    ev.val('maximum |plotted - printed| tap x', psrc, 'cases', index, 'maximum_tap_x_error', unit='c'),
                    ev.val('maximum Cp difference between repeated endpoints', psrc, 'cases', index, 'duplicate_endpoint_max_cp_error'),
                    ev.val('pressure_uncertainty', psrc, 'pressure_uncertainty'),
                    ev.quote(RCA+'FINDINGS.md', 'There is no recovered case-specific Cp uncertainty budget, so these are descriptive errors, not statistical significance tests.')]
                   + control_quotes, [GEO+'pressure/experiment.csv', RCA+'excluded_taps.csv'])
    mark('cl104', 'experiment_pressure', 'exists (Figure 12(b) copy)', f'{len(markers)} markers; {len(supported)} supported taps',
         [GEO+'pressure/experiment.csv'], where, row_id)
    appendix = dict(reynolds='1000000', configuration='free', CL='1.041')
    arow = ev.row(TRIP+'experimental_tables.csv', **appendix)
    where = (f'Appendix table PDF p. {arow["pdf_page"]} / printed p. {arow["printed_page"]}; {page13}; '
             f'caption of Figure 7(k), PDF p. {item["snapshot_page"]}')
    row_id = M.add('C', EXP, control_case, 'drag (wake-survey section CD)', 'admitted',
                   'Regression control after a proposed remedy (already inspected); matched-lift drag reference',
                   'Blind validation; per-element drag', where, cond_free,
                   f'appendix row printed as "{arow["raw_tokens"]}" (alpha, CL, CD, CM)', drag_uncertainty, ['U19', 'U23'],
                   [ev.cell('appendix row (alpha CL CD CM as printed)', TRIP+'experimental_tables.csv', 'raw_tokens', **appendix),
                    ev.cell('caption CD', DRAG+'target_references.csv', 'caption_CD', case='cl104'),
                    ev.cell('experiment CD interpolated from Figure 13(c) markers', DRAG+'target_references.csv', 'experiment_CD_interpolated', case='cl104'),
                    ev.cell('experiment CD from the appendix table at CL 1.041', TRIP+'matched_lift_targets.csv', 'experiment_free_CD', reynolds=1000000, CL='1.041'),
                    ev.val('caption of Figure 7(k)', TRIP+'report_page_inventory.json', index_k, 'caption_text'), fig13]
                   + article_drag+control_quotes, [DRAG+'experiment_drag.csv'])
    mark('cl104', 'experiment_drag', 'exists', f'CD {arow["CD"]} (appendix row at CL {arow["CL"]})', [TRIP+'experimental_tables.csv', DRAG+'target_references.csv'], where, row_id)
    vertices = sum(r['case'] == 'cl104' for r in ev.rows(GEO+'pressure/mses_2010_envelope.csv'))
    where = f'Figure {item["figure"]}, PDF p. {item["page"]}'
    row_id = M.add('C', MSES10, control_case, 'pressure (published curve)', 'conditional',
                   'Regression control: published-MSES pressure at nominal CL 1.04 beside the CL 1.041 measurement',
                   'Experiment; an exact-lift match; the 2014 full-e^N series', where, cond_mses,
                   f'{vertices} published curve vertices over four surface branches', 'not reported', ['U03', 'U09'],
                   [ev.val('nominal theoretical CL of the published curve', psrc, 'cases', index, 'nominal_cl'),
                    ev.val('associated experimental CL', psrc, 'cases', index, 'experimental_cl'),
                    ev.derived('theory curve vertices', vertices, GEO+'pressure/mses_2010_envelope.csv', 'count(case=cl104)')]+control_quotes)
    mark('cl104', 'mses2010_pressure', 'exists (nominal CL 1.04)', f'{vertices} vertices', [GEO+'pressure/mses_2010_envelope.csv'], where, row_id)
    low, high = neighbours(mses_free, 1.041)
    row_id = M.add('C', MSES10, control_case, 'drag (published curve)', 'admitted',
                   'Regression control: published-MSES drag at the same lift', 'Experiment; the 2014 full-e^N series',
                   page13, cond_mses,
                   f'piecewise-linear published curve; nearest vertices at CL {low:.4f} and {high:.4f}', 'not reported', ['U03'],
                   [ev.cell('published MSES CD interpolated at CL 1.041', DRAG+'target_references.csv', 'mses_CD_interpolated', case='cl104'),
                    ev.cell('MSES minus experiment (Figure 13(c) experimental markers)', DRAG+'target_references.csv', 'mses_minus_experiment_counts', unit='counts', case='cl104'),
                    ev.cell('MSES minus experiment (appendix table)', TRIP+'matched_lift_targets.csv', 'mses_minus_experiment_free_counts', unit='counts', reynolds=1000000, CL='1.041'),
                    fig13]
                   + control_quotes, [TRIP+'published_mses_drag.csv', DRAG+'mses_drag.csv'],
                   notes='Two archived bias values differ only through the experimental reference convention.')
    mark('cl104', 'mses2010_drag', 'exists', 'interpolated at CL 1.041', [DRAG+'target_references.csv', TRIP+'published_mses_drag.csv'], page13, row_id)
    row_id = M.add('C', EXP, control_case, 'transition / separation / reattachment location', 'missing',
                   'none as a quantitative reference', 'Any numerical transition or reattachment target',
                   f'Figure 7(k), PDF p. {item["snapshot_page"]}; report discussion printed pp. 9-11 (PDF pp. 11-13)', cond_free, 'no measured coordinates', 'not reported',
                   ['U04', 'U24'],
                   [ev.val('measured_transition_endpoints_available', psrc, 'measured_transition_endpoints_available'),
                    ev.val('experimental_transition_endpoints_recovered', GEO+'verification.json', 'scientific_limits', 'experimental_transition_endpoints_recovered'),
                    ev.quote(GEO+'FINDINGS.md', 'It does not provide exact experimental transition or reattachment coordinates'), discussion])
    mark('cl104', 'experiment_transition_reattachment_location', 'missing', 'no measured coordinates', [psrc], 'none', row_id)
    row_id = M.add('C', MSES10, control_case, 'transition / separation / reattachment location', 'missing',
                   'none: no MSES transition or bubble location is recorded in the existing files; qualitative author interpretation only (see quotation)',
                   'Inferring an MSES transition location from a Cp kink', 'none recorded', cond_mses, 'no data', 'not reported',
                   ['U09', 'U03'],
                   [ev.derived('columns of the archived MSES pressure file', list(ev.rows(GEO+'pressure/mses_2010_envelope.csv')[0].keys()),
                               GEO+'pressure/mses_2010_envelope.csv', 'header (no transition field)'), bubble_quote])
    mark('cl104', 'mses_transition', 'missing', 'not recorded in any archived file', [GEO+'pressure/mses_2010_envelope.csv'], 'none', row_id)

    # ------------------------------------------------------------------ C: archived MultiFoil states
    case_of = {'-0.072': 'clm0072', '0.056': 'cl0056', '0.175': 'cl0175', '0.284': 'cl028', '1.041': 'cl104'}
    states = {}
    explicit = {}   # archived state -> later jobs with a recorded Jacobian mode that reproduce it exactly
    for analysis, study in STUDIES:
        recorded = ev.json(analysis+'summary.json').get('controls')
        for job, control in (recorded.get('jobs', {}) if isinstance(recorded, dict) else {}).items():
            mode = ev.json(study+job+'/result.json').get('amplification_jacobian')
            if mode and control['passed']:
                explicit.setdefault(control['reference'].rstrip('/'), []).append((study+job, mode, analysis+'summary.json'))
    for analysis, study in STUDIES:
        s = ev.json(analysis+'summary.json')
        audit = (ev.val('full_state_audit_pass', analysis+'summary.json', 'full_state_audit_pass') if 'full_state_audit_pass' in s
                 else ev.val('selected_states_freshly_audited', analysis+'summary.json', 'selected_states_freshly_audited'))
        controls = {}
        if (ROOT/analysis/'controls.json').is_file():
            c = ev.json(analysis+'controls.json')
            controls[c.get('repeat') or c['control']] = (analysis+'controls.json', c['reference'])
        for job, c in s.get('controls', {}).get('jobs', {}).items() if isinstance(s.get('controls'), dict) else []:
            controls[job] = (analysis+'summary.json', c['reference'])
        for t in ev.rows(analysis+'targets.csv'):
            job = t['job']
            rel = study+job+'/result.json'
            result = ev.json(rel)
            feasible = result['feasible']
            basis = result.get('comparison_basis', 'fixed_lift')
            target = result.get('target_cl', result.get('reference_experimental_cl'))
            key = case_of[show(target)]
            mode = result.get('amplification_jacobian', 'not recorded in result.json')
            evidence = [
                ev.val('feasible', rel, 'feasible'), ev.val('CL', rel, 'CL'), ev.val('CD', rel, 'CD'),
                ev.val('alpha', rel, 'alpha_deg', unit='deg'), ev.val('scenario', rel, 'scenario'),
                ev.val('ncrit', rel, 'ncrit'), ev.val('wake_length', rel, 'wake_length', unit='c'),
                ev.val('amplification_jacobian', rel, 'amplification_jacobian') if 'amplification_jacobian' in result else
                ev.derived('amplification_jacobian', mode, rel, 'key absent'),
                ev.val('selected_path', rel, 'selected_path') if feasible else ev.val('error', rel, 'error'),
                audit,
                ev.val('numerical_convergence_established', analysis+'summary.json', 'numerical_convergence_established'),
            ]
            for field in ('recomputed_residual', 'residual'):
                if t.get(field):
                    evidence.append(ev.cell(field, analysis+'targets.csv', field, job=job))
                    break
            note = ''
            for path, later_mode, source in explicit.get(study+job, []):
                evidence.append(ev.derived('reproduced exactly by a later job with a recorded Jacobian mode', f'{path} ({later_mode})', source, 'controls/jobs reference'))
            if job in controls:
                evidence.append(ev.derived('exact-repeat control of', controls[job][1], controls[job][0], 'controls reference'))
                note = 'Exact repeat of an earlier archived state. '
            cite = [rel]
            if feasible:
                for name in ('fore_state.npz', 'aft_state.npz'):
                    require((ROOT/study/job/name).is_file(), 'Missing saved state: '+study+job+'/'+name)
                    cite.append(study+job+'/'+name)
            label = (f'CL {show(target)} target, {result["scenario"]}' if basis == 'fixed_lift' else
                     f'fixed alpha {result["target_alpha"]} deg (experimental CL {show(target)}), {result["scenario"]}')
            setup = study_of[analysis]
            condition = (f'Re {setup["reynolds"]:,.0f}; Mach {setup["mach"]}; Ncrit {result["ncrit"]}; panels '
                         f'{setup["panels"][0]}/{setup["panels"][1]} (fore/aft); {result["wake_length"]}c wake; '
                         f'comparison basis {basis}; amplification_jacobian {mode}')
            evidence += [ev.val('reynolds', study+'study.json', 'reynolds'), ev.val('mach', study+'study.json', 'mach'),
                         ev.val('panels', study+'study.json', 'panels'), ev.val('cl_tolerance', study+'study.json', 'cl_tolerance'),
                         ev.val('residual_tolerance', study+'study.json', 'residual_tolerance')]
            if not feasible:
                status, use, not_for = 'missing', 'none: failed target; forces are null and no nearby evaluation substitutes for it', 'Any finite force value or stall inference'
                note += 'Failed target preserved; '+result['error']
                unknown = ['U18', 'U17']
            elif basis != 'fixed_lift':
                status, use = 'conditional', 'Regression control on the identical-incidence (mechanism) basis only'
                not_for = 'Matched-lift drag validation; blind validation'
                unknown = ['U17', 'U01']
            else:
                status, use = 'admitted', 'Regression control after a proposed remedy (already inspected, fully audited saved state)'
                not_for = 'Blind validation; certification of small design gains; proof of grid convergence or uniqueness'
                unknown = ['U17', 'U01']
            row_id = M.add('C', MODEL, f'{job} ({label})',
                           'archived audited solver state' if feasible else 'archived solver target (failed; no state)', status, use, not_for,
                           'not a source document; archived run directory '+study+job, condition,
                           f'{setup["panels"][0]} fore / {setup["panels"][1]} aft panels; '
                           + (f'fixed-lift tolerance {setup["cl_tolerance"]} CL' if basis == 'fixed_lift' else 'requested incidence imposed (lift not matched)')
                           + f'; residual tolerance {setup["residual_tolerance"]}',
                           'not applicable (numerical state); grid convergence and uniqueness not established', unknown,
                           evidence+control_quotes, cite, notes=note.strip())
            states.setdefault(key, []).append(dict(job=job, path=study+job, basis=basis, feasible=feasible,
                                                   mode=mode, row_id=row_id, files=cite, analysis=analysis,
                                                   scenario=result['scenario'], repeat=job in controls))
    for key, items in states.items():
        good = [i for i in items if i['feasible']]
        lift = [i for i in good if i['basis'] == 'fixed_lift']
        repeats, variants = sum(i['repeat'] for i in lift), sum(i['scenario'] != 'nominal' for i in lift)
        text = f'{len(good)} audited state director{"y" if len(good) == 1 else "ies"}: {len(lift)} matched-lift'
        if repeats or variants:
            text += f' ({repeats} exact-repeat controls, {variants} other-geometry scenarios)'
        if len(good) > len(lift):
            text += f', {len(good)-len(lift)} fixed-angle'
        if len(items) > len(good):
            text += f'; {len(items)-len(good)} failed target without a state'
        detail = text+' | '+'; '.join(
            f'{i["path"]} [{i["basis"]}; scenario {i["scenario"]}; jacobian {"not recorded" if i["mode"].startswith("not") else i["mode"]}'
            + ('; exact-repeat control' if i['repeat'] else '')+('' if i['feasible'] else '; FAILED TARGET, no state')
            + f'; audit {i["analysis"]}]' for i in items)
        mark(key, 'multifoil_archived_audited_state', 'exists' if good else 'missing', detail,
             [f for i in items for f in i['files']], 'not a source document', ','.join(i['row_id'] for i in items))

    # ------------------------------------------------------------------ D: 2014 material
    article_entry = next(i for i, e in enumerate(ev.json(REVIEW+'catalog.json')['entries']) if e['id'] == 'N02')
    manifest_entry = next(i for i, e in enumerate(ev.json(GEO+'source_manifest.json')['existing_sources'])
                          if e['path'] == ARTICLE_PDF)
    cond14 = 'Low-speed experiments Re 0.5-1.5 million; MSES full e^N, Ncrit 9; per-figure Mach and transition condition not recorded in the existing files'
    basis14 = [ev.quote(REGISTER, 'Low-speed experiments span Re 0.5–1.5 million; the paper compares lift, drag and pressure with MSES and OVERFLOW.'),
               ev.quote(REGISTER, 'MSES uses full e^N at Ncrit 9; its earlier report uses an approximate envelope method.')]
    version14 = ev.val('version', GEO+'source_manifest.json', 'existing_sources', manifest_entry, 'version')
    separation = [ev.quote(PLAN, 'Separate, condition-qualified supporting series'),
                  ev.quote(PLAN, 'They must not be merged with the 2010 envelope series or counted as independent experiments where overlapping')]
    M.add('D', MSES14, 'source document', 'archived article (PDF and text rendition)', 'admitted',
          'Bibliographic source for the transition-method distinction and the stated measurement uncertainty',
          'A digitized data series (none has been extracted)', version14['value'],
          cond14, 'whole document', 'not reported for the document as a whole; the stated force-coefficient uncertainty is in the experimental row of this group', [],
          [ev.val('recorded sha256 of the archived PDF', GEO+'source_manifest.json', 'existing_sources', manifest_entry, 'sha256'),
           version14, ev.val('register id', REVIEW+'catalog.json', 'entries', article_entry, 'id')]+basis14,
          [ARTICLE_PDF, ARTICLE_TXT], strict=NA,
          notes='The text rendition is not pinned by an earlier manifest; its hash is recorded here for the first time.')
    M.add('D', MSES14, 'all conditions in the article', 'lift and drag curves (full e^N)', 'not_extracted',
          'Separate, condition-qualified supporting series once digitized and condition-qualified; nothing is usable numerically yet',
          'Merging with, or substituting for, the 2010 envelope-MSES series', 'article figures (not digitized)', cond14,
          'no values archived', 'not reported', ['U15', 'U03'],
          basis14+[ev.quote(DRAG+'FINDINGS.md', 'Do not substitute the 2014 full-e^N formulation.')]+separation, [ARTICLE_PDF])
    M.add('D', MSES14, 'Re 0.7/1.5 million, CL 0.493/0.518 (article Figures 11-12)', 'pressure (full e^N curves and measurements)', 'not_extracted',
          'Separate supporting pressure reference at other Reynolds numbers once digitized; nothing is usable numerically yet',
          'Merging into the 2010 Re 1e6 envelope series; substituting design Mach numbers', 'article Figures 11-12 (not digitized)',
          'Re 0.7 and 1.5 million; CL 0.493 and 0.518; Mach not recorded in the existing files', 'no values archived', 'not reported',
          ['U15'],
          [ev.quote(GEO+'FINDINGS.md', 'The 2014 Figures 11–12 offer another pressure reference at Re=0.7/1.5 million and CL=0.493/0.518; their full-e^N curves have not been digitized here.'),
           ev.srcline(ARTICLE_TXT, 'is 0.493, whereas for R')]+separation, [ARTICLE_PDF])
    M.add('D', 'overflow_2014 (Coder, Maughmer and Somers; not digitized)', 'all conditions in the article',
          'lift, drag and pressure (OVERFLOW)', 'not_extracted',
          'none yet: mentioned in the register only', 'Any numerical comparison', 'article figures (not digitized)', cond14,
          'no values archived', 'not reported', ['U16'],
          basis14+separation, [ARTICLE_PDF])
    M.add('D', EXP14, 'all conditions in the article', 'experimental lift, drag, moment and their stated uncertainty', 'conditional',
          'Supporting context and the only archived statement of force-coefficient uncertainty for this campaign',
          'Counting as an experiment independent of the 2010 report; a Cp uncertainty',
          'article section on experimental methods (text rendition lines cited)', cond14,
          'no measurement values archived from the article',
          'drag 0.00005 (0.5 drag counts) in the low-drag range, increasing to 0.00015 (1.5 drag counts) approaching stall; '
          "lift 0.005 'at higher lift'; pitching moment 0.002 [archived article text; relation/sign glyphs unreadable in the rendition]",
          ['U05'],
          article_drag+article_lift+[
                        ev.srcline(ARTICLE_TXT, '0.13 mm (0.005 in.) of the prescribed shape. The model as installed in'),
                        ev.quote(GEO+'FINDINGS.md', 'These are overlapping experiments, not independent 2010 and 2014 validation datasets.'),
                        ev.quote(REGISTER, "This report overlaps N02's experiment rather than adding another independent validation dataset.")]
          + basis14+separation, [ARTICLE_PDF])
    M.add('D', EXP14, 'lower drag-bucket boundary', 'authors\' interpretation of the drag dip (fore-wake / aft-bubble interaction)', 'conditional',
          'Qualitative hypothesis to keep in view; an interpretation',
          'Measured transported turbulence; a closure-fitting target', 'article section VI', cond14, 'qualitative',
          'not reported', ['U10', 'U24'],
          [ev.quote(GEO+'FINDINGS.md', 'The proposed wake mechanism remains an interpretation, not measured transported turbulence.'),
           ev.val('note', WAKE+'source_trace.json', 'notes', 1)]+basis14, [ARTICLE_PDF], strict=NA)
    M.add('D', 'twiss_2019_thesis (separate campaign; unqualified lead)', 'not qualified', 'S414 high-lift measurements and oil-flow observations', 'not_admitted',
          'none: discovery lead only, not acquired', 'Any scored comparison', 'not acquired',
          'conditions and geometry not qualified', 'no values archived', 'not reported', ['U22'],
          [ev.val('included_in_scored_data', TRIP+'research_search.json', 'additional_lead', 'included_in_scored_data'),
           ev.val('local_download_status', TRIP+'research_search.json', 'additional_lead', 'local_download_status'),
           ev.val('qualification', TRIP+'research_search.json', 'additional_lead', 'qualification'),
           ev.quote(TRIP+'FINDINGS.md', 'That is a separate experimental campaign and an unqualified follow-up lead; none of its measurements are mixed into this comparison.')],
          strict=NA)

    # ------------------------------------------------------------------ A/B/E: complete polars and trip series
    admission = ev.json(quali)['admission']
    trip_quote = ev.quote(TRIP+'FINDINGS.md', 'Admit the data as qualified global-force context; reject admission as a clean causal wake-transfer benchmark or a target for fitting a closure.')
    secondary_quote = ev.quote(TRIP+'FINDINGS.md', 'the secondary condition is Re 700,000, Mach 0.07, where the author describes a stronger low-lift drag dip')
    labels = {'free': 'transition free', 'fore_fixed': 'fore element tripped, aft free', 'both_fixed': 'fore and aft elements tripped'}
    for reynolds, mach in ((primary['reynolds'], primary['mach']), (secondary['reynolds'], secondary['mach'])):
        for configuration in ('free', 'fore_fixed', 'both_fixed'):
            table = [r for r in experiments if r['reynolds'] == str(reynolds) and r['configuration'] == configuration]
            curve = [r for r in mses_curves if r['reynolds'] == str(reynolds) and r['configuration'] == configuration]
            lead = reynolds == primary['reynolds'] and configuration == 'free'
            condition = f'Re {reynolds:,}; Mach {mach}; {labels[configuration]}; reference chord {primary["reference_chord_m"]} m'
            case = f'Re {reynolds:,}, {labels[configuration]} (complete polar)'
            tripped = configuration != 'free'
            M.add('A' if lead else 'E', EXP, case, 'lift/drag/moment table (printed appendix)',
                  'admitted' if not tripped or admission['matched_global_force_context'] else 'conditional',
                  ('Source table for the matched-lift drag references of the four cases and the CL 1.041 control' if lead else
                   'Qualified global-force context at matched lift' if tripped else 'Secondary-condition context (stronger low-lift drag dip described by the author)'),
                  ('Extrapolation; filling missing drag; post-stall rows in interpolation' if not tripped else
                   'Clean wake-disturbance intervention; matched tripped aft-pressure benchmark; closure fitting'),
                  f'Appendix table PDF p. {table[0]["pdf_page"]} / printed p. {table[0]["printed_page"]}', condition,
                  f'{len(table)} printed rows ({sum(r["prestall_branch"] == "True" for r in table)} on the pre-stall branch; '
                  f'{sum(r["drag_available"] != "True" for r in table)} without drag)',
                  drag_uncertainty, ['U19', 'U23']+(['U12', 'U11'] if tripped else []),
                  [ev.derived('printed rows', len(table), TRIP+'experimental_tables.csv', f'count(reynolds={reynolds},configuration={configuration})'),
                   ev.derived('rows on the pre-stall branch', sum(r['prestall_branch'] == 'True' for r in table), TRIP+'experimental_tables.csv', 'count(prestall_branch=True)'),
                   ev.derived('rows without drag', sum(r['drag_available'] != 'True' for r in table), TRIP+'experimental_tables.csv', 'count(drag_available!=True)'),
                   ev.val('maximum plotted-marker vs appendix CL difference', tsum, 'maximum_plot_table_CL_difference'),
                   ev.val('maximum plotted-marker vs appendix CD difference', tsum, 'maximum_plot_table_CD_difference')]
                  + ([trip_quote, ev.val('clean_wake_disturbance_intervention', quali, 'admission', 'clean_wake_disturbance_intervention')] if tripped else [])
                  + ([secondary_quote] if not tripped and not lead else []) + article_drag, [TRIP+'experimental_plot_markers.csv', TRIP+'table_plot_crosschecks.csv'])
            cls = [float(r['CL']) for r in curve]
            M.add('B' if lead else 'E', MSES10, case, 'drag curve (published)',
                  'admitted' if not tripped or admission['matched_global_force_context'] else 'conditional',
                  ('Source curve for the published-MSES drag references' if lead else
                   'Qualified global-force context at matched lift' if tripped else 'Secondary-condition context'),
                  ('Extrapolation below the first vertex; the 2014 full-e^N series' if not tripped else
                   'Equivalence of physical trips and theoretical transition forcing; closure fitting'),
                  f'Figure {"13" if configuration == "free" else "14" if configuration == "fore_fixed" else "15"}, PDF p. {curve[0]["pdf_page"]}',
                  condition+f'; Ncrit {ev.json(quali)["mses"]["ncrit"]} for free transition; {ev.json(quali)["mses"]["transition_method"]}',
                  f'{len(curve)} published vertices, CL {min(cls):.5f} to {max(cls):.5f}', 'not reported',
                  ['U03']+(['U12'] if tripped else []),
                  [ev.derived('published vertices', len(curve), TRIP+'published_mses_drag.csv', f'count(reynolds={reynolds},configuration={configuration})'),
                   ev.derived('CL range of the published curve', [min(cls), max(cls)], TRIP+'published_mses_drag.csv', 'min/max(CL)'),
                   ev.val('actual_trip_displacement_equivalence_established', quali, 'mses', 'actual_trip_displacement_equivalence_established')]
                  + ([trip_quote] if tripped else []))
    targets = TRIP+'matched_lift_targets.csv'
    pages_table = sorted({r['pdf_page'] for r in experiments if r['reynolds'] == str(primary['reynolds'])})
    pages_curve = sorted({r['pdf_page'] for r in mses_curves if r['reynolds'] == str(primary['reynolds'])})
    for reference, prefix, label, pages, uncertainty in (
            (EXP, 'experiment', 'experiment', f'Appendix tables PDF pp. {", ".join(pages_table)}', drag_uncertainty),
            (MSES10, 'mses', 'published MSES', f'Figures 13(c), 14(c), 15(c), PDF pp. {", ".join(pages_curve)}', 'not reported')):
        M.add('E', reference, f'Re {primary["reynolds"]:,}, CL 0.284 and 1.041', 'matched-lift trip drag contrasts (tripped minus free)', 'conditional',
              'Total-force context: the fore-trip drag penalty at matched lift, to be read beside the other series, never merged with it',
              'Clean causal wake-transfer benchmark; isolated bubble-drag measurement; closure fitting',
              pages, f'Re {primary["reynolds"]:,}; Mach {primary["mach"]}; matched lift; incidences differ between configurations',
              'piecewise-linear interpolation at the same CL', uncertainty,
              ['U12', 'U11', 'U10']+(['U03'] if prefix == 'mses' else ['U19']),
              [ev.cell(f'{label} fore trips minus free, CL 0.284', targets, prefix+'_fore_minus_free_counts', unit='counts', reynolds=1000000, CL='0.284'),
               ev.cell(f'{label} both trips minus fore trips, CL 0.284', targets, prefix+'_both_minus_fore_counts', unit='counts', reynolds=1000000, CL='0.284'),
               ev.cell(f'{label} fore trips minus free, CL 1.041', targets, prefix+'_fore_minus_free_counts', unit='counts', reynolds=1000000, CL='1.041'),
               ev.cell(f'{label} both trips minus fore trips, CL 1.041', targets, prefix+'_both_minus_fore_counts', unit='counts', reynolds=1000000, CL='1.041'),
               ev.val('matched_global_force_context', quali, 'admission', 'matched_global_force_context'),
               ev.val('clean_wake_disturbance_intervention', quali, 'admission', 'clean_wake_disturbance_intervention'),
               ev.val('closure_fit_justified', quali, 'admission', 'closure_fit_justified'), trip_quote])
    M.add('E', EXP, 'tripped configurations', 'trip locations and sizes (Table III)', 'conditional',
          'Nominal trip description; fore positions use the total chord, aft positions the aft-element chord',
          'Proof of transition onset at the trip; a controlled disturbance intervention',
          f'Table III, PDF p. {ev.json(quali)["trip_table_pdf_page"]}; trip confounds PDF p. {ev.json(quali)["trip_confounds_pdf_page"]}',
          'Re 0.7 and 1.0 million tripped configurations', 'nominal positions and sizes as printed', 'not reported', ['U12'],
          [ev.val('trips', quali, 'trips'), ev.val('confounds', quali, 'confounds')], [TRIP+'source_page_21.png'], strict=NA)
    M.add('E', EXP, 'tripped configurations', 'pressure with trips (aft or fore Cp)', 'missing',
          'none: no tripped pressure distributions were found in the report', 'A matched tripped aft-pressure benchmark',
          f'none (page inventory of all {len(pages_report)} report pages)', 'Re 0.7 and 1.0 million tripped configurations', 'no data', 'not reported', ['U11'],
          [ev.val('fore_only_fixed_Cp_found_in_report', quali, 'pressure_coverage', 'fore_only_fixed_Cp_found_in_report'),
           ev.val('both_fixed_Cp_found_in_report', quali, 'pressure_coverage', 'both_fixed_Cp_found_in_report'),
           ev.val('matched_tripped_aft_pressure_benchmark', quali, 'admission', 'matched_tripped_aft_pressure_benchmark'),
           ev.val('outcome', TRIP+'research_search.json', 'outcome')], [TRIP+'report_page_inventory.json'])
    M.add('E', EXP, 'all configurations', 'wake turbulence / velocity profiles and oil-flow bubble extents', 'missing',
          'none: qualitative oil-flow description only', 'Any quantitative wake-disturbance or bubble-length reference',
          'report printed pp. 9-11 (PDF pp. 11-13), qualitative', 'all archived conditions', 'no data', 'not reported', ['U10', 'U24'],
          [ev.val('wake_turbulence_profiles_found', quali, 'pressure_coverage', 'wake_turbulence_profiles_found'),
           ev.val('oil_flow', quali, 'pressure_coverage', 'oil_flow'), discussion], strict=NA)
    for i, excluded in enumerate(ev.json(quali)['excluded_conditions']):
        mismatch = 'not' in excluded['reason'].lower() and 'matched' in excluded['reason'].lower()
        tag = f'R = {excluded["reynolds"]/1e6:.2f} × 106'
        hits = [(j, q) for j, q in enumerate(pages_report) if 'Figure 13' in q['caption_text'] and tag in q['caption_text']]
        require(len(hits) == 1, 'Expected one Figure 13 caption for '+tag)
        M.add('E', 'report_2010 (experimental and published-MSES panels; nothing extracted)', f'Re {excluded["reynolds"]:,}', 'polars at an excluded Reynolds number',
              'not_admitted' if mismatch else 'not_extracted',
              'none: '+excluded['reason'], 'Any matched intervention comparison at this Reynolds number',
              f'free-transition panel of Figure 13, PDF p. {hits[0][1]["pdf_page"]} (tripped panels and appendix tables: see page inventory)',
              f'Re {excluded["reynolds"]:,}; Mach as in the archived caption', 'not extracted', 'not reported', ['U27'],
              [ev.val('reynolds', quali, 'excluded_conditions', i, 'reynolds'), ev.val('reason', quali, 'excluded_conditions', i, 'reason'),
               ev.val('caption of the free-transition panel', TRIP+'report_page_inventory.json', hits[0][0], 'caption_text')])

    # ------------------------------------------------------------------ uncertainty classes
    def klass(name, item, record, quantifies, not_as, status='quantified'):
        # For an unknown item the record is the archived statement or flag showing that it is unknown.
        unknown = status == 'unknown'
        classes.append(dict(uncertainty_class=name, item=item,
                            value='unknown' if unknown else show(record['value']), unit='' if unknown else record['unit'],
                            status=status, what_the_files_quantify=quantifies, not_to_be_read_as=not_as,
                            basis_for_unknown=f'{record["name"]}={show(record["value"])}' if unknown else '',
                            kind=record['kind'], source_path=record['source'], locator=record['locator'],
                            source_sha256=ev.digest(record['source'])))

    c1, c2, c3, c4, c5, c6 = ('extraction repeatability', 'reported manufacturing tolerance', 'reconstruction spread',
                              'actual as-built uncertainty', 'printed-table versus plotted tap coordinates', 'unmeasured between-tap nose peak')
    not_physical = 'measurement uncertainty, manufactured accuracy or model error'
    klass(c1, f'outline: max fore coordinate difference over {len(outlines)} repeated report drawings', geometry_evidence[10],
          'agreement of repeated vector drawings of the same published outline', not_physical)
    klass(c1, f'outline: max aft vertex-to-segment distance over {len(outlines)} repeated report drawings', geometry_evidence[11],
          'agreement of repeated vector drawings of the same published outline', not_physical)
    klass(c1, 'outline: Figure 1 vector geometry identical in two independent PDF copies', ev.val('source_copies_geometry_identical', summary, 'source_copies_geometry_identical'),
          'identity of the drawing objects in the mirror and archive copies', not_physical)
    klass(c1, 'pressure: max Cp difference between the Figure 7(d) and 12(a) copies (strict 3e-5 check failed)', ev.val('max Cp difference', tsum, 'maximum_repeated_pressure_Cp_difference'),
          'difference between two plotted copies of one measurement', not_physical)
    klass(c1, 'pressure: max plotted-x difference between the Figure 7(d) and 12(a) copies (strict 1e-5c check failed)', ev.val('max x difference', tsum, 'maximum_repeated_pressure_x_difference', unit='c'),
          'difference between two plotted copies of one measurement', not_physical)
    klass(c1, 'pressure: max Cp difference between repeated endpoint markers, Figure 12(a)', ev.val('duplicate endpoint Cp', psrc, 'cases', 0, 'duplicate_endpoint_max_cp_error'),
          'difference between markers plotted twice in one figure', not_physical)
    klass(c1, 'pressure: max Cp difference between repeated endpoint markers, Figure 12(b)', ev.val('duplicate endpoint Cp', psrc, 'cases', 1, 'duplicate_endpoint_max_cp_error'),
          'difference between markers plotted twice in one figure', not_physical)
    klass(c1, f'drag: max CD difference between {ev.json(tsum)["independent_plot_table_checks"]} plotted markers and their appendix rows', ev.val('max CD difference', tsum, 'maximum_plot_table_CD_difference'),
          'agreement of plotted and printed values', not_physical)
    klass(c1, 'lift: max CL difference between plotted markers and appendix rows', ev.val('max CL difference', tsum, 'maximum_plot_table_CL_difference'),
          'agreement of plotted and printed values', not_physical)
    klass(c1, 'lift: CL offset between corresponding markers of the drag and lift plots of Figure 13(c)', ev.val('max CL offset', DRAG+'source.json', 'maximum_drag_lift_CL_offset'),
          'plotting offset between two panels of one figure', not_physical)
    klass(c2, 'model contour: measured contour within this distance of the prescribed shape', ev.val('reported_model_contour_tolerance_mm', summary, 'reported_model_contour_tolerance_mm', unit='mm'),
          "the authors' stated contour tolerance of the machined model", 'as-built contour error distribution, rigging uncertainty or a confidence bound')
    klass(c2, 'same tolerance in chord units as used for the aft z-shift sensitivity', ev.val('aft_shift_c', STUDIES[0][1]+'study.json', 'aft_shift_c', unit='c'),
          f'unit conversion with the {study_of[RCA]["model_chord_mm"]} mm model chord', 'measured rigging uncertainty')
    klass(c2, 'pressure orifice diameter', ev.val('reported_tap_diameter_mm', summary, 'reported_tap_diameter_mm', unit='mm'),
          "the authors' stated orifice diameter", 'a justification for a fitted pressure-averaging correction')
    klass(c3, 'max interpolation displacement PCHIP versus natural cubic, fore / aft', ev.quote(DRAG+'FINDINGS.md', 'Maximum interpolation displacement is 0.05875 mm on the fore element and 0.04814 mm on the aft element.'),
          'difference between two interpolants of the same published vertices', 'manufacturing or as-built uncertainty')
    klass(c3, 'nominal solver nose: max vertex-to-source-polyline distance', ev.val('distance', NOSE+'summary.json', 'geometry', 2, 'same_distance_mm', unit='mm'),
          'departure of the paneled nose from the published polyline', 'measured geometry tolerance or Hausdorff bound')
    klass(c3, 'cubic solver nose: max vertex-to-source-polyline distance', ev.val('distance', NOSE+'summary.json', 'geometry', 4, 'same_distance_mm', unit='mm'),
          'departure of the paneled nose from the published polyline', 'measured geometry tolerance or Hausdorff bound')
    klass(c3, 'nominal solver nose: negative local turns on the lower nose (source polyline has 0)', ev.val('negative turns', NOSE+'summary.json', 'geometry', 2, 'negative_lower_nose_turns'),
          'a sampling-induced turning reversal absent from the source polyline', 'a feature of the manufactured contour')
    klass(c3, 'fore-only cubic versus nominal: max nose vertex-to-other-polyline distance', ev.val('distance', FRA+'summary.json', 'low_lift_comparison', 'nose_vertex_to_other_polyline_max_mm', unit='mm'),
          'difference between two reconstructions of the fore nose', 'manufacturing uncertainty bound')
    klass(c3, 'drag response to interpolant / aft +-0.13 mm z-shift at CL 0.284 (max |change|)', ev.val('max change', RCA+'summary.json', 'scenario_changes', 'cl028', 'maximum_absolute_CD_change_counts', unit='counts'),
          'model drag sensitivity to four declared reconstructions at matched lift', 'a statistical uncertainty envelope or certified design gain')
    klass(c3, 'drag response to interpolant / aft +-0.13 mm z-shift at CL 1.041 (max |change|)', ev.val('max change', RCA+'summary.json', 'scenario_changes', 'cl104', 'maximum_absolute_CD_change_counts', unit='counts'),
          'model drag sensitivity to four declared reconstructions at matched lift', 'a statistical uncertainty envelope or certified design gain')
    klass(c3, 'drag response to the fore-only cubic reconstruction at CL 0.056', ev.val('change', FRA+'summary.json', 'low_lift_comparison', 'delta_CD_counts', unit='counts'),
          'model drag sensitivity to one alternative fore interpolant at matched lift', 'a geometry correction')
    klass(c3, 'drag response to the fore-only cubic reconstruction at CL 0.175',
          ev.val('archived paired drag change of the fore-cubic scenario', FRA+'summary.json', 'paired_drag_changes_counts', 'fore_source_cubic'),
          'nothing: the target failed and the archived value is null', 'a finite value inferred from nearby evaluations', status='unknown')
    klass(c3, 'same-angle (-2.56888 deg) cubic minus nominal drag: not matched lift', ev.val('change', FRA+'summary.json', 'same_angle_endpoints', 0, 'delta_CD_counts', unit='counts'),
          'response near the transition-sensitive interval at identical incidence but different lift', 'a matched-lift drag sensitivity')
    for r in public:
        klass(c3, f'public {r["file"]} {r["element"]}: max vertex distance to the test outline (file not admitted)',
              ev.cell('distance', GEO+'public_geometry_comparison.csv', 'maximum_vertex_distance_mm', unit='mm', file=r['file'], element=r['element']),
              'distance after a declared scale normalisation of an unqualified sharp-TE file', 'experimental or as-built error')
    klass(c4, 'as-built contour deviation of the tested model',
          ev.val('original_author_coordinates_recovered', summary, 'original_author_coordinates_recovered'),
          'nothing: original coordinates and a contour survey were not recovered', 'the 0.13 mm reported tolerance or any reconstruction spread', status='unknown')
    klass(c4, 'as-built element placement / rigging (gap, overlap, deflection)',
          ev.quote(GEO+'FINDINGS.md', 'The gap and overlap are our explicit geometric diagnostics, not a recovered rigging table.'),
          'nothing: no rigging table recovered; gap and overlap are diagnostics of the drawing', 'the +-0.13 mm aft-shift sensitivity', status='unknown')
    klass(c4, 'statement that the +-0.13 mm scale is not a measured rigging uncertainty', ev.quote(RC_PROTOCOL, 'it is not measured rigging uncertainty, a probability distribution, or a confidence bound'),
          'a protocol limitation', 'a quantified as-built uncertainty', status='limitation recorded')
    klass(c4, 'three printed tap centres outside the reconstructed branch extents', tap_records[7],
          'an unresolved inconsistency between Table II and the reconstructed outline', 'a measured contour error', status='unresolved inconsistency')
    klass(c5, 'max |plotted - printed| tap x, Figure 12(a)', plotted_evidence[0], 'offset between plotted marker x and printed Table II x', 'PDF rounding or measurement uncertainty', status='unresolved inconsistency')
    klass(c5, 'max |plotted - printed| tap x, Figure 12(b)', plotted_evidence[1], 'offset between plotted marker x and printed Table II x', 'PDF rounding or measurement uncertainty', status='unresolved inconsistency')
    klass(c5, 'max |plotted - printed| tap x, Figure 7(a-d)', plotted_evidence[5], 'offset between plotted marker x and printed Table II x', 'PDF rounding or measurement uncertainty', status='unresolved inconsistency')
    klass(c5, 'taps failing the 2e-5c exact-table check, Figures 12(a)/12(b)', plotted_evidence[4], 'number of taps whose plotted x differs from the table by more than 2e-5c', 'a failed tap association (ordered association passes)', status='unresolved inconsistency')
    for t in lb_targets:
        klass(c5, f'max change of fore-lower/aft-upper regional Cp RMSE between conventions, CL {t["target_cl"]}',
              ev.val('RMSE change', LBA+'diagnosis.json', 'cases', ev.find(LBA+'diagnosis.json', ['cases'], case=t['case']), 'maximum_fore_lower_aft_upper_RMSE_coordinate_change'),
              'sensitivity of the model-experiment score to the coordinate convention', 'experimental uncertainty')
    klass(c6, 'first two supported fore-lower taps x/c', tap_records[2], 'the only measured nose stations on the fore lower surface', 'coverage of the suction peak', status='unmeasured between taps')
    klass(c6, 'spacing of those taps', tap_records[5], 'the unmeasured interval', 'tap size', status='unmeasured between taps')
    for record in nose_evidence:
        if record['kind'] == 'derived' and record['name'].startswith('MODEL'):
            klass(c6, record['name'], record, 'location and level of the archived MODEL minimum inside the unmeasured interval', 'a measurement or a pressure correction', status='unmeasured between taps')
    klass(c6, 'measured peak Cp between the first two fore-lower taps',
          ev.quote(PLAN, "The model's unresolved between-tap nose peak must remain visibly unmeasured."),
          'nothing: no tap between them', 'the model minimum', status='unknown')
    c7 = 'reported measurement uncertainty (kept separate from the four classes)'
    stated = "the authors' statement in the 2014 article, which overlaps the 2010 campaign"
    klass(c7, 'drag coefficient, low-drag range', article_drag[1], stated, 'a 2010-report statement or a Cp uncertainty')
    klass(c7, 'drag coefficient, approaching stall', article_drag[2], stated, 'a 2010-report statement or a Cp uncertainty')
    klass(c7, "lift coefficient, 'at higher lift'", article_lift[1], stated, 'a 2010-report statement or a Cp uncertainty')
    klass(c7, 'pitching-moment coefficient', article_lift[2], stated, 'a 2010-report statement or a Cp uncertainty')
    klass(c7, 'force-coefficient uncertainty stated by the 2010 report itself', unknowns['U19']['basis'][0],
          'nothing: the archived 2010 report text contains no uncertainty statement', 'the 2014 statement', status='unknown')
    klass(c7, 'surface pressure coefficient', ev.val('pressure_uncertainty', psrc, 'pressure_uncertainty'),
          'nothing: no Cp uncertainty budget recovered', 'extraction offsets', status='unknown')

    # ------------------------------------------------------------------ reconciliation of earlier manifests
    reconciliation = reconcile(ev)
    changed = sorted({name for r in reconciliation for name in r['changed_files'].split('; ') if name})

    def number(rel, field, **where):
        return float(ev.row(rel, **where)[field])

    def note(topic, values, explanation):
        """values: list of (description, value, source file)."""
        for _, _, rel in values:
            ev.path(rel)
        notes.append(dict(topic=topic, explanation=explanation,
                          values=[dict(description=d, value=v, source=rel, sha256=ev.digest(rel)) for d, v, rel in values]))

    note('Published-MSES minus experiment drag at CL 0.284 [counts]',
         [('experiment from Figure 13(c) markers', number(DRAG+'target_references.csv', 'mses_minus_experiment_counts', case='cl028'), DRAG+'target_references.csv'),
          ('experiment from the appendix table', number(LBA+'targets.csv', 'mses_error_CD_counts', case='cl028'), LBA+'targets.csv')],
         'Same published-MSES value; only the experimental reference convention differs. Both are retained.')
    note('Published-MSES minus experiment drag at CL 1.041 [counts]',
         [('experiment from Figure 13(c) markers', number(DRAG+'target_references.csv', 'mses_minus_experiment_counts', case='cl104'), DRAG+'target_references.csv'),
          ('experiment from the appendix table', number(TRIP+'matched_lift_targets.csv', 'mses_minus_experiment_free_counts', reynolds=1000000, CL='1.041'), TRIP+'matched_lift_targets.csv')],
         'Same cause: Figure 13(c) markers versus appendix table as the experimental reference. Both are retained.')
    note('MultiFoil minus experiment drag at CL 0.284, identical saved model state [counts]',
         [('experiment from Figure 13(c) markers', number(RCA+'targets.csv', 'error_CD_counts', job='cl028_nominal'), RCA+'targets.csv'),
          ('experiment from the appendix table', number(LBA+'targets.csv', 'error_CD_counts', case='cl028'), LBA+'targets.csv')],
         'The model state is an exact repeat; the difference is again the experimental reference convention.')
    region = dict(case='cl028', surface='aft_upper', coordinate='table', region='aft_bubble_recovery')
    note('Aft-upper recovery Cp RMSE (x/c 0.86-0.94) at CL 0.284',
         [('MultiFoil, scored against the Figure 7(d) copy', number(LBA+'pressure_metrics.csv', 'rmse', method='multifoil', **region), LBA+'pressure_metrics.csv'),
          ('MultiFoil, scored against the Figure 12(a) copy', number(RCA+'pressure_metrics.csv', 'rmse', method='nominal', **region), RCA+'pressure_metrics.csv'),
          ('published MSES, scored against the Figure 7(d) copy', number(LBA+'pressure_metrics.csv', 'rmse', method='mses', **region), LBA+'pressure_metrics.csv'),
          ('published MSES, scored against the Figure 12(a) copy', number(RCA+'pressure_metrics.csv', 'rmse', method='mses', **region), RCA+'pressure_metrics.csv')],
         f'Two plotted copies of the same measurement differ by up to {ev.json(tsum)["maximum_repeated_pressure_Cp_difference"]:.6f} Cp '
         '(strict copy check failed; both copies retained, no alignment).')
    note('Figure 12 pressure panels',
         [('panels listed as published-MSES pressure comparisons', ev.json(quali)['pressure_coverage']['published_MSES_free_figures'], quali),
          ('panels actually extracted', [c['figure'] for c in ev.json(psrc)['cases']], psrc)],
         'The report has three Figure 12 panels; only (a) and (b) were extracted. Panel (c) is recorded in the matrix as not extracted.')
    note('Superseded next-step text in the geometry qualification',
         [('next_step', ev.json(summary)['next_step'], summary),
          ('Figure 13(c) extraction that was completed afterwards', ev.json(DRAG+'source.json')['product'], DRAG+'source.json')],
         'The dated summary still lists the Figure 13 drag extraction as remaining; the dated file is left unchanged.')
    current = dict(ev.reconciled, **ev.living)
    notes.append(dict(
        topic='Files whose current sha256 differs from a hash recorded in an earlier S414 receipt',
        explanation='Every recorded OUTPUT hash of the evidence directories still matches. The differing files are living '
                    'documents and production/analysis code edited after those receipts were issued; the earlier receipts are left intact.',
        values=[dict(description=f'differs from the hash recorded in {sum(name in r["changed_files"].split("; ") for r in reconciliation)} earlier sha256 maps',
                     value=name, source=name, sha256=current[name]) for name in changed]))

    # ------------------------------------------------------------------ coverage table
    kinds = ['experiment_pressure', 'experiment_drag', 'mses2010_drag', 'mses2010_pressure', 'mses_transition',
             'experiment_transition_reattachment_location', 'multifoil_archived_audited_state']
    for cl, key in case_of.items():
        for kind in kinds:
            entry = cover.get(key, {}).get(kind)
            require(entry is not None, f'No coverage entry for {key} {kind}')
            for rel in entry['sources']:
                ev.path(rel)
            coverage.append(dict(case='CL '+cl, evidence=kind, status=entry['status'], detail=entry['detail'],
                                 source_page_figure=entry['page'], matrix_row_ids=entry['row_id'],
                                 source_paths='; '.join(entry['sources']),
                                 source_sha256='; '.join(ev.digest(s) for s in entry['sources'])))

    # ------------------------------------------------------------------ outputs
    rows = M.ordered()
    for item in unknowns.values():
        item['rows'].sort()
    fields = ['row_id', 'reference_group', 'reference', 'case', 'quantity', 'admission_status', 'admitted_use',
              'not_admitted_for', 'strict_matched_geometry_validation', 'conditions', 'source_page_figure',
              'measurement_resolution', 'reported_uncertainty', 'unknown_quantities', 'values', 'evidence_quotes',
              'source_paths', 'source_sha256', 'notes']
    write_csv(out/'evidence_matrix.csv', M.flat(), fields)
    unknown_rows = [dict(id=k, unknown_quantity=v['text'], matrix_row_ids=','.join(v['rows']),
                         basis=' | '.join(fmt(b) for b in v['basis']),
                         source_paths='; '.join(dict.fromkeys(b['source'] for b in v['basis'])),
                         source_sha256='; '.join(ev.digest(x) for x in dict.fromkeys(b['source'] for b in v['basis'])))
                    for k, v in unknowns.items()]
    write_csv(out/'unknown_quantities.csv', unknown_rows)
    write_csv(out/'case_coverage.csv', coverage)
    write_csv(out/'uncertainty_classes.csv', classes)
    write_csv(out/'manifest_reconciliation.csv', reconciliation)
    write_json(out/'evidence_matrix.json', dict(
        gate=gate['value'], gate_source=PLAN, statuses=dict(
            admitted='qualified for the stated limited use', conditional='usable only under the stated condition',
            not_admitted='explicitly rejected for use by an existing finding',
            not_extracted='present in an archived source document but never digitized',
            missing='not available in the archived evidence'),
        reference_groups=GROUPS,
        value_kinds=dict(data='read directly from the cited CSV/JSON field',
                         derived='count, extremum or difference computed here from the cited archived column',
                         quote='literal quotation checked to be present in the cited file',
                         source_text='literal fragment located on one line of an archived source-text rendition'),
        rows=rows, unknown_registry=unknowns, reconciliation_notes=notes))
    reconciliation_only = dict(sorted((k, v) for k, v in ev.reconciled.items() if k not in ev.files and k not in ev.living))
    write_json(out/'input_manifest.json', dict(
        files=dict(sorted(ev.files.items())), living_documents=dict(sorted(ev.living.items())),
        reconciliation_only_files=reconciliation_only,
        code_sha256={'validation/report_s414_evidence_matrix.py': sha(Path(__file__).resolve())},
        categories=dict(files='evidence read for content or cited by the matrix; any change invalidates the matrix',
                        living_documents='plan and reference register, quoted here and expected to be edited later',
                        reconciliation_only_files='hashed only to compare with hashes recorded in earlier S414 receipts')))
    status_count = {x: sum(r['admission_status'] == x for r in rows) for x in STATUSES}
    stated = [r['row_id'] for r in rows if '2014 article' in r['reported_uncertainty'] or 'reported manufacturing tolerance' in r['reported_uncertainty']
              or r['reported_uncertainty'].startswith('drag 0.00005')]
    not_applicable = [r['row_id'] for r in rows if r['reported_uncertainty'].startswith('not applicable')]
    not_reported = [r['row_id'] for r in rows if r['row_id'] not in stated and r['row_id'] not in not_applicable]
    require(all('not reported' in r['reported_uncertainty'] for r in rows if r['row_id'] in not_reported), 'Unclassified uncertainty entry')
    facts = dict(taps=len([r for r in ev.rows(LBA+'pressure_at_taps.csv') if r['case'] == 'cl0056' and r['coordinate'] == 'table'
                           and r['method'] == 'multifoil' and r['included'] == 'True']),
                 lower=lower_taps[:2], tolerance=ev.json(summary)['reported_model_contour_tolerance_mm'],
                 reynolds=primary['reynolds'], mach=primary['mach'], not_reported=not_reported, stated=stated,
                 not_applicable=not_applicable)
    write_findings(out, rows, coverage, classes, unknowns, notes, reconciliation, changed, kinds, case_of, gate, status_count, facts)
    write_json(out/'summary.json', dict(
        complete=True, stage='V0', generated_utc=datetime.now(timezone.utc).isoformat(),
        gate=gate['value'], strict_matched_geometry_validation='conditional',
        exact_tested_geometry_established=ev.json(summary)['tested_geometry_benchmark_admitted'],
        rows=len(rows), rows_by_group={GROUPS[g]: n for g, n in M.count.items()},
        admitted=status_count['admitted'], conditional=status_count['conditional'],
        not_admitted=status_count['not_admitted'], not_extracted=status_count['not_extracted'], missing=status_count['missing'],
        unknown_items=len(unknowns), unknown_item_ids=list(unknowns),
        rows_with_uncertainty_not_reported=len(not_reported), row_ids_with_uncertainty_not_reported=not_reported,
        rows_with_a_stated_uncertainty_or_tolerance=len(stated), row_ids_with_a_stated_uncertainty_or_tolerance=stated,
        rows_where_uncertainty_is_not_applicable=len(not_applicable),
        coverage_entries=len(coverage), coverage_missing=sum(c['status'].startswith('missing') for c in coverage),
        uncertainty_class_items=len(classes), uncertainty_items_unknown=sum(c['status'] == 'unknown' for c in classes),
        evidence_files=len(ev.files), living_documents=len(ev.living), reconciliation_only_files=len(reconciliation_only),
        earlier_receipt_hash_maps_checked=len(reconciliation),
        earlier_receipt_hashes_checked=sum(r['entries'] for r in reconciliation),
        earlier_receipt_hashes_changed=sum(r['changed'] for r in reconciliation),
        earlier_receipt_hashes_unresolved=sum(r['unresolved'] for r in reconciliation),
        files_changed_since_earlier_receipts=changed, reconciliation_notes=len(notes),
        new_flow_solves=0, new_mses_runs=0, new_downloads=0, new_pdf_extractions=0,
        runtime=dict(python=platform.python_version())))
    print(f'{len(rows)} rows: '+', '.join(f'{k} {v}' for k, v in status_count.items())
          + f'; {len(unknowns)} unknown items; {len(coverage)} coverage entries; {len(classes)} uncertainty items')


def write_findings(out, rows, coverage, classes, unknowns, notes, reconciliation, changed, kinds, case_of, gate, status_count, facts):
    def link(rel, label=None):
        return f'[{label or rel.split("/")[-1]}]({os.path.relpath(ROOT/rel, out)})'

    def ids(items, limit=6):
        return ', '.join(items[:limit])+(f' and {len(items)-limit} more' if len(items) > limit else '') if items else 'registry only'

    short = {'experiment_pressure': 'Experiment pressure', 'experiment_drag': 'Experiment drag',
             'mses2010_drag': 'MSES-2010 drag', 'mses2010_pressure': 'MSES-2010 pressure',
             'mses_transition': 'MSES transition', 'experiment_transition_reattachment_location': 'Experiment transition / reattachment',
             'multifoil_archived_audited_state': 'Archived audited MultiFoil state'}
    lines = [
        '# S414 evidence and coverage matrix (stage V0) — 4 October 2026', '',
        f'**Gate: {gate["value"]}** Every row that supports a model-versus-reference comparison therefore carries '
        '`strict_matched_geometry_validation = conditional (geometry not qualified)`.', '',
        'This is a consolidation of the existing manifests only. No flow solve, MSES run, download, literature search '
        'or PDF extraction was performed. Every value is read by the script from an archived CSV/JSON file, is a '
        'labelled count/extremum/difference of archived columns, or is a literal quotation checked against the cited '
        'file. Quantities recorded nowhere are written as `unknown`.', '',
        f'The matrix has **{len(rows)} rows**: {status_count["admitted"]} admitted for a stated limited use, '
        f'{status_count["conditional"]} conditional, {status_count["not_admitted"]} not admitted, '
        f'{status_count["not_extracted"]} present in an archived source but not extracted, and '
        f'{status_count["missing"]} missing. It registers **{len(unknowns)} explicitly unknown quantities**. '
        f'No uncertainty is reported for {len(facts["not_reported"])} rows. {len(facts["stated"])} rows carry only the '
        'force-coefficient uncertainty stated in the overlapping 2014 article or the reported contour tolerance; for the '
        f'{len(facts["not_applicable"])} archived model-state rows a measurement uncertainty does not apply.', '',
        '- [evidence_matrix.csv](evidence_matrix.csv) / [evidence_matrix.json](evidence_matrix.json): one row per reference series and case/quantity',
        '- [case_coverage.csv](case_coverage.csv): which evidence exists for each of the five cases',
        '- [uncertainty_classes.csv](uncertainty_classes.csv): the four uncertainty classes kept distinct, plus the separate items',
        '- [unknown_quantities.csv](unknown_quantities.csv): the registry of unknown quantities and the rows they affect',
        '- [manifest_reconciliation.csv](manifest_reconciliation.csv): every sha256 map of the earlier S414 receipts compared with the files as they are now',
        '- [input_manifest.json](input_manifest.json), [summary.json](summary.json), [verification.json](verification.json)', '',
        '## Coverage by case', '',
        f'Re {facts["reynolds"]:,}, Mach {facts["mach"]}, free transition. "exists" means an archived file with a recorded hash; '
        'it does not mean admitted without condition.', '',
        '| Evidence | '+' | '.join('CL '+cl for cl in case_of)+' |', '|---|'+'---|'*len(case_of)]
    for kind in kinds:
        cells = []
        for cl in case_of:
            entry = next(c for c in coverage if c['case'] == 'CL '+cl and c['evidence'] == kind)
            first = entry['matrix_row_ids'].split(',')
            text = '**missing**' if entry['status'].startswith('missing') else entry['status']
            if kind == 'multifoil_archived_audited_state':
                text += ': '+entry['detail'].split(' | ')[0]
            cells.append(f'{text} ({first[0]}{" …" if len(first) > 1 else ""})')
        lines.append(f'| {short[kind]} | '+' | '.join(cells)+' |')
    lines += ['', 'Row identifiers in parentheses point into the matrix. File paths and hashes for every cell are in '
              '[case_coverage.csv](case_coverage.csv). The CL 1.041 measurement comes from its Figure 12(b) copy. '
              'The published-MSES pressure curves are at nominal CL 0.28 and 1.04, not at the measured 0.284 and 1.041. '
              'Several archived model states are exact repeats of one another, so the number of state directories is not '
              'a number of independent solutions. One archived target (fore-only cubic reconstruction at CL 0.175) failed and has no state.', '',
              '## What each reference may and may not be used for', '',
              '| Reference group | May be used for | May not be used for | Rows |', '|---|---|---|---|']
    uses = [
        ('A', f'Qualified pressure markers at {facts["taps"]} supported taps per case (original masks), matched-lift drag from the printed appendix table, and the stated lift/normal-force reduction, on the nominal reconstructed outline',
         'Exact matched-geometry validation; measured transition, bubble or reattachment locations; the nose-peak magnitude between taps; a Cp uncertainty band; a fitted angle, alignment or pressure-averaging correction'),
        ('B', 'Published drag at CL 0.056, 0.175 and 0.284; published pressure at nominal CL 0.28 as the 0.284 control; the lift curve as an initialisation seed only',
         'Experiment; drag at CL -0.072; pressure at -0.072, 0.056 or 0.175; any transition location; a run of an installed MSES; assumed identical geometry or settings; mixing with the 2014 full-e^N series'),
        ('C', 'Regression controls after a proposed remedy: the CL 1.041 pressure/drag records and the audited archived MultiFoil states (fixed-angle states only on the identical-incidence basis)',
         'Blind validation (already inspected); grid-convergence or uniqueness claims; the failed CL 0.175 cubic target as a finite value'),
        ('D', 'Separate, condition-qualified supporting series; the only archived statement of force-coefficient uncertainty; the transition-method distinction',
         'Merging with the 2010 envelope series; counting the overlapping experiment as independent; any numerical comparison (no curve is digitized)'),
        ('E', 'Total-force context at matched lift for free, fore-tripped and both-tripped configurations at Re 0.7 and 1.0 million',
         'A clean wake-disturbance benchmark; tripped pressure comparisons; closure fitting; the Re 1.5 million set (Mach mismatch)'),
    ]
    for group, may, may_not in uses:
        members = [r['row_id'] for r in rows if r['row_id'].startswith(group)]
        lines.append(f'| {GROUPS[group]} | {may} | {may_not} | {members[0]}–{members[-1]} |')
    lines += ['', 'Published-MSES and model values are never listed as experiment. The 2010 envelope and 2014 full-e^N series '
              'carry separate reference labels throughout.', '', '## Explicitly unknown quantities and missing evidence', '',
              '| Id | Unknown quantity | Matrix rows |', '|---|---|---|']
    for key, item in unknowns.items():
        lines.append(f'| {key} | {item["text"]} | {ids(item["rows"])} |')
    absent = [r for r in rows if r['admission_status'] in ('missing', 'not_extracted')]
    lines += ['', 'Rows with no usable evidence:', '', '| Row | Reference | Case | Quantity | Status |', '|---|---|---|---|---|']
    lines += [f'| {r["row_id"]} | {r["reference"].split(" (")[0]} | {r["case"]} | {r["quantity"]} | {r["admission_status"].replace("_", " ")} |' for r in absent]
    lines += ['', '## Uncertainty classes kept distinct', '',
              '| Class | Items with a value | Items unknown | What it is not |', '|---|---:|---:|---|']
    for name in dict.fromkeys(c['uncertainty_class'] for c in classes):
        items = [c for c in classes if c['uncertainty_class'] == name]
        lines.append(f'| {name} | {sum(c["status"] != "unknown" for c in items)} | {sum(c["status"] == "unknown" for c in items)} | {items[0]["not_to_be_read_as"]} |')
    lines += ['', f'The actual as-built uncertainty is unknown. The reported {facts["tolerance"]} mm contour tolerance, the extraction '
              'repeatability and the reconstruction spread are different things and are not combined into one band. '
              'The printed-table versus plotted tap-coordinate difference is an unresolved source inconsistency. '
              f'There is no tap between x/c {facts["lower"][0]} and {facts["lower"][1]} on the fore lower surface, and the archived model '
              'pressure minimum lies inside that interval in the three states with an archived nose-pressure file: the peak is unmeasured.', '',
              '## Reconciliation of the existing manifests', '',
              f'{len(reconciliation)} sha256 maps from the earlier S414 receipts were compared with the current files '
              f'({sum(r["entries"] for r in reconciliation)} recorded hashes; {sum(r["changed"] for r in reconciliation)} differ; '
              f'{sum(r["unresolved"] for r in reconciliation)} could not be located). Every recorded output hash of the evidence '
              f'directories still matches. The differing entries all point to these {len(changed)} files:', '']
    lines += [f'- `{name}`' for name in changed]
    lines += ['', 'These are living documents and code edited after the receipts were issued. The earlier receipts are '
              'left intact; a guard that pins one of these files will not pass against the current tree.', '',
              'Archived values that describe the same quantity but differ:', '']
    for item in notes[:-1]:
        values = '; '.join(f'{show(v["value"])} ({v["description"]}, {link(v["source"])})' for v in item['values'])
        lines.append(f'- **{item["topic"]}**: {values}. {item["explanation"]}')
    lines += ['', '## Sources and limits', '',
              'Primary evidence directories: '+', '.join(link(d+'FINDINGS.md', d.split('/')[-2]) for d in (GEO, DRAG, TRIP, LBA, RCA, FRA, NOSE, PB))
              + '; protocol '+link(LB_PROTOCOL)+'; register '+link(REGISTER)+'; plan '+link(PLAN)+'.', '',
              'Limits of this matrix: it adds no evidence and resolves no inconsistency. Values taken from prose are '
              'literal quotations, marked `quote` or `source_text` in the matrix, and were not re-derived from the source PDFs. '
              'The 2014 uncertainty statement is read from an archived text rendition in which the relation and sign glyphs '
              'are unreadable. Exact tested geometry, grid convergence and solution uniqueness remain unestablished.', '',
              '## Reproduction', '', '```bash', 'python3 validation/report_s414_evidence_matrix.py   # refuses to overwrite a completed output',
              'python3 validation/verify_s414_evidence_matrix.py', '```', '',
              'Both scripts use only the Python standard library. Use `--out` for a new output directory.', '']
    (out/'FINDINGS.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
