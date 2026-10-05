#!/usr/bin/env python3
"""Archive bounded anonymous source checks for the S702 condition audit."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PDF = 'validation/data/s702/NASA-CR-20220014964.pdf'
PDF_HASH = '01b8bb6fc9b4fb2bc3246611106d64f980c03c41b18dcd1740a59a769e66cf03'
SOURCES = [
    ('author_home', 'https://airfoils.com/'),
    ('author_publications', 'https://airfoils.com/publications/'),
    ('uiuc_coordinates', 'https://m-selig.ae.illinois.edu/ads/coord_database.html'),
    ('psu_article', 'https://pure.psu.edu/en/publications/exploration-of-a-low-speed-flapped-slotted-natural-laminar-flow-a-2/'),
    ('publisher_article', 'https://arc.aiaa.org/doi/10.2514/1.C037802'),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_conditions_audit_20260922')
    out = ap.parse_args().out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest = out/'source_manifest.json'
    if manifest.exists():
        raise ValueError('Preserve the original audit; use a new output directory')
    digest = hashlib.sha256((ROOT/PDF).read_bytes()).hexdigest()
    if digest != PDF_HASH:
        raise ValueError('Archived NASA report changed')
    def fetch(item):
        label, url = item
        row = dict(id=label, url=url, access_mode='public anonymous HTTPS GET',
                   retrieved_utc=datetime.now(timezone.utc).isoformat(),
                   scope='One landing/index page; no bulk crawl or private coordinate request',
                   supplied_checksum=None)
        try:
            with urlopen(Request(url, headers={'User-Agent': 'MultiFoil-source-audit/1.0'}), timeout=20) as response:
                data = response.read()
                path = out/(label+'.html'); path.write_bytes(data)
                row.update(status='downloaded', http_status=response.status, final_url=response.url,
                           file=path.name, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            row.update(status='access_failed', http_status=getattr(exc, 'code', None), error=str(exc),
                       interpretation='Access failure alone does not establish an account, licence or subscription requirement.')
        return row
    with ThreadPoolExecutor(max_workers=5) as pool:
        sources = list(pool.map(fetch, SOURCES))
    sources.insert(0, dict(id='nasa_report', product='NASA/CR-20220014964, December 2022',
        url='https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf',
        access_mode='Previously archived public report, read locally; no new download',
        status='local_checksum_verified', file=PDF, sha256=digest, supplied_checksum=None,
        scope='Printed pages 7–9, 13, 15; Figure 1 and fixed/free-transition experiment descriptions'))
    sources.extend(dict(id=Path(name).stem, file=name, status='local_input_verified',
                        sha256=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()) for name in (
        'validation/data/s702/source.json', 'validation/data/s702/geometry_extraction.json',
        'validation/data/s702/audit.json', 'validation/data/s702/pressure_cl1.04/source.json',
        'validation/results/s702_ncrit_analysis/pressure_errors.csv'))
    manifest.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        audit_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        qualification='Bounded public-source search, not proof that coordinates are unavailable everywhere.',
        sources=sources), indent=2)+'\n')
    print(json.dumps([dict(id=r['id'], status=r['status'], http_status=r.get('http_status')) for r in sources], indent=2))


if __name__ == '__main__':
    main()
