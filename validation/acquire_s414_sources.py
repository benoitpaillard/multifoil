"""Acquire the bounded S414 geometry sources, retaining hashes and failed access.

Run from any directory. Successful downloads are checked and reused on restart.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'validation/results/s414_geometry_qualification_20261002'
SOURCES = {
    'dtic_access_test.html': 'https://apps.dtic.mil/sti/pdfs/ADA545821.pdf',
    'metu_catalog.html': 'https://www.ae.metu.edu.tr/tuncer/ae546/multi/',
    's414_report_2010.pdf': 'https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf',
    's414.pts': 'https://www.ae.metu.edu.tr/tuncer/ae546/multi/s414.pts',
    's414-1.pts': 'https://www.ae.metu.edu.tr/tuncer/ae546/multi/s414-1.pts',
    'archive_metadata.json': 'https://archive.org/metadata/DTIC_ADA545821',
    'archive_s414_report_2010.pdf': 'https://archive.org/download/DTIC_ADA545821/DTIC_ADA545821.pdf',
    'airfoils_publications.html': 'https://airfoils.com/publications/',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(item):
    name, url = item
    record = dict(name=name, url=url, access_mode='anonymous HTTPS GET',
                  request_parameters={}, supplied_checksum=None,
                  retrieved_utc=datetime.now(timezone.utc).isoformat(),
                  scope='S414 2010 tested contour and 2014 comparison qualification')
    try:
        r = requests.get(url, timeout=30)
        record.update(http_status=r.status_code, content_type=r.headers.get('Content-Type'),
                      final_url=r.url, bytes=len(r.content))
        if len(r.content) > 20_000_000:
            raise ValueError('20 MB bounded acquisition limit exceeded')
        path = OUT / 'sources' / name
        path.write_bytes(r.content)
        record.update(path=str(path.relative_to(ROOT)), sha256=digest(path))
        r.raise_for_status()
        if name.endswith('.pdf') and not r.content.startswith(b'%PDF'):
            raise ValueError('Response is not PDF')
        record['download_status'] = 'saved'
    except Exception as error:
        record.update(download_status='failed', error=str(error),
                      blocker='Access failure alone does not establish an account requirement')
    return record


def main():
    (OUT / 'sources').mkdir(parents=True, exist_ok=True)
    path = OUT / 'source_manifest.json'
    manifest = json.loads(path.read_text()) if path.exists() else {'sources': []}
    for record in manifest['sources']:
        if record.get('path'):
            assert digest(ROOT / record['path']) == record['sha256'], record['path']
    done = {r['url'] for r in manifest['sources']}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for record in pool.map(fetch, [(k, v) for k, v in SOURCES.items() if v not in done]):
            manifest['sources'].append(record)
            path.write_text(json.dumps(manifest, indent=2) + '\n')
            print(record['download_status'], record['name'], record.get('bytes'), flush=True)
    manifest.update(product='S414 tested-geometry and comparison-source qualification',
                    versions='RDECOM TR 10-D-112 (August 2010); Coder et al., Journal of Aircraft 51(6), 2014; public .pts versions unspecified',
                    acquisition_script_sha256=digest(Path(__file__)))
    metadata_path=OUT/'sources/archive_metadata.json'
    if metadata_path.exists():
        item=next(f for f in json.loads(metadata_path.read_text())['files'] if f['name']=='DTIC_ADA545821.pdf')
        for record in manifest['sources']:
            if record['name']=='archive_s414_report_2010.pdf':
                record['supplied_checksum']={k:item[k] for k in ('sha1','md5')}
                record['checksum_source']='https://archive.org/metadata/DTIC_ADA545821'
    old=json.loads((ROOT/'validation/results/reference_review_20260924/source_manifest.json').read_text())
    manifest['existing_sources']=[r for r in old['sources'] if r.get('path','').endswith('/s414_article.pdf')]
    for record in manifest['existing_sources']:
        assert digest(ROOT/record['path'])==record['sha256']
    path.write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__':
    main()
