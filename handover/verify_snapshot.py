#!/usr/bin/env python3
"""Verify all files in the delivered MultiFoil snapshot, without running CFD.

Uses only Python's standard library. Extra files, such as a local environment
or new study outputs, are permitted. Intentional edits to manifest-listed files
are reported as changes; retain this manifest as the original release record.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / 'handover/archive_manifest.json').read_text())
    if manifest.get('schema_version') != 1:
        raise ValueError('Unsupported archive manifest')
    failures, checked, checked_bytes, seen = [], 0, 0, set()
    for item in manifest['files']:
        relative = PurePosixPath(item['path'])
        if relative.is_absolute() or '..' in relative.parts or item['path'] in seen:
            raise ValueError('Unsafe or duplicate manifest path')
        seen.add(item['path'])
        path = root.joinpath(*relative.parts)
        if not path.resolve().is_relative_to(root) or path.is_symlink():
            failures.append(dict(path=item['path'], reason='Path escapes root or is a symlink'))
        elif not path.is_file():
            failures.append(dict(path=item['path'], reason='Missing'))
        elif path.stat().st_size != item['bytes'] or digest(path) != item['sha256']:
            failures.append(dict(path=item['path'], reason='Content changed'))
        else:
            checked += 1
            checked_bytes += item['bytes']
    print(json.dumps(dict(status='failed' if failures else 'passed',
                          checked_files=checked, checked_bytes=checked_bytes,
                          failures=failures, scope='Snapshot integrity, not aerodynamic validation'), indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
