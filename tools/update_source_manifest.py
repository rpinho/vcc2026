"""Fingerprint the source export using the original publication procedure."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]


def main():
    paths = sorted(
        p for p in ROOT.rglob('*')
        if p.is_file() and '.git' not in p.relative_to(ROOT).parts
        and p.name not in {'PUBLISH_CHECK.md', 'SOURCE_MANIFEST.sha256'}
    )
    (ROOT / 'SOURCE_MANIFEST.sha256').write_text(''.join(
        hashlib.sha256(p.read_bytes()).hexdigest() + '  '
        + str(p.relative_to(ROOT)) + '\n' for p in paths
    ))
    print(f'Manifest: {len(paths)} files')


if __name__ == '__main__':
    main()
