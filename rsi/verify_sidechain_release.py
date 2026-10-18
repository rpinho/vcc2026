"""Run the author's pinned synthetic kernel fixtures with network access disabled.

This verifies the released kernels in this environment, not SER17's missing
submission inputs or hidden-context score. No biological benchmark is opened.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
RELEASE = ROOT / 'vendor/sidechain_release'


def main():
    receipt = json.loads((ROOT / 'artifacts/sidechain_release_receipt.json').read_text())
    assert receipt['status'].startswith('PASS')
    for name, record in receipt['files'].items():
        data = (RELEASE / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == record['sha256'], name
    # Refuse accidental public I/O by test code or dependencies, including DNS.
    def offline(*args, **kwargs):
        raise RuntimeError('Network disabled during public-kernel verification')
    socket.socket.connect = offline
    socket.socket.connect_ex = offline
    socket.create_connection = offline
    socket.getaddrinfo = offline
    os.environ['HF_HUB_OFFLINE'] = '1'
    os.environ['TRANSFORMERS_OFFLINE'] = '1'
    sys.path.insert(0, str(RELEASE / 'src'))
    import pytest
    import numpy
    import scipy
    import anndata
    import cell_eval2
    # The neighbour suite imports torch at collection time. Do not download a
    # large unrelated runtime merely to deserialize its absent STRING table.
    deferred = {'tests/test_neighbour_arm.py': 'torch is absent; public STRING geometry/pool also absent'}
    tests = sorted(str(RELEASE / name) for name in receipt['files']
                   if name.startswith('tests/') and name not in deferred)
    output = ROOT / 'artifacts/sidechain_release_pytest.xml'
    rc = pytest.main(['-q', '-p', 'no:cacheprovider', '--basetemp',
                      str(ROOT / 'artifacts/sidechain_fixture_tmp'),
                      '--junitxml', str(output), *tests])
    suites = ET.parse(output).getroot()
    counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites.iter('testsuite'))
              for key in ['tests', 'failures', 'errors', 'skipped']}
    record = {
        'utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS' if rc == 0 else 'FAIL', 'pytest_rc': int(rc), **counts,
        'commit': receipt['commit'], 'network_disabled': True,
        'deferred_suites': deferred,
        'release_receipt_sha256': hashlib.sha256((ROOT / 'artifacts/sidechain_release_receipt.json').read_bytes()).hexdigest(),
        'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'junit_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'environment': {'python': sys.version, 'numpy': numpy.__version__,
                        'scipy': scipy.__version__, 'anndata': anndata.__version__,
                        'cell_eval2': cell_eval2.__version__},
        'benchmark_or_holdout_scored': False,
        'limitation': 'Author synthetic fixtures only; missing exact SER17 release command, scatter table, neighbour geometry/pool and source moments remain unresolved. No external official score can be paired with a local hybrid.'}
    (ROOT / 'artifacts/sidechain_kernel_audit.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2), flush=True)
    return int(rc)


if __name__ == '__main__':
    raise SystemExit(main())
