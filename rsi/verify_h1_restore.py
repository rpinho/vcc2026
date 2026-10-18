"""Network-free failure/restart checks for H1 payload reuse and transfer charges."""
import runtime as rt
import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import numpy as np

from download_checkpoint import DownloadCheckpoint
from transfer_budget import TransferBudget

specs = {'data': ((8,), np.float32), 'indices': ((8,), np.int32)}
identity = {'source_etag': 'fixed', 'plan_sha256': 'fixed'}
with tempfile.TemporaryDirectory(dir=rt.ROOT/'artifacts/tmp') as tmp:
    root = Path(tmp)/'payload'
    journal = Path(tmp)/'charges.jsonl'
    subprocess.run([sys.executable, '-c', '''
import os, sys, numpy as np
from download_checkpoint import DownloadCheckpoint
from transfer_budget import TransferBudget
c=DownloadCheckpoint(sys.argv[1], {'data':((8,),np.float32),'indices':((8,),np.int32)}, {'source_etag':'fixed','plan_sha256':'fixed'})
b=TransferBudget(sys.argv[2],100)
b.charge(16,'first committed range')
for a in c.arrays.values(): a[:4]=np.arange(4)
c.commit(1)
b.charge(16,'attempt interrupted before checkpoint')
for a in c.arrays.values(): a[4:]=99
os._exit(0)
''', str(root), str(journal)], check=True)
    c = DownloadCheckpoint(root, specs, identity)
    b = TransferBudget(journal, 100)
    assert c.next_plan == 1 and b.total == 32
    for array in c.arrays.values():
        np.testing.assert_array_equal(array[:4], np.arange(4))
    b.charge(16, 'replay of uncommitted range')
    for array in c.arrays.values():
        array[4:] = np.arange(4, 8)
    # A failed payload flush must never publish a later cursor.
    with patch('download_checkpoint.os.fsync', side_effect=OSError('injected failure')):
        try:
            c.commit(2)
        except OSError:
            pass
        else:
            raise AssertionError('Failure injection did not run')
    assert json.loads(c.path.read_text())['next_plan'] == 1
    c.commit(2)
    for array in c.arrays.values():
        np.testing.assert_array_equal(array, np.arange(8))
        array._mmap.close()
    del c
    c = DownloadCheckpoint(root, specs, identity)
    assert c.next_plan == 2 and b.total == 48
    for changed_specs, changed_identity in [
        (specs, dict(identity, source_etag='changed')),
        (specs, dict(identity, plan_sha256='changed')),
        (dict(specs, data=((9,), np.float32)), identity),
    ]:
        try:
            DownloadCheckpoint(root, changed_specs, changed_identity)
        except AssertionError:
            pass
        else:
            raise AssertionError('Changed checkpoint identity accepted')
    c.remove()
    b.close()
    assert not root.exists()

rt.save('artifacts/h1_restore_fixture.json', dict(
    status='PASS', utc=datetime.now(timezone.utc).isoformat(),
    checks=['committed payload and cursor survive abrupt process exit',
            'interrupted and replayed reads both remain charged',
            'failed fsync cannot advance cursor',
            'resumed arrays equal uninterrupted arrays exactly',
            'changed source, plan or array layout refused',
            'checkpoint removed after successful publication'],
    reader_sha256=rt.sha('acquire_h1_rich.py'),
    checkpoint_sha256=rt.sha('download_checkpoint.py'),
    budget_sha256=rt.sha('transfer_budget.py'),
    limitation='Network-free failure fixture; live H1 source and output identities require separate verification'))
print('PASS H1 payload checkpoint and durable transfer charges (6 checks)')
