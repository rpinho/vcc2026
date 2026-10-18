"""Network-free crash, concurrency, ownership and identity checks for reader charges."""
import runtime as rt
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from transfer_budget import TransferBudget


with tempfile.TemporaryDirectory(dir=rt.ROOT/'artifacts/tmp') as tmp:
    path = Path(tmp)/'charges.jsonl'
    # A process can die after I/O reservation without a statistics checkpoint.
    subprocess.run([sys.executable, '-c',
                    'from transfer_budget import TransferBudget; import os,sys; '
                    'b=TransferBudget(sys.argv[1],100); b.charge(37,"failed request"); os._exit(0)',
                    str(path)], check=True)
    b = TransferBudget(path, 100)
    assert b.total == 37
    b.bind('source', {'etag': 'frozen-file'})
    try:
        duplicate = TransferBudget(path, 100)
    except BlockingIOError:
        pass
    else:
        raise AssertionError('Duplicate reader acquired journal')

    def attempt(_):
        try:
            b.charge(7, 'concurrent request')
            return True
        except AssertionError:
            return False

    with ThreadPoolExecutor(max_workers=16) as pool:
        results = list(pool.map(attempt, range(32)))
    assert sum(results) == 9 and b.total == 100
    b.close()
    b = TransferBudget(path, 100)
    assert b.total == 100
    b.bind('source', {'etag': 'frozen-file'})
    try:
        b.bind('source', {'etag': 'changed-file'})
    except AssertionError:
        pass
    else:
        raise AssertionError('Changed source accepted')
    b.close()
    try:
        changed_cap = TransferBudget(path, 101)
    except AssertionError:
        pass
    else:
        raise AssertionError('Changed allocation accepted')
    with path.open('a') as out:
        out.write('{"event":')
    try:
        torn = TransferBudget(path, 100)
    except json.JSONDecodeError:
        pass
    else:
        raise AssertionError('Torn journal accepted')

rt.save('artifacts/transfer_budget_fixture.json', dict(
    status='PASS', utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
    checks=['abrupt process exit retains charge', 'concurrent cap enforcement',
            'duplicate reader refused', 'source change refused',
            'allocation change refused', 'torn journal fails closed'],
    budget_sha256=rt.sha('transfer_budget.py'), reader_sha256=rt.sha('acquire_cd4_source.py'),
    limitation='Network-free accounting fixture; live bank numerical reproduction remains separate'))
print('PASS durable transfer budget (6 checks)')
