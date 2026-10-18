"""Exercise the actual CD4 reader on local HDF5 bytes, including failed restarts."""
import runtime as rt
import gc
import io
import json
import os
import runpy
import tempfile
from pathlib import Path
from unittest.mock import patch
from datetime import datetime, timezone

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import tree_prefetch

reader = rt.ROOT/'acquire_cd4_destination.py'
original_root = rt.ROOT
original_write = ad.AnnData.write_h5ad
checks = []
with tempfile.TemporaryDirectory(dir=rt.ROOT/'artifacts/tmp') as temporary:
    root = Path(temporary)
    job = root/'job'
    task = job/'artifacts/tasks'
    task.mkdir(parents=True)
    (root/'data').mkdir()
    pd.DataFrame({'gene': ['b', 'a']}).to_csv(root/'data/gene_names.csv', index=False)
    counts = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [3, 2, 1]], np.float32)
    source = ad.AnnData(sparse.csr_matrix(counts), var=pd.DataFrame(
        {'gene_name': ['a', 'b', 'a']}, index=['id1', 'id2', 'id3']))
    original_write(source, root/'source.h5ad')
    payload = (root/'source.h5ad').read_bytes()
    selection = pd.DataFrame(dict(source_row=[3, 0, 2, 1], target_gene=['non-targeting', 'x', 'non-targeting', 'x'],
        batch=['lane']*4, guide_id=['g']*4,
        evaluation_role=['reference_control', 'perturbation', 'input_control', 'perturbation']))
    selection.to_parquet(task/'CD4_matched_search_selection.parquet', index=False)
    rt.save(task/'cd4_wide_blueprint.json', dict(n_real=96, selection_sha256=rt.sha(task/'CD4_matched_search_selection.parquet')))
    rt.save(task/'manifest.json', {})
    rt.save(job/'artifacts/tree_prefetch_fixture.json', dict(status='PASS bitwise index addresses and data', implementation_sha256=rt.sha(original_root/'tree_prefetch.py')))
    rt.save(job/'artifacts/CD4_cell_feasibility.json', dict(etag='fixed', source_size=len(payload)))
    (job/'artifacts/CD4_cell_census.parquet').write_bytes(b'fixture identity')
    (job/'tree_prefetch.py').write_bytes((original_root/'tree_prefetch.py').read_bytes())
    calls = []
    fail_payload = [True]

    class FakeRemote(io.BytesIO):
        def __init__(self, url, cap, block):
            super().__init__(payload)
            self.etag, self.size, self.bytes = 'fixed', len(payload), 0
            self.http = self
        def fetch(self, start, n):
            self.bytes += n
            return payload[start:start+n]
        def enable_tree_prefetch(self): pass
        def disable_tree_prefetch(self): pass

    class FakeSession:
        def get(self, url, headers, **kwargs):
            rows = [json.loads(line) for line in (job/'artifacts/CD4_destination_charges.jsonl').read_text().splitlines()]
            a, b = map(int, headers['Range'].removeprefix('bytes=').split('-'))
            assert rows[-1]['event']=='charge' and rows[-1]['bytes']==b-a+1
            calls.append((a, b))
            if fail_payload[0]: raise OSError('injected request failure')
            class Response:
                status_code = 206
                raw = io.BytesIO(payload[a:b+1])
                def __enter__(self): return self
                def __exit__(self, *args): pass
            return Response()

    def run_failure(message):
        try: runpy.run_path(str(reader), run_name='__main__')
        except OSError as error: assert message in str(error)
        else: raise AssertionError('Expected injected failure')
        gc.collect()

    previous_cwd = Path.cwd()
    try:
        os.chdir(job)
        with patch.object(rt, 'ROOT', job), patch.object(rt, 'guard'), patch.object(tree_prefetch, 'PrefetchRemote', FakeRemote), patch('requests.Session', FakeSession):
            run_failure('injected request failure')
            assert json.loads((job/'artifacts/CD4_destination_payload_checkpoint/state.json').read_text())['next_plan']==0
            assert len(calls)>=3
            checks.append('Every failed HTTP attempt is durably charged; failed batch cannot advance checkpoint')
            fail_payload[0] = False
            with patch.object(ad.AnnData, 'write_h5ad', side_effect=OSError('injected publication failure')):
                run_failure('injected publication failure')
            gc.collect()  # Release the mock's exception traceback, as process exit would.
            saved = json.loads((job/'artifacts/CD4_destination_payload_checkpoint/state.json').read_text())
            assert saved['next_plan']>0
            number_calls = len(calls)
            runpy.run_path(str(reader), run_name='__main__')
            assert len(calls)==number_calls
            checks.append('Restart after publication failure reuses complete payload without another payload request')
            result = ad.read_h5ad(task/'CD4_matched_search.h5ad')
            np.testing.assert_array_equal(result.X.toarray(), np.column_stack((counts[:, 1], counts[:, 0]+counts[:, 2])))
            assert list(result.var_names)==['b', 'a'] and list(result.obs.source_row)==[0, 1, 2, 3]
            checks.append('Actual reader output preserves sorted frozen rows, official gene order and duplicate-gene sums exactly')
            assert not (job/'artifacts/CD4_destination_payload_checkpoint').exists()
            rows = [json.loads(line) for line in (job/'artifacts/CD4_destination_charges.jsonl').read_text().splitlines()]
            charged = sum(row['bytes'] for row in rows if row['event']=='charge')
            assert charged==3*251_000_000+sum(b-a+1 for a,b in calls)
            assert json.loads((task/'CD4_matched_search.json').read_text())['bytes']==charged
            checks.append('All invocation reserves and retries survive restarts; checkpoint removed only after task/receipt/manifest publication')
    finally:
        os.chdir(previous_cwd)

rt.save('artifacts/cd4_restore_fixture.json', dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
    checks=checks, reader_sha256=rt.sha(reader), checkpoint_sha256=rt.sha('download_checkpoint.py'),
    budget_sha256=rt.sha('transfer_budget.py'), limitation='Local HDF5 bytes and mocked HTTP; live source/output verification remains separate'))
print('PASS actual CD4 reader restart/output fixture', len(checks), 'checks')
