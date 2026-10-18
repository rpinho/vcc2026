"""Run the actual incremental reader on synthetic HDF5, with failed restarts."""
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

reader = rt.ROOT/'acquire_cd4_depth192.py'
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
    counts = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [3, 2, 1], [4, 1, 2], [8, 2, 5]], np.float32)
    source = ad.AnnData(sparse.csr_matrix(counts), var=pd.DataFrame(
        {'gene_name': ['a', 'b', 'a']}, index=['id1', 'id2', 'id3']))
    original_write(source, root/'source.h5ad')
    payload = (root/'source.h5ad').read_bytes()
    projected = np.column_stack((counts[:, 1], counts[:, 0]+counts[:, 2])).astype(np.int32)
    selection = pd.DataFrame(dict(source_row=list(range(6)),
        target_gene=['x', 'x', 'non-targeting', 'non-targeting', 'x', 'x'],
        batch=['lane']*6, guide_id=['g']*6,
        evaluation_role=['perturbation', 'perturbation', 'input_control', 'reference_control', 'perturbation', 'perturbation']))
    unused = selection.iloc[[0]].copy(); unused['source_row'] = 999; unused['target_gene'] = 'ignored'
    parent_obs = pd.concat([selection.iloc[:4], unused], ignore_index=True)
    base = ad.AnnData(sparse.csr_matrix(np.vstack([projected[:4], [9, 9]])), obs=parent_obs, var=pd.DataFrame(index=['b', 'a']))
    original_write(base, task/'CD4_depth120_search.h5ad')
    name = 'CD4_depth192_search'
    selection.to_parquet(task/f'{name}_selection.parquet', index=False)
    selection.iloc[4:].to_parquet(task/f'{name}_increment.parquet', index=False)
    journal = job/'artifacts/CD4_destination_charges.jsonl'
    journal.write_text(json.dumps(dict(event='charge', bytes=100, total=100))+'\n')
    second_journal = job/'artifacts/CD4_depth120_charges.jsonl'
    second_journal.write_text(json.dumps(dict(event='charge', bytes=200, total=200))+'\n')
    rt.save(job/'artifacts/download_allocations.json', dict(new_allocations={'CD4_depth192': 1_200_000_000}, limit_bytes=15_000_000_000))
    blueprint = task/f'{name}_blueprint.json'
    rt.save(blueprint, dict(selection_sha256=rt.sha(task/f'{name}_selection.parquet'),
        increment_sha256=rt.sha(task/f'{name}_increment.parquet'), parent_task_sha256=rt.sha(task/'CD4_depth120_search.h5ad'),
        previous_charge_journals={str(journal.relative_to(job)): rt.sha(journal), str(second_journal.relative_to(job)): rt.sha(second_journal)},
        increment_ceiling_bytes=1_200_000_000))
    rt.save(job/'artifacts/CD4_depth192_design_audit.json', dict(
        status='PASS independent count-only depth192 design audit',
        blueprint_sha256=rt.sha(blueprint), additional_cells=2))
    rt.save(task/'manifest.json', {'base_preserved': True})
    rt.save(job/'artifacts/tree_prefetch_fixture.json', dict(status='PASS bitwise index addresses and data', implementation_sha256=rt.sha(original_root/'tree_prefetch.py')))
    rt.save(job/'artifacts/CD4_cell_feasibility.json', dict(etag='fixed', source_size=len(payload)))
    (job/'artifacts/CD4_cell_census.parquet').write_bytes(b'fixture identity')
    (job/'tree_prefetch.py').write_bytes((original_root/'tree_prefetch.py').read_bytes())
    rt.save(job/'artifacts/cd4_depth192_fixture.json', dict(status='PASS', reader_sha256=rt.sha(reader)))
    calls, metadata_calls = [], []
    fail_payload = [True]

    class FakeRemote(io.BytesIO):
        def __init__(self, url, cap, block):
            super().__init__(payload)
            metadata_calls.append(url)
            self.etag, self.size, self.bytes = 'fixed', len(payload), 0
            self.http = self
        def fetch(self, start, n):
            self.bytes += n
            return payload[start:start+n]
        def enable_tree_prefetch(self): pass

    class FakeSession:
        def get(self, url, headers, **kwargs):
            rows = [json.loads(line) for line in (job/'artifacts/CD4_depth192_charges.jsonl').read_text().splitlines()]
            a, b = map(int, headers['Range'].removeprefix('bytes=').split('-'))
            assert rows[-1]['event'] == 'charge' and rows[-1]['bytes'] == b-a+1
            assert headers['If-Match'] == 'fixed'
            calls.append((a, b))
            if fail_payload[0]:
                raise OSError('injected request failure')
            class Response:
                status_code = 206
                raw = io.BytesIO(payload[a:b+1])
                def __enter__(self): return self
                def __exit__(self, *args): pass
            return Response()

    def run_failure(message, kind=OSError):
        try:
            runpy.run_path(str(reader), run_name='__main__')
        except kind as error:
            assert message in str(error)
        else:
            raise AssertionError('Expected injected failure')
        gc.collect()

    previous_cwd = Path.cwd()
    try:
        os.chdir(job)
        with patch.object(rt, 'ROOT', job), patch.object(rt, 'guard'), patch.object(tree_prefetch, 'PrefetchRemote', FakeRemote), patch('requests.Session', FakeSession):
            run_failure('injected request failure')
            state = job/'artifacts/CD4_depth192_payload_checkpoint/state.json'
            assert json.loads(state.read_text())['next_plan'] == 0
            assert len(calls) >= 3 and len(metadata_calls) == 1
            checks.append('Failed HTTP attempts charged before GET; failed batch cannot advance checkpoint')
            number_calls = len(calls)
            cached = job/'artifacts/CD4_depth192_index.npz'
            original_cache = cached.read_bytes()
            cached.write_bytes(original_cache+b'corrupt')
            run_failure('', AssertionError)
            assert len(calls) == number_calls and len(metadata_calls) == 1
            cached.write_bytes(original_cache)
            checks.append('Corrupted committed index rejected before any request')
            original_journal = journal.read_bytes()
            journal.write_bytes(original_journal+b' ')
            run_failure('', AssertionError)
            assert len(calls) == number_calls
            journal.write_bytes(original_journal)
            checks.append('Changed original CD4 accounting refused before any request')
            fail_payload[0] = False
            with patch.object(ad.AnnData, 'write_h5ad', side_effect=OSError('injected publication failure')):
                run_failure('injected publication failure')
            gc.collect()
            assert json.loads(state.read_text())['next_plan'] > 0
            number_calls = len(calls)
            runpy.run_path(str(reader), run_name='__main__')
            assert len(calls) == number_calls and len(metadata_calls) == 1
            checks.append('Publication failure restarts with no payload or metadata refetch; only one metadata reserve')
            result = ad.read_h5ad(task/f'{name}.h5ad')
            np.testing.assert_array_equal(result.X.toarray(), projected)
            np.testing.assert_array_equal(result.X[:4].toarray(), base.X[:4].toarray())
            assert list(result.var_names) == ['b', 'a'] and list(result.obs.source_row) == list(range(6))
            checks.append('All inherited counts/control identities unchanged; increment projected and duplicate genes summed exactly')
            assert not state.exists()
            charges = [json.loads(line) for line in (job/'artifacts/CD4_depth192_charges.jsonl').read_text().splitlines()]
            charged = sum(row['bytes'] for row in charges if row['event'] == 'charge')
            assert charged == 251_000_000+sum(b-a+1 for a, b in calls)
            result_receipt = json.loads((task/f'{name}.json').read_text())
            assert result_receipt['cumulative_cd4_bytes'] == charged+300 <= 2_600_000_000+1_200_000_000
            assert json.loads((task/'manifest.json').read_text())['base_preserved']
            assert name in json.loads((task/f'{name}_manifest_extension.json').read_text())
            checks.append('Original plus incremental retry charges match receipt; frozen manifest unchanged; separate extension published before checkpoint removal')
    finally:
        os.chdir(previous_cwd)

rt.save('artifacts/cd4_depth192_fixture.json', dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
    checks=checks, reader_sha256=rt.sha(reader), checkpoint_sha256=rt.sha('download_checkpoint.py'),
    budget_sha256=rt.sha('transfer_budget.py'), tree_prefetch_sha256=rt.sha('tree_prefetch.py'),
    limitation='Local HDF5 and mocked HTTP. Live count/identity/quota verification is separate; more cells do not guarantee valid scales.'))
print('PASS actual CD4 incremental reader', len(checks), 'checks')
