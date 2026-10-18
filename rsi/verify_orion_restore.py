"""Network-free tests of the actual Orion bridge and reader failure/restart paths."""
import runtime as rt
from datetime import datetime, timezone
import io
import gzip
import os
from pathlib import Path
import tempfile
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import requests
from urllib3.response import HTTPResponse

import acquire_orion as reader
from durable_public_store import PublicStore
from transfer_budget import TransferBudget


class Response:
    def __init__(self, body=b'', headers=None, status=200, url='https://public.invalid/fixed'):
        self.headers = requests.structures.CaseInsensitiveDict(headers or {})
        self.raw = HTTPResponse(body=io.BytesIO(body), headers=self.headers, preload_content=False)
        self.status_code, self.url = status, url

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError('fixture error')


class Remote:
    def __init__(self, budget):
        self.budget = budget
        self.payload = b'0123456789abcdefghij'
        self.etag = '"fixture-v1"'
        self.calls = []
        self.fail_next = False
        self.truncate_next = False

    def get(self, url, **kwargs):
        if '/api/datasets/' in url:
            assert self.budget.total >= 1_000_000
            return Response(rt.json.dumps(dict(sha='a' * 40)).encode())
        start, end, partial = PublicStore.bounds(kwargs['headers'].get('Range'), len(self.payload))
        amount = end - start + 1
        # This must already be on disk when the remote request begins.
        rows = [rt.json.loads(line) for line in self.budget.path.read_text().splitlines()]
        assert rows[-1]['event'] == 'charge' and rows[-1]['bytes'] == amount
        self.calls.append((start, end))
        if self.fail_next:
            self.fail_next = False
            raise requests.ConnectionError('fixture failed request')
        body = self.payload[start:end + 1]
        if self.truncate_next:
            self.truncate_next = False
            body = body[:2]
        headers = {'Content-Length': str(amount), 'ETag': self.etag}
        if partial:
            headers['Content-Range'] = f'bytes {start}-{end}/{len(self.payload)}'
        return Response(body, headers, 206 if partial else 200)

    def head(self, url, **kwargs):
        assert '/resolve/' + 'a' * 40 + '/' in url
        assert self.budget.total >= 64_000
        return Response(headers={'Content-Length': str(len(self.payload)), 'ETag': self.etag})

    def close(self):
        pass


checks = []
with tempfile.TemporaryDirectory(dir=rt.ROOT / 'artifacts/tmp') as tmp:
    # Exercise actual urllib3 decoding, including a server that ignores identity.
    for encoding in ['identity', 'gzip']:
        budget = TransferBudget(Path(tmp) / (encoding + '.jsonl'), 2_000_000)
        remote = Remote(budget)
        body = rt.json.dumps(dict(sha='a' * 40)).encode()
        payload = gzip.compress(body) if encoding == 'gzip' else body
        def metadata(url, **kwargs):
            assert kwargs['headers']['Accept-Encoding'] == 'identity'
            return Response(payload, {'Content-Encoding': encoding, 'Content-Length': str(len(payload))})
        remote.get = metadata
        bridge = PublicStore(budget, session=remote)
        assert bridge.revision == 'a' * 40 and budget.total == 1_000_000
        bridge.close()
        budget.close()
    budget = TransferBudget(Path(tmp) / 'oversize.jsonl', 2_000_000)
    remote = Remote(budget)
    remote.get = lambda *args, **kwargs: Response(gzip.compress(b' ' * 1_000_001), {'Content-Encoding': 'gzip'})
    try:
        PublicStore(budget, session=remote)
    except AssertionError:
        pass
    else:
        raise AssertionError('Oversized decoded metadata accepted')
    assert budget.total == 1_000_000 and 'repository_revision' not in budget.bindings
    budget.close()
    checks.append('Actual urllib3 identity/gzip metadata decoding; decompressed size bound; failed attempt stays charged')

    old_inputs = dict(selection='selection', genes='genes', offsets='offsets', version=908,
        implementation='be119c8701c38551f9621a7969a84d0037ca56b3ee757413045f2dc0bcc3faf2',
        bridge='a5d48a092335c7198461c73a5940d7cdc5f7d3a343497cd13d0cf376ce010396')
    inputs = dict(old_inputs, implementation=rt.sha('acquire_orion.py'), bridge=rt.sha('durable_public_store.py'))
    for case in ['valid', 'changed_selection', 'payload', 'source_bound', 'extra_charge']:
        journal = Path(tmp) / (case + '.jsonl')
        budget = TransferBudget(journal, 3_000_000_000)
        budget.bind('frozen_inputs', old_inputs)
        budget.charge(1_000_000, 'repository revision metadata upper bound')
        prefix = journal.read_bytes()
        checkpoint = Path(tmp) / (case + '_payload')
        requested = dict(inputs)
        if case == 'changed_selection': requested['selection'] = 'changed'
        if case == 'payload': checkpoint.mkdir()
        if case == 'source_bound': budget.bind('repository_revision', 'a' * 40)
        if case == 'extra_charge': budget.charge(1, 'payload attempt')
        try:
            reader.bind_inputs(budget, requested, checkpoint)
        except AssertionError:
            assert case != 'valid'
            assert 'metadata_decode_fix_v1' not in budget.bindings
        else:
            assert case == 'valid'
            assert journal.read_bytes().startswith(prefix) and budget.total == 1_000_000
            budget.bind('repository_revision', 'a' * 40)
            budget.close()
            budget = TransferBudget(journal, 3_000_000_000)
            reader.bind_inputs(budget, requested, checkpoint)
            assert budget.total == 1_000_000 and budget.bindings['frozen_inputs'] == old_inputs
            try:
                reader.bind_inputs(budget, dict(requested, bridge='another edit'), checkpoint)
            except AssertionError:
                pass
            else:
                raise AssertionError('Further untested migration accepted')
        budget.close()
    checks.append('Append-only pre-payload repair migration and restart; refuses changed inputs, existing payload/source, extra charges and further edits')

with tempfile.TemporaryDirectory(dir=rt.ROOT / 'artifacts/tmp') as tmp:
    journal = Path(tmp) / 'bridge.jsonl'
    budget = TransferBudget(journal, 2_000_000)
    remote = Remote(budget)
    bridge = PublicStore(budget, session=remote)
    client = requests.Session()
    client.trust_env = False
    url = bridge.options()['endpoint'] + '/public/data/HEK293T/expression.lance/fixture.lance'
    try:
        response = client.get(url, headers={'Range': 'bytes=2-6'}, timeout=5)
        assert response.status_code == 206 and response.content == b'23456'
        assert budget.total == 1_064_005
        response = client.get(url, headers={'Range': 'bytes=-3'}, timeout=5)
        assert response.content == b'hij'
        remote.fail_next = True
        response = client.get(url, headers={'Range': 'bytes=0-4'}, timeout=5)
        assert response.status_code == 503 and budget.total == 1_064_013
        response = client.get(url, headers={'Range': 'bytes=0-4'}, timeout=5)
        assert response.content == b'01234' and budget.total == 1_064_018
        remote.truncate_next = True
        try:
            client.get(url, headers={'Range': 'bytes=0-4'}, timeout=5)
        except requests.RequestException:
            pass
        else:
            raise AssertionError('Truncated payload accepted')
        assert budget.total == 1_064_023
        checks.append('Actual HTTP bridge: exact/suffix ranges; pre-I/O disk charges; failed and truncated requests retain full charge')
        budget.charge(budget.limit - budget.total - 2, 'fixture remaining budget')
        calls = len(remote.calls)
        response = client.get(url, headers={'Range': 'bytes=0-4'}, timeout=5)
        assert response.status_code == 503 and len(remote.calls) == calls
        checks.append('Cumulative cap refuses GET before remote connection')
    finally:
        bridge.close()
        budget.close()
        client.close()
    budget = TransferBudget(journal, 2_000_000)
    assert budget.bindings['repository_revision'] == 'a' * 40
    # Change identity and give HEAD enough budget in this independent journal.
    budget.close()
    identity_journal = Path(tmp) / 'identity.jsonl'
    budget = TransferBudget(identity_journal, 2_000_000)
    budget.bind('repository_revision', 'a' * 40)
    key = 'data/HEK293T/expression.lance/fixture.lance'
    budget.bind('file:' + key, dict(size=20, etag='"earlier-version"'))
    remote = Remote(budget)
    bridge = PublicStore(budget, session=remote)
    client = requests.Session()
    client.trust_env = False
    try:
        response = client.get(bridge.options()['endpoint'] + '/public/' + key, timeout=5)
        assert response.status_code == 503 and not remote.calls
        assert budget.total == 64_000
        checks.append('Restart reuses immutable revision and refuses changed file ETag before payload')
    finally:
        bridge.close()
        budget.close()
        client.close()


class FakeStore:
    current = None

    def __init__(self, budget):
        self.budget = budget
        self.revision = 'b' * 40
        budget.bind('repository_revision', self.revision)
        FakeStore.current = self

    def options(self):
        return {}

    def close(self):
        pass


class Expression:
    version = 908
    fail_at = None
    starts = []

    def take(self, ix, columns):
        lo = int(ix[0] // 3)
        self.starts.append(lo)
        FakeStore.current.budget.charge(len(ix) * 12, 'fixture Lance payload attempt')
        if lo == self.fail_at:
            raise OSError('fixture interrupted Lance read')
        return pa.table(dict(cell_integer_id=ix // 3,
                             gene_integer_id=np.tile([2, 0, 1], len(ix) // 3),
                             value=np.tile([1, 2, 3], len(ix) // 3)))


with tempfile.TemporaryDirectory(dir=rt.ROOT / 'artifacts/tmp') as tmp:
    prior_cwd = Path.cwd()
    try:
        os.chdir(tmp)
        folder = Path('artifacts/tasks')
        folder.mkdir(parents=True)
        ids = np.arange(260, dtype=np.int64)
        sel = pd.DataFrame(dict(cell_integer_id=ids, source_row=ids, n_genes_by_counts=3,
            cell_start_index=ids * 3, total_counts=6, target_gene=np.where(ids < 4, 'non-targeting', 'T'),
            batch='batch', evaluation_role=np.where(ids < 2, 'input_control',
                np.where(ids < 4, 'reference_control', 'perturbation'))))
        selection = folder / 'HEK293T_holdout_selection.parquet'
        sel.to_parquet(selection)
        pd.DataFrame(dict(gene_id=['A', 'A', 'B'])).to_parquet('artifacts/HEK293T_genes.parquet')
        np.savez('artifacts/HEK293T_offset_mapping.npz', version=908, original_start=[0], logical_start=[0])
        rt.save(folder / 'manifest.json', dict(HEK293T_holdout=dict(context='HEK293T', selection_sha256=rt.sha(selection))))
        rt.save('artifacts/download_allocations.json', dict(new_allocations=dict(HEK293T=3_000_000_000)))
        Expression.fail_at, Expression.starts = 128, []
        with patch.object(reader, 'PublicStore', FakeStore), patch.object(reader.lance, 'dataset', return_value=Expression()), patch.object(rt, 'guard'):
            try:
                reader.main('HEK293T')
            except OSError:
                pass
            else:
                raise AssertionError('Interrupted read not exercised')
            root = Path('artifacts/HEK293T_payload_checkpoint')
            assert rt.json.loads((root / 'state.json').read_text())['next_plan'] == 128
            budget = TransferBudget('artifacts/HEK293T_charges.jsonl', 3_000_000_000)
            assert budget.total == 256 * 3 * 12
            budget.close()
            checks.append('Actual reader: interrupted second batch keeps first 128 cells and failed-attempt charges')
            Expression.fail_at = None
            original_write = ad.AnnData.write_h5ad

            def fail_publication(self, *args, **kwargs):
                original_write(self, *args, **kwargs)
                raise OSError('fixture interrupted publication')

            with patch.object(ad.AnnData, 'write_h5ad', fail_publication):
                try:
                    reader.main('HEK293T')
                except OSError:
                    pass
                else:
                    raise AssertionError('Failed publication not exercised')
            assert Expression.starts == [0, 128, 128, 256]
            assert rt.json.loads((root / 'state.json').read_text())['next_plan'] == 260
            np.testing.assert_array_equal(np.load(root / 'data.npy'), np.tile([1, 2, 3], 260))
            np.testing.assert_array_equal(np.load(root / 'indices.npy'), np.tile([1, 0, 0], 260))
            checks.append('Restart skips committed batch; duplicate-gene publication failure leaves raw checkpoint intact')
            with patch.object(reader.lance, 'dataset', side_effect=AssertionError('Payload refetched')):
                reader.main('HEK293T')
            assert not root.exists()
            a = ad.read_h5ad(folder / 'HEK293T_holdout.h5ad')
            np.testing.assert_array_equal(a.var_names, ['A', 'B'])
            np.testing.assert_array_equal(a.X.toarray(), np.tile([5, 1], (260, 1)))
            for column in ['target_gene', 'batch', 'evaluation_role', 'source_row', 'cell_integer_id']:
                np.testing.assert_array_equal(a.obs[column].to_numpy(), sel[column].to_numpy())
            receipt = rt.json.loads((folder / 'HEK293T_holdout.json').read_text())
            assert receipt['sha256'] == rt.sha(folder / 'HEK293T_holdout.h5ad')
            assert receipt['bytes'] == (260 + 128) * 3 * 12
            checks.append('Completed checkpoint publishes without Lance I/O; exact counts/rows/gene aggregation and receipt hash')
    finally:
        os.chdir(prior_cwd)

rt.save('artifacts/orion_restore_fixture.json', dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
    checks=checks, network='No public network; real localhost HTTP bridge, mocked upstream and Lance payloads',
    implementation_sha256={file: rt.sha(file) for file in ['acquire_orion.py', 'durable_public_store.py',
        'download_checkpoint.py', 'transfer_budget.py']}, fixture_sha256=rt.sha(__file__),
    limitations='Actual public Lance compatibility, source revision, task identity and production transfer use remain UNVERIFIED'))
print('PASS', len(checks), 'Orion bridge/reader failure-restart checks')
