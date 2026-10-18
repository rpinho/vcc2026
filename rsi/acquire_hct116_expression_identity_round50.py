"""Metadata-only immutable expression manifest, complete inventory and gene axis.

No expression data files are fetched. Encoded page metadata is a later step.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import urllib.parse

import runtime as rt
import requests
import lance
import numpy as np
import pandas as pd

from acquire_hct116_full_metadata_round48 import fetch_file, verify
from durable_public_store import PublicStore
from transfer_budget import TransferBudget

ROOT = Path('artifacts/HCT116_expression_identity_round50')
KEY = 'HCT116_expression_identity_round50'
LIMIT = 15_000_000
REVISION = '598aa5442a4ee9a8f1a383a843ab5ce142904d04'
REPO = 'slaf-project/X-Atlas-Orion'
PREFIX = 'data/HCT116/'
API = f'https://huggingface.co/api/datasets/{REPO}/tree/{REVISION}/'
JOURNAL = Path('artifacts/'+KEY+'_charges.jsonl')
PLAN = ROOT/'plan_retry.json'
CAP = 768_000


def load(p):
    return rt.json.loads(Path(p).read_text())


def inventory(session, budget, relative, recursive=False):
    """Each page is separately charged and cached, including its next link."""
    stem = relative.replace('/', '_')
    url = API+PREFIX+relative+f'?recursive={str(recursive).lower()}&limit=1000'
    initial = url
    entries = []
    seen = set()
    for page in range(10):
        assert url not in seen, 'Inventory pagination cycle'
        seen.add(url)
        dest = ROOT/f'{stem}_page{page}.json'
        if dest.exists():
            saved = load(dest)
            assert saved['url'] == url
            body = bytes.fromhex(saved['body_hex'])
            assert saved['sha256'] == rt.hashlib.sha256(body).hexdigest()
            next_url = saved['next_url']
        else:
            rt.guard()
            budget.charge(CAP, 'Bounded pinned inventory page attempt '+url)
            with session.get(url, headers={'Accept-Encoding':'identity'}, stream=True,
                             allow_redirects=False, timeout=(10,60)) as response:
                response.raise_for_status()
                assert response.status_code == 200
                body = response.raw.read(CAP, decode_content=True)
                assert len(body) < CAP, 'Inventory page exceeds charged cap'
                links = response.links
                assert set(links) <= {'next'}, 'Unexpected inventory link relation'
                next_url = links.get('next', {}).get('url')
                if next_url:
                    next_url = urllib.parse.urljoin(url, next_url)
                    before, after = urllib.parse.urlsplit(initial), urllib.parse.urlsplit(next_url)
                    assert (after.scheme, after.netloc, urllib.parse.unquote(after.path)) == (before.scheme, before.netloc, urllib.parse.unquote(before.path)), next_url
            value = rt.json.loads(body)
            assert isinstance(value, list) and len(value) <= 1000
            rt.save(dest, dict(url=url, next_url=next_url, body_hex=body.hex(),
                              sha256=rt.hashlib.sha256(body).hexdigest()))
        rows = rt.json.loads(body)
        budget.bind('inventory_page:'+str(dest), rt.sha(dest))
        assert isinstance(rows, list) and len(rows) <= 1000
        if next_url:
            before, after = urllib.parse.urlsplit(initial), urllib.parse.urlsplit(next_url)
            assert (after.scheme, after.netloc, urllib.parse.unquote(after.path)) == (before.scheme, before.netloc, urllib.parse.unquote(before.path)), next_url
        for e in rows:
            assert e['path'].startswith(PREFIX+relative+'/') and '..' not in e['path']
        entries.extend(rows)
        if not next_url:
            assert len({e['path'] for e in entries}) == len(entries), 'Duplicate inventory path'
            return dict(revision=REVISION, complete=True, pages=page+1, entries=entries)
        url = next_url
    raise AssertionError('Inventory pagination exceeds frozen ten-page ceiling')


def freeze():
    ROOT.mkdir(exist_ok=True)
    assert not PLAN.exists()
    assert load('artifacts/HCT116_extension_preflight_round49/audit.json')['status'].startswith('PASS')
    paths = [Path(__file__), Path('acquire_hct116_full_metadata_round48.py'),
        Path('verify_hct116_expression_identity_reader_round50.py'), Path('durable_public_store.py'),
        Path('transfer_budget.py'), Path('artifacts/download_allocations.json'),
        Path('artifacts/HCT116_extension_preflight_round49/preflight.json'),
        Path('artifacts/HCT116_extension_preflight_round49/audit.json'),
        Path('artifacts/HCT116_extension_preflight_round49/logical_record_ranges.parquet'),
        Path('../model_v4/artifacts/HCT116_fragments.json'),
        Path('../model_v4/artifacts/HCT116_genes.parquet'),
        Path('artifacts/HCT116_offset_mapping.npz'), Path('artifacts/frozen_metadata_index.json'),
        Path('artifacts/tasks/manifest.json'), Path('artifacts/benchmark_roles.json'),
        Path('artifacts/calibration_gate.json')]
    rt.save(PLAN, dict(scope='Metadata identity only; no expression data or model/response scores',
        revision=REVISION, expression_version=682, allocation_bytes=LIMIT,
        sha256={str(p):rt.sha(p) for p in paths}, expression_payloads_allowed=False,
        inventory_page_cap=CAP, max_pages_per_directory=10))


def main():
    plan = load(PLAN)
    for p, h in plan['sha256'].items():
        assert rt.sha(p) == h, p
    fixture = load('artifacts/HCT116_expression_identity_reader_fixture_round50.json')
    assert fixture['status'].startswith('PASS') and fixture['reader_sha256'] == rt.sha(__file__)
    assert fixture['verifier_sha256'] == rt.sha('verify_hct116_expression_identity_reader_round50.py')
    assert load('artifacts/download_allocations.json')['new_allocations'][KEY] == LIMIT
    budget = TransferBudget(JOURNAL, LIMIT)
    budget.bind('repository_revision', REVISION)
    budget.bind('plan_retry_sha256', rt.sha(PLAN))
    session = requests.Session(); session.trust_env = False
    store = None
    try:
        invs = {}
        for folder, recursive in [('expression.lance/_versions',False),
                                  ('expression.lance/data',False), ('genes.lance',True)]:
            invs[folder] = inventory(session,budget,folder,recursive)
        rt.save(ROOT/'inventory.json', invs)
        version_path = PREFIX+'expression.lance/_versions/'+str(2**64-1-682)+'.manifest'
        manifests = [e for e in invs['expression.lance/_versions']['entries'] if e['path'] == version_path]
        assert len(manifests) == 1
        genes = [e for e in invs['genes.lance']['entries'] if e['type']=='file' and '/_transactions/' not in e['path']]
        assert genes
        files = manifests+genes
        assert budget.total+sum(e['size']+64_000 for e in files) <= LIMIT
        store = PublicStore(budget)
        receipts = {}
        for e in files:
            assert '/expression.lance/data/' not in e['path'], 'Expression payload forbidden'
            dest = ROOT/'release'/e['path'][len(PREFIX):]
            h = fetch_file(e,store,budget,dest)
            receipts[e['path']] = dict(path=str(dest),size=e['size'],sha256=h,lfs=e.get('lfs'),git_oid=e['oid'])
            rt.save(ROOT/'progress.json', dict(files=receipts,charged_bytes=budget.total))
        # All Lance operations below are local. Dataset metadata does not need data-file bodies.
        ds = lance.dataset(str(ROOT/'release/expression.lance'), version=682)
        frags = []
        for f in ds.get_fragments():
            assert len(f.metadata.files) == 1 and not f.metadata.deletion_file
            frags.append(dict(id=f.fragment_id, physical_rows=f.metadata.physical_rows,
                              path=f.metadata.files[0].path, metadata=f.metadata.to_json()))
        old = load('../model_v4/artifacts/HCT116_fragments.json')
        assert [{k:f[k] for k in ['id','physical_rows','path']} for f in frags] == old
        data = {e['path']:e for e in invs['expression.lance/data']['entries'] if e['type']=='file'}
        for f in frags:
            key = PREFIX+'expression.lance/'+f['path']
            if key not in data:
                key = PREFIX+'expression.lance/data/'+f['path']
            assert key in data
            f['inventory_entry'] = data[key]
        gene_ds = lance.dataset(str(ROOT/'release/genes.lance'))
        axis = gene_ds.to_table().to_pandas()
        historical = pd.read_parquet('../model_v4/artifacts/HCT116_genes.parquet')
        pd.testing.assert_frame_equal(axis, historical, check_dtype=True)
        axis.to_parquet(ROOT/'genes.parquet',index=False)
        offsets = np.load('artifacts/HCT116_offset_mapping.npz')
        lengths = np.array([f['physical_rows'] for f in frags],np.int64)
        groups = offsets['first_cell']//5000
        prefix = np.r_[0,np.cumsum(lengths)]
        for g in range(len(offsets['logical_start'])):
            ix = np.flatnonzero(groups == g)
            if not len(ix):
                assert offsets['logical_start'][g] == -1 and offsets['block_length'][g] == 0
            else:
                assert np.array_equal(ix,np.arange(ix[0],ix[-1]+1))
                assert prefix[ix[0]] == offsets['logical_start'][g]
                assert lengths[ix].sum() == offsets['block_length'][g]
        rt.save(ROOT/'fragments.json',frags)
        rt.save(ROOT/'receipt.json',dict(status='PASS fresh immutable manifest/inventory/gene identity; independent audit pending',
            utc=datetime.now(timezone.utc).isoformat(), revision=REVISION, expression_version=ds.version,
            expression_schema=str(ds.schema), fragments=len(frags), expression_rows=int(lengths.sum()),
            inventory_data_files=len(data), referenced_data_files=len(frags), genes=len(axis),
            files=receipts, charged_bytes=budget.total, journal=str(JOURNAL), plan_sha256=rt.sha(PLAN),
            artifacts_sha256={str(p):rt.sha(p) for p in [ROOT/'inventory.json',ROOT/'fragments.json',ROOT/'genes.parquet']},
            historical_fragments_exact=True, block_lengths_logical_offsets_exact=True, gene_axis_exact=True,
            boundary_cell_ids_freshly_read=False, expression_payloads_read=False,
            encoded_request_plan_verified=False, source_activated=False, holdout_scored=False))
        print('PASS fresh expression metadata identity',len(frags),'fragments; charged',budget.total,flush=True)
    finally:
        session.close()
        if store: store.close()
        budget.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--freeze',action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()
