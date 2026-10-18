"""Exact-range HCT116 staging only; no source-bank or scoring activation.

Reuse independently audited metadata. Fetch only frozen miniblocks, decode with
the pinned official Lance reader, and verify every selected record's cell ID.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import os
import threading
import tempfile

import runtime as rt
import numpy as np
import pandas as pd
from lance.file import LanceFileReader
from durable_public_store import PublicStore
from transfer_budget import TransferBudget
import lance_encoded_plan_round51 as lp

ROOT = Path('artifacts/HCT116_payload_round52')
META = Path('artifacts/HCT116_encoded_metadata_round51')
PREFLIGHT = Path('artifacts/HCT116_extension_preflight_round49')
PLAN = ROOT/'plan.json'
KEY = 'HCT116_payload_round52'
LIMIT = 1_200_000_000
JOURNAL = Path('artifacts/'+KEY+'_charges.jsonl')
REV = '598aa5442a4ee9a8f1a383a843ab5ce142904d04'
ALLOC = Path('artifacts/download_allocations.json')
BEFORE = Path('artifacts/download_allocations_before_payload_round52.json')


def load(p):
    return rt.json.loads(Path(p).read_text())


def journal(path):
    total = 0
    bindings = {}
    for s in Path(path).read_text().splitlines():
        e = rt.json.loads(s)
        if e['event'] == 'charge':
            assert e['bytes'] > 0
            total += e['bytes']
            assert total == e['total']
        else:
            assert e['event'] == 'bind' and e['key'] not in bindings
            bindings[e['key']] = e['value']
    return total, bindings


def durable_bytes(path, body):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    with tmp.open('wb') as out:
        out.write(body)
        out.flush()
        os.fsync(out.fileno())
    return tmp


def sync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class PayloadCache:
    def __init__(self, root, budget, store):
        self.root = Path(root)
        self.budget = budget
        self.store = store
        self.allowed = set()

    def permit(self, fid, ranges):
        self.allowed.update((int(fid), int(a), int(b)) for a, b in ranges)

    def read(self, fid, entry, a, b):
        assert (fid, a, b) in self.allowed, 'Unplanned payload range'
        assert 0 <= a < b <= entry['size']
        dest = self.root/str(fid)/f'{a}-{b}.bin'
        key = 'payload_cache:'+str(dest)
        if dest.exists():
            assert dest.stat().st_size == b-a
            assert self.budget.bindings.get(key) == rt.sha(dest), 'Corrupt cache'
            return dest.read_bytes()
        self.budget.bind('inventory:'+entry['path'], entry)
        url, source = self.store.resolve(entry['path'])
        assert source['size'] == entry['size'] and url.startswith('https://')
        # resolve() must match the metadata-stage ETag prebound by the caller.
        assert self.budget.bindings['file:'+entry['path']] == source
        rt.guard()
        self.budget.charge(b-a, 'Exact payload range attempt '+entry['path']+f' [{a},{b})')
        with self.store.http.get(url, headers={'Accept-Encoding':'identity',
                'If-Match':source['etag'], 'Range':f'bytes={a}-{b-1}'},
                stream=True, allow_redirects=False, timeout=(10,90)) as response:
            assert response.status_code == 206
            assert int(response.headers['Content-Length']) == b-a
            assert response.headers['Content-Range'] == f'bytes {a}-{b-1}/{entry["size"]}'
            assert response.headers['ETag'] == source['etag']
            assert response.headers.get('Content-Encoding','identity') == 'identity'
            body = response.raw.read(b-a+1)
            assert len(body) == b-a, 'Truncated or oversized payload'
        tmp = durable_bytes(dest, body)
        self.budget.bind(key, rt.hashlib.sha256(body).hexdigest())
        tmp.replace(dest)
        sync_dir(dest.parent)
        return body


def metadata_bytes(fid, a, b, bindings):
    path = META/str(fid)/f'{a}-{b}.bin'
    assert path.stat().st_size == b-a
    assert bindings['cache:'+str(path)] == rt.sha(path)
    return path.read_bytes()


def decode(sparse, spans, first_cell, last_cell, physical_rows, n_genes):
    official = LanceFileReader(str(sparse))
    assert official.metadata().num_rows == physical_rows
    # Boundary reads request ONLY column 0: other boundary columns were not paid.
    cell_reader = LanceFileReader(str(sparse), columns=['cell_integer_id'])
    for position, expected in [(0,first_cell), (physical_rows-1,last_cell)]:
        assert cell_reader.read_range(position,1).to_table()['cell_integer_id'][0].as_py() == expected
    parts = {k:[] for k in ['cell_integer_id','gene_integer_id','value']}
    for r in spans.itertuples(index=False):
        n = int(r.record_stop-r.record_start)
        table = official.read_range(int(r.record_start),n).to_table()
        arrays = {k:table[k].to_numpy() for k in parts}
        assert len(table) == n
        assert np.all(arrays['cell_integer_id'] == r.cell_integer_id), 'Cell offset identity mismatch'
        assert arrays['cell_integer_id'].dtype == np.dtype('uint32')
        assert arrays['gene_integer_id'].dtype == arrays['value'].dtype == np.dtype('uint16')
        assert np.all(arrays['gene_integer_id'] < n_genes)
        assert np.all(arrays['value'] > 0), 'Sparse expression must have positive integer counts'
        for k in parts:
            parts[k].append(arrays[k])
    return {k:np.concatenate(v) for k,v in parts.items()}


def build_file(fid, fileplan, entry, spans, cache, metadata_bindings, first, last, root=ROOT):
    root = Path(root)
    dest = root/'decoded'/f'{fid}.npz'
    receipt = root/'decoded'/f'{fid}.json'
    key = 'decoded:'+str(dest)
    if receipt.exists():
        row = load(receipt)
        assert row['sha256'] == cache.budget.bindings.get(key) == rt.sha(dest)
        assert row['fileplan_sha256'] == rt.hashlib.sha256(rt.json.dumps(fileplan, sort_keys=True).encode()).hexdigest()
        return row
    cache.permit(fid, fileplan['payload_ranges'])
    with tempfile.TemporaryDirectory(prefix=f'payload_{fid}_', dir='artifacts/tmp') as tmpdir:
        sparse = Path(tmpdir)/'sparse.lance'
        with sparse.open('wb') as out:
            out.truncate(entry['size'])
            for p in sorted((META/str(fid)).glob('*.bin')):
                a,b = map(int,p.stem.split('-'))
                out.seek(a)
                out.write(metadata_bytes(fid,a,b,metadata_bindings))
            for a,b in fileplan['payload_ranges']:
                out.seek(a)
                out.write(cache.read(fid,entry,a,b))
        arrays = decode(sparse,spans,first,last,fileplan['physical_rows'],38606)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **arrays,
        record_start=spans.record_start.to_numpy(np.int64),
        record_stop=spans.record_stop.to_numpy(np.int64),
        span_cell_id=spans.cell_integer_id.to_numpy(np.int64))
    with tmp.open('rb') as f:
        os.fsync(f.fileno())
    digest = rt.sha(tmp)
    cache.budget.bind(key,digest)
    tmp.replace(dest)
    sync_dir(dest.parent)
    row = dict(fragment_id=fid,sha256=digest,records=len(arrays['value']),
        spans=len(spans),first_cell_verified=int(first),last_cell_verified=int(last),
        all_selected_cell_ids_verified=True,
        fileplan_sha256=rt.hashlib.sha256(rt.json.dumps(fileplan,sort_keys=True).encode()).hexdigest())
    rt.save(receipt,row)
    return row


def freeze():
    ROOT.mkdir(exist_ok=True)
    assert not PLAN.exists() and not BEFORE.exists()
    old = load(META/'plan.json')
    for p,h in old['sha256'].items():
        assert rt.sha(p) == h,p
    for name in ['audit.json','scope_audit.json']:
        assert load(META/name)['status'].startswith('PASS')
        assert load(META/name)['receipt_sha256'] == rt.sha(META/'receipt.json')
    steps = [rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    end = [e for e in steps if e['name']=='acquire_HCT116_encoded_metadata_round51'][-1]
    assert end['event']=='end' and end['rc']==0
    _, bindings = journal('artifacts/HCT116_encoded_metadata_round51_charges.jsonl')
    receipt = load(META/'receipt.json')
    files = {}
    extra = 0
    for path,h in receipt['file_plan_sha256'].items():
        assert rt.sha(path)==h
        saved = load(path)
        fid = saved['fragment_id']
        edges = []
        for p in saved['pages']:
            assert p['layout']==1
            if p['column']!=0:
                continue
            for rowlo,rowhi,a,b in lp.chunks(p,metadata_bytes(fid,*p['index_range'],bindings)):
                if rowlo==0 or rowhi==saved['physical_rows']:
                    edges.append([a,b])
        ranges = lp.merge(saved['payload']['payload_ranges']+edges)
        nbytes = sum(b-a for a,b in ranges)
        extra += nbytes-saved['payload']['payload_bytes']
        files[str(fid)] = dict(payload_ranges=ranges,payload_bytes=nbytes,
            boundary_cell_id_ranges=lp.merge(edges),physical_rows=saved['physical_rows'],
            metadata_plan=path,metadata_plan_sha256=h)
    payload = sum(f['payload_bytes'] for f in files.values())
    head = 64000*len(files)
    assert payload+head < LIMIT
    rt.save(ROOT/'files.json',files)
    # Only completed-reader unused reserves move. Their original journals stay byte exact.
    a = load(ALLOC)
    assert sum(a['new_allocations'].values())==a['total_allocated_upper_bound_bytes']<=15_000_000_000
    donors = [('H1_rich','acquire_h1_rich','artifacts/H1_rich_charges.jsonl'),
        ('HEK293T','acquire_orion_HEK293T','artifacts/HEK293T_charges.jsonl'),
        ('CD4_depth192','acquire_CD4_depth192_round17','artifacts/CD4_depth192_charges.jsonl')]
    reconciled = []
    for donor,step,path in donors:
        last = [e for e in steps if e['name']==step][-1]
        assert last['event']=='end' and last['rc']==0
        total,_ = journal(path)
        assert total <= a['new_allocations'][donor]-400_000_000
        reconciled.append(dict(reader=donor,journal=path,journal_sha256=rt.sha(path),
            charged_bytes=total,before_ceiling=a['new_allocations'][donor],
            after_ceiling=a['new_allocations'][donor]-400_000_000))
    BEFORE.write_bytes(ALLOC.read_bytes())
    for d in reconciled:
        a['new_allocations'][d['reader']]=d['after_ceiling']
        a['reader_accounting'][d['reader']]['effective_ceiling_bytes_after_round52']=d['after_ceiling']
        a['reader_accounting'][d['reader']]['reallocation_policy']='Completed reader. Historical journal and hard-coded limit are retained as history only; no restart under an obsolete ceiling. Effective current allocation is new_allocations.'
    a['new_allocations'][KEY]=LIMIT
    a['reader_accounting'][KEY]=dict(ceiling_bytes=LIMIT,journal=str(JOURNAL),
        policy='Frozen exact miniblocks including cell-ID boundaries; fresh64KB HEAD resolutions, original-stage ETags, pre-I/O charges and hash-bound durable per-range/decoded checkpoints; no implicit retries or coalescing gaps. Three lanes. No source activation.')
    a['reallocations'].append(dict(round=52,utc=datetime.now(timezone.utc).isoformat(),
        from_readers=reconciled,to_reader=KEY,bytes=LIMIT,before_snapshot=str(BEFORE),
        before_sha256=rt.sha(BEFORE),reason='Staging unchanged CD4 panel72 missing HCT116 sources plus original controls; no bank activation. One complete exact request plan and remaining explicitly capped failed-attempt/restart envelope.'))
    assert sum(a['new_allocations'].values())==a['total_allocated_upper_bound_bytes']
    rt.save(ALLOC,a)
    paths = [Path(__file__),Path('verify_hct116_payload_round52.py'),Path('audit_hct116_payload_round52.py'),
        ROOT/'files.json',META/'plan.json',META/'receipt.json',META/'audit.json',META/'scope_audit.json',
        Path('artifacts/HCT116_encoded_metadata_round51_charges.jsonl'),BEFORE,ALLOC,
        Path('artifacts/HCT116_expression_identity_round50/fragments.json'),
        PREFLIGHT/'selection.parquet',PREFLIGHT/'logical_record_ranges.parquet',PREFLIGHT/'audit.json',
        Path('artifacts/HCT116_offset_mapping.npz'),Path('artifacts/cd4_banks/HCT116.npz'),
        Path('artifacts/frozen_metadata_index.json'),Path('artifacts/tasks/manifest.json'),
        Path('artifacts/benchmark_roles.json'),Path('artifacts/calibration_gate.json'),
        Path('durable_public_store.py'),Path('transfer_budget.py'),Path(lp.__file__),
        Path('vendor/deps/lance/lance.abi3.so'),Path('vendor/deps/lance/file.py')]
    paths += [Path(d['journal']) for d in reconciled]
    paths += list(Path('artifacts/scores/H1_rich_search').glob('*_s42.json'))
    paths += list(Path('artifacts/scores/CD4_depth192_search').glob('*_s42.json'))
    rt.save(PLAN,dict(utc=datetime.now(timezone.utc).isoformat(),revision=REV,files=len(files),
        allocation_bytes=LIMIT,encoded_payload_bytes=payload,boundary_extra_bytes=extra,
        fresh_head_reserve_bytes=head,first_attempt_bytes=payload+head,
        additional_failure_restart_envelope_bytes=LIMIT-payload-head,
        no_implicit_retries=True,lanes=3,coalescing_gap_bytes=0,
        scope='Staged raw records only; no bank build, centering, recipe/protocol change, model scoring or HOLDOUT access',
        partial_file_full_digest_verified=False,
        sha256={str(p):rt.sha(p) for p in paths}))
    print('FROZEN',payload,'payload bytes;',extra,'extra boundary bytes;',LIMIT-payload-head,'failure reserve',flush=True)


def main():
    plan = load(PLAN)
    for p,h in plan['sha256'].items():
        assert rt.sha(p)==h,p
    for path in [ROOT/'fixture.json',ROOT/'preflight_audit.json']:
        evidence=load(path)
        assert evidence['status'].startswith('PASS') and evidence['plan_sha256']==rt.sha(PLAN)
    assert not (ROOT/'receipt.json').exists(), 'Preserve completed acquisition'
    files = load(ROOT/'files.json')
    frags = {f['id']:f for f in load('artifacts/HCT116_expression_identity_round50/fragments.json')}
    spans = pd.read_parquet(PREFLIGHT/'logical_record_ranges.parquet')
    by_file = {int(fid):rows for fid,rows in spans.groupby('fragment_id',sort=True)}
    offsets = np.load('artifacts/HCT116_offset_mapping.npz')
    _, metadata_bindings = journal('artifacts/HCT116_encoded_metadata_round51_charges.jsonl')
    budget = TransferBudget(JOURNAL,LIMIT)
    budget.bind('repository_revision',REV)
    budget.bind('plan_sha256',rt.sha(PLAN))
    # Every fresh file resolution must replay the immutable metadata-stage ETag.
    for fid in by_file:
        path = frags[fid]['inventory_entry']['path']
        budget.bind('file:'+path,metadata_bindings['file:'+path])
    stop = threading.Event()
    completed = {}
    lock = threading.Lock()
    def worker(ids):
        store = PublicStore(budget)
        cache = PayloadCache(ROOT/'ranges',budget,store)
        try:
            for fid in ids:
                if stop.is_set():
                    return
                row=build_file(fid,files[str(fid)],frags[fid]['inventory_entry'],by_file[fid],
                    cache,metadata_bindings,offsets['first_cell'][fid],offsets['last_cell'][fid])
                with lock:
                    completed[str(fid)]=row
                    if len(completed)%25==0:
                        print('staged',len(completed),'/2018 charged',budget.total,flush=True)
        except BaseException:
            stop.set()
            raise
        finally:
            store.close()
    try:
        ids=sorted(by_file)
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures=[executor.submit(worker,ids[i::3]) for i in range(3)]
            for future in futures:
                future.result()
        assert len(completed)==2018
        assert sum(r['records'] for r in completed.values())==215_800_779
        rt.save(ROOT/'receipt.json',dict(status='STAGED; independent final audit pending',
            utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(PLAN),
            files=len(completed),logical_records=215_800_779,charged_bytes=budget.total,
            journal=str(JOURNAL),decoded=completed,source_activated=False,
            moments_verified=False,holdout_scored=False,full_partial_file_digest_verified=False))
    finally:
        budget.close()
    import audit_hct116_payload_round52 as audit
    audit.complete()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    freeze() if args.freeze else main()
