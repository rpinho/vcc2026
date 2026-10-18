"""Frozen payload-free metadata acquisition; three bounded network lanes.

Every fetched range is one footer, column metadata region or certified mini-
block index. Compressed values, scalar/dictionary buffers and broad tail reads
are forbidden. ETags bind range requests; partial-file SHA256 is not asserted.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import argparse,os,threading
import runtime as rt
import pandas as pd
from durable_public_store import PublicStore
from transfer_budget import TransferBudget
import lance_encoded_plan_round51 as lp

ROOT=Path('artifacts/HCT116_encoded_metadata_round51')
PLAN=ROOT/'plan.json'
KEY='HCT116_encoded_metadata_round51'
LIMIT=250_000_000
REV='598aa5442a4ee9a8f1a383a843ab5ce142904d04'
FRAGMENTS=Path('artifacts/HCT116_expression_identity_round50/fragments.json')
LOGICAL=Path('artifacts/HCT116_extension_preflight_round49/logical_record_ranges.parquet')
JOURNAL=Path('artifacts/'+KEY+'_charges.jsonl')

def load(p):return rt.json.loads(Path(p).read_text())

class RangeCache:
    def __init__(self,root,budget,store):
        self.root=Path(root);self.budget=budget;self.store=store;self.allowed=set()
    def permit(self,fragment,ranges):
        self.allowed.update((fragment,int(a),int(b)) for a,b in ranges)
    def read(self,fragment,entry,start,stop):
        assert (fragment,start,stop) in self.allowed, 'Range outside certified metadata whitelist'
        size=entry['size'];assert 0<=start<stop<=size
        dest=self.root/str(fragment)/f'{start}-{stop}.bin'
        cachekey='cache:'+str(dest)
        if dest.exists():
            assert cachekey in self.budget.bindings and self.budget.bindings[cachekey]==rt.sha(dest), 'Cache identity mismatch'
            assert dest.stat().st_size==stop-start
            return dest.read_bytes()
        key=entry['path']
        self.budget.bind('inventory:'+key,entry)
        url,source=self.store.resolve(key)
        assert source['size']==size
        assert url.startswith('https://'), 'Public origin must retain TLS'
        self.budget.charge(stop-start,'Exact metadata range attempt '+key+f' [{start},{stop})')
        rt.guard()
        with self.store.http.get(url,headers={'Accept-Encoding':'identity','If-Match':source['etag'],
            'Range':f'bytes={start}-{stop-1}'},stream=True,allow_redirects=False,timeout=(10,90)) as response:
            assert response.status_code==206, 'Range ignored or public range failed'
            assert int(response.headers['Content-Length'])==stop-start
            assert response.headers['Content-Range']==f'bytes {start}-{stop-1}/{size}'
            assert response.headers['ETag']==source['etag']
            assert response.headers.get('Content-Encoding','identity')=='identity'
            body=response.raw.read(stop-start+1)
            assert len(body)==stop-start, 'Truncated/extra encoded metadata'
        dest.parent.mkdir(parents=True,exist_ok=True)
        temp=dest.with_suffix('.tmp')
        with temp.open('wb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        # Commit the durable binding first. A crash before rename spends budget
        # but cannot publish a checkpoint without its verified byte identity.
        self.budget.bind(cachekey,rt.hashlib.sha256(body).hexdigest());temp.replace(dest)
        return body


def describe_file(fragment,logical,cache,root):
    fid=fragment['id'];entry=fragment['inventory_entry'];size=entry['size']
    cache.permit(fid,[[size-40,size]])
    foot=lp.footer(cache.read(fid,entry,size-40,size),size)
    cache.permit(fid,[[foot['column_start'],size-40]])
    meta=cache.read(fid,entry,foot['column_start'],size-40)
    pages=lp.pages_from_metadata(meta,foot,fragment['physical_rows'])
    # Global buffer 0 is certified by the GBO table as the file schema. Read
    # exactly that buffer to permit independent official metadata replay.
    assert foot['num_global_buffers']==1, 'Extra global buffers require a new plan'
    schema_start,schema_size=lp.struct.unpack_from('<QQ',meta,foot['gbo_start']-foot['column_start'])
    assert 0<schema_size<=64_000 and 0<=schema_start<schema_start+schema_size<=foot['column_start']
    assert all(b<=schema_start for p in pages for a,b in p['buffers'])
    cache.permit(fid,[[schema_start,schema_start+schema_size]])
    cache.read(fid,entry,schema_start,schema_start+schema_size)
    indexes=[p['index_range'] for p in pages if p['layout']==1 and any(a<p['row_stop'] and b>p['row_start'] for a,b in logical)]
    cache.permit(fid,indexes)
    payload=lp.plan_payload(pages,logical,lambda p:cache.read(fid,entry,*p['index_range']))
    index_bytes=sum(b-a for a,b in lp.merge(indexes))
    result=dict(fragment_id=fid,path=entry['path'],size=size,physical_rows=fragment['physical_rows'],
        footer=foot,pages=pages,schema_range=[schema_start,schema_start+schema_size],
        metadata_ranges=lp.merge([[schema_start,schema_start+schema_size],[size-40,size],[foot['column_start'],size-40]]+indexes),
        metadata_unique_bytes=size-foot['column_start']+schema_size+index_bytes,
        index_bytes=index_bytes,payload=payload,expression_values_fetched=False,
        partial_file_full_digest_verified=False,range_identity='Immutable repository path/inventory size plus consistent resolved ETag; full-file raw LFS/Git digest not checked on partial reads')
    rt.save(Path(root)/str(fid)/'plan.json',result)
    return result


def allocate():
    ap=Path('artifacts/download_allocations.json');a=load(ap)
    if KEY in a['new_allocations']:
        assert a['new_allocations'][KEY]==LIMIT;return
    events=[rt.json.loads(s) for s in Path('artifacts/H1_rich_charges.jsonl').read_text().splitlines()]
    total=0
    for e in events:
        if e['event']=='charge':total+=e['bytes'];assert e['total']==total
    assert total==4_700_011_240
    steps=[rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert [s for s in steps if s['name']=='acquire_h1_rich'][-1]['event']=='end'
    assert [s for s in steps if s['name']=='acquire_h1_rich'][-1]['rc']==0
    assert a['new_allocations']['H1_rich']==5_500_000_000
    before=Path('artifacts/download_allocations_before_encoded_metadata_round51.json');assert not before.exists();before.write_bytes(ap.read_bytes())
    a['new_allocations']['H1_rich']-=LIMIT;a['new_allocations'][KEY]=LIMIT
    assert a['new_allocations']['H1_rich']>=total
    a['reader_accounting']['H1_rich']['effective_ceiling_bytes_after_round51']=a['new_allocations']['H1_rich']
    a['reader_accounting']['H1_rich']['reallocation_policy']='Completed reader; original journal historical5.5GB binding preserved. Do not rerun old hard-coded ceiling without reconciling new5.25GB effective ceiling.'
    a['reader_accounting'][KEY]=dict(ceiling_bytes=LIMIT,journal=str(JOURNAL),policy='Immutable v682 metadata only:64KB HEAD resolution per file, exact40-byte footer, certified column metadata and miniblock index ranges only. Durable pre-I/O attempts, strict206/size/ETag/Content-Range, hash-bound cache reuse. No payload fallback or implicit retry. Three bounded lanes.')
    a['reallocations'].append(dict(round=51,utc=datetime.now(timezone.utc).isoformat(),from_reader='H1_rich',to_reader=KEY,bytes=LIMIT,before_snapshot=str(before),before_sha256=rt.sha(before),from_charged_bytes=total,from_journal_sha256=rt.sha('artifacts/H1_rich_charges.jsonl'),reason='Completed H1 reader has799,988,760 unused bytes. Reassign250MB to exact metadata/index plan for2018 files, including129.152MB initial HEAD charges; no encoded values.'))
    assert sum(a['new_allocations'].values())==a['total_allocated_upper_bound_bytes']<=a['limit_bytes']
    rt.save(ap,a)


def freeze():
    ROOT.mkdir(exist_ok=True);assert not PLAN.exists()
    assert load('artifacts/HCT116_expression_identity_round50/audit.json')['status'].startswith('PASS')
    allocate()
    paths=[Path(__file__),Path('lance_encoded_plan_round51.py'),Path('verify_lance_encoded_reader_round51.py'),
        Path('audit_hct116_encoded_plan_round51.py'),Path('durable_public_store.py'),Path('transfer_budget.py'),FRAGMENTS,LOGICAL,
        Path('artifacts/HCT116_expression_identity_round50/audit.json'),Path('artifacts/HCT116_expression_identity_round50/receipt.json'),
        Path('artifacts/HCT116_extension_preflight_round49/selection.parquet'),Path('artifacts/download_allocations.json'),
        Path('artifacts/download_allocations_before_encoded_metadata_round51.json'),Path('artifacts/H1_rich_charges.jsonl'),
        Path('artifacts/frozen_metadata_index.json'),Path('artifacts/tasks/manifest.json'),Path('artifacts/benchmark_roles.json'),
        Path('artifacts/calibration_gate.json'),Path('artifacts/cd4_banks/HCT116.npz'),
        Path('artifacts/Lance_format_round51/receipt.json'),Path('pin_lance_format_round51.py'),
        Path('vendor/deps/lance/lance.abi3.so'),Path('vendor/deps/lance/file.py'),
        Path('requirements-frozen.txt')]
    paths+=sorted(Path('artifacts/scores/H1_rich_search').glob('*_s42.json'))
    paths+=sorted(Path('artifacts/scores/CD4_depth192_search').glob('*_s42.json'))
    paths+=sorted(Path('artifacts/Lance_format_round51').glob('*.rs'))+sorted(Path('artifacts/Lance_format_round51').glob('*.proto'))
    rt.save(PLAN,dict(utc=datetime.now(timezone.utc).isoformat(),scope='Metadata/index only for frozen round49 logical spans; no expression values, source activation, model scores or target changes',revision=REV,allocation_bytes=LIMIT,lanes=3,coalescing_gap_bytes=0,retries=0,metadata_per_file_max=256_000,head_reserve_per_file=64_000,
        accounting='Pre-I/O amount per attempt; failures spend budget. Completion hashes fsynced before atomic rename. Reuse only journal-bound exact ranges. Do not reset journal or expand scope on restart.',sha256={str(p):rt.sha(p) for p in paths}))
    print('FROZEN',PLAN,flush=True)


def main():
    plan=load(PLAN)
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p
    fixture=load(ROOT/'reader_fixture.json');assert fixture['status'].startswith('PASS')
    for p,h in fixture['sha256'].items():assert rt.sha(p)==h,p
    assert not (ROOT/'receipt.json').exists(), 'Preserve completed evidence'
    frame=pd.read_parquet(LOGICAL)
    jobs=[(f,frame.loc[frame.fragment_id==f['id'],['record_start','record_stop']].to_numpy().tolist()) for f in load(FRAGMENTS) if (frame.fragment_id==f['id']).any()]
    assert len(jobs)==2018
    budget=TransferBudget(JOURNAL,LIMIT);budget.bind('repository_revision',REV);budget.bind('plan_sha256',rt.sha(PLAN))
    stop=threading.Event();results={};lock=threading.Lock()
    lanes=[jobs[i::3] for i in range(3)]
    def worker(jobs):
        store=PublicStore(budget);cache=RangeCache(ROOT,budget,store)
        try:
            for fragment,logical in jobs:
                if stop.is_set():return
                result=describe_file(fragment,logical,cache,ROOT)
                with lock:
                    results[fragment['id']]=result
                    if len(results)%25==0:print('metadata complete',len(results),'/',len(frame.fragment_id.unique()),'charged',budget.total,flush=True)
        except BaseException:stop.set();raise
        finally:store.close()
    try:
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures=[executor.submit(worker,j) for j in lanes]
            for future in futures:future.result()
        assert len(results)==2018
        ordered=[results[k] for k in sorted(results)]
        summary=dict(status='PASS metadata acquisition; independent cost audit pending',utc=datetime.now(timezone.utc).isoformat(),files=len(results),logical_spans=len(frame),logical_rows=sum(x['payload']['logical_rows'] for x in ordered),
            metadata_unique_bytes=sum(x['metadata_unique_bytes'] for x in ordered),encoded_value_payload_bytes=sum(x['payload']['payload_bytes'] for x in ordered),schema_bytes_cached=sum(x['schema_range'][1]-x['schema_range'][0] for x in ordered),
            encoded_value_requests=sum(len(x['payload']['payload_ranges']) for x in ordered),future_payload_head_reserve_bytes=64_000*len(ordered),future_retry_policy='No implicit retries. Each full attempt must be charged again. A later isolated reader must allocate value+schema+HEAD cost and its explicit retry envelope before I/O.',
            charged_bytes=budget.total,journal=str(JOURNAL),plan_sha256=rt.sha(PLAN),
            file_plan_sha256={str(ROOT/str(k)/'plan.json'):rt.sha(ROOT/str(k)/'plan.json') for k in sorted(results)},
            expression_values_fetched=False,source_activated=False,new_biological_scores=False,holdout_scored=False,full_partial_file_digest_verified=False)
        rt.save(ROOT/'receipt.json',summary)
        print('METADATA COMPLETE',summary['encoded_value_payload_bytes'],'future encoded value bytes',flush=True)
    finally:budget.close()
    # Independent frozen auditor uses only saved metadata, no network or biology.
    import audit_hct116_encoded_plan_round51 as audit
    audit.main()

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');args=parser.parse_args()
    freeze() if args.freeze else main()
