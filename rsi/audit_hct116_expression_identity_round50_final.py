"""Independent local metadata/inventory, accounting and conditional-range audit."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
from pathlib import Path
import urllib.parse

import runtime as rt
import lance
import numpy as np
import pandas as pd

ROOT=Path('artifacts/HCT116_expression_identity_round50')
PLAN=ROOT/'audit_plan_retry.json'


def load(p): return rt.json.loads(Path(p).read_text())


def freeze():
    assert not PLAN.exists()
    paths=[Path(__file__),ROOT/'plan_retry.json',ROOT/'receipt.json',ROOT/'inventory.json',
        ROOT/'fragments.json',ROOT/'genes.parquet',
        Path('artifacts/HCT116_expression_identity_reader_fixture_round50.json')]
    rt.save(PLAN,dict(scope='Independent metadata-only identity/conditional range/accounting audit',
        sha256={str(p):rt.sha(p) for p in paths}))


def main():
    assert not (ROOT/'audit.json').exists()
    for p,h in load(PLAN)['sha256'].items(): assert rt.sha(p)==h,p
    producer=load(ROOT/'plan_retry.json')
    for p,h in producer['sha256'].items(): assert rt.sha(p)==h,p
    receipt=load(ROOT/'receipt.json'); inv=load(ROOT/'inventory.json')
    for p,h in receipt['artifacts_sha256'].items(): assert rt.sha(p)==h,p
    assert receipt['plan_sha256']==rt.sha(ROOT/'plan_retry.json')
    steps=[rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    ends=[s for s in steps if s.get('name')=='acquire_HCT116_expression_identity_round50_retry' and s['event']=='end']
    assert ends[-1]['rc']==0
    events=[rt.json.loads(s) for s in Path(receipt['journal']).read_text().splitlines()]
    total=0;bindings={}
    for e in events:
        if e['event']=='bind':
            assert e['key'] not in bindings; bindings[e['key']]=e['value']
        else:
            assert e['event']=='charge' and e['bytes']>0
            total+=e['bytes'];assert e['total']==total<=15_000_000
            assert '/expression.lance/data/' not in e['reason'], 'Expression data GET/HEAD forbidden'
    assert total==receipt['charged_bytes']==6_125_800
    revision='598aa5442a4ee9a8f1a383a843ab5ce142904d04'
    assert bindings['repository_revision']==receipt['revision']==revision
    entries={};pages=0
    for folder,inventory in inv.items():
        assert inventory['complete'] and inventory['revision']==revision
        found=[];seen=set();previous=None
        for page in range(inventory['pages']):
            p=ROOT/(folder.replace('/','_')+f'_page{page}.json')
            assert bindings['inventory_page:'+str(p)]==rt.sha(p)
            x=load(p);body=bytes.fromhex(x['body_hex'])
            assert hashlib.sha256(body).hexdigest()==x['sha256']
            url=urllib.parse.urlsplit(x['url'])
            assert url.scheme=='https' and url.netloc=='huggingface.co'
            assert urllib.parse.unquote(url.path)==f'/api/datasets/slaf-project/X-Atlas-Orion/tree/{revision}/data/HCT116/{folder}'
            assert x['url'] not in seen;seen.add(x['url'])
            if previous is not None: assert x['url']==previous
            previous=x['next_url']
            found.extend(rt.json.loads(body));pages+=1
        assert previous is None and found==inventory['entries']
        assert len({e['path'] for e in found})==len(found)
        for e in found:
            assert e['path'].startswith('data/HCT116/'+folder+'/') and '..' not in e['path']
            assert e['path'] not in entries;entries[e['path']]=e
    assert pages==5
    for key,r in receipt['files'].items():
        assert '/expression.lance/data/' not in key
        e=entries[key];data=Path(r['path']).read_bytes()
        assert len(data)==e['size']==r['size']
        assert hashlib.sha256(data).hexdigest()==r['sha256']
        if e.get('lfs'): assert r['sha256']==e['lfs']['oid']
        else: assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['oid']
    ds=lance.dataset(str(ROOT/'release/expression.lance'),version=682)
    assert ds.version==receipt['expression_version']==682
    assert ds.schema.names==['cell_integer_id','gene_integer_id','value']
    old=load('../model_v4/artifacts/HCT116_fragments.json');fresh=load(ROOT/'fragments.json')
    actual=ds.get_fragments();assert len(actual)==len(old)==len(fresh)==2078
    prefix=[0];file_paths=set();versions=Counter()
    for i,(f,o,saved) in enumerate(zip(actual,old,fresh)):
        m=f.metadata
        assert f.fragment_id==o['id']==saved['id']==i
        assert m.physical_rows==o['physical_rows']==saved['physical_rows']
        assert len(m.files)==1 and not m.deletion_file and not m.overlays
        file=m.files[0]
        assert file.path==o['path']==saved['path']
        assert saved['metadata']==m.to_json()
        key='data/HCT116/expression.lance/data/'+file.path
        assert key not in file_paths;file_paths.add(key)
        assert entries[key]==saved['inventory_entry']
        assert file.file_size_bytes==entries[key]['size']
        if entries[key].get('lfs'):
            assert entries[key]['lfs']['size']==entries[key]['size']
            assert len(entries[key]['lfs']['oid'])==64
        else:
            assert len(entries[key]['oid'])==40
        versions[file.file_major_version,file.file_minor_version]+=1
        prefix.append(prefix[-1]+m.physical_rows)
    assert prefix[-1]==receipt['expression_rows']==17_395_402_897
    offsets=np.load('artifacts/HCT116_offset_mapping.npz')
    # Rebuild group offsets with scalar counters, rather than producer vector sums.
    starts={};counts=Counter();previous=-1
    for i,f in enumerate(fresh):
        first=int(offsets['first_cell'][i]);last=int(offsets['last_cell'][i]);group=first//5000
        assert last//5000==group and first<=last
        if group!=previous:
            assert group not in starts
            starts[group]=prefix[i];previous=group
        counts[group]+=f['physical_rows']
    assert int(offsets['version'])==682
    for group in range(682):
        assert int(offsets['logical_start'][group])==starts.get(group,-1)
        assert int(offsets['block_length'][group])==counts[group]
    assert set(range(682))-set(starts)=={374}
    spans=pd.read_parquet('artifacts/HCT116_extension_preflight_round49/logical_record_ranges.parquet')
    assert len(spans)==41_791
    record_sum=0;touched=set();by_cell=Counter()
    selection=pd.read_parquet('artifacts/HCT116_extension_preflight_round49/selection.parquet')
    expected_records={int(r.cell_integer_id):int(r.n_genes_by_counts) for r in selection.itertuples(index=False)}
    expected_start={int(r.cell_integer_id):int(r.cell_start_index)-int(offsets['original_start'][int(r.cell_integer_id)//5000])+int(offsets['logical_start'][int(r.cell_integer_id)//5000]) for r in selection.itertuples(index=False)}
    for r in spans.itertuples(index=False):
        f=fresh[r.fragment_id]
        assert r.path==f['path'] and 0<=r.record_start<r.record_stop<=f['physical_rows']
        assert prefix[r.fragment_id]+r.record_start==expected_start[r.cell_integer_id]+by_cell[r.cell_integer_id]
        n=r.record_stop-r.record_start;by_cell[r.cell_integer_id]+=n;record_sum+=n;touched.add(r.fragment_id)
    assert dict(by_cell)==expected_records and record_sum==215_800_779 and len(touched)==2018
    axis=lance.dataset(str(ROOT/'release/genes.lance')).to_table().to_pandas()
    oldaxis=pd.read_parquet('../model_v4/artifacts/HCT116_genes.parquet')
    pd.testing.assert_frame_equal(axis,oldaxis)
    pd.testing.assert_frame_equal(axis,pd.read_parquet(ROOT/'genes.parquet'))
    assert axis.gene_integer_id.tolist()==list(range(len(axis)))
    frozen=load('artifacts/frozen_metadata_index.json')
    assert rt.sha(frozen['archive'])==frozen['archive_sha256']
    for p,e in frozen['entries'].items(): assert rt.sha(p)==e['sha256'],p
    gate=load('artifacts/calibration_gate.json'); assert gate['status']=='CLOSED' and not gate['phase2_allowed']
    all_data={e['path'] for e in inv['expression.lance/data']['entries'] if e['type']=='file'}
    rt.save(ROOT/'audit.json',dict(status='PASS independent immutable metadata/inventory/accounting audit',
        utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=rt.sha(__file__),
        audit_plan_sha256=rt.sha(PLAN),receipt_sha256=rt.sha(ROOT/'receipt.json'),
        revision=revision,expression_version=682,fragments=2078,expression_rows=prefix[-1],
        exact_historical_fragment_paths_lengths=True,manifest_file_sizes_equal_inventory=True,
        block_offsets_conditionally_replayed=True,boundary_cell_ids_freshly_read=False,
        logical_spans_verified=len(spans),logical_records_verified=record_sum,touched_fragments=len(touched),
        data_inventory_files=len(all_data),unreferenced_files=sorted(all_data-file_paths),
        manifest_referenced_file_bytes=sum(entries[k]['size'] for k in file_paths),
        lance_file_versions={str(k):v for k,v in versions.items()},
        gene_rows=len(axis),unique_gene_ids=int(axis.gene_id.nunique()),
        frozen_members=len(frozen['entries']),charged_bytes=total,inventory_pages=pages,
        expression_payloads_read=False,encoded_request_plan_verified=False,source_activated=False,
        holdout_scored=False,phase2_allowed=False))
    print('PASS independent manifest/inventory/gene/accounting audit',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');a=p.parse_args()
    freeze() if a.freeze else main()
