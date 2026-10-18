"""Offline actual-payload-reader failure, identity and exact decoding fixtures."""
from pathlib import Path
import tempfile
import runtime as rt
import numpy as np
import pandas as pd
import pyarrow as pa
from lance.file import LanceFileWriter
from transfer_budget import TransferBudget
import acquire_hct116_payload_round52 as reader
import acquire_hct116_encoded_metadata_round51 as metadata
import lance_encoded_plan_round51 as lp
from verify_lance_encoded_reader_round51 import Store


def main():
    assert not (reader.ROOT/'fixture.json').exists()
    checks=[]
    with tempfile.TemporaryDirectory(prefix='payload_fixture_r52_',dir='artifacts/tmp') as tmp:
        root=Path(tmp)
        rng=np.random.default_rng(20261018)
        n=42007
        arrays=dict(cell_integer_id=(np.arange(n)//17).astype('uint32'),
            gene_integer_id=(np.arange(n)%17).astype('uint16'),
            value=rng.integers(1,200,n,dtype=np.uint16))
        source=root/'source.lance'
        with LanceFileWriter(str(source),version='2.1',max_page_bytes=32768) as writer:
            writer.write_batch(pa.table(arrays))
        body=source.read_bytes()
        entry=dict(path='data/HCT116/expression.lance/data/synthetic.lance',size=len(body),oid='0'*40)
        frag=dict(id=0,inventory_entry=entry,physical_rows=n)
        spans=pd.DataFrame([dict(record_start=a,record_stop=b,cell_integer_id=a//17)
            for a,b in [(17,34),(4080,4097),(n-34,n-17)]])
        logical=spans[['record_start','record_stop']].to_numpy().tolist()
        mb=TransferBudget(root/'metadata.jsonl',10_000_000)
        store=Store(body,mb)
        meta_root=root/'metadata'
        mc=metadata.RangeCache(meta_root,mb,store)
        saved=metadata.describe_file(frag,logical,mc,root/'metaplans')
        _,bindings=reader.journal(mb.path)
        mb.close()
        original_meta=reader.META
        reader.META=meta_root
        try:
            boundary=[]
            for p in saved['pages']:
                if p['column']==0:
                    for lo,hi,a,b in lp.chunks(p,reader.metadata_bytes(0,*p['index_range'],bindings)):
                        if lo==0 or hi==n:
                            boundary.append([a,b])
            ranges=lp.merge(saved['payload']['payload_ranges']+boundary)
            fp=dict(payload_ranges=ranges,physical_rows=n)
            journal=root/'payload.jsonl'
            budget=TransferBudget(journal,10_000_000)
            budget.bind('file:'+entry['path'],dict(size=len(body),etag='"fixture"'))
            store=Store(body,budget)
            cache=reader.PayloadCache(root/'ranges',budget,store)
            result=reader.build_file(0,fp,entry,spans,cache,bindings,0,int(arrays['cell_integer_id'][-1]),root/'stage')
            expected=np.concatenate([np.arange(a,b) for a,b in logical])
            path=root/'stage/decoded/0.npz'
            with np.load(path) as d:
                for key in arrays:
                    np.testing.assert_array_equal(d[key],arrays[key][expected])
                np.testing.assert_array_equal(d['record_start'],spans.record_start)
            assert store.ranges==ranges
            assert budget.total==sum(b-a for a,b in ranges)
            assert result['records']==51
            checks+=['exact official selected-record decode and cell-ID boundaries from sparse file',
                'exact frozen requests with zero gap overfetch and durable hashes']
            calls=store.calls
            total=budget.total
            budget.close()
            budget=TransferBudget(journal,10_000_000)
            store.budget=budget
            cache=reader.PayloadCache(root/'ranges',budget,store)
            assert result==reader.build_file(0,fp,entry,spans,cache,bindings,0,int(arrays['cell_integer_id'][-1]),root/'stage')
            assert store.calls==calls and budget.total==total
            checks+=['completed decoded checkpoint no-I/O restart reuse']
            # Refuse corrupted committed arrays before network access.
            with path.open('ab') as out:
                out.write(b'corrupt')
            try:
                reader.build_file(0,fp,entry,spans,cache,bindings,0,int(arrays['cell_integer_id'][-1]),root/'stage')
            except AssertionError:
                pass
            else:
                raise AssertionError('Corrupt decoded checkpoint accepted')
            assert store.calls==calls
            checks+=['corrupt decoded checkpoint refusal']
            # Reconstruct from already downloaded ranges under a new output path.
            bad=spans.copy()
            bad.loc[0,'cell_integer_id']+=1
            for stage,selected,first in [('bad_cell',bad,0),('bad_boundary',spans,999)]:
                try:
                    reader.build_file(0,fp,entry,selected,cache,bindings,first,int(arrays['cell_integer_id'][-1]),root/stage)
                except AssertionError:
                    pass
                else:
                    raise AssertionError('Cell identity violation accepted')
                assert not list((root/stage).glob('decoded/*.npz'))
            assert store.calls==calls and budget.total==total
            checks+=['wrong selected or boundary cell ID refused without decoded publication']
            # Incomplete file after a committed first range: resume only misses.
            budget.close()
            budget=TransferBudget(root/'partial.jsonl',10_000_000)
            budget.bind('file:'+entry['path'],dict(size=len(body),etag='"fixture"'))
            store=Store(body,budget)
            cache=reader.PayloadCache(root/'partial_ranges',budget,store)
            cache.permit(0,ranges)
            cache.read(0,entry,*ranges[0])
            store.mode='network'
            try:
                cache.read(0,entry,*ranges[1])
            except OSError:
                pass
            else:
                raise AssertionError('Failure injection ignored')
            spent=budget.total
            budget.close()
            budget=TransferBudget(root/'partial.jsonl',10_000_000)
            store.budget=budget
            store.mode=None
            cache=reader.PayloadCache(root/'partial_ranges',budget,store)
            reader.build_file(0,fp,entry,spans,cache,bindings,0,int(arrays['cell_integer_id'][-1]),root/'resumed')
            assert store.ranges.count(ranges[0])==1
            assert budget.total==spent+sum(b-a for a,b in ranges[1:])
            checks+=['partial committed range reuse; failed ranges recharge after restart']
            budget.close()
        finally:
            reader.META=original_meta
        data=b'A'*100
        ent=dict(path='data/HCT116/expression.lance/data/failure.lance',size=100,oid='0'*40)
        for mode in ['network','truncated','ignored','range','etag','encoding']:
            path=root/(mode+'.jsonl')
            budget=TransferBudget(path,1000)
            budget.bind('file:'+ent['path'],dict(size=100,etag='"fixture"'))
            store=Store(data,budget)
            store.mode=mode
            cache=reader.PayloadCache(root/mode,budget,store)
            cache.permit(0,[[60,100]])
            try:
                cache.read(0,ent,60,100)
            except (OSError,AssertionError):
                pass
            else:
                raise AssertionError(mode)
            assert budget.total==40 and not list((root/mode).glob('**/*.bin'))
            budget.close()
            budget=TransferBudget(path,1000)
            store.budget=budget
            store.mode=None
            cache=reader.PayloadCache(root/mode,budget,store)
            cache.permit(0,[[60,100]])
            assert cache.read(0,ent,60,100)==b'A'*40 and budget.total==80
            dest=root/mode/'0/60-100.bin'
            dest.write_bytes(b'X'*40)
            before=(store.calls,budget.total)
            try:
                cache.read(0,ent,60,100)
            except AssertionError:
                pass
            else:
                raise AssertionError('Corrupt range accepted')
            assert before==(store.calls,budget.total)
            budget.close()
        checks+=['pre-I/O fsynced failure charges survive all six HTTP failure modes',
            'truncation/status/content-range/ETag/encoding refusal and corrupt cache refusal']
        budget=TransferBudget(root/'exhausted.jsonl',39)
        budget.bind('file:'+ent['path'],dict(size=100,etag='"fixture"'))
        store=Store(data,budget)
        cache=reader.PayloadCache(root/'exhausted',budget,store)
        cache.permit(0,[[60,100]])
        for a,b in [(60,100),(0,1)]:
            before=(store.calls,store.resolves,budget.total)
            try:
                cache.read(0,ent,a,b)
            except AssertionError:
                pass
            else:
                raise AssertionError('Budget or whitelist violation accepted')
            assert store.calls==before[0] and budget.total==before[2]
            if (a,b)==(0,1):
                assert store.resolves==before[1]
        budget.close()
        checks+=['budget exhausted before GET and unplanned range refused before resolution']
        # A fresh PublicStore resolution uses TransferBudget.bind; changed sources
        # cannot replace the metadata-stage prebound size/ETag even on restart.
        budget=TransferBudget(root/'source.jsonl',1000)
        budget.bind('file:'+ent['path'],dict(size=100,etag='"fixture"'))
        try:
            budget.bind('file:'+ent['path'],dict(size=100,etag='"changed"'))
        except AssertionError:
            pass
        else:
            raise AssertionError('Source rebinding accepted')
        budget.close()
        checks+=['metadata-stage file identity cannot change during fresh resolution']
    rt.save(reader.ROOT/'fixture.json',dict(status='PASS actual payload reader offline fixtures',
        utc=reader.datetime.now(reader.timezone.utc).isoformat(),plan_sha256=rt.sha(reader.PLAN),
        verifier_sha256=rt.sha(__file__),checks=checks,network_io=False,
        exact_selected_records=51,source_activated=False,holdout_scored=False))
    print('PASS',len(checks),'actual payload fixture groups',flush=True)


if __name__=='__main__':
    main()
