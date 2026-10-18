"""Actual-reader failure/reuse tests and official Lance sparse reconstruction."""
import io,tempfile
from pathlib import Path
import runtime as rt
import numpy as np
import pyarrow as pa
from lance.file import LanceFileReader,LanceFileWriter
from transfer_budget import TransferBudget
import acquire_hct116_encoded_metadata_round51 as reader
import lance_encoded_plan_round51 as lp

class Response:
    def __init__(self,data,start,stop,size,mode):
        self.status_code=200 if mode=='ignored' else 206
        self.headers={'Content-Length':str(stop-start),'Content-Range':f'bytes {start}-{stop-1}/{size}','ETag':'"fixture"'}
        if mode=='range':self.headers['Content-Range']='bytes 0-1/2'
        if mode=='etag':self.headers['ETag']='"changed"'
        if mode=='encoding':self.headers['Content-Encoding']='gzip'
        self.raw=io.BytesIO(data[start:stop-(1 if mode=='truncated' else 0)])
    def __enter__(self):return self
    def __exit__(self,*args):pass

class Store:
    def __init__(self,data,budget):self.data=data;self.budget=budget;self.mode=None;self.calls=0;self.resolves=0;self.ranges=[];self.http=self
    def resolve(self,key):self.resolves+=1;return 'https://fixture.invalid/pinned',dict(size=len(self.data),etag='"fixture"')
    def get(self,url,**kwargs):
        # Inspect the fsynced journal from a separate file handle before I/O.
        events=[rt.json.loads(s) for s in self.budget.path.read_text().splitlines()]
        assert events[-1]['event']=='charge' and events[-1]['total']==self.budget.total
        assert kwargs['allow_redirects'] is False and kwargs['stream']
        assert kwargs['headers']['If-Match']=='"fixture"'
        start,stop=map(int,kwargs['headers']['Range'][6:].split('-'));stop+=1
        assert events[-1]['bytes']==stop-start
        self.calls+=1;self.ranges.append([start,stop])
        if self.mode=='network':raise OSError('Injected failed request')
        return Response(self.data,start,stop,len(self.data),self.mode)


def main():
    root=reader.ROOT;root.mkdir(exist_ok=True)
    assert not (root/'reader_fixture.json').exists(), 'Preserve fixture evidence'
    checks=[];fixture_records=[]
    with tempfile.TemporaryDirectory(dir='artifacts/tmp',prefix='lance_encoded_r51_') as tmp:
        base=Path(tmp);rng=np.random.default_rng(20261018)
        n=42000
        cases=[dict(cell_integer_id=np.arange(n,dtype=np.uint32)//31,
                    gene_integer_id=rng.integers(0,38606,n,dtype=np.uint16),
                    value=rng.integers(0,200,n,dtype=np.uint16)),
               dict(cell_integer_id=np.full(9000,9,np.uint32),gene_integer_id=np.full(9000,7,np.uint16),value=np.full(9000,3,np.uint16))]
        for ci,arrays in enumerate(cases):
            table=pa.table(arrays);path=base/f'full{ci}.lance'
            with LanceFileWriter(str(path),version='2.1',max_page_bytes=32768) as writer:writer.write_batch(table)
            data=path.read_bytes();size=len(data);rows=len(table)
            entry=dict(path=f'data/HCT116/expression.lance/data/fixture{ci}.lance',size=size,oid='0'*40)
            fragment=dict(id=ci,inventory_entry=entry,physical_rows=rows)
            journal=base/f'fixture{ci}.jsonl';budget=TransferBudget(journal,10_000_000);store=Store(data,budget)
            cache=reader.RangeCache(base/f'cache{ci}',budget,store)
            logical=[[0,23],[4095,4102],[7300,7380],[rows-10,rows]]
            result=reader.describe_file(fragment,logical,cache,base/f'plans{ci}')
            requested=list(store.ranges)
            true_meta=LanceFileReader(str(path)).metadata()
            assert true_meta.num_rows==rows
            actual=[[[b.position,b.position+b.size] for b in p.buffers] for c in true_meta.columns for p in c.pages]
            assert actual==[p['buffers'] for p in result['pages']]
            values=[r for p in result['pages'] for r in (p['buffers'][1:] if p['layout']==1 else p['buffers'])]
            assert all(b<=x or y<=a for a,b in requested for x,y in values), 'Reader fetched synthetic value data'
            assert sum(b-a for a,b in requested)==result['metadata_unique_bytes']
            # Populate a sparse file with exactly the later planned requests;
            # official decoding must recover every selected value exactly.
            sparse=base/f'sparse{ci}.lance'
            ranges=lp.merge(requested+[result['schema_range']]+result['payload']['payload_ranges'])
            with sparse.open('wb') as out:
                out.truncate(size)
                for a,b in ranges:out.seek(a);out.write(data[a:b])
            official=LanceFileReader(str(sparse))
            for a,b in logical:
                assert official.read_range(a,b-a).to_table().equals(table.slice(a,b-a))
            # Exact row coverage and chunk framing are challenged at boundaries.
            assert result['payload']['logical_rows']==sum(b-a for a,b in logical)
            calls=store.calls;total=budget.total;budget.close()
            budget=TransferBudget(journal,10_000_000);store.budget=budget
            resumed=reader.RangeCache(base/f'cache{ci}',budget,store)
            replay=reader.describe_file(fragment,logical,resumed,base/f'replay{ci}')
            assert replay==result and store.calls==calls and budget.total==total
            # Values not in the metadata whitelist must fail before resolution/charge.
            forbidden=values[0] if values else [1,2]
            before=(store.calls,store.resolves,budget.total)
            try:resumed.read(ci,entry,*forbidden)
            except AssertionError:pass
            else:raise AssertionError('Unplanned value read accepted')
            assert before==(store.calls,store.resolves,budget.total)
            first=base/f'cache{ci}'/str(ci)/f'{size-40}-{size}.bin';first.write_bytes(b'X'*40)
            try:resumed.read(ci,entry,size-40,size)
            except AssertionError:pass
            else:raise AssertionError('Corrupt committed cache accepted')
            assert store.calls==calls
            fixture_records.append(dict(rows=rows,pages=len(result['pages']),layouts=sorted(set(p['layout'] for p in result['pages'])),metadata_bytes=result['metadata_unique_bytes'],planned_payload_bytes=result['payload']['payload_bytes'],selected_chunks=result['payload']['selected_chunks'],official_exact=True))
            budget.close()
        checks+=['official page buffer offset equality','no encoded values fetched','official sparse reconstruction exact at boundaries','constant and miniblock layouts','durable no-I/O restart reuse','forbidden value read before charge or I/O','corrupt cache refusal']
        # Failure charges persist across restart; no failed attempt becomes a cache.
        data=b'A'*100;entry=dict(path='data/HCT116/expression.lance/data/failure.lance',size=100,oid='0'*40)
        for mode in ['network','truncated','ignored','range','etag','encoding']:
            journal=base/f'{mode}.jsonl';budget=TransferBudget(journal,1000);store=Store(data,budget);store.mode=mode
            cache=reader.RangeCache(base/mode,budget,store);cache.permit(0,[[60,100]])
            try:cache.read(0,entry,60,100)
            except (AssertionError,OSError):pass
            else:raise AssertionError('Injected failure accepted '+mode)
            assert budget.total==40 and not list((base/mode).glob('**/*.bin'))
            budget.close();budget=TransferBudget(journal,1000);store.budget=budget;store.mode=None
            cache=reader.RangeCache(base/mode,budget,store);cache.permit(0,[[60,100]])
            assert cache.read(0,entry,60,100)==b'A'*40 and budget.total==80
            budget.close()
        checks+=['pre-I/O fsynced failed-attempt charges','network/truncation/full-response/Content-Range/ETag/encoding refusal','restart retains failed charges and retries cost again']
        budget=TransferBudget(base/'exhausted.jsonl',39);store=Store(data,budget);cache=reader.RangeCache(base/'exhausted',budget,store);cache.permit(0,[[60,100]])
        try:cache.read(0,entry,60,100)
        except AssertionError:pass
        else:raise AssertionError('Budget exhaustion ignored')
        assert store.calls==0;budget.close();checks+=['budget exhaustion before GET']
    paths=[Path(__file__),Path(reader.__file__),Path(lp.__file__),Path('durable_public_store.py'),Path('transfer_budget.py')]
    rt.save(root/'reader_fixture.json',dict(status='PASS actual bounded range reader and official sparse decode fixtures',checks=checks,fixtures=fixture_records,network_io=False,sha256={str(p):rt.sha(p) for p in paths}))
    print('PASS actual-reader fixtures',len(checks),'groups',flush=True)
if __name__=='__main__':main()
