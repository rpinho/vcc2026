"""Network-free actual-reader failure, retry accounting and reuse checks."""
import hashlib
import io
from pathlib import Path
import tempfile
from types import SimpleNamespace

import acquire_hct116_full_metadata_round48 as reader
from transfer_budget import TransferBudget


class Response:
    status_code=200

    def __init__(self,data,length,etag='"fixture"'):
        self.raw=io.BytesIO(data)
        self.headers={'Content-Length':str(length),'ETag':etag}

    def raise_for_status(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self,*args):
        pass


def main():
    payload=b'offline HCT116 metadata reader fixture'
    entry=dict(path=reader.PREFIX+'data/fixture.lance',size=len(payload),
               oid=hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest())
    with tempfile.TemporaryDirectory(prefix='hct_metadata_fixture_',dir='artifacts/tmp') as folder:
        folder=Path(folder);journal=folder/'charges.jsonl';destination=folder/'fixture.lance'
        budget=TransferBudget(journal,1000)
        state={'body':payload[:-1],'etag':'"fixture"','calls':0}
        def get(*args,**kwargs):
            assert budget.total==(state['calls']+1)*len(payload), 'GET before durable charge'
            state['calls']+=1
            return Response(state['body'],len(payload),state['etag'])
        def resolve(path):
            assert path==entry['path']
            return 'http://offline.fixture',dict(size=len(payload),etag='"fixture"')
        store=SimpleNamespace(resolve=resolve,http=SimpleNamespace(get=get))
        try:
            reader.fetch_file(entry,store,budget,destination)
            raise RuntimeError('Truncated response accepted')
        except AssertionError as e:
            assert 'Truncated' in str(e)
        assert not destination.exists() and budget.total==len(payload)
        budget.close();budget=TransferBudget(journal,1000)
        state.update(body=payload,etag='"wrong"')
        try:
            reader.fetch_file(entry,store,budget,destination)
            raise RuntimeError('Changed ETag accepted')
        except AssertionError:
            pass
        assert not destination.exists() and budget.total==2*len(payload)
        state.update(body=b'x'*len(payload),etag='"fixture"')
        try:
            reader.fetch_file(entry,store,budget,destination)
            raise RuntimeError('Changed payload accepted')
        except AssertionError as e:
            assert 'Git blob mismatch' in str(e)
        assert not destination.exists() and budget.total==3*len(payload)
        state['body']=payload
        h=reader.fetch_file(entry,store,budget,destination)
        assert h==hashlib.sha256(payload).hexdigest() and budget.total==4*len(payload)
        assert reader.fetch_file(entry,store,budget,destination)==h and state['calls']==4
        lfs=dict(entry,lfs=dict(oid=h))
        assert reader.verify(destination,lfs)==h
        destination.write_bytes(b'corrupt')
        try:
            reader.fetch_file(entry,store,budget,destination)
            raise RuntimeError('Corrupt completed file accepted')
        except AssertionError:
            pass
        assert state['calls']==4
        budget.close()
    reader.rt.save('artifacts/HCT116_full_metadata_reader_fixture_round48.json',dict(
        status='PASS network-free actual downloader fixture',reader_sha256=reader.rt.sha(reader.__file__),
        verifier_sha256=reader.rt.sha(__file__),
        checks=['pre-I/O durable charge','truncation refusal','restart charge retention','ETag refusal',
                'git-blob mismatch refusal','successful atomic retry','no-I/O completed reuse',
                'LFS SHA256 branch','corrupt completed file refusal'],network_io=False))
    print('PASS nine grouped actual-reader checks')


if __name__=='__main__':
    main()
