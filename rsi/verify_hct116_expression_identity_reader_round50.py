"""Network-free tests of the actual inventory reader and inherited file reader."""
import io
from pathlib import Path
import tempfile

import acquire_hct116_expression_identity_round50 as reader
import verify_hct116_full_metadata_reader_round48 as file_fixture
from transfer_budget import TransferBudget


class Response:
    status_code = 200
    def __init__(self,body,links=None):
        class Raw(io.BytesIO):
            def read(self,size=-1,decode_content=None): return super().read(size)
        self.raw=Raw(body); self.links=links or {}
    def raise_for_status(self): pass
    def __enter__(self): return self
    def __exit__(self,*args): pass


def main():
    # Runs the same full-file function imported by this reader, including durable
    # failed-attempt charges, immutable identities, atomic completion and reuse.
    prior=Path('artifacts/HCT116_full_metadata_reader_fixture_round48.json').read_bytes()
    file_fixture.main()
    # Preserve the historical receipt's byte identity (the fixture is deterministic).
    assert Path('artifacts/HCT116_full_metadata_reader_fixture_round48.json').read_bytes() == prior
    original_root=reader.ROOT
    with tempfile.TemporaryDirectory(prefix='hct_expr_reader_fixture_',dir='artifacts/tmp') as d:
        reader.ROOT=Path(d)
        journal=reader.ROOT/'charges.jsonl'
        budget=TransferBudget(journal,reader.CAP*10)
        folder='expression.lance/data'
        base=reader.API+reader.PREFIX+folder+'?recursive=false&limit=1000'
        row=lambda name: dict(type='file',path=reader.PREFIX+folder+'/'+name,size=1,oid='0'*40)
        state=dict(calls=0,fail=True)
        class Session:
            def get(self,url,**kwargs):
                assert budget.total == (state['calls']+1)*reader.CAP
                state['calls']+=1
                assert kwargs['allow_redirects'] is False
                if state['fail']:
                    return Response(b'[')  # Failed JSON must never become a checkpoint.
                if url==base:
                    return Response(reader.rt.json.dumps([row('a')]).encode(),{'next':{'url':base+'&cursor=second'}})
                assert url==base+'&cursor=second'
                return Response(reader.rt.json.dumps([row('b')]).encode())
        session=Session()
        try:
            reader.inventory(session,budget,folder)
            raise RuntimeError('Malformed response accepted')
        except reader.rt.json.JSONDecodeError: pass
        assert not list(reader.ROOT.glob('*page*.json'))
        budget.close(); budget=TransferBudget(journal,reader.CAP*10)
        state['fail']=False
        result=reader.inventory(session,budget,folder)
        assert result['complete'] and result['pages']==2 and len(result['entries'])==2
        assert budget.total==reader.CAP*3
        # Reopen the real journal, reuse pages with no requests and no new charge.
        budget.close(); budget=TransferBudget(journal,reader.CAP*10)
        assert reader.inventory(session,budget,folder)==result and state['calls']==3
        path=reader.ROOT/'expression.lance_data_page0.json'
        p=reader.load(path);p['next_url']='https://invalid.example/';reader.rt.save(path,p)
        try:
            reader.inventory(session,budget,folder)
            raise RuntimeError('Changed completed inventory accepted')
        except AssertionError: pass
        assert state['calls']==3
        budget.close()
    reader.ROOT=original_root
    reader.rt.save('artifacts/HCT116_expression_identity_reader_fixture_round50.json',dict(
        status='PASS actual inventory reader and inherited full-file reader, network-free',
        reader_sha256=reader.rt.sha(reader.__file__),verifier_sha256=reader.rt.sha(__file__),
        inherited_file_reader_sha256=reader.rt.sha(file_fixture.reader.__file__),
        checks=['pre-I/O durable per-page charge','malformed JSON not cached','restart retains failed charge',
            'complete pagination','no-I/O cached inventory reuse','changed cache refusal',
            'inherited nine full-file accounting/identity/reuse checks'],network_io=False))
    print('PASS inventory/file identity accounting fixtures')


if __name__=='__main__': main()
