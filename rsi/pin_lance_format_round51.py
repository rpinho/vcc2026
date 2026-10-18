"""Replay/restore the pinned primary Lance format sources with bounded reads.

Original discovery errors and charges remain in the two journals. Fresh replay
requires the same allocation entries or separately versioned journals, never
reset these ledgers. This reader fetches code only, no expression data.
"""
from pathlib import Path
import gzip
import runtime as rt
import requests
from transfer_budget import TransferBudget
ROOT=Path('artifacts/Lance_format_round51')
COMMIT='0923349baa19fd9af776b7d4626dc758dc610c80'
BASE='https://raw.githubusercontent.com/lance-format/lance/'+COMMIT+'/'

def main():
    a=rt.json.loads(Path('artifacts/download_allocations.json').read_text())
    assert a['new_allocations']['Lance_format_code_round51']==2_000_000
    assert a['new_allocations']['Lance_primitive_code_round51']==600_000
    b=TransferBudget(ROOT/'charges.jsonl',2_000_000)
    extra=TransferBudget(ROOT/'primitive_charges.jsonl',600_000)
    session=requests.Session();session.trust_env=False
    files=[]
    try:
        assert b.bindings['commit']==COMMIT
        directories=[('directory.json','rust/lance-file/src/', ['reader.rs','format.rs']),
                     ('primitive_directory.json','rust/lance-encoding/src/encodings/logical/primitive/', ['miniblock.rs','chunk_index.rs','constant.rs','fullzip.rs']),
                     ('logical_directory.json','rust/lance-encoding/src/encodings/logical/', ['primitive.rs']),
                     ('protos_directory.json','protos/', ['encodings_v2_1.proto','file2.proto'])]
        for listing,prefix,names in directories:
            path=ROOT/listing
            assert b.bindings['cache:'+listing]==rt.sha(path)
            entries=rt.json.loads(path.read_text())
            for name in names:
                e=next(x for x in entries if x['name']==name)
                dest=ROOT/('file2_decoded.proto' if name=='file2.proto' else name)
                budget=extra if name=='primitive.rs' else b
                if dest.exists():
                    assert budget.bindings['cache:'+dest.name]==rt.sha(dest)
                    data=dest.read_bytes()
                else:
                    budget.charge(e['size']+1,'Pinned exact-size primary source '+prefix+name)
                    with session.get(BASE+prefix+name,stream=True,allow_redirects=False,timeout=(10,60),headers={'Accept-Encoding':'identity'}) as r:
                        r.raise_for_status();assert r.status_code==200
                        data=r.raw.read(e['size']+1,decode_content=True)
                    assert len(data)==e['size']
                    assert rt.hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['sha']
                    dest.write_bytes(data);budget.bind('cache:'+dest.name,rt.sha(dest))
                assert len(data)==e['size']
                assert rt.hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['sha']
                files.append(dict(path=str(dest),source=BASE+prefix+name,size=e['size'],git_blob=e['sha'],sha256=rt.sha(dest)))
        assert gzip.decompress((ROOT/'file2.proto').read_bytes())==(ROOT/'file2_decoded.proto').read_bytes()
        rt.save(ROOT/'receipt.json',dict(status='PASS pinned primary v13.0.0 sources; raw source git-blob identities verified',commit=COMMIT,source='https://github.com/lance-format/lance/tree/'+COMMIT,
            files=files,charged_bytes=b.total+extra.total,journals={str(x.path):dict(sha256=rt.sha(x.path),charged_bytes=x.total,limit_bytes=x.limit) for x in [b,extra]},
            notes='Discovery404 attempts and symlink response retained; compressed raw file2 original retained separately and decoded within its charged100KB cap. Original failed script retained. No payload reads.'))
        print('PASS primary sources',len(files),'total charged',b.total+extra.total,flush=True)
    finally:session.close();b.close();extra.close()
if __name__=='__main__':main()
