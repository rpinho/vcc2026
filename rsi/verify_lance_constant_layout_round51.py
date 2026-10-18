"""Version-scope refusal fixture and separate completed-live layout gate.

Does not modify any frozen producer/auditor, transfer journal or cached range.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,struct,tempfile
import runtime as rt
import numpy as np
import pyarrow as pa
from google.protobuf import any_pb2
from lance.file import LanceFileReader,LanceFileWriter
from audit_hct116_encoded_plan_round51 import messages
import lance_encoded_plan_round51 as lp
ROOT=Path('artifacts/HCT116_encoded_metadata_round51')

def load(p):return rt.json.loads(Path(p).read_text())

def fixture():
    dest=ROOT/'constant_scope_fixture.json';assert not dest.exists()
    with tempfile.TemporaryDirectory(prefix='lance_constant_scope_r51_',dir='artifacts/tmp') as tmp:
        root=Path(tmp);n=9000;values=[9,7,3]
        table=pa.table(dict(cell_integer_id=np.full(n,9,np.uint32),gene_integer_id=np.full(n,7,np.uint16),value=np.full(n,3,np.uint16)))
        full=root/'writer.lance'
        with LanceFileWriter(str(full),version='2.1') as writer:writer.write_batch(table)
        original=full.read_bytes();f=lp.footer(original[-40:],len(original));a,b,c=f['column_start'],f['cmo_start'],f['gbo_start']
        types=messages();columns=[];offsets=[];cursor=a
        for i,v in enumerate(values):
            cm=types['Column']();page=cm.pages.add();page.length=n
            layout=types['Layout']();layout.constant.layers.append(1);layout.constant.inline=v.to_bytes(4 if i==0 else 2,'little')
            page.encoding.direct.encoding=any_pb2.Any(type_url='/lance.encodings21.PageLayout',value=layout.SerializeToString()).SerializeToString()
            encoded=cm.SerializeToString();columns.append(encoded);offsets.append([cursor,len(encoded)]);cursor+=len(encoded)
        data=original[:a]+b''.join(columns)+b''.join(struct.pack('<QQ',*r) for r in offsets)+original[c:-40]+struct.pack('<QQQIIHH4s',a,cursor,cursor+48,1,3,2,1,b'LANC')
        path=root/'unsupported.lance';path.write_bytes(data)
        try:LanceFileReader(str(path))
        except OSError as e:
            assert 'v2.1 only accepts the all-null form of constant page layout' in str(e)
        else:raise AssertionError('Expected installed v2.1 refusal')
        foot=lp.footer(data[-40:],len(data));pages=lp.pages_from_metadata(data[foot['column_start']:-40],foot,n)
        assert [p['layout'] for p in pages]==[2,2,2]
        rt.save(dest,dict(status='PASS official version-scope refusal fixture; planner layout2 branch NOT certified',
            utc=datetime.now(timezone.utc).isoformat(),official_rejects_nonnull_constant_layout2_in_v21=True,
            generic_planner_accepts_but_not_supported_by_file_version=True,
            clarification='The writer-generated constant-valued fixture emitted layout1, so its caption constant and miniblock layouts is imprecise. It covers random and constant-valued data through miniblocks. Non-null layout2 is outside Lance2.1; the generic planner branch alone is not sufficient format validation. Require the separate all-file miniblock scope audit before accepting live costs or starting a payload reader.',
            no_live_metadata_or_producer_mutated=True,network_io=False,
            sha256={str(p):rt.sha(p) for p in [Path(__file__),Path(lp.__file__),ROOT/'reader_fixture.json']}))
    print('PASS official version refusal; layout2 unsupported in this planner scope',flush=True)

def live():
    assert not (ROOT/'scope_audit.json').exists()
    plan=load(ROOT/'plan.json');receipt=load(ROOT/'receipt.json');audit=load(ROOT/'audit.json');scope=load(ROOT/'constant_scope_fixture.json')
    assert audit['status'].startswith('PASS') and audit['receipt_sha256']==rt.sha(ROOT/'receipt.json')
    assert scope['status'].startswith('PASS')
    for p,h in scope['sha256'].items():assert rt.sha(p)==h,p
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p
    assert audit['files']==receipt['files']==2018
    assert audit['layouts'].get('constant',0)==0
    pages=0
    for p,h in receipt['file_plan_sha256'].items():
        assert rt.sha(p)==h
        row=load(p);assert row['pages']
        assert all(x['layout']==1 for x in row['pages']), 'Any non-miniblock requires a new format contract before accepting cost'
        pages+=len(row['pages'])
    assert audit['layouts']=={'miniblock':pages}
    rt.save(ROOT/'scope_audit.json',dict(status='PASS all live files in official-tested miniblock scope',files=2018,pages=pages,
        utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=rt.sha(__file__),main_audit_sha256=rt.sha(ROOT/'audit.json'),receipt_sha256=rt.sha(ROOT/'receipt.json'),scope_fixture_sha256=rt.sha(ROOT/'constant_scope_fixture.json'),expression_values_fetched=False,source_activated=False))
    print('PASS full live layout scope',2018,'files',pages,'pages',flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--live',action='store_true');args=parser.parse_args()
    live() if args.live else fixture()
