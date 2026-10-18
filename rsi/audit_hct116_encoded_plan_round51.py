"""Independent dynamic-protobuf, scalar chunk/cost and official metadata audit.

No planner helpers, network I/O, measured responses or scorer outputs are read.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
import hashlib,tempfile,re
import runtime as rt
import pandas as pd
from google.protobuf import descriptor_pb2,descriptor_pool,message_factory,any_pb2
from lance.file import LanceFileReader

ROOT=Path('artifacts/HCT116_encoded_metadata_round51')

def load(p):return rt.json.loads(Path(p).read_text())

def messages():
    """Independent official-schema subset, parsed by the protobuf library."""
    fd=descriptor_pb2.FileDescriptorProto(name='r51_independent.proto',package='r51',syntax='proto3')
    def msg(name):m=fd.message_type.add();m.name=name;return m
    def field(m,name,number,typ,repeated=False,message=None):
        f=m.field.add();f.name=name;f.number=number;f.type=typ;f.label=3 if repeated else 1
        if message:f.type_name='.r51.'+message
    T=descriptor_pb2.FieldDescriptorProto
    direct=msg('Direct');field(direct,'encoding',1,T.TYPE_BYTES)
    enc=msg('Encoding');field(enc,'direct',2,T.TYPE_MESSAGE,message='Direct')
    page=msg('Page')
    for name,num in [('offsets',1),('sizes',2)]:field(page,name,num,T.TYPE_UINT64,True)
    field(page,'length',3,T.TYPE_UINT64);field(page,'encoding',4,T.TYPE_MESSAGE,message='Encoding');field(page,'priority',5,T.TYPE_UINT64)
    col=msg('Column');field(col,'pages',2,T.TYPE_MESSAGE,True,'Page');field(col,'offsets',3,T.TYPE_UINT64,True);field(col,'sizes',4,T.TYPE_UINT64,True)
    mini=msg('Mini')
    for name,num in [('rep',1),('defn',2),('values',3),('dictionary',4)]:field(mini,name,num,T.TYPE_BYTES)
    for name,num in [('num_dict',5),('num_buffers',7),('rep_index_depth',8),('num_items',9)]:field(mini,name,num,T.TYPE_UINT64)
    field(mini,'layers',6,T.TYPE_UINT64,True);field(mini,'large',10,T.TYPE_BOOL)
    constant=msg('Constant');field(constant,'layers',5,T.TYPE_UINT64,True);field(constant,'inline',6,T.TYPE_BYTES)
    layout=msg('Layout');field(layout,'mini',1,T.TYPE_MESSAGE,message='Mini');field(layout,'constant',2,T.TYPE_MESSAGE,message='Constant')
    pool=descriptor_pool.DescriptorPool();pool.Add(fd)
    return {name:message_factory.GetMessageClass(pool.FindMessageTypeByName('r51.'+name)) for name in ['Column','Layout']}


def union(ranges):
    events=defaultdict(int)
    for start,end in ranges:
        assert start<=end
        if start==end:continue
        events[start]+=1;events[end]-=1
    count=0;out=[];opening=None
    for at,change in sorted(events.items()):
        before=count;count+=change;assert count>=0
        if before==0 and count>0:opening=at
        if before>0 and count==0:out.append([opening,at])
    assert count==0
    return out


def main():
    assert not (ROOT/'audit.json').exists(), 'Preserve completed audit'
    plan=load(ROOT/'plan.json');receipt=load(ROOT/'receipt.json')
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p
    assert receipt['plan_sha256']==rt.sha(ROOT/'plan.json')
    fixture=load(ROOT/'reader_fixture.json');assert fixture['status'].startswith('PASS')
    for p,h in fixture['sha256'].items():assert rt.sha(p)==h,p
    events=[rt.json.loads(s) for s in Path(receipt['journal']).read_text().splitlines()]
    total=0;bindings={};attempts=[]
    for e in events:
        if e['event']=='bind':assert e['key'] not in bindings;bindings[e['key']]=e['value']
        else:
            assert e['event']=='charge' and e['bytes']>0
            total+=e['bytes'];assert total==e['total']<=250_000_000;attempts.append(e)
    assert total==receipt['charged_bytes'] and bindings['plan_sha256']==rt.sha(ROOT/'plan.json')
    assert bindings['repository_revision']==plan['revision']=='598aa5442a4ee9a8f1a383a843ab5ce142904d04'
    frame=pd.read_parquet('artifacts/HCT116_extension_preflight_round49/logical_record_ranges.parquet')
    frags={f['id']:f for f in load('artifacts/HCT116_expression_identity_round50/fragments.json')}
    touched=sorted(set(map(int,frame.fragment_id)));assert len(touched)==receipt['files']==2018
    schema=messages();counts=Counter();all_legal={};metadata_total=0;payload_total=0;requests_total=0;logical_total=0;schema_total=0
    cache_hashes={};official_verified=[]
    # Every saved range is independently addressed and checked against durable
    # hashes before any offset/row/affordability conclusion is considered.
    def cached(fid,a,b):
        p=ROOT/str(fid)/f'{a}-{b}.bin';assert p.is_file()
        data=p.read_bytes();assert len(data)==b-a
        h=hashlib.sha256(data).hexdigest();assert bindings['cache:'+str(p)]==h
        cache_hashes[str(p)]=h;return data
    with tempfile.TemporaryDirectory(prefix='lance_r51_official_',dir='artifacts/tmp') as tmp:
        for fid in touched:
            f=frags[fid];entry=f['inventory_entry'];size=entry['size'];key=entry['path'];rows=f['physical_rows']
            savedpath=ROOT/str(fid)/'plan.json';assert receipt['file_plan_sha256'][str(savedpath)]==rt.sha(savedpath)
            saved=load(savedpath)
            assert bindings['inventory:'+key]==entry
            assert bindings['file:'+key]['size']==size
            foot=cached(fid,size-40,size)
            # Scalar LE footer decoder, independent of planner struct.unpack.
            values=[int.from_bytes(foot[i:j],'little') for i,j in [(0,8),(8,16),(16,24),(24,28),(28,32),(32,34),(34,36)]]
            a,b,c,ng,nc,major,minor=values
            assert foot[36:]==b'LANC' and (major,minor)==(2,1)
            assert nc==3 and ng==1 and 0<a<=b<c<size-40
            assert c-b==48 and size-40-c==16
            meta=cached(fid,a,size-40);legal=[[a,size-40],[size-40,size]]
            schema_start=int.from_bytes(meta[c-a:c-a+8],'little');schema_size=int.from_bytes(meta[c-a+8:c-a+16],'little')
            assert 0<schema_size<=64000 and 0<=schema_start<schema_start+schema_size<=a
            schema_range=[schema_start,schema_start+schema_size];schema_bytes=cached(fid,*schema_range);legal.append(schema_range);schema_total+=schema_size
            selected=frame.loc[frame.fragment_id==fid,['record_start','record_stop']].to_numpy().tolist()
            logical=union(selected);assert all(0<=x<y<=rows for x,y in logical)
            logical_total+=sum(y-x for x,y in logical)
            expected=[];page_buffers=[];saved_pages=[];overfetch=[0,0,0];chunks_count=0;used_pages=[]
            for col in range(3):
                cp=int.from_bytes(meta[b-a+col*16:b-a+col*16+8],'little');cl=int.from_bytes(meta[b-a+col*16+8:b-a+col*16+16],'little')
                assert a<=cp<=cp+cl<=b
                cm=schema['Column']();cm.ParseFromString(meta[cp-a:cp+cl-a]);assert not cm.offsets and not cm.sizes
                row=0
                for pi,p in enumerate(cm.pages):
                    assert p.length>0 and p.priority==row
                    start,end=row,row+p.length;row=end
                    assert len(p.offsets)==len(p.sizes)
                    buffers=[[x,x+y] for x,y in zip(p.offsets,p.sizes)];page_buffers.append(buffers)
                    assert all(0<=x<=y<=schema_start for x,y in buffers)
                    anymsg=any_pb2.Any();anymsg.ParseFromString(p.encoding.direct.encoding)
                    assert anymsg.type_url=='/lance.encodings21.PageLayout'
                    layout=schema['Layout']();layout.ParseFromString(anymsg.value)
                    assert layout.HasField('mini') != layout.HasField('constant')
                    mini=layout.HasField('mini');counts['miniblock' if mini else 'constant']+=1
                    useful=[(max(start,x),min(end,y)) for x,y in logical if x<end and y>start]
                    wanted=sum(y-x for x,y in useful)
                    if useful:used_pages.append([col,pi])
                    s=saved['pages'][len(saved_pages)];saved_pages.append(s)
                    assert [s['column'],s['page'],s['row_start'],s['row_stop'],s['buffers'],s['layout']]==[col,pi,start,end,buffers,1 if mini else 2]
                    if mini:
                        m=layout.mini;assert list(m.layers)==[1] and not m.rep_index_depth and not m.large and m.num_items==p.length
                        assert not m.rep and not m.defn
                        assert len(buffers)==(3 if m.dictionary else 2)
                        if useful:
                            ir=buffers[0];legal.append(ir);index=cached(fid,*ir);assert len(index)%2==0
                            words=[int.from_bytes(index[i:i+2],'little') for i in range(0,len(index),2)]
                            assert words
                            rr=start;bb=buffers[1][0]
                            for wi,word in enumerate(words):
                                exponent=word%16;span=(2**exponent) if wi<len(words)-1 else end-rr
                                assert span>0 and (wi==len(words)-1 or exponent>0)
                                if wi==len(words)-1 and exponent:assert span==2**exponent
                                nbytes=(word//16+1)*8
                                if any(rr<y and rr+span>x for x,y in useful):
                                    expected.append([bb,bb+nbytes]);overfetch[col]+=span;chunks_count+=1
                                rr+=span;bb+=nbytes
                            assert rr==end and bb==buffers[1][1]
                            overfetch[col]-=wanted
                            expected.extend(buffers[2:])
                    else:
                        assert list(layout.constant.layers)==[1]
                        if useful:expected.extend(buffers)
                assert row==rows
            expected=union(expected);legal=union(legal)
            assert saved['payload']['payload_ranges']==expected
            assert saved['payload']['payload_bytes']==sum(y-x for x,y in expected)
            assert saved['payload']['selected_chunks']==chunks_count
            assert saved['payload']['encoded_row_overfetch_by_column']==overfetch
            assert saved['payload']['used_pages']==used_pages
            assert saved['metadata_ranges']==legal and saved['schema_range']==schema_range
            assert saved['payload']['logical_rows']==sum(y-x for x,y in logical)
            assert saved['metadata_unique_bytes']==sum(y-x for x,y in legal)
            assert not saved['expression_values_fetched'] and not saved['partial_file_full_digest_verified']
            # Every metadata byte interval must be outside each encoded value,
            # dictionary or scalar buffer; adjacency is allowed, overlap is not.
            value_ranges=[r for s in saved_pages for r in (s['buffers'][1:] if s['layout']==1 else s['buffers'])]
            assert all(y<=u or v<=x for x,y in legal for u,v in value_ranges)
            all_legal[key]=legal
            metadata_total+=sum(y-x for x,y in legal);payload_total+=sum(y-x for x,y in expected);requests_total+=len(expected)
            # Deterministic id stride, plus both endpoint files; no outcome choice.
            if fid==touched[0] or fid==touched[-1] or fid%100==0:
                local=Path(tmp)/f'{fid}.lance'
                with local.open('wb') as out:
                    out.truncate(size)
                    # Copy individual cached ranges, not their coalesced names.
                    for p in (ROOT/str(fid)).glob('*.bin'):
                        x,y=map(int,p.stem.split('-'));out.seek(x);out.write(cached(fid,x,y))
                official=LanceFileReader(str(local)).metadata()
                assert official.num_rows==rows and official.schema.names==['cell_integer_id','gene_integer_id','value']
                assert str(official.schema.field(0).type)=='uint32' and str(official.schema.field(1).type)=='uint16' and str(official.schema.field(2).type)=='uint16'
                actual=[[[q.position,q.position+q.size] for q in p.buffers] for col in official.columns for p in col.pages]
                assert actual==page_buffers
                official_verified.append(fid);local.unlink()
            if len(all_legal)%100==0:print('Independent audit',len(all_legal),'/2018',flush=True)
    # Account for *every* public attempt, including failures/retries. Requests
    # must lie entirely within a proved metadata interval for a touched file.
    heads=0;ranges=0
    for e in attempts:
        reason=e['reason']
        if reason.startswith('Orion HEAD metadata upper bound '):
            key=reason.removeprefix('Orion HEAD metadata upper bound ');assert key in all_legal and e['bytes']==64000;heads+=1
        else:
            m=re.fullmatch(r'Exact metadata range attempt (.+) \[(\d+),(\d+)\)',reason);assert m,reason
            key,x,y=m[1],int(m[2]),int(m[3]);assert e['bytes']==y-x
            assert any(a<=x<y<=b for a,b in all_legal[key]);ranges+=1
    assert metadata_total==receipt['metadata_unique_bytes']
    assert payload_total==receipt['encoded_value_payload_bytes'] and requests_total==receipt['encoded_value_requests']
    assert schema_total==receipt['schema_bytes_cached']
    assert logical_total==receipt['logical_rows']==215_800_779
    assert receipt['future_payload_head_reserve_bytes']==2018*64000
    assert not any(receipt[k] for k in ['expression_values_fetched','source_activated','new_biological_scores','holdout_scored','full_partial_file_digest_verified'])
    # Existing frozen inputs, old bank/score/gate identities are in the plan.
    # The historical metadata index also independently checks archive members.
    frozen=load('artifacts/frozen_metadata_index.json')
    assert rt.sha(frozen['archive'])==frozen['archive_sha256']
    for p,e in frozen['entries'].items():assert rt.sha(p)==e['sha256'],p
    assert len(frozen['entries'])==38
    gate=load('artifacts/calibration_gate.json')
    assert gate['status']=='CLOSED' and not gate['phase2_allowed']
    a=load('artifacts/download_allocations.json')
    assert sum(a['new_allocations'].values())==a['total_allocated_upper_bound_bytes']<=15_000_000_000
    assert a['new_allocations']['H1_rich']==5_250_000_000
    rt.save(ROOT/'audit.json',dict(status='PASS independent dynamic-protobuf metadata/request-cost audit; expression extraction unverified',utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=rt.sha(__file__),receipt_sha256=rt.sha(ROOT/'receipt.json'),plan_sha256=rt.sha(ROOT/'plan.json'),files=2018,logical_rows=logical_total,metadata_unique_bytes=metadata_total,encoded_value_payload_bytes=payload_total,encoded_value_requests=requests_total,head_attempts=heads,range_attempts=ranges,charged_bytes=total,layouts=dict(counts),official_metadata_file_ids=official_verified,
        cache_range_hashes_checked=len(cache_hashes),all38_frozen_members_unchanged=True,
        first_future_attempt_minimum_bytes=payload_total+2018*64000,
        future_double_payload_attempt_reserve_bytes=2*payload_total+2*2018*64000,
        future_reader_still_requires_allocation_freeze_identity_and_failure_audit=True,
        all_cell_ids_and_moments_unverified=True,source_activated=False,holdout_scored=False,phase2_allowed=False))
    print('PASS independent request cost',payload_total,'encoded values bytes; no value reads',flush=True)
if __name__=='__main__':main()
