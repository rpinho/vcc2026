"""Payload-free Lance2.1 footer/page/miniblock index and exact chunk-range planner.

Primary Apache2 sources are pinned in artifacts/Lance_format_round51. This
module deliberately does not decode or retrieve any encoded expression values.
"""
import struct


def varint(data, pos=0):
    value=0
    for shift in range(0,70,7):
        assert pos<len(data), 'Truncated protobuf varint'
        x=data[pos];pos+=1;value|=(x&127)<<shift
        if x<128:return value,pos
    raise AssertionError('Overlong protobuf varint')


def wire(data):
    out={};pos=0
    while pos<len(data):
        key,pos=varint(data,pos);field,kind=key>>3,key&7
        assert field>0
        if kind==0:value,pos=varint(data,pos)
        elif kind==2:
            n,pos=varint(data,pos);assert pos+n<=len(data);value=data[pos:pos+n];pos+=n
        elif kind in [1,5]:
            n={1:8,5:4}[kind];assert pos+n<=len(data);value=data[pos:pos+n];pos+=n
        else:raise AssertionError('Unsupported protobuf wire type')
        out.setdefault(field,[]).append((kind,value))
    return out


def one(msg,field,default=None):
    values=msg.get(field,[]);assert len(values)<=1
    return values[0][1] if values else default


def packed(msg,field):
    out=[]
    for kind,value in msg.get(field,[]):
        if kind==0:out.append(value)
        else:
            assert kind==2;pos=0
            while pos<len(value):x,pos=varint(value,pos);out.append(x)
    return out


def merge(ranges,gap=0):
    result=[]
    for a,b in sorted((int(a),int(b)) for a,b in ranges if a<b):
        if result and a<=result[-1][1]+gap:result[-1][1]=max(result[-1][1],b)
        else:result.append([a,b])
    return result


def footer(data,size):
    assert len(data)==40
    a,b,c,ng,nc,major,minor,magic=struct.unpack('<QQQIIHH4s',data)
    assert magic==b'LANC' and (major,minor)==(2,1)
    assert nc==3 and 0<a<=b<=c<=size-40
    assert c-b==16*nc and size-40-c==16*ng
    assert size-40-a<=256_000, 'Column metadata ceiling'
    return dict(column_start=a,cmo_start=b,gbo_start=c,num_global_buffers=ng,num_columns=nc,size=size)


def pages_from_metadata(data,foot,rows):
    """data covers only [column_start, size-40), excludes schema/value buffers."""
    a,b,c=foot['column_start'],foot['cmo_start'],foot['gbo_start']
    assert len(data)==foot['size']-40-a
    result=[]
    for col in range(3):
        pos,n=struct.unpack_from('<QQ',data,b-a+16*col)
        assert a<=pos<=pos+n<=b
        cm=wire(data[pos-a:pos+n-a]);assert not packed(cm,3) and not packed(cm,4), 'Column buffers unsupported'
        page_rows=0
        for pi,(kind,pb) in enumerate(cm.get(2,[])):
            assert kind==2;p=wire(pb);length=one(p,3,0);priority=one(p,5,0)
            assert length>0 and priority==page_rows
            offsets=packed(p,1);sizes=packed(p,2);assert len(offsets)==len(sizes)
            buffers=[[int(o),int(o+s)] for o,s in zip(offsets,sizes)]
            assert all(0<=x<=y<=a for x,y in buffers)
            enc=wire(one(p,4));assert set(enc)=={2}, 'Deferred/unknown encoding refused'
            direct=wire(one(enc,2));anymsg=wire(one(direct,1))
            assert one(anymsg,1)==b'/lance.encodings21.PageLayout'
            layout=wire(one(anymsg,2));assert len(layout)==1
            tag=next(iter(layout));details=wire(one(layout,tag))
            page=dict(column=col,page=pi,row_start=page_rows,row_stop=page_rows+length,buffers=buffers,layout=tag)
            if tag==1:
                assert one(details,8,0)==0 and not one(details,10,0), 'Nested/large miniblock refused'
                assert packed(details,6)==[1] and 1 not in details and 2 not in details
                assert one(details,9)==length
                hasdict=4 in details
                assert len(buffers)==(3 if hasdict else 2)
                assert buffers[0][1]>buffers[0][0] and (buffers[0][1]-buffers[0][0])%2==0
                assert buffers[0][1]-buffers[0][0]<=256_000, 'Index buffer memory ceiling'
                page['index_range']=buffers[0];page['value_range']=buffers[1]
                page['dictionary_ranges']=buffers[2:] if hasdict else []
            elif tag==2:
                assert packed(details,5)==[1] and 7 not in details and 8 not in details
                # An older constant layout may store its scalar in a data buffer.
                # Those bytes belong to the later payload plan, never this reader.
                page['scalar_inline']=6 in details
            else:raise AssertionError(f'Unsupported structural layout {tag}; no payload fallback')
            result.append(page);page_rows+=length
        assert page_rows==rows, (col,page_rows,rows)
    # Page index ranges must never overlap any compressed value/scalar/dict range.
    index=[p['index_range'] for p in result if p['layout']==1]
    payload=[]
    for p in result:payload.extend(p['buffers'][1:] if p['layout']==1 else p['buffers'])
    assert all(y<=u or v<=x for x,y in index for u,v in payload)
    return result


def chunks(page,data):
    assert page['layout']==1
    assert len(data)==page['index_range'][1]-page['index_range'][0]
    words=struct.unpack('<'+'H'*(len(data)//2),data);assert words
    row=page['row_start'];byte=page['value_range'][0];out=[]
    for i,w in enumerate(words):
        log=w&15
        n=(1<<log) if i<len(words)-1 else page['row_stop']-row
        assert (i==len(words)-1 or log>0) and n>0
        if i==len(words)-1 and log:assert n==1<<log
        size=((w>>4)+1)*8
        assert row+n<=page['row_stop'] and byte+size<=page['value_range'][1]
        out.append([row,row+n,byte,byte+size]);row+=n;byte+=size
    assert row==page['row_stop'] and byte==page['value_range'][1]
    return out


def plan_payload(pages,logical,get_index):
    logical=merge(logical);payload=[];selected_chunks=0;overfetch_rows=[0,0,0]
    used_pages=[]
    for p in pages:
        lo,hi=p['row_start'],p['row_stop']
        wanted=merge([[max(lo,a),min(hi,b)] for a,b in logical if a<hi and b>lo])
        if not wanted:continue
        count=sum(b-a for a,b in wanted);used_pages.append([p['column'],p['page']])
        if p['layout']==1:
            cs=chunks(p,get_index(p))
            selected=[x for x in cs if any(a<x[1] and b>x[0] for a,b in wanted)]
            payload.extend([[x[2],x[3]] for x in selected]);payload.extend(p['dictionary_ranges'])
            selected_chunks+=len(selected);overfetch_rows[p['column']]+=sum(x[1]-x[0] for x in selected)-count
        else:payload.extend(p['buffers'])
    return dict(payload_ranges=merge(payload),used_pages=used_pages,selected_chunks=selected_chunks,
                encoded_row_overfetch_by_column=overfetch_rows,
                payload_bytes=sum(b-a for a,b in merge(payload)),
                logical_rows=sum(b-a for a,b in logical),coalescing_gap_bytes=0)
