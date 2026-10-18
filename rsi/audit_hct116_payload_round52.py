"""Independent scalar request-plan audit and bounded staged-record audit.

Does not import the payload producer or its planner. Final full-cell counts and
depths are compared with frozen metadata, with no source activation or scoring.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import argparse
import re
import tempfile
import runtime as rt
import numpy as np
import pandas as pd
from lance.file import LanceFileReader

ROOT=Path('artifacts/HCT116_payload_round52')
META=Path('artifacts/HCT116_encoded_metadata_round51')
PREFLIGHT=Path('artifacts/HCT116_extension_preflight_round49')


def load(path):
    return rt.json.loads(Path(path).read_text())


def accounting(path,limit):
    total=0
    bindings={}
    charges=[]
    for line in Path(path).read_text().splitlines():
        e=rt.json.loads(line)
        if e['event']=='bind':
            assert e['key'] not in bindings
            bindings[e['key']]=e['value']
        else:
            assert e['event']=='charge' and e['bytes']>0
            total+=e['bytes']
            assert total==e['total']<=limit
            charges.append(e)
    return total,bindings,charges


def union(intervals):
    # Independent event sweep (no producer merge/chunk functions).
    events=Counter()
    for a,b in intervals:
        assert a<b
        events[a]+=1
        events[b]-=1
    level=0
    result=[]
    for position,change in sorted(events.items()):
        old=level
        level+=change
        if old==0 and level>0:
            start=position
        if old>0 and level==0:
            result.append([start,position])
    assert level==0
    return result


def preflight():
    assert not (ROOT/'preflight_audit.json').exists()
    plan=load(ROOT/'plan.json')
    for path,h in plan['sha256'].items():
        assert rt.sha(path)==h,path
    before=ROOT.parent/'download_allocations_before_payload_round52.json'
    old=load(META/'plan.json')
    for path,h in old['sha256'].items():
        current=before if path=='artifacts/download_allocations.json' else Path(path)
        assert rt.sha(current)==h,path
    main=load(META/'audit.json')
    scope=load(META/'scope_audit.json')
    assert main['status'].startswith('PASS') and scope['status'].startswith('PASS')
    assert main['receipt_sha256']==scope['receipt_sha256']==rt.sha(META/'receipt.json')
    assert scope['main_audit_sha256']==rt.sha(META/'audit.json')
    assert rt.sha('audit_hct116_encoded_plan_round51.py')==main['verifier_sha256']
    total,bindings,_=accounting('artifacts/HCT116_encoded_metadata_round51_charges.jsonl',250_000_000)
    assert total==main['charged_bytes']
    files=load(ROOT/'files.json')
    r=load(META/'receipt.json')
    assert len(files)==len(r['file_plan_sha256'])==2018
    payload=0
    extra=0
    requests=0
    for path,h in r['file_plan_sha256'].items():
        assert rt.sha(path)==h
        saved=load(path)
        fid=saved['fragment_id']
        actual=files[str(fid)]
        assert actual['metadata_plan']==path and actual['metadata_plan_sha256']==h
        boundary=[]
        for p in saved['pages']:
            assert p['layout']==1
            if p['column']!=0:
                continue
            a,b=p['index_range']
            source=META/str(fid)/f'{a}-{b}.bin'
            assert bindings['cache:'+str(source)]==rt.sha(source)
            data=source.read_bytes()
            words=[int.from_bytes(data[i:i+2],'little') for i in range(0,len(data),2)]
            pos=p['value_range'][0]
            row=p['row_start']
            for wi,word in enumerate(words):
                span=2**(word%16) if wi<len(words)-1 else p['row_stop']-row
                size=(word//16+1)*8
                if row==0 or row+span==saved['physical_rows']:
                    boundary.append([pos,pos+size])
                pos+=size
                row+=span
            assert pos==p['value_range'][1] and row==p['row_stop']
        expected=union(saved['payload']['payload_ranges']+boundary)
        assert actual['boundary_cell_id_ranges']==union(boundary)
        assert actual['payload_ranges']==expected
        assert actual['physical_rows']==saved['physical_rows']
        size=sum(b-a for a,b in expected)
        assert actual['payload_bytes']==size
        payload+=size
        extra+=size-saved['payload']['payload_bytes']
        requests+=len(expected)
    assert payload==plan['encoded_payload_bytes']
    assert extra==plan['boundary_extra_bytes']
    assert plan['first_attempt_bytes']==payload+2018*64000
    assert plan['additional_failure_restart_envelope_bytes']==plan['allocation_bytes']-plan['first_attempt_bytes']>0
    a=load('artifacts/download_allocations.json')
    b=load(before)
    assert sum(a['new_allocations'].values())==a['total_allocated_upper_bound_bytes']<=15_000_000_000
    assert a['total_allocated_upper_bound_bytes']==b['total_allocated_upper_bound_bytes']
    reallocation=a['reallocations'][-1]
    assert reallocation['before_sha256']==rt.sha(before)
    assert reallocation['bytes']==a['new_allocations']['HCT116_payload_round52']==1_200_000_000
    for d in reallocation['from_readers']:
        assert rt.sha(d['journal'])==d['journal_sha256']
        used,_,_=accounting(d['journal'],d['before_ceiling'])
        assert used==d['charged_bytes']<=d['after_ceiling']
        assert b['new_allocations'][d['reader']]==d['before_ceiling']
        assert a['new_allocations'][d['reader']]==d['after_ceiling']==d['before_ceiling']-400_000_000
    fixture=load(ROOT/'fixture.json')
    assert fixture['status'].startswith('PASS') and fixture['plan_sha256']==rt.sha(ROOT/'plan.json')
    frozen=load('artifacts/frozen_metadata_index.json')
    assert rt.sha(frozen['archive'])==frozen['archive_sha256']
    for p,e in frozen['entries'].items():
        assert rt.sha(p)==e['sha256']
    assert len(frozen['entries'])==38
    assert load('artifacts/calibration_gate.json')['status']=='CLOSED'
    rt.save(ROOT/'preflight_audit.json',dict(status='PASS independent scalar payload-cost/allocation/identity audit',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(ROOT/'plan.json'),
        verifier_sha256=rt.sha(__file__),files=2018,exact_encoded_requests=requests,
        encoded_payload_bytes=payload,boundary_extra_bytes=extra,first_attempt_bytes=payload+2018*64000,
        additional_failure_restart_envelope_bytes=plan['additional_failure_restart_envelope_bytes'],
        fixture_sha256=rt.sha(ROOT/'fixture.json'),all38_frozen_members_unchanged=True,
        source_activated=False,holdout_scored=False))
    print('PASS independent preflight',payload,'payload bytes',requests,'requests',flush=True)


def complete():
    assert not (ROOT/'audit.json').exists()
    plan=load(ROOT/'plan.json')
    for p,h in plan['sha256'].items():
        assert rt.sha(p)==h,p
    receipt=load(ROOT/'receipt.json')
    assert receipt['plan_sha256']==rt.sha(ROOT/'plan.json')
    total,bindings,charges=accounting(receipt['journal'],plan['allocation_bytes'])
    assert total==receipt['charged_bytes']
    assert bindings['plan_sha256']==rt.sha(ROOT/'plan.json')
    files=load(ROOT/'files.json')
    fragments={f['id']:f for f in load('artifacts/HCT116_expression_identity_round50/fragments.json')}
    paths={fragments[int(fid)]['inventory_entry']['path']:f for fid,f in files.items()}
    _,mb,_=accounting('artifacts/HCT116_encoded_metadata_round51_charges.jsonl',250_000_000)
    heads=ranges=0
    for e in charges:
        if e['reason'].startswith('Orion HEAD metadata upper bound '):
            path=e['reason'].removeprefix('Orion HEAD metadata upper bound ')
            assert path in paths and e['bytes']==64000
            heads+=1
        else:
            m=re.fullmatch(r'Exact payload range attempt (.+) \[(\d+),(\d+)\)',e['reason'])
            assert m
            path,a,b=m[1],int(m[2]),int(m[3])
            assert [a,b] in paths[path]['payload_ranges'] and e['bytes']==b-a
            ranges+=1
    spans=pd.read_parquet(PREFLIGHT/'logical_record_ranges.parquet')
    by_file={int(fid):rows for fid,rows in spans.groupby('fragment_id',sort=True)}
    counts=Counter()
    depths=Counter()
    records=0
    offsets=np.load('artifacts/HCT116_offset_mapping.npz')
    for fid,rows in by_file.items():
        path=ROOT/'decoded'/f'{fid}.npz'
        local=load(path.with_suffix('.json'))
        assert local==receipt['decoded'][str(fid)]
        assert rt.sha(path)==local['sha256']==bindings['decoded:'+str(path)]
        fileplan=files[str(fid)]
        assert local['fileplan_sha256']==rt.hashlib.sha256(rt.json.dumps(fileplan,sort_keys=True).encode()).hexdigest()
        entry=fragments[fid]['inventory_entry']
        assert bindings['file:'+entry['path']]==mb['file:'+entry['path']]
        for a,b in fileplan['payload_ranges']:
            cache=ROOT/'ranges'/str(fid)/f'{a}-{b}.bin'
            assert cache.stat().st_size==b-a
            assert rt.sha(cache)==bindings['payload_cache:'+str(cache)]
        with np.load(path) as d:
            np.testing.assert_array_equal(d['record_start'],rows.record_start)
            np.testing.assert_array_equal(d['record_stop'],rows.record_stop)
            np.testing.assert_array_equal(d['span_cell_id'],rows.cell_integer_id)
            cells=d['cell_integer_id']
            genes=d['gene_integer_id']
            values=d['value']
            assert cells.dtype==np.dtype('uint32')
            assert genes.dtype==values.dtype==np.dtype('uint16')
            n=int((rows.record_stop-rows.record_start).sum())
            assert len(values)==len(genes)==len(cells)==n==local['records']
            assert np.all(genes<38606) and np.all(values>0)
            # Independent full replay from journal-bound encoded bytes, using
            # the pinned official decoder without importing producer helpers.
            with tempfile.TemporaryDirectory(prefix='independent_payload_r52_',dir='artifacts/tmp') as temp:
                sparse=Path(temp)/'sparse.lance'
                with sparse.open('wb') as out:
                    out.truncate(entry['size'])
                    for p in (META/str(fid)).glob('*.bin'):
                        a,b=map(int,p.stem.split('-'))
                        assert mb['cache:'+str(p)]==rt.sha(p) and p.stat().st_size==b-a
                        out.seek(a)
                        out.write(p.read_bytes())
                    for a,b in fileplan['payload_ranges']:
                        out.seek(a)
                        out.write((ROOT/'ranges'/str(fid)/f'{a}-{b}.bin').read_bytes())
                official=LanceFileReader(str(sparse))
                boundary=LanceFileReader(str(sparse),columns=['cell_integer_id'])
                for position,expected in [(0,offsets['first_cell'][fid]),
                        (fileplan['physical_rows']-1,offsets['last_cell'][fid])]:
                    assert boundary.read_range(position,1).to_table()['cell_integer_id'][0].as_py()==expected
                cursor=0
                for row in rows.itertuples(index=False):
                    length=int(row.record_stop-row.record_start)
                    stop=cursor+length
                    table=official.read_range(int(row.record_start),length).to_table()
                    for key,array in [('cell_integer_id',cells),('gene_integer_id',genes),('value',values)]:
                        np.testing.assert_array_equal(table[key].to_numpy(),array[cursor:stop])
                    assert np.all(cells[cursor:stop]==row.cell_integer_id)
                    counts[int(row.cell_integer_id)]+=length
                    depths[int(row.cell_integer_id)]+=int(values[cursor:stop].sum(dtype=np.uint64))
                    cursor=stop
                assert cursor==n
            records+=n
        if fid%100==0:
            print('independent staged-cell audit',fid,'records',records,flush=True)
    selected=pd.read_parquet(PREFLIGHT/'selection.parquet')
    assert set(counts)==set(depths)==set(map(int,selected.cell_integer_id))
    for row in selected.itertuples(index=False):
        assert counts[row.cell_integer_id]==row.n_genes_by_counts
        assert depths[row.cell_integer_id]==row.total_counts
    assert records==215_800_779 and len(by_file)==receipt['files']==2018
    for fid in by_file:
        local=receipt['decoded'][str(fid)]
        assert local['first_cell_verified']==offsets['first_cell'][fid]
        assert local['last_cell_verified']==offsets['last_cell'][fid]
    rt.save(ROOT/'audit.json',dict(status='PASS independent staged record/count/depth/accounting audit; moments not built',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(ROOT/'plan.json'),
        receipt_sha256=rt.sha(ROOT/'receipt.json'),verifier_sha256=rt.sha(__file__),
        files=len(by_file),selected_cells=len(selected),logical_records=records,
        selected_cell_ids_counts_and_full_depths_verified=True,boundary_receipts_verified=True,
        independent_official_redecode_all_selected_records=True,
        charged_bytes=total,head_attempts=heads,payload_range_attempts=ranges,
        source_activated=False,moments_verified=False,holdout_scored=False,
        full_partial_file_digest_verified=False,phase2_allowed=False))
    print('PASS independent staged-cell audit; moments unbuilt',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--complete',action='store_true')
    args=parser.parse_args()
    complete() if args.complete else preflight()
