"""Independent streamed counters and scalar fragment-range replay; no model access."""
from bisect import bisect_right
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import lance
import numpy as np
import pandas as pd

ROOT=Path('artifacts/HCT116_extension_preflight_round49')
META=Path('artifacts/HCT116_full_cells_metadata_round48')


def load(p):
    return rt.json.loads(Path(p).read_text())


def main():
    out=ROOT/'audit.json';assert not out.exists()
    audit_plan=load(ROOT/'audit_plan.json')
    for p,h in audit_plan['sha256'].items():assert rt.sha(p)==h,p
    plan=load(ROOT/'plan_retry.json');result=load(ROOT/'preflight.json')
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p
    assert result['plan_sha256']==rt.sha(ROOT/'plan_retry.json')
    for p,h in result['artifact_sha256'].items():assert rt.sha(ROOT/p)==h,p
    ends=[rt.json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert [s for s in ends if s.get('name')=='preflight_HCT116_extension_round49_retry' and s['event']=='end'][-1]['rc']==0
    manifest=load('artifacts/tasks/manifest.json');panel=manifest['CD4_depth192_search']['targets']
    with np.load('artifacts/cd4_banks/HCT116.npz') as d:
        old=dict(zip(d['targets'].tolist(),d['n_cells'].tolist()))
    with np.load('artifacts/HCT116_offset_mapping.npz') as d:
        logical=d['logical_start'];original=d['original_start'];lengths=d['block_length']
    ds=lance.dataset(str((META/'release/cells.lance').resolve()))
    count=Counter();records=Counter();sample_count=Counter();ctrl=Counter()
    selected_ids={};last_end=0;last_id=-1;full=0;block_nn=Counter()
    panel_set=set(panel);eligible_ids={t:[] for t in panel}
    # Deliberately scalar streaming implementation, distinct from builder pandas selection.
    for batch in ds.to_batches(columns=plan['columns'],batch_size=65536):
        cols=batch.to_pydict()
        for vals in zip(*(cols[k] for k in plan['columns'])):
            cid,t,s,passed,start,nn,depth=vals
            assert cid==last_id+1 and start==last_end and nn>0 and depth>0
            last_id=cid;last_end=start+nn;full+=1;block_nn[cid//5000]+=nn
            if cid%5000==0:assert start==int(original[cid//5000])
            if not passed:continue
            present=int(logical[cid//5000])>=0
            role=None
            if t=='Non-Targeting':
                ctrl[s]+=1
                if ctrl[s]<=250 and present:role='source_control'
            if present:
                count[t]+=1;records[t]+=nn;sample_count[t,s]+=1
                if t in panel_set:eligible_ids[t].append(cid)
            if role or (t in panel_set and old.get(t,0)<20 and present):
                selected_ids[cid]=(t,s,start,nn,depth,role or 'append_target')
    assert full==result['full_metadata_cells']==3409169
    for i in range(len(logical)):
        if int(logical[i])>=0:assert block_nn[i]==int(lengths[i])
    assert set(i for i in range(len(logical)) if int(logical[i])<0)=={374}
    for t,n in old.items():assert count[t]==n
    actions={t:('retain_existing' if old.get(t,0)>=20 else 'append' if count[t]>=20 else 'abstain_below20') for t in panel}
    assert [t for t in panel if actions[t]=='append']==result['append_targets']
    assert not set(result['append_targets']) & set(manifest['H1_holdout']['targets'])
    assert [t for t in panel if actions[t]=='abstain_below20']==['HMGXB4'] and count['HMGXB4']==2
    selected_ids={cid:v for cid,v in selected_ids.items() if v[-1]=='source_control' or actions[v[0]]=='append'}
    sel=pd.read_parquet(ROOT/'selection.parquet')
    assert sel.cell_integer_id.is_unique and set(sel.cell_integer_id)==set(selected_ids)
    by_sample=Counter();by_target=Counter();total_nnz=Counter();metadata_by_id={}
    for r in sel.itertuples(index=False):
        assert selected_ids[r.cell_integer_id]==(r.gene_target,r.sample,r.cell_start_index,r.n_genes_by_counts,r.total_counts,r.source_role)
        by_sample[r.source_role,r.sample]+=1;by_target[r.gene_target]+=1
        total_nnz[r.source_role]+=int(r.n_genes_by_counts);metadata_by_id[r.cell_integer_id]=r
    assert len([s for role,s in by_sample if role=='source_control'])==109
    expected_control_counts=Counter(v[1] for v in selected_ids.values() if v[-1]=='source_control')
    assert {s:n for (role,s),n in by_sample.items() if role=='source_control'}==dict(expected_control_counts)
    assert all(0<n<=250 for n in expected_control_counts.values())
    assert 250*109-sum(expected_control_counts.values())==result['omitted_control_cells_without_expression']==82
    for t in result['append_targets']:
        assert by_target[t]==count[t] and count[t]>=20
        assert all(by_sample['source_control',s]>0 for tt,s in sample_count if tt==t)
    targets=pd.read_csv(ROOT/'targets.tsv',sep='\t')
    assert targets.target.tolist()==panel
    for r in targets.itertuples(index=False):
        assert r.old_bank_cells==old.get(r.target,0) and r.block_present_cells==count[r.target]
        assert r.expression_records==records[r.target] and r.action==actions[r.target]
        assert r.source_samples==sum(1 for t,s in sample_count if t==r.target)
    cross=pd.read_csv(ROOT/'target_sample_counts.tsv',sep='\t',index_col=0)
    for t in cross.index:
        for s in cross.columns:assert int(cross.loc[t,s])==sample_count[t,s]
    frags=load('../model_v4/artifacts/HCT116_fragments.json')
    prefix=[0]
    for f in frags:prefix.append(prefix[-1]+f['physical_rows'])
    expected=[]
    for cid,r in metadata_by_id.items():
        group=cid//5000;pos=int(r.cell_start_index-original[group]+logical[group]);end=pos+int(r.n_genes_by_counts)
        assert int(logical[group])<=pos<end<=int(logical[group]+lengths[group])
        while pos<end:
            k=bisect_right(prefix,pos)-1;stop=min(end,prefix[k+1])
            expected.append((cid,frags[k]['id'],frags[k]['path'],pos-prefix[k],stop-prefix[k]));pos=stop
    ranges=pd.read_parquet(ROOT/'logical_record_ranges.parquet')
    actual=list(ranges.itertuples(index=False,name=None))
    assert sorted(actual)==sorted(expected)
    assert sum(e[-1]-e[-2] for e in expected)==sum(total_nnz.values())==result['total_expression_records']
    assert total_nnz['append_target']==result['target_expression_records']
    assert total_nnz['source_control']==result['control_expression_records']
    assert by_target['Non-Targeting']==result['control_cells']
    assert len(sel)-by_target['Non-Targeting']==result['target_cells']
    assert len(set(r[1] for r in expected))==result['touched_historical_fragments']
    assert result['exact_http_payload_bytes'] is None and not result['expression_identity_freshly_verified']
    frozen=load('artifacts/frozen_metadata_index.json')
    for p,e in frozen['entries'].items():assert rt.sha(p)==e['sha256'],p
    assert rt.sha(frozen['archive'])==frozen['archive_sha256']
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p
    assert load('artifacts/calibration_gate.json')['status']=='CLOSED'
    obj=dict(status='PASS metadata selection and conditional logical ranges; HTTP cost UNVERIFIED',
        utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=rt.sha(__file__),
        audit_plan_sha256=rt.sha(ROOT/'audit_plan.json'),
        plan_sha256=rt.sha(ROOT/'plan_retry.json'),preflight_sha256=rt.sha(ROOT/'preflight.json'),
        panel_targets=79,retained_supported_targets=6,append_targets=len(result['append_targets']),
        source_supported_targets=6+len(result['append_targets']),abstention=dict(target='HMGXB4',cells=2,minimum=20),
        selected_cells=len(sel),target_cells=result['target_cells'],control_cells=result['control_cells'],
        logical_expression_records=sum(total_nnz.values()),touched_historical_fragments=result['touched_historical_fragments'],
        all_original_bank_counts_replayed=437,first250_before_block_filter_all109_samples=True,omitted_control_cells_without_expression=82,
        independent_scalar_fragment_spans=len(expected),frozen_metadata_members_verified=len(frozen['entries']),
        expression_identity_freshly_verified=False,exact_http_payload_bytes=None,
        network_io=False,expression_payloads_read=False,source_activated=False,holdout_scored=False,phase2_allowed=False)
    text=rt.json.dumps(obj,indent=2,allow_nan=False)+'\n'
    with out.open('x') as f:f.write(text)
    print(obj['status'],flush=True)


if __name__=='__main__':main()
