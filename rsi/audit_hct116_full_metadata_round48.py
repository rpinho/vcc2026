"""Independent immutable full-release census audit; counts and identities only."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import runtime as rt
import lance
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT=Path('artifacts/HCT116_full_cells_metadata_round48')
PLAN=ROOT/'audit_plan_retry.json'
OUT=ROOT/'audit_retry.json'
TABLE=ROOT/'target_counts_retry.tsv'
REVISION='598aa5442a4ee9a8f1a383a843ab5ce142904d04'
COLS=['cell_integer_id','gene_target','sample','pass_guide_filter']


def load(p):
    return rt.json.loads(Path(p).read_text())


def save_new(p,v):
    text=rt.json.dumps(v,indent=2,allow_nan=False)+'\n'
    with p.open('x') as f:f.write(text)


def bindings(plan):
    for p,h in plan['sha256'].items():assert rt.sha(p)==h,p


def freeze():
    assert not PLAN.exists() and not OUT.exists() and not TABLE.exists()
    receipt=load(ROOT/'receipt.json')
    assert receipt['revision']==REVISION
    paths=[Path(__file__),Path('acquire_hct116_full_metadata_round48.py'),
        Path('verify_hct116_full_metadata_reader_round48.py'),ROOT/'inventory.json',ROOT/'receipt.json',
        Path('artifacts/HCT116_full_cells_metadata_round48_charges.jsonl'),
        Path('artifacts/HCT116_full_metadata_reader_fixture_round48.json'),
        Path('artifacts/HCT116_census_scope_audit_round48.json'),
        Path('../model_v4/artifacts/HCT116_cells.parquet'),
        Path('../model_v4/artifacts/HCT116_offset_mapping.npz'),
        Path('../model_v4/artifacts/targets.json'),Path('artifacts/GWPS_raw_obs.parquet'),
        Path('artifacts/CD4_cell_census.parquet'),Path('artifacts/cd4_sources/K562_GWPS_bulk.npz'),
        Path('artifacts/cd4_banks/K562_GWPS_raw.npz'),
        Path('artifacts/calibration_gate.json'),Path('artifacts/benchmark_roles.json'),
        Path('artifacts/tasks/manifest.json'),Path('artifacts/download_allocations.json'),
        Path('artifacts/frozen_metadata_index.json')]
    paths += [Path(v['path']) for v in receipt['files'].values()]
    paths += sorted(Path('artifacts/scores').glob('*/*_s42.json'))
    save_new(PLAN,dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Complete pinned HCT116 cell metadata only; no expression or response metric reads. Independently verify release files and row counters, reconcile original selected cache, quantify omitted target coverage and necessary precision bounds. No target selection or source/task activation.',
        columns=COLS, source_members=['targets','n_cells'],
        fixed_populations='HCT116 pass-guide cells: all observed and full cached raw GWPS>=200 minus sealed H1 target names. k79/300 common-sample necessary depth ceilings, unlimited controls. CD4: existing79 depth192 targets plus original round44 symmetric prepared GWPS200 superset, with full HCT>=20. Retain all24 lanes and H1 target-name exclusion in conservative scope; separately report relaxed superset. No proposed panel or lower-depth retry.',
        limits='File identities prove present immutable release, not its original biological completeness. Cell metadata existence does not prove expression-block support; separately flag absent original682 expression blocks without asserting new expression identity. Bounds are necessary and not joint feasibility. Target support is not source activation.',
        sha256={str(p):rt.sha(p) for p in paths}))
    print('FROZEN',PLAN,flush=True)


def upper_bound(table,k):
    if len(table)<k:return dict(status='INSUFFICIENT TARGETS',targets=len(table))
    a=table.to_numpy(np.int64)
    values=np.partition(a,len(a)-k,axis=0)[len(a)-k]
    independent=[sorted(table[c].tolist(),reverse=True)[k-1] for c in table.columns]
    assert values.tolist()==independent
    return dict(targets=len(table),panel_size=k,unlimited_control_depth_upper_bound=int(values.sum()),
        per_sample_bounds=dict(zip(map(str,table.columns),map(int,values))),
        excludes192=bool(values.sum()<192),excludes400=bool(values.sum()<400))


def main():
    assert not OUT.exists() and not TABLE.exists()
    plan=load(PLAN);bindings(plan)
    receipt=load(ROOT/'receipt.json');inv=load(ROOT/'inventory.json')
    assert receipt['revision']==inv['revision']==REVISION and inv['complete']
    entries={e['path']:e for e in inv['entries'] if e['type']=='file' and '/_transactions/' not in e['path']}
    assert set(entries)==set(receipt['files'])
    for name,meta in receipt['files'].items():
        p=Path(meta['path']);e=entries[name]
        assert p.stat().st_size==e['size']==meta['size'] and rt.sha(p)==meta['sha256']
        if e.get('lfs'):assert rt.sha(p)==e['lfs']['oid']
        else:
            data=p.read_bytes()
            assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['oid']
    fixture=load('artifacts/HCT116_full_metadata_reader_fixture_round48.json')
    assert fixture['status'].startswith('PASS') and fixture['reader_sha256']==rt.sha('acquire_hct116_full_metadata_round48.py')==receipt['reader_sha256']
    events=[rt.json.loads(x) for x in Path(receipt['journal']).read_text().splitlines()]
    charged=sum(e['bytes'] for e in events if e['event']=='charge')
    assert charged==receipt['charged_bytes']<=180_000_000
    ends=[rt.json.loads(x) for x in Path('ledger/steps.jsonl').read_text().splitlines()]
    assert [e for e in ends if e['event']=='end' and e['name']=='acquire_HCT116_full_metadata_round48'][-1]['rc']==0

    dataset=lance.dataset(str((ROOT/'release/cells.lance').resolve()))
    frame=dataset.to_table(columns=COLS).to_pandas()
    assert len(frame)==dataset.count_rows() and not frame.isna().any().any()
    ids=frame.cell_integer_id.to_numpy(np.int64)
    assert np.array_equal(ids,np.arange(len(frame))), 'Cell row/index convention changed'
    frame.gene_target=frame.gene_target.astype(str);frame['sample']=frame['sample'].astype(str)
    counter=Counter()
    for batch in dataset.to_batches(columns=['gene_target','sample','pass_guide_filter'],batch_size=65536):
        counter.update(zip(*(batch.column(i).to_pylist() for i in range(3))))
    assert sum(counter.values())==len(frame)
    tab=frame.loc[frame.pass_guide_filter].groupby(['gene_target','sample'],observed=True).size().unstack(fill_value=0).sort_index().sort_index(axis=1)
    for t,row in tab.iterrows():assert row.tolist()==[counter[(t,b,True)] for b in tab.columns]
    assert int(tab.to_numpy().sum())==int(frame.pass_guide_filter.sum())
    old=pd.read_parquet('../model_v4/artifacts/HCT116_cells.parquet',columns=COLS)
    match=frame.iloc[old.cell_integer_id.to_numpy(np.int64)]
    for col in COLS:
        assert match[col].astype(str).tolist()==old[col].astype(str).tolist(),col
    labels=set(load('../model_v4/artifacts/targets.json'))|{'Non-Targeting'}
    query=frame.loc[frame.gene_target.isin(labels)]
    assert set(query.cell_integer_id)==set(old.cell_integer_id), 'Original selected cache differs at pinned release'
    controls=tab.loc['Non-Targeting'];pert=tab.drop(index='Non-Targeting')
    by_target=pert.sum(axis=1)
    with np.load('../model_v4/artifacts/HCT116_offset_mapping.npz',allow_pickle=False) as d:
        block=ids//5000
        assert block.max()<len(d['logical_start'])
        present=d['logical_start'][block]>=0
        expression_version=int(d['version'])
    supported=frame.loc[present & frame.pass_guide_filter].groupby('gene_target',observed=True).size()
    old_counts=old.groupby('gene_target',observed=True).size()
    full_gwps=pd.read_parquet('artifacts/GWPS_raw_obs.parquet',columns=['gene']).gene.value_counts()
    gwps_ge200=set(full_gwps[full_gwps>=200].index)-{'non-targeting'}
    manifest=load('artifacts/tasks/manifest.json')
    sealed=set(manifest['H1_holdout']['targets'])
    hct_scopes={}
    for name,ts in dict(all_pass=set(pert.index),full_raw_gwps200_unsealed=set(pert.index)&gwps_ge200-sealed).items():
        subset=pert.loc[sorted(ts)]
        hct_scopes[name]=dict(targets=len(ts),targets_ge400=int((subset.sum(axis=1)>=400).sum()),
            panels={str(k):upper_bound(subset,k) for k in [79,300]})
    active=manifest['CD4_depth192_search']['targets']
    cd4_active=dict(targets=len(active),cached_hct_ge20=int(sum(old_counts.get(t,0)>=20 for t in active)),
        full_hct_ge20=int(sum(by_target.get(t,0)>=20 for t in active)),
        full_hct_original_expression_blocks_ge20=int(sum(supported.get(t,0)>=20 for t in active)))
    cd4=pd.read_parquet('artifacts/CD4_cell_census.parquet',columns=['target_gene','batch'])
    cd4tab=pd.crosstab(cd4.target_gene,cd4.batch).drop(index='non-targeting')
    assert len(cd4tab.columns)==24
    prepared=[]
    for path in ['artifacts/cd4_sources/K562_GWPS_bulk.npz','artifacts/cd4_banks/K562_GWPS_raw.npz']:
        with np.load(path,allow_pickle=False) as d:
            prepared.append({str(t) for t,n in zip(d['targets'],d['n_cells']) if n>=200})
    symmetric=set.intersection(*prepared)&set(cd4tab.index)
    hct_ge20={t for t,n in by_target.items() if n>=20}
    cd4_scopes={}
    for name,pool in dict(conservative=symmetric-sealed-set(manifest['H1_rich_search']['targets']),relaxed=symmetric-sealed).items():
        ts=sorted(pool&hct_ge20)
        subset=cd4tab.loc[ts]
        cd4_scopes[name]=dict(targets=len(ts),targets_all24_lanes_ge8=int((subset>=8).all(axis=1).sum()),
            panel79=upper_bound(subset,79))
    rows=pd.DataFrame([dict(target=t,full_guide_pass_cells=int(by_target[t]),
        original_cache_cells=int(old_counts.get(t,0)),queried=t in labels,
        original_expression_block_cells=int(supported.get(t,0)),full_raw_gwps_ge200=t in gwps_ge200,
        cd4_active_target=t in active) for t in sorted(pert.index)])
    with TABLE.open('x') as f:rows.to_csv(f,sep='\t',index=False)
    pd.testing.assert_frame_equal(rows,pd.read_csv(TABLE,sep='\t'))
    index=load('artifacts/frozen_metadata_index.json')
    for p,m in index['entries'].items():assert rt.sha(p)==m['sha256'],p
    assert rt.sha(index['archive'])==index['archive_sha256']
    bindings(plan)
    save_new(OUT,dict(status='PASS immutable full cells-metadata certificate only',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(PLAN),revision=REVISION,
        local_lance_version=int(dataset.version),schema=str(dataset.schema),rows=len(frame),
        guide_pass_rows=int(frame.pass_guide_filter.sum()),guide_fail_rows=int((~frame.pass_guide_filter).sum()),
        perturbation_targets=len(pert),samples=len(tab.columns),controls=int(controls.sum()),
        cache_exact_query_replay=True,cache_rows=len(old),added_rows=len(frame)-len(old),
        queried_label_scope=len(labels),novel_observed_target_names=len(set(pert.index)-labels),
        full_metadata_sample_names_missing_from_cache=sorted(set(frame['sample'])-set(old['sample'])),
        original_expression_block_scope=dict(version=expression_version,missing_metadata_rows=int((~present).sum()),
            status='Historical offset-block membership only; no new expression identity or payload verification'),
        fixed_cd4_panel_support=cd4_active,hct_precision_bounds=hct_scopes,cd4_source_comparability=cd4_scopes,
        charged_bytes=charged,allocation_bytes=180_000_000,table_sha256=rt.sha(TABLE),
        frozen_metadata_members_verified=len(index['entries']),network_io=False,expression_payloads_read=False,
        biological_scores_read=False,source_activated=False,task_activated=False,holdout_scored=False,phase2_allowed=False))
    print('PASS full metadata audit',len(frame),'cells',len(pert),'targets',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');args=parser.parse_args()
    freeze() if args.freeze else main()
