"""Merge independent fetch shards and emit single logical screen response banks."""
import common as job
import io_guard;io_guard.install()
import argparse,numpy as np,pandas as pd,shutil,gc
from prepare import commit

ROW_KEYS=['target_count_sums','target_mean_cpm','matched_control_probability','matched_control_mean_cpm',
          'delta_cpm_variance','target_mean_cpm_variance','matched_control_mean_cpm_variance','n_cells']

def merge_bank(paths,dest):
    # Shared rows first, then new shards; baseline target rows always win.
    first=None;seen=set();pieces=[]
    for path in paths:
        with np.load(path) as z:d={k:z[k] for k in z.files}
        if first is None:first=d
        assert np.array_equal(first['genes'],d['genes']),str(path)
        if 'conditions' in first:assert np.array_equal(first['conditions'],d['conditions']),str(path)
        if 'measured_genes' in first:assert np.array_equal(first['measured_genes'],d['measured_genes']),str(path)
        keep=np.array([str(t) not in seen for t in d['targets']]);seen.update(map(str,d['targets']))
        pieces.append((path,keep,int(keep.sum())))
        del d
    nt=sum(n for _,_,n in pieces);tmp=job.ROOT/'artifacts/tmp/merge';tmp.mkdir(parents=True,exist_ok=True)
    merged={};ng=len(first['genes'])
    if 'conditions' in first:
        keys=['targets','log2fc','lfcSE','adjusted_p','available','quality_pass','n_cells']
        for k in first:
            if k not in keys:merged[k]=first[k]
        for k in keys:
            vals=[]
            for p,keep,_ in pieces:
                with np.load(p) as z:vals.append(z[k][keep] if k=='targets' else z[k][:,keep])
            merged[k]=np.concatenate(vals,axis=0 if k=='targets' else 1)
    else:
        for k,v in first.items():
            if k not in ROW_KEYS+['targets']:merged[k]=v
        for k in ['targets']+[k for k in ROW_KEYS if k in first]:
            sample=first[k];shape=(nt,)+sample.shape[1:]
            if sample.dtype.kind in 'US':arr=np.empty(shape,dtype='<U100')
            else:arr=np.lib.format.open_memmap(tmp/f'{k}.npy',mode='w+',shape=shape,dtype=sample.dtype)
            off=0
            for p,keep,n in pieces:
                with np.load(p) as z:
                    if k not in z:
                        # Upstream older shared banks may omit component variance, so omit that field.
                        break
                    arr[off:off+n]=z[k][keep]
                off+=n
            if off==nt:merged[k]=arr
    partial=dest.with_name(dest.stem+'.partial.npz');np.savez_compressed(partial,**merged);partial.replace(dest)
    del merged,first;gc.collect()
    for p in tmp.glob('*.npy'):p.unlink()
    print('MERGED MOMENT',dest.name,nt,flush=True)

def convert(path,dest):
    with np.load(path) as z:
        name=path.stem;family='CD4' if name.startswith('CD4_') else 'K562' if name.startswith('K562') else name
        oldroot=job.ROOT
        # commit uses common.ROOT. Keep the exact formula, only direct its output to dest.
        dest.mkdir(parents=True,exist_ok=True)
        original_data=oldroot/'data'
        staging=dest.parent/('convert_'+dest.name);(staging/'data').mkdir(parents=True,exist_ok=True)
        job.ROOT=staging
        try:
            commit(name,family,z['targets'],z['genes'],z['target_count_sums'],z['matched_control_probability'],z['n_cells'],
                dict(path=str(path),sha256=job.sha(path),transform='unchanged model_billm.prepare.commit',
                     control='upstream matched_control_probability; CD4 pseudobulk approximation retained'))
        finally:job.ROOT=oldroot
        for p in (staging/'data').iterdir():p.replace(dest/p.name)
        # Native source gene list contains extra zero-filled coordinates for Orion. Preserve genuine measured masks.
        out=dest/(name+'.npz')
        with np.load(out) as r:d={k:r[k] for k in r.files}
        if 'measured_genes' in z:
            lookup=pd.Index(z['genes'].astype(str)).get_indexer(d['genes']);ok=lookup>=0
            mask=np.zeros(len(lookup),bool);mask[ok]=z['measured_genes'][lookup[ok]]
            d['measured'] &= mask;d['response'][:,~d['measured']]=0;d['control'][:,~d['measured']]=0
        np.savez_compressed(out,**d)
        job.save(out.with_suffix('.json'),dict(rows=len(d['targets']),family=family,sha256=job.sha(out),
            bytes=out.stat().st_size,measured_genes=int(d['measured'].sum()),moment_bank_sha256=job.sha(path)))

def consolidate_native(prefix):
    paths=sorted((job.ROOT/'data').glob(prefix+'_full_*.npz'))
    if not paths:return
    out=job.ROOT/f'data/{prefix}_full.npz';meta=out.with_suffix('.json')
    if out.exists() and meta.exists():
        done=job.json.loads(meta.read_text())
        if done.get('consolidation_complete') and done['sha256']==job.sha(out):
            for p in paths:p.unlink();p.with_suffix('.json').unlink(missing_ok=True)
            print('REUSE CONSOLIDATED',prefix,done['rows'],flush=True);return
    seen=set();pieces=[];base={};source_receipts=[]
    if prefix=='K562_GWPS':
        # Preserve raw matched GWPS and essential K562 rows over globally controlled full-screen duplicates.
        for name in ['K562_GWPS_raw','K562']:
            with np.load(job.ROOT/f'data/{name}.npz') as z:seen.update(z['targets'].astype(str))
        from response_union import raw_targets
        seen.update(raw_targets(job.ROOT/'data'))
    else:
        cached=job.ROOT/'data/H1_cached.npz'
        if cached.exists():
            with np.load(cached) as z:seen.update(z['targets'].astype(str))
    for p in paths:
        if p.with_suffix('.json').exists():source_receipts.append(job.json.loads(p.with_suffix('.json').read_text()))
        with np.load(p) as z:
            keep=~np.isin(z['targets'],list(seen));pieces.append((p,keep,int(keep.sum())))
            seen.update(z['targets'].astype(str));base={k:z[k] for k in ['genes','measured','family']}
    nt=sum(n for _,_,n in pieces);arrays={**base};tmp=job.ROOT/'artifacts/tmp'
    for key,dt,shape in [('response','f4',(nt,len(base['genes']))),('control','f4',(nt,len(base['genes']))),('n_cells','i8',(nt,)),('targets','<U100',(nt,))]:
        arr=np.lib.format.open_memmap(tmp/f'consolidate_{key}.npy',mode='w+',dtype=dt,shape=shape);off=0
        for p,keep,n in pieces:
            with np.load(p) as z:arr[off:off+n]=z[key][keep]
            off+=n
        arrays[key]=arr
    partial=out.with_suffix('.partial.npz');np.savez_compressed(partial,**arrays);partial.replace(out);del arrays;gc.collect()
    job.save(out.with_suffix('.json'),dict(name=out.stem,rows=nt,family=str(base['family']),sha256=job.sha(out),bytes=out.stat().st_size,
        response='log1p(50000 * native target pooled probability) minus global non-targeting pooled control; NOT batch matched',
        provenance='Original cached matched rows retained separately; disjoint newly summarized source targets',source_shard_receipts=source_receipts,
        consolidation_complete=True))
    for p in tmp.glob('consolidate_*.npy'):p.unlink()
    for p in paths:p.unlink();p.with_suffix('.json').unlink(missing_ok=True)
    print('CONSOLIDATED',prefix,nt,flush=True)

def complete_local_targets():
    # One source-label discrepancy exists: native ZNF735 versus bulk ZNF735_NM.
    # Preserve labels and fill any eligible bulk-only target through the explicitly permitted fallback.
    wanted=set(pd.read_csv(job.ROOT/'targets_local.csv').target_gene.astype(str));have=set()
    for p in (job.ROOT/'data').glob('*.npz'):
        with np.load(p) as z:
            if str(z['family'])=='K562':have.update(z['targets'].astype(str))
    missing=sorted(wanted-have)
    if missing:
        from source_runtime import readers
        m=readers();moment=job.ROOT/'artifacts/local_bulk_fallback/K562_GWPS_bulk.npz';moment.parent.mkdir(exist_ok=True)
        bank=m.local_bank('K562_GWPS_bulk',missing);m.savez(moment,**bank)
        fresh=moment.parent/'responses';convert(moment,fresh)
        from response_union import append
        append(job.ROOT/'data','K562_GWPS_bulk_missing',fresh/moment.name)
        have.update(bank['targets'][bank['n_cells']>=30].astype(str))
        job.save(moment.with_suffix('.json'),dict(targets=missing,sha256=job.sha(moment),
            source='public_data/K562_GWPS_bulk.h5ad; pooled mean UMI per guide weighted by num_cells_filtered',
            control='global core non-targeting guide pool; not batch matched; bulk mean reconstruction, not raw-cell integer counts'))
    assert wanted<=have,sorted(wanted-have)
    fallback=[];fallback_path=job.ROOT/'data/K562_GWPS_bulk_missing.npz'
    if fallback_path.exists():
        with np.load(fallback_path) as z:fallback=z['targets'].astype(str).tolist()
    job.save(job.ROOT/'artifacts/local_target_completeness.json',dict(status='PASS every frozen eligible bulk target represented in K562 family',
        targets=len(wanted),bulk_fallback_targets=fallback,bulk_fallback_added_this_attempt=missing,extra_native_label='ZNF735 retained separately from bulk ZNF735_NM',
        target_file_sha256=job.sha(job.ROOT/'targets_local.csv')))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--moments',action='store_true');ap.add_argument('--responses',action='store_true');a=ap.parse_args()
    if a.moments:
        from source_selection import select
        selected=select()
        for name in ['HCT116','HEK293T','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','CD4_DE_statistics']:
            shared=job.REPO/f'public_data/model_v4_signature_banks/{name}.npz'
            source='CD4' if name.startswith('CD4_') else name
            paths=([shared] if shared.exists() else [])+[folder/f'{name}.npz' for folder in selected[source]]
            assert len(paths)>1,f'No expansion shards: {name}'
            merge_bank(paths,job.ROOT/f'banks/{name}.npz')
    if a.responses:
        check=job.ROOT/'artifacts/response_union_check.json'
        if not check.exists() or job.json.loads(check.read_text()).get('helper_sha256')!=job.sha(job.ROOT/'response_union.py'):
            from check_response_union import main as check_union
            check_union()
        for name in ['HCT116','HEK293T','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']:convert(job.ROOT/f'banks/{name}.npz',job.ROOT/'data')
        # The mandatory mock fetch also supplies new matched raw GWPS rows.
        # Retain billmA's cached rows, then prefer these fresh matched rows over global-pool duplicates.
        from response_union import append,raw_priority
        fresh=job.ROOT/'artifacts/raw_gwps_upgrade';convert(job.ROOT/'banks/K562_GWPS_raw.npz',fresh)
        added=append(job.ROOT/'data','K562_GWPS_raw',fresh/'K562_GWPS_raw.npz')
        from response_union import raw_targets
        with np.load(job.REPO/'model_billm/data/K562_GWPS_raw.npz') as z:baseline_raw=set(z['targets'].astype(str))
        job.save(job.ROOT/'artifacts/raw_gwps_upgrade.json',dict(new_matched_rows=len(raw_targets(job.ROOT/'data')-baseline_raw),
            rows_added_this_attempt=added,baseline_sha256=job.sha(job.REPO/'model_billm/data/K562_GWPS_raw.npz')))
        consolidate_native('K562_GWPS');consolidate_native('H1')
        # Retain only one row per K562-family target; cached raw GWPS takes priority.
        raw=job.ROOT/'data/K562_GWPS_raw.npz'
        with np.load(raw) as z:rawtargets=z['targets']
        p=job.ROOT/'data/K562.npz'
        with np.load(p) as z:d={k:z[k] for k in z.files}
        keep=~np.isin(d['targets'],rawtargets)
        for k in ['targets','response','control','n_cells']:d[k]=d[k][keep]
        np.savez_compressed(p,**d)
        raw_priority(job.ROOT/'data')
        complete_local_targets()
        manifest={};families={};dups={}
        for p in (job.ROOT/'data').glob('*.npz'):
            with np.load(p) as z:
                family=str(z['family']);tt=z['targets'].astype(str);n=len(tt)
                assert (z['n_cells']>=30).all() and np.isfinite(z['response']).all() and z['response'].shape[1]==18533
                assert len(set(tt))==n
                if family!='CD4':
                    overlap=set(tt)&dups.setdefault(family,set());assert not overlap,(family,overlap)
                    dups[family].update(tt)
                families[family]=families.get(family,0)+n
                manifest[p.stem]=dict(rows=n,family=family,sha256=job.sha(p),bytes=p.stat().st_size)
                meta=p.with_suffix('.json')
                if meta.exists():previous=job.json.loads(meta.read_text())
                else:
                    original=job.REPO/f'model_billm/data/{"H1" if p.stem=="H1_cached" else p.stem}.json'
                    previous=job.json.loads(original.read_text()) if original.exists() else {}
                previous.update(manifest[p.stem],measured_genes=int(z['measured'].sum()))
                previous['note']='Final post-deduplication summary hash/row count; prior provenance retained'
                job.save(meta,previous)
        job.save(job.ROOT/'data/manifest.json',manifest);job.save(job.ROOT/'artifacts/source_family_counts.json',families)
        print('RESPONSE MANIFEST',families,flush=True)
if __name__=='__main__':main()
