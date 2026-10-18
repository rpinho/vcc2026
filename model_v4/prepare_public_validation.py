"""Response-blind Orion whole-lineage holdouts with aggregate batch-matched controls.
Actual 400-cell groups and disjoint 1600-cell input/reference controls; unknown readouts excluded.
"""
import runtime,sys,os,hashlib
os.environ['RAYON_NUM_THREADS']='1';os.environ['TOKIO_WORKER_THREADS']='1';sys.path.insert(0,'vendor/deps')
import lance,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from public_store import PublicStore
from orion_offsets import offset_mapping
rt=runtime.ROOT;requested=set(pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.astype(str));official=pd.read_csv(rt.parent/'data/gene_names.csv').gene_name.to_numpy(str)
for name in ['HCT116','HEK293T']:
 path=rt/f'artifacts/validation_{name}.h5ad';marker=rt/f'artifacts/validation_{name}.json'
 if path.exists() and marker.exists():assert runtime.sha(path)==runtime.json.loads(marker.read_text())['sha256'];continue
 c=pd.read_parquet(rt/f'artifacts/{name}_cells.parquet');v=pd.read_parquet(rt/f'artifacts/{name}_genes.parquet');genes=official[np.isin(official,v.gene_id.astype(str))];mapping=pd.Index(genes).get_indexer(v.gene_id.to_numpy(str));p=PublicStore();meta=lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/expression.lance',storage_options={'hf_token':''});e=lance.dataset(f's3://public/data/{name}/expression.lance',version=meta.version,storage_options=p.options(),block_size=4096);original,logical=offset_mapping(name,e,p)
 c=c.loc[logical[c.cell_integer_id.to_numpy(np.int64)//5000]>=0].copy();samples=sorted(c['sample'].astype(str).unique());nt=c[c.gene_target=='Non-Targeting'];nc=nt.groupby('sample',observed=True).size().reindex(samples,fill_value=0).to_numpy();capacities=nc//8
 counts=c.groupby(['gene_target','sample'],observed=True).size().unstack(fill_value=0).reindex(columns=samples,fill_value=0);eligible=[t for t in requested&set(genes)&set(counts.index) if counts.loc[t].sum()>=400];eligible=sorted(eligible,key=lambda t:(-int(counts.loc[t].sum()),hashlib.sha256(('20261016:'+t).encode()).digest()))
 targets=sorted(eligible,key=lambda t:hashlib.sha256(('20261016:'+t).encode()).digest())[:12]
 if len(targets)<4:
  runtime.save(marker,dict(status='UNAVAILABLE',reason='Fewer than four actual 400-cell groups',targets=targets));p.close();print('Unavailable validation',name,targets,flush=True);continue
 rng=np.random.default_rng(20261016);chosen=[];roles=[];batch_total=np.zeros(len(samples),int)
 for target in targets:
  rows=c[c.gene_target==target];pick=rng.choice(len(rows),400,replace=False);rr=rows.iloc[pick];chosen.extend(rr.index);roles.extend(['perturbation']*400);batch_total+=rr.groupby('sample',observed=True).size().reindex(samples,fill_value=0).to_numpy()
 # Each input/reference arm matches the pooled perturbation batch mixture.
 # Individual targets retain their actual source batch mix; exact 400-cell shared quotas are infeasible.
 capacity=nc//2;ideal=1600*batch_total/batch_total.sum();quota=np.minimum(np.floor(ideal).astype(int),capacity)
 while quota.sum()<1600:
  allowed=np.flatnonzero(quota<capacity);assert len(allowed)>0;i=max(allowed,key=lambda i:(ideal[i]-quota[i],-i));quota[i]+=1
 for sample,q in zip(samples,quota):
  if not q:continue
  rows=nt[nt['sample'].astype(str)==sample];pick=rng.choice(len(rows),2*int(q),replace=False);chosen.extend(rows.iloc[pick].index);roles.extend(['input_control']*int(q)+['reference_control']*int(q))

 selected=c.loc[chosen].copy();selected['evaluation_role']=roles;selected=selected.sort_values('cell_integer_id').reset_index(drop=True);assert selected.cell_integer_id.is_unique;selected.to_parquet(rt/f'artifacts/validation_{name}_selection.parquet');blocks=[]
 print(name,'validation',len(targets),'targets',len(selected),'cells',len(genes),'genes',flush=True)
 try:
  for lo in range(0,len(selected),128):
   runtime.guard();rr=selected.iloc[lo:lo+128];length=rr.n_genes_by_counts.to_numpy(np.int64);g=rr.cell_integer_id.to_numpy(np.int64)//5000;starts=rr.cell_start_index.to_numpy(np.int64)-original[g]+logical[g];ix=np.concatenate([np.arange(s,s+k,dtype=np.uint64) for s,k in zip(starts,length)]);tab=e.take(ix,columns=['cell_integer_id','gene_integer_id','value']);local=np.repeat(np.arange(len(rr)),length);assert np.array_equal(tab['cell_integer_id'].to_numpy(),np.repeat(rr.cell_integer_id.to_numpy(),length));value=tab['value'].to_numpy().astype(np.int32);actual=np.bincount(local,weights=value,minlength=len(rr));np.testing.assert_array_equal(actual,rr.total_counts.to_numpy());col=mapping[tab['gene_integer_id'].to_numpy()];keep=col>=0;x=sparse.csr_matrix((value[keep],(local[keep],col[keep])),shape=(len(rr),len(genes)));x.sum_duplicates();x.sort_indices();assert np.all(np.asarray(x.sum(1)).ravel()>0);blocks.append(x)
   if lo%1280==0:print(name,'validation',lo,len(selected),'MB',round(p.bytes/1e6),flush=True)
 finally:p.close()
 obs=pd.DataFrame(dict(target_gene=selected.gene_target.replace({'Non-Targeting':'non-targeting'}).to_numpy(),evaluation_role=selected.evaluation_role.to_numpy(),batch=selected['sample'].astype(str).to_numpy(),source_cell_id=selected.cell_integer_id.to_numpy()),index=[f'{name}_{i}' for i in selected.cell_integer_id]);a=ad.AnnData(sparse.vstack(blocks,format='csr'),obs=obs,var=pd.DataFrame(index=genes));a.uns['evaluation_only']=True;a.write_h5ad(path,compression='gzip',compression_opts=1)
 runtime.save(marker,dict(sha256=runtime.sha(path),read_bytes=p.bytes,targets=targets,source=name,source_version=meta.version,source_family_excluded=name,batch_quotas=dict(zip(samples,map(int,quota))),genes=len(genes),official_unmeasured_readouts_excluded=int(len(official)-len(genes)),cells_per_target=400,input_controls=1600,reference_controls=1600,selection='Fixed hash order among targets with at least 400 cells; no expression/response scores used',batch_matching='Input and reference control quotas match the pooled perturbation sample mixture. Individual target mixtures differ; this is not exact per-target batch balance.'))
 print('DONE validation',name,flush=True)
