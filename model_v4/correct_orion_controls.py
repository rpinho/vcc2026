"""First-250-per-sample NTCs and compact, conservative CPM sampling variances.
Legacy per-batch arrays are reduced in bounded chunks, never loaded in full.
"""
import runtime,sys,os,zipfile
os.environ['RAYON_NUM_THREADS']='1';os.environ['TOKIO_WORKER_THREADS']='1';sys.path.insert(0,'vendor/deps')
import lance,numpy as np,pandas as pd
from public_store import PublicStore
from orion_offsets import offset_mapping
rt=runtime.ROOT

def compact_blocks(name,bank):
 path=rt/f'artifacts/banks/{name}_blocks.npz'
 with np.load(path) as d:
  if 'target_batch_n' in d:return {k:d[k] for k in d.files}
  keys=d['keys'];n=d['n'];depth=d['depth_sum']
 targets=bank['targets'];samples=np.array(sorted(set(keys[:,0])));ti=pd.Index(targets).get_indexer(keys[:,1]);bi=pd.Index(samples).get_indexer(keys[:,0]);bn=np.zeros((len(targets),len(samples)));bl=np.zeros_like(bn);valid=ti>=0;np.add.at(bn,(ti[valid],bi[valid]),n[valid]);np.add.at(bl,(ti[valid],bi[valid]),depth[valid]);q2=np.zeros_like(bank['target_mean_cpm'],dtype=float)
 with zipfile.ZipFile(path) as z,z.open('sum_squared_proportion.npy') as f:
  version=np.lib.format.read_magic(f);shape,fortran,dtype=np.lib.format.read_array_header_1_0(f) if version==(1,0) else np.lib.format.read_array_header_2_0(f);assert not fortran and shape==(len(keys),len(bank['genes']))
  for lo in range(0,len(keys),64):
   runtime.guard();hi=min(lo+64,len(keys));size=(hi-lo)*shape[1]*dtype.itemsize;buf=f.read(size);assert len(buf)==size;array=np.frombuffer(buf,dtype=dtype).reshape(hi-lo,shape[1])
   for i,row in enumerate(array,start=lo):
    if ti[i]>=0:q2[ti[i]]+=row
 n0=bank['n_cells'];var=np.maximum(np.divide(q2,n0[:,None],out=np.zeros_like(q2),where=n0[:,None]>0)-(bank['target_mean_cpm'].astype(float)/1e6)**2,0)/np.maximum(n0[:,None]-1,1)*1e12
 compact=dict(targets=targets,samples=samples,target_batch_n=bn,target_batch_library=bl,target_mean_cpm_variance=var.astype('f4'));tmp=path.with_suffix('.tmp.npz');runtime.guard(128*1024**2);np.savez_compressed(tmp,**compact);tmp.replace(path);print('Compacted legacy blocks',name,flush=True);return compact

for name in (sys.argv[1:] or ['HCT116','HEK293T']):
 marker=rt/f'artifacts/{name}_control_correction.json'
 if marker.exists():continue
 path=rt/f'artifacts/banks/{name}.npz'
 with np.load(path) as d:bank={k:d[k] for k in d.files}
 compact=compact_blocks(name,bank);samples=compact['samples'];bn=compact['target_batch_n'];bl=compact['target_batch_library'];pv=compact['target_mean_cpm_variance'].astype(float)
 c=pd.read_parquet(rt/f'artifacts/{name}_cells.parquet');c=c[c.gene_target=='Non-Targeting'].sort_values(['sample','cell_integer_id']).groupby('sample',observed=True).head(250).sort_values('cell_integer_id').reset_index(drop=True);v=pd.read_parquet(rt/f'artifacts/{name}_genes.parquet');genes=bank['genes'];mapping=pd.Index(genes).get_indexer(v.gene_id.to_numpy(str));p=PublicStore();meta=lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/expression.lance',storage_options={'hf_token':''});e=lance.dataset(f's3://public/data/{name}/expression.lance',version=meta.version,storage_options=p.options(),block_size=4096);original,logical=offset_mapping(name,e,p);c=c.loc[logical[c.cell_integer_id.to_numpy(np.int64)//5000]>=0].reset_index(drop=True);codes=pd.Index(samples).get_indexer(c['sample'].astype(str));assert np.all(codes>=0);n=np.bincount(codes,minlength=len(samples));shape=(len(samples),len(genes));counts=np.zeros(shape);qp=np.zeros(shape);q2=np.zeros(shape)
 try:
  for lo in range(0,len(c),128):
   runtime.guard();rr=c.iloc[lo:lo+128];length=rr.n_genes_by_counts.to_numpy(np.int64);g=rr.cell_integer_id.to_numpy(np.int64)//5000;starts=rr.cell_start_index.to_numpy(np.int64)-original[g]+logical[g];ix=np.concatenate([np.arange(s,s+k,dtype=np.uint64) for s,k in zip(starts,length)]);tab=e.take(ix,columns=['cell_integer_id','gene_integer_id','value']);local=np.repeat(np.arange(len(rr)),length);assert np.array_equal(tab['cell_integer_id'].to_numpy(),np.repeat(rr.cell_integer_id.to_numpy(),length));val=tab['value'].to_numpy().astype(float);np.testing.assert_array_equal(np.bincount(local,weights=val,minlength=len(rr)),rr.total_counts.to_numpy());col=mapping[tab['gene_integer_id'].to_numpy()];keep=col>=0;val=val[keep];local=local[keep];col=col[keep];group=codes[lo:lo+128][local];q=val/rr.total_counts.to_numpy()[local];np.add.at(counts,(group,col),val);np.add.at(qp,(group,col),q);np.add.at(q2,(group,col),q*q)
   if lo%2560==0:print(name,'control refresh',lo,len(c),'MB',round(p.bytes/1e6),flush=True)
 finally:p.close()
 gp=counts.sum(0)/counts.sum();gm=qp.sum(0)/n.sum();gv=np.maximum(q2.sum(0)/n.sum()-gm*gm,0)/max(n.sum()-1,1);cp=np.divide(counts,counts.sum(1,keepdims=True),out=np.tile(gp,(len(n),1)),where=counts.sum(1,keepdims=True)>0);cm=np.divide(qp,n[:,None],out=np.tile(gm,(len(n),1)),where=n[:,None]>0);cv=np.maximum(np.divide(q2,n[:,None],out=np.zeros_like(q2),where=n[:,None]>0)-cm*cm,0)/np.maximum(n[:,None]-1,1);cv[n==0]=gv;tn=bank['n_cells'];cw=np.divide(bn,tn[:,None],out=np.zeros_like(bn),where=tn[:,None]>0);lw=np.divide(bl,bl.sum(1,keepdims=True),out=np.zeros_like(bl),where=bl.sum(1,keepdims=True)>0);mp=lw@cp;mm=cw@cm*1e6;vc=(cw*cw)@cv*1e12
 bank.update(matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=(gm*1e6).astype('f4'),delta_cpm_variance=(pv+vc).astype('f4'),target_mean_cpm_variance=pv.astype('f4'),matched_control_mean_cpm_variance=vc.astype('f4'));tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**bank);tmp.replace(path)
 runtime.save(marker,dict(read_bytes=p.bytes,controls=len(c),samples=len(samples),orion_lance_version=meta.version,control_selection='Released first 250 per sample, excluding unavailable expression blocks',target_variance='Conservative CPM variance pooled across source batches',bank_sha256=runtime.sha(path)))
 print('Corrected controls',name,flush=True)
