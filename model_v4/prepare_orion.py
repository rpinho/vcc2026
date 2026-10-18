"""Stream Orion cells into compact target moments and batch weight tables.
No target-by-batch-by-gene dense allocation. All data fetched with anonymous public ranges.
"""
import runtime,sys,os,time,hashlib
os.environ['RAYON_NUM_THREADS']='1';os.environ['TOKIO_WORKER_THREADS']='1';sys.path.insert(0,'vendor/deps')
import lance,numpy as np,pandas as pd,anndata as ad
from public_store import PublicStore
from orion_offsets import offset_mapping
rt=runtime.ROOT;name=sys.argv[1];out=rt/f'artifacts/banks/{name}.npz';out.parent.mkdir(exist_ok=True)
if out.exists():raise SystemExit('Existing bank')
c=pd.read_parquet(rt/f'artifacts/{name}_cells.parquet');v=pd.read_parquet(rt/f'artifacts/{name}_genes.parquet');ntc=c[c.gene_target=='Non-Targeting'].sort_values(['sample','cell_integer_id']).groupby('sample',observed=True).head(250);c=pd.concat([c[c.gene_target!='Non-Targeting'],ntc]).sort_values('cell_integer_id').reset_index(drop=True)
h1=ad.read_h5ad(rt.parent/'public_data/H1.h5ad',backed='r');h1genes=h1.var_names.to_numpy(str);h1.file.close();genes=np.array(list(dict.fromkeys(pd.read_csv(rt.parent/'data/gene_names.csv').gene_name.tolist()+list(h1genes))));mapping=pd.Index(genes).get_indexer(v.gene_id.to_numpy(str));targets=np.array(runtime.json.loads((rt/'artifacts/targets.json').read_text()));nt=len(targets);ng=len(genes)
p=PublicStore();meta=lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/expression.lance',storage_options={'hf_token':''});e=lance.dataset(f's3://public/data/{name}/expression.lance',version=meta.version,storage_options=p.options(),block_size=4096);original,logical=offset_mapping(name,e,p);missing=logical[c.cell_integer_id.to_numpy(np.int64)//5000]<0;c=c.loc[~missing].reset_index(drop=True);c.to_parquet(rt/f'artifacts/{name}_selected.parquet');samples=sorted(c['sample'].astype(str).unique());bi=pd.Index(samples).get_indexer(c['sample'].astype(str));ti=pd.Index(targets).get_indexer(c.gene_target.astype(str));group=np.where(ti>=0,ti,nt+bi);n=np.bincount(group,minlength=nt+len(samples));shape=(len(n),ng);sums=np.zeros(shape);qp=np.zeros(shape);q2=np.zeros(shape);dep=np.zeros(len(n));bn=np.zeros((nt,len(samples)));bl=np.zeros_like(bn);valid=ti>=0;np.add.at(bn,(ti[valid],bi[valid]),1)
checkpoint=rt/f'artifacts/{name}_checkpoint.npz';resume=0;records=0;prior_bytes=0;start=time.time();digest=hashlib.sha256(c.cell_integer_id.to_numpy().tobytes()).hexdigest()
if checkpoint.exists():
 with np.load(checkpoint) as d:
  assert str(d['format'])=='compact_v2' and str(d['selection_sha256'])==digest;sums=d['sums'];qp=d['qp'];q2=d['q2'];dep=d['dep'];bl=d['bl'];resume=int(d['next_cell']);records=int(d['records']);prior_bytes=int(d['read_bytes'])
print(name,'compact arrays MB',round((sums.nbytes+qp.nbytes+q2.nbytes+bn.nbytes+bl.nbytes)/1e6),'selected',len(c),'omitted missing cells',int(missing.sum()),flush=True)
try:
 for lo in range(resume,len(c),128):
  runtime.guard();hi=min(lo+128,len(c));rr=c.iloc[lo:hi];length=rr.n_genes_by_counts.to_numpy(np.int64);gg=rr.cell_integer_id.to_numpy(np.int64)//5000;starts=rr.cell_start_index.to_numpy(np.int64)-original[gg]+logical[gg];ix=np.concatenate([np.arange(s,s+k,dtype=np.uint64) for s,k in zip(starts,length)])
  for attempt in range(5):
   try:tab=e.take(ix,columns=['cell_integer_id','gene_integer_id','value']);break
   except OSError:
    if attempt==4:raise
    time.sleep(15)
  local=np.repeat(np.arange(len(rr)),length);assert np.array_equal(tab['cell_integer_id'].to_numpy(),np.repeat(rr.cell_integer_id.to_numpy(),length));val=tab['value'].to_numpy().astype(float);assert np.all(val==np.floor(val));depth=rr.total_counts.to_numpy(float);np.testing.assert_array_equal(np.bincount(local,weights=val,minlength=len(rr)),depth);col=mapping[tab['gene_integer_id'].to_numpy()];keep=col>=0;val=val[keep];col=col[keep];local=local[keep];g=group[lo:hi][local];q=val/depth[local];np.add.at(sums,(g,col),val);np.add.at(qp,(g,col),q);np.add.at(q2,(g,col),q*q);mapped_depth=np.bincount(local,weights=val,minlength=len(rr));np.add.at(dep,group[lo:hi],mapped_depth);ok=ti[lo:hi]>=0;np.add.at(bl,(ti[lo:hi][ok],bi[lo:hi][ok]),mapped_depth[ok]);records+=len(val)
  if lo%1280==0 or hi==len(c):
   read_bytes=prior_bytes+p.bytes;print(name,hi,len(c),'GB',round(read_bytes/1e9,3),'seconds',round(time.time()-start),flush=True);runtime.save(rt/f'artifacts/{name}_active_transfer.json',dict(read_bytes=read_bytes,cells=hi,selected=len(c)))
  if hi%8192==0:
   # Bound the complete replacement before any bytes are written; compressed checkpoint stays small.
   runtime.guard(sums.nbytes+qp.nbytes+q2.nbytes+bl.nbytes+2**28);tmp=checkpoint.with_suffix('.tmp.npz');np.savez_compressed(tmp,format=np.array('compact_v2'),selection_sha256=np.array(digest),sums=sums,qp=qp,q2=q2,dep=dep,bl=bl,next_cell=hi,records=records,read_bytes=prior_bytes+p.bytes);tmp.replace(checkpoint)
finally:p.close()
cn=n[nt:];control=sums[nt:];control_q=qp[nt:];control_q2=q2[nt:];gp=control.sum(0)/control.sum();gm=control_q.sum(0)/cn.sum();gv=np.maximum(control_q2.sum(0)/cn.sum()-gm*gm,0)/max(cn.sum()-1,1);cp=np.divide(control,dep[nt:,None],out=np.tile(gp,(len(cn),1)),where=dep[nt:,None]>0);cm=np.divide(control_q,cn[:,None],out=np.tile(gm,(len(cn),1)),where=cn[:,None]>0);cv=np.divide(control_q2,cn[:,None],out=np.zeros_like(control_q2),where=cn[:,None]>0)-cm*cm;cv=np.maximum(cv,0)/np.maximum(cn[:,None]-1,1);cv[cn==0]=gv
tn=n[:nt].astype(np.int64);tc=sums[:nt];tm=np.divide(qp[:nt]*1e6,tn[:,None],out=np.zeros_like(tc),where=tn[:,None]>0);pv=np.maximum(np.divide(q2[:nt],tn[:,None],out=np.zeros_like(tc),where=tn[:,None]>0)-(tm/1e6)**2,0)/np.maximum(tn[:,None]-1,1)*1e12;cw=np.divide(bn,tn[:,None],out=np.zeros_like(bn),where=tn[:,None]>0);lw=np.divide(bl,dep[:nt,None],out=np.zeros_like(bl),where=dep[:nt,None]>0);mp=lw@cp;mm=cw@cm*1e6;vc=(cw*cw)@cv*1e12
np.savez_compressed(out,source=np.array(name),targets=targets,genes=genes,n_cells=tn,target_count_sums=tc.astype('f4'),target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=(gm*1e6).astype('f4'),measured_genes=np.isin(np.arange(ng),mapping[mapping>=0]),delta_cpm_variance=(pv+vc).astype('f4'),target_mean_cpm_variance=pv.astype('f4'),matched_control_mean_cpm_variance=vc.astype('f4'),arithmetic_cpm_observed=np.array(True))
np.savez_compressed(rt/f'artifacts/banks/{name}_blocks.npz',targets=targets,samples=np.array(samples),target_batch_n=bn,target_batch_library=bl,target_mean_cpm_variance=pv.astype('f4'))
runtime.save(rt/f'artifacts/{name}_control_correction.json',dict(read_bytes=0,control_selection='Released first 250 per sample, excluding unavailable expression blocks',controls=int(cn.sum()),samples=len(samples),orion_lance_version=meta.version,bank_sha256=runtime.sha(out)))
runtime.save(rt/f'artifacts/{name}_transfer.json',dict(read_bytes=prior_bytes+p.bytes,cells=len(c),targets=int((tn>=20).sum()),records=records,bank_bytes=out.stat().st_size,seconds=time.time()-start,variance='Conservative target CPM variance across batches plus matched control variance',omitted_cells_without_expression=int(missing.sum())))
checkpoint.unlink(missing_ok=True);(rt/f'artifacts/{name}_active_transfer.json').unlink(missing_ok=True);print('DONE',name,flush=True)
