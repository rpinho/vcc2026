"""Stream CD4 public pseudobulk and DE summaries; preserve donor/condition matching.
Pseudobulk-derived arithmetic CPM is explicitly NOT a mean-per-cell CPM observation.
"""
import runtime
import numpy as np,pandas as pd,h5py,anndata as ad,time
from remote_h5 import Remote
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor
rt=runtime.ROOT;out=rt/'artifacts/banks';out.mkdir(exist_ok=True)
targets=np.array(runtime.json.loads((rt/'artifacts/targets.json').read_text()));tset=set(targets)
obs=pd.read_parquet(rt/'artifacts/GWCD4i.pseudobulk_merged.h5ad.obs.parquet');var=pd.read_parquet(rt/'artifacts/GWCD4i.pseudobulk_merged.h5ad.var.parquet')
sel=obs.perturbed_gene_name.astype(str).isin(tset)&obs.keep_for_DE.astype(bool)
# NTC guide counts are capped equally within each matched donor/condition/run.
ctrl=obs[obs.guide_type.astype(str)=='non-targeting'].sort_index().groupby(['10xrun_id','donor_id','culture_condition'],observed=True).head(20)
sel|=obs.index.isin(ctrl.index);rows=np.flatnonzero(sel);chosen=obs.iloc[rows].copy();chosen.to_parquet(rt/'artifacts/CD4_selected.parquet')
r=Remote('https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.pseudobulk_merged.h5ad',cap=4_000_000_000,block=65536)
start=time.time();print('CD4 selected',len(rows),'target rows',int(sel.sum()-len(ctrl)),'NTC',len(ctrl),flush=True)
with h5py.File(r,'r') as h:
 g=h['X'];ptr=g['indptr'][:].astype(np.int64);lengths=ptr[rows+1]-ptr[rows];op=np.r_[0,np.cumsum(lengths)];arrays={}
 for kind in ['data','indices']:
  d=g[kind];assert d.compression is None and not d.shuffle;chunk=d.chunks[0];needed=set()
  for row in rows:
   a,b=ptr[row:row+2]
   if b>a:needed.update(range(int(a)//chunk,(int(b)-1)//chunk+1))
  locations={}
  def visit(info):
   number=int(info.chunk_offset[0])//chunk
   if number in needed:locations[number]=(info.byte_offset,info.size)
  d.id.chunk_iter(visit);assert needed<=locations.keys();plans=[]
  for k,row in enumerate(rows):
   a=int(ptr[row]);b=int(ptr[row+1]);dest=int(op[k])
   while a<b:
    number=a//chunk;take=min(b,(number+1)*chunk)-a;off=locations[number][0]+(a-number*chunk)*d.dtype.itemsize
    plans.append((off,take*d.dtype.itemsize,dest,take));a+=take;dest+=take
  groups=[]
  for plan in sorted(plans):
   off,size,*_=plan
   if groups and off<=groups[-1][1]+16384 and max(groups[-1][1],off+size)-groups[-1][0]<=4*1024**2:groups[-1][1]=max(groups[-1][1],off+size);groups[-1][2].append(plan)
   else:groups.append([off,off+size,[plan]])
  amount=sum(b-a for a,b,_ in groups);assert r.bytes+amount<r.cap
  arr=np.empty(int(op[-1]),dtype=d.dtype);print(kind,'dtype',d.dtype,'entries',len(arr),'groups',len(groups),'MB',amount//10**6,flush=True)
  # Read exact public byte ranges, with two network workers and explicit budget.
  import requests,threading
  lock=threading.Lock();done=[0]
  def fetch(group):
   a,b,members=group
   s=requests.Session();s.trust_env=False
   for attempt in range(5):
    with lock:
     assert r.bytes+b-a<r.cap;r.bytes+=b-a
    try:
     with s.get(r.url,headers={'Range':f'bytes={a}-{b-1}','If-Match':r.etag,'Accept-Encoding':'identity'},timeout=(10,90),stream=True) as resp:
      assert resp.status_code==206;buf=resp.raw.read(b-a+1);assert len(buf)==b-a
     for off,size,dest,count in members:arr[dest:dest+count]=np.frombuffer(buf,dtype=d.dtype,count=count,offset=off-a)
     break
    except Exception:
     if attempt==4:raise
     time.sleep(2)
   with lock:
    done[0]+=1
    if done[0]%100==0:print(kind,done[0],len(groups),'GB',round(r.bytes/1e9,3),'sec',round(time.time()-start),flush=True)
  with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(fetch,groups))
  arrays[kind]=arr
 x=sparse.csr_matrix((arrays['data'],arrays['indices'],op),shape=(len(rows),len(var)))
assert np.all(x.data>=0) and np.all(x.data==np.floor(x.data));actual=np.asarray(x.sum(1)).ravel();np.testing.assert_allclose(actual,chosen.total_counts.to_numpy(),rtol=1e-6)
genes=var.gene_name.to_numpy(str);assert len(set(genes))==len(genes)
# Save compact actual pseudobulk rows temporarily for bank assembly and exact spot checks.
sparse.save_npz(rt/'artifacts/CD4_pseudobulk_selected.npz',x)
conditions=sorted(chosen.culture_condition.astype(str).unique());blocks=[]
for cond in conditions:
 c=chosen[chosen.culture_condition.astype(str)==cond];keys=list(dict.fromkeys(zip(c['10xrun_id'].astype(str),c.donor_id.astype(str))))
 nt=len(targets);ng=len(genes);counts=np.zeros((nt,ng));mp=np.zeros_like(counts);cm=np.zeros_like(counts);tm=np.zeros_like(counts);nc=np.zeros(nt,np.int64);dv=np.zeros_like(counts)
 indices=np.flatnonzero(chosen.culture_condition.astype(str).to_numpy()==cond);ctrl_rows=indices[chosen.iloc[indices].guide_type.astype(str).to_numpy()=='non-targeting'];gc=np.asarray(x[ctrl_rows].sum(0)).ravel();gp=gc/gc.sum();gm=gp*1e6
 for ti,t in enumerate(targets):
  rr=indices[chosen.iloc[indices].perturbed_gene_name.astype(str).to_numpy()==t]
  if not len(rr):continue
  n=chosen.iloc[rr].n_cells.to_numpy(float);nc[ti]=int(n.sum());counts[ti]=np.asarray(x[rr].sum(0)).ravel();depth=actual[rr];tm[ti]=counts[ti]/max(counts[ti].sum(),1)*1e6
  for run,donor in keys:
   pr=rr[(chosen.iloc[rr]['10xrun_id'].astype(str).to_numpy()==run)&(chosen.iloc[rr].donor_id.astype(str).to_numpy()==donor)]
   if not len(pr):continue
   cr=ctrl_rows[(chosen.iloc[ctrl_rows]['10xrun_id'].astype(str).to_numpy()==run)&(chosen.iloc[ctrl_rows].donor_id.astype(str).to_numpy()==donor)]
   if not len(cr):raise ValueError('Missing matched CD4 controls')
   cp=np.asarray(x[cr].sum(0)).ravel();cp/=cp.sum();lw=actual[pr].sum()/depth.sum();cw=chosen.iloc[pr].n_cells.sum()/n.sum();mp[ti]+=lw*cp;cm[ti]+=cw*cp*1e6
   # Guide-level replicate variance of CPM profiles, not invented single-cell variance.
   pp=x[pr].multiply((1e6/actual[pr])[:,None]).toarray();cc=x[cr].multiply((1e6/actual[cr])[:,None]).toarray()
   if len(pr)>1:dv[ti]+=cw*cw*np.var(pp,axis=0,ddof=1)/len(pr)
   if len(cr)>1:dv[ti]+=cw*cw*np.var(cc,axis=0,ddof=1)/len(cr)
 np.savez_compressed(out/f'CD4_{cond}.npz',source=np.array('CD4_'+cond),targets=targets,genes=genes,n_cells=nc,target_count_sums=counts.astype('f4'),target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=cm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=gm.astype('f4'),measured_genes=np.ones(ng,bool),delta_cpm_variance=dv.astype('f4'),arithmetic_cpm_observed=np.array(False))
# Read the publisher DE estimates, SEs and QC without the 16.8-GB full download.
dobs=pd.read_parquet(rt/'artifacts/GWCD4i.DE_stats.h5ad.obs.parquet');dvar=pd.read_parquet(rt/'artifacts/GWCD4i.DE_stats.h5ad.var.parquet');dr=np.flatnonzero(dobs.target_contrast_gene_name.astype(str).isin(tset));de=r.__class__('https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad',cap=1_000_000_000,block=65536)
with h5py.File(de,'r') as h:
 chosen_de=dobs.iloc[dr];layers={}
 for k in ['log_fc','lfcSE','adj_p_value']:
  ds=h['layers'][k];assert ds.chunks is None
  layers[k]=np.stack([ds[int(row),:] for row in dr]);print('DE',k,'MB',de.bytes//10**6,flush=True)
 conditions=sorted(chosen_de.culture_condition.astype(str).unique());shape=(len(conditions),len(targets),len(dvar));effects=np.zeros(shape,'f4');ses=np.zeros(shape,'f4');q=np.ones(shape,'f4');avail=np.zeros(shape[:2],bool);qc=np.zeros_like(avail);n=np.zeros(shape[:2],np.int64)
 for ci,cond in enumerate(conditions):
  for ti,t in enumerate(targets):
   rr=np.flatnonzero((chosen_de.culture_condition.astype(str).to_numpy()==cond)&(chosen_de.target_contrast_gene_name.astype(str).to_numpy()==t))
   if not len(rr):continue
   m=chosen_de.iloc[rr];effects[ci,ti]=layers['log_fc'][rr].mean(0);ses[ci,ti]=layers['lfcSE'][rr].mean(0);q[ci,ti]=layers['adj_p_value'][rr].max(0);avail[ci,ti]=True;n[ci,ti]=int(m.n_cells_target.min());qc[ci,ti]=bool(((m.n_guides>=2)&~m.single_guide_estimate.astype(bool)&m.ontarget_significant.astype(bool)&~m.distal_offtarget_flag.astype(bool)&~m.low_target_gex.astype(bool)).all())
 np.savez_compressed(out/'CD4_DE_statistics.npz',targets=targets,genes=dvar.gene_name.to_numpy(str) if 'gene_name' in dvar else dvar.index.to_numpy(str),conditions=np.array(conditions),log2fc=effects,lfcSE=ses,adjusted_p=q,available=avail,quality_pass=qc,n_cells=n)
runtime.save(rt/'artifacts/CD4_transfer.json',dict(read_bytes=r.bytes+de.bytes,pseudobulk_rows=len(rows),controls=len(ctrl),seconds=time.time()-start,arithmetic_CPM='Not measured: guide pseudobulk proportions only; published DE LFC used on covered coordinates'))
print('DONE CD4',flush=True)
