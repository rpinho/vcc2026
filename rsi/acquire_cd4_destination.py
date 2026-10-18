"""Restore frozen CD4 rows with durable charges and resumable raw payloads."""
import runtime as rt
from transfer_budget import TransferBudget
from download_checkpoint import DownloadCheckpoint
import hashlib
from tree_prefetch import PrefetchRemote as Remote,read_uncompressed_vector
from pathlib import Path
import h5py,anndata as ad,numpy as np,pandas as pd,requests,time,threading
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor
name='CD4_matched_search';folder=Path('artifacts/tasks');path=folder/(name+'.h5ad');checkpoint_root=Path('artifacts/CD4_destination_payload_checkpoint')
if path.exists() and not (checkpoint_root/'state.json').exists():raise SystemExit('CD4 destination exists')
m=rt.json.loads((folder/'cd4_wide_blueprint.json').read_text());selection=folder/(name+'_selection.parquet');sel=pd.read_parquet(selection).sort_values('source_row').reset_index(drop=True);url='https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/D4_Stim8hr.assigned_guide.h5ad';fixture=rt.json.loads(Path('artifacts/tree_prefetch_fixture.json').read_text());assert fixture['status']=='PASS bitwise index addresses and data' and fixture['implementation_sha256']==rt.sha('tree_prefetch.py')
assert m['n_real']==96 and rt.sha(selection)==m['selection_sha256'],'Restore must preserve the frozen selection and quota'
budget=TransferBudget('artifacts/CD4_destination_charges.jsonl',2_600_000_000)
budget.bind('frozen_inputs',dict(source=url,blueprint=rt.sha(folder/'cd4_wide_blueprint.json'),selection=rt.sha(selection),census=rt.sha('artifacts/CD4_cell_census.parquet'),implementation=rt.sha(__file__)))
# PrefetchRemote charges each attempted metadata range inside its 250 MB cap.
# The extra 1 MB covers the initial one-byte size probe and conservative overhead.
budget.charge(251_000_000,'metadata invocation upper bound including retries and probe')
r=Remote(url,cap=250_000_000,block=65536);expected=rt.json.loads(Path('artifacts/CD4_cell_feasibility.json').read_text());assert r.etag==expected['etag'] and r.size==expected['source_size'];budget.bind('source_identity',dict(url=url,etag=r.etag,size=r.size));start=time.time()
receipt=Path('artifacts/CD4_destination_transfer.json')
def row_layout(ptr):
 rows=np.sort(sel.source_row.to_numpy(np.int64));length=ptr[rows+1]-ptr[rows];return rows,np.r_[0,np.cumsum(length)]
locations={};dtypes={};chunks={}
with h5py.File(r,'r') as h:
 var=ad.io.read_elem(h['var']);g=h['X'];ptr=read_uncompressed_vector(g['indptr'],r).astype(np.int64);print('CD4 row pointers, metadata MB',round(r.bytes/1e6,2),flush=True);r.enable_tree_prefetch();rows,outptr=row_layout(ptr)
 for kind in ['data','indices']:
  ds=g[kind];assert ds.compression is None and not ds.shuffle;chunk=ds.chunks[0];dtypes[kind]=ds.dtype;chunks[kind]=chunk;needed=set()
  for row in rows:
   a,b=map(int,ptr[row:row+2])
   if b>a:needed.update(range(a//chunk,(b-1)//chunk+1))
  loc={}
  def visit(info):
   number=int(info.chunk_offset[0])//chunk
   if number in needed:loc[number]=(info.byte_offset,info.size)
  ds.id.chunk_iter(visit);assert needed<=loc.keys();locations[kind]=loc;print('CD4 indexed',kind,'metadata MB',round(r.bytes/1e6,2),flush=True)
r.disable_tree_prefetch()
def plans_for(rows,outptr):
 plans=[]
 for kind in ['data','indices']:
  pieces=[];chunk=chunks[kind];dtype=dtypes[kind]
  for i,row in enumerate(rows):
   a,b=map(int,ptr[row:row+2]);dest=int(outptr[i])
   while a<b:
    num=a//chunk;take=min(b,(num+1)*chunk)-a;off=locations[kind][num][0]+(a-num*chunk)*dtype.itemsize;pieces.append((off,take*dtype.itemsize,dest,take,kind));a+=take;dest+=take
  for piece in sorted(pieces):
   off,n,*_=piece
   if plans and plans[-1][3]==kind and off<=plans[-1][1]+8192 and max(off+n,plans[-1][1])-plans[-1][0]<=2*1024**2:plans[-1][1]=max(off+n,plans[-1][1]);plans[-1][2].append(piece)
   else:plans.append([off,off+n,[piece],kind])
 return plans
plans=plans_for(rows,outptr);amount=sum(b-a for a,b,*_ in plans)
identity=dict(ranges_sha256=hashlib.sha256(np.asarray([(a,b) for a,b,*_ in plans],dtype=np.int64).tobytes()).hexdigest(),rows_sha256=hashlib.sha256(rows.tobytes()).hexdigest(),outptr_sha256=hashlib.sha256(outptr.tobytes()).hexdigest(),source=budget.bindings['source_identity'],inputs=budget.bindings['frozen_inputs'])
budget.bind('selection_and_plan',identity)
rt.guard(int(outptr[-1])*sum(v.itemsize for v in dtypes.values()))
checkpoint=DownloadCheckpoint(checkpoint_root,{k:((int(outptr[-1]),),v) for k,v in dtypes.items()},identity);arrays=checkpoint.arrays;resume=checkpoint.next_plan;assert resume<=len(plans)
remaining=sum(b-a for a,b,*_ in plans[resume:]);assert budget.total+remaining<=budget.limit,'Insufficient cumulative budget to finish from checkpoint'
local=threading.local();print('CD4 destination',len(rows),'rows',len(plans),'ranges','resume',resume,'planned GB',round(amount/1e9,3),'charged GB',round(budget.total/1e9,3),flush=True)
def fetch(plan):
 a,b,pieces,kind=plan
 if not hasattr(local,'http'):local.http=requests.Session();local.http.trust_env=False
 for attempt in range(3):
  rt.guard()
  budget.charge(b-a,'payload range attempt')
  try:
   with local.http.get(url,headers={'Range':f'bytes={a}-{b-1}','Accept-Encoding':'identity','If-Match':r.etag},stream=True,timeout=(10,60)) as response:
    assert response.status_code==206;buf=response.raw.read(b-a+1);assert len(buf)==b-a
   for off,n,dest,take,k in pieces:arrays[k][dest:dest+take]=np.frombuffer(buf,dtype=arrays[k].dtype,count=take,offset=off-a)
   return
  except Exception:
   if attempt==2:raise
try:
 with ThreadPoolExecutor(max_workers=3) as executor:
  for batch_start in range(resume,len(plans),64):
   batch_end=min(batch_start+64,len(plans))
   for _ in executor.map(fetch,plans[batch_start:batch_end]):pass
   checkpoint.commit(batch_end)
   if batch_start//64%8==0 or batch_end==len(plans):
    print('CD4 ranges',batch_end,len(plans),'GB',round(budget.total/1e9,3),'seconds',round(time.time()-start),flush=True);rt.save(receipt,dict(bytes=budget.total,ranges=batch_end,total_ranges=len(plans),status='in progress',charge_journal=str(budget.path)))
finally:rt.save(receipt,dict(bytes=budget.total,metadata_read_bytes=r.bytes,metadata_invocation_reserved_bytes=251_000_000,planned_bytes=amount,source=url,etag=r.etag,source_size=r.size,seconds=time.time()-start,status='ranges ended',charge_journal=str(budget.path),resumed_at_plan=resume));r.http.close()
x=sparse.csr_matrix((arrays['data'],arrays['indices'],outptr),shape=(len(rows),len(var)));del arrays
for lo in range(0,len(x.data),1_000_000):
 v=x.data[lo:lo+1_000_000];assert np.isfinite(v).all() and np.all(v>=0) and np.equal(v,np.floor(v)).all()
assert np.array_equal(rows,sel.source_row.to_numpy());x.data=x.data.astype(np.int32);gene_names=var.gene_name.astype(str).to_numpy() if 'gene_name' in var else var.index.to_numpy(str);official=pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:,0].astype(str).to_numpy();measured=set(gene_names);genes=np.array([s for s in official if s in measured])
if len(measured)==len(gene_names):x=x[:,pd.Index(gene_names).get_indexer(genes)]
else:
 mapped=pd.Index(genes).get_indexer(gene_names);valid=mapped>=0;project=sparse.csr_matrix((np.ones(valid.sum(),np.int32),(np.flatnonzero(valid),mapped[valid])),shape=(len(gene_names),len(genes)));x=x@project
x.sum_duplicates();x.sort_indices();totals=np.asarray(x.sum(1)).ravel();assert np.all((totals>0)&(totals<=1_000_000));a=ad.AnnData(x,obs=sel[['target_gene','batch','guide_id','evaluation_role','source_row']].copy(),var=pd.DataFrame(index=genes));a.uns['evaluation_only']=True;rt.guard(2*1024**3);temp=path.with_suffix('.partial.h5ad');a.write_h5ad(temp,compression='gzip',compression_opts=1);temp.replace(path);m.update(readouts=len(genes),readout_policy='measured official intersection');rt.save(path.with_suffix('.json'),dict(sha256=rt.sha(path),source=url,etag=r.etag,source_size=r.size,bytes=budget.total,selection_sha256=rt.sha(selection),readout_policy='measured official intersection',source_census_sha256=rt.sha('artifacts/CD4_cell_census.parquet')));manifest=rt.json.loads((folder/'manifest.json').read_text());manifest[name]=m;rt.save(folder/'manifest.json',manifest);print('DONE CD4 matched destination',a.shape,'n_real',m['n_real'],flush=True)
checkpoint.remove();budget.close()
