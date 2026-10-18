"""Anonymous sparse HDF5 range extraction, no full 15.5 GB download."""
import runtime as rt
import fcntl,os,subprocess
_owner_lock=open('artifacts/H1_rich_owner.lock','a')
try:fcntl.flock(_owner_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit('H1 rich acquisition already owns the lock')
_owner_path=rt.ROOT/'artifacts/H1_rich_owner.json'
if _owner_path.exists():
 _old=rt.json.loads(_owner_path.read_text())
 _check=subprocess.run(['ps','-p',str(_old['pid']),'-o','lstart='],capture_output=True,text=True)
 if _check.returncode==0 and _check.stdout.strip()==_old['started']:raise SystemExit('H1 rich acquisition already active; leave its owner running')
rt.save(_owner_path,dict(pid=os.getpid(),started=subprocess.check_output(['ps','-p',str(os.getpid()),'-o','lstart='],text=True).strip()))
from remote_h5 import Remote
from transfer_budget import TransferBudget
from download_checkpoint import DownloadCheckpoint
import h5py,anndata as ad,numpy as np,pandas as pd,requests,threading,time
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
manifest=rt.json.loads(Path('artifacts/tasks/h1_rich_blueprint.json').read_text());url='https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad?generation=1765904883947296'
selections={k:pd.read_parquet('artifacts/tasks/'+k+'_selection.parquet') for k,m in manifest.items() if m['context']=='H1' and m['source']=='H1_full_rich'}
if all(Path('artifacts/tasks/'+k+'.h5ad').exists() for k in selections) and not Path('artifacts/H1_rich_payload_checkpoint/state.json').exists():raise SystemExit('H1 exists')
rows=np.sort(np.unique(np.concatenate([s.source_row.to_numpy() for s in selections.values()])))
budget=TransferBudget('artifacts/H1_rich_charges.jsonl',5_500_000_000)
budget.bind('frozen_inputs',dict(source=url,blueprint=rt.sha('artifacts/tasks/h1_rich_blueprint.json'),selections={k:rt.sha('artifacts/tasks/'+k+'_selection.parquet') for k in selections},observations=rt.sha('artifacts/H1_full_obs.parquet')))
# Remote's 100-MB metadata cap can overshoot by four final 4-KiB retries,
# plus its one-byte fallback probe. Charge this upper bound before any I/O.
budget.charge(101_000_000,'metadata invocation upper bound including retries and probe')
r=Remote(url,cap=100_000_000,block=4096)
assert r.etag and r.size>0,'Missing public source identity'
budget.bind('source_identity',dict(url=url,etag=r.etag,size=r.size))
specs={};plans=[]
with h5py.File(r,'r') as h:
 var=ad.io.read_elem(h['var']);check_labels=ad.io.read_elem(h['obs/target_gene']);np.testing.assert_array_equal(np.asarray(check_labels).astype(str),pd.read_parquet('artifacts/H1_full_obs.parquet').target_gene.astype(str).to_numpy());g=h['X'];ptr=g['indptr'][:].astype(np.int64);length=ptr[rows+1]-ptr[rows];outptr=np.r_[0,np.cumsum(length)]
 for kind in ['data','indices']:
  ds=g[kind];assert ds.compression is None and not ds.shuffle;chunk=ds.chunks[0];loc={}
  def visit(info):loc[int(info.chunk_offset[0])//chunk]=(info.byte_offset,info.size)
  ds.id.chunk_iter(visit);pieces=[];specs[kind]=((int(outptr[-1]),),ds.dtype)
  for i,row in enumerate(rows):
   a,b=map(int,ptr[row:row+2]);dest=int(outptr[i])
   while a<b:
    num=a//chunk;take=min(b,(num+1)*chunk)-a;off=loc[num][0]+(a-num*chunk)*ds.dtype.itemsize;pieces.append((off,take*ds.dtype.itemsize,dest,take,kind));a+=take;dest+=take
  for piece in sorted(pieces):
   off,n,*_=piece
   if plans and plans[-1][3]==kind and off<=plans[-1][1]+65536 and max(off+n,plans[-1][1])-plans[-1][0]<=8*1024**2:plans[-1][1]=max(off+n,plans[-1][1]);plans[-1][2].append(piece)
   else:plans.append([off,off+n,[piece],kind])
assert r.bytes<101_000_000
import hashlib
plan_identity=dict(ranges_sha256=hashlib.sha256(np.asarray([(a,b) for a,b,*_ in plans],dtype=np.int64).tobytes()).hexdigest(),rows_sha256=hashlib.sha256(rows.tobytes()).hexdigest(),outptr_sha256=hashlib.sha256(outptr.tobytes()).hexdigest(),source=budget.bindings['source_identity'],inputs=budget.bindings['frozen_inputs'])
budget.bind('selection_and_plan',plan_identity)
checkpoint=DownloadCheckpoint('artifacts/H1_rich_payload_checkpoint',specs,plan_identity);arrays=checkpoint.arrays
resume=checkpoint.next_plan;assert resume<=len(plans)
amount=sum(b-a for a,b,*_ in plans);remaining=sum(b-a for a,b,*_ in plans[resume:]);print('H1 rows',len(rows),'plans',len(plans),'resume',resume,'planned bytes',amount,'charged bytes',budget.total,flush=True);assert budget.total+remaining<=budget.limit,'Insufficient cumulative budget to finish from checkpoint'
local=threading.local();start=time.time()
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
with ThreadPoolExecutor(max_workers=3) as executor:
 for batch_start in range(resume,len(plans),64):
  batch_end=min(batch_start+64,len(plans))
  for _ in executor.map(fetch,plans[batch_start:batch_end]):pass
  checkpoint.commit(batch_end)
  if batch_start//64%8==0 or batch_end==len(plans):
   print(batch_end,len(plans),round(budget.total/1e9,3),round(time.time()-start),flush=True);rt.save('artifacts/H1_rich_acquisition_progress.json',dict(bytes=budget.total,plans_done=batch_end,plans=len(plans),metadata_read_bytes=r.bytes,charge_journal='artifacts/H1_rich_charges.jsonl'))
x=sparse.csr_matrix((arrays['data'],arrays['indices'],outptr),shape=(len(rows),len(var)));assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
for name,sel in selections.items():
 ix=np.searchsorted(rows,sel.source_row.to_numpy());assert np.array_equal(rows[ix],sel.source_row.to_numpy());official=pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:,0].astype(str).to_numpy();cols=pd.Index(var.index).get_indexer(official);genes=official[cols>=0];cols=cols[cols>=0];a=ad.AnnData(x[ix][:,cols].astype(np.int32),obs=sel[['target_gene','batch','guide_id','evaluation_role','source_row']].copy(),var=var.iloc[cols].copy());a.var_names=genes;a.uns['evaluation_only']=True;path=Path('artifacts/tasks')/(name+'.h5ad');temporary=path.with_suffix('.partial.h5ad');a.write_h5ad(temporary,compression='gzip',compression_opts=1);temporary.replace(path);rt.save(path.with_suffix('.json'),dict(sha256=rt.sha(path),source=url,selected_rows_sha256=rt.sha(Path('artifacts/tasks')/(name+'_selection.parquet')),bytes_total=budget.total,readout_policy='measured official intersection',source_generation='1765904883947296',source_etag=r.etag,source_size=r.size));print('Saved',name,a.shape,flush=True)
rt.save('artifacts/H1_rich_acquisition.json',dict(bytes=budget.total,planned_bytes=amount,seconds=time.time()-start,source=url,source_etag=r.etag,source_size=r.size,metadata_read_bytes=r.bytes,accounting='Durable pre-I/O charges, including retries and 101 MB metadata reserve per invocation',charge_journal='artifacts/H1_rich_charges.jsonl',resumed_at_plan=resume));r.http.close()

current=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());rt.save('artifacts/tasks/manifest_before_h1_rich.json',current);current.update(manifest);rt.save('artifacts/tasks/manifest.json',current)
checkpoint.remove();budget.close()
