"""Anonymous sparse HDF5 range extraction, no full 15.5 GB download."""
import runtime as rt
from remote_h5 import Remote
import h5py,anndata as ad,numpy as np,pandas as pd,requests,threading,time
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());url='https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad'
selections={k:pd.read_parquet('artifacts/tasks/'+k+'_selection.parquet') for k,m in manifest.items() if m['context']=='H1' and m['source']=='H1_full'}
if all(Path('artifacts/tasks/'+k+'.h5ad').exists() for k in selections):raise SystemExit('H1 exists')
rows=np.sort(np.unique(np.concatenate([s.source_row.to_numpy() for s in selections.values()])));r=Remote(url,cap=100_000_000,block=4096)
arrays={};plans=[]
with h5py.File(r,'r') as h:
 var=ad.io.read_elem(h['var']);g=h['X'];ptr=g['indptr'][:].astype(np.int64);length=ptr[rows+1]-ptr[rows];outptr=np.r_[0,np.cumsum(length)]
 for kind in ['data','indices']:
  ds=g[kind];assert ds.compression is None and not ds.shuffle;chunk=ds.chunks[0];loc={}
  def visit(info):loc[int(info.chunk_offset[0])//chunk]=(info.byte_offset,info.size)
  ds.id.chunk_iter(visit);pieces=[];arrays[kind]=np.empty(outptr[-1],ds.dtype)
  for i,row in enumerate(rows):
   a,b=map(int,ptr[row:row+2]);dest=int(outptr[i])
   while a<b:
    num=a//chunk;take=min(b,(num+1)*chunk)-a;off=loc[num][0]+(a-num*chunk)*ds.dtype.itemsize;pieces.append((off,take*ds.dtype.itemsize,dest,take,kind));a+=take;dest+=take
  for piece in sorted(pieces):
   off,n,*_=piece
   if plans and plans[-1][3]==kind and off<=plans[-1][1]+65536 and max(off+n,plans[-1][1])-plans[-1][0]<=8*1024**2:plans[-1][1]=max(off+n,plans[-1][1]);plans[-1][2].append(piece)
   else:plans.append([off,off+n,[piece],kind])
amount=sum(b-a for a,b,*_ in plans);print('H1 rows',len(rows),'plans',len(plans),'planned bytes',amount,flush=True);assert amount+r.bytes<6_000_000_000
lock=threading.Lock();downloaded=[r.bytes];local=threading.local();start=time.time()
def fetch(plan):
 a,b,pieces,kind=plan
 if not hasattr(local,'http'):local.http=requests.Session();local.http.trust_env=False
 for attempt in range(3):
  rt.guard()
  with lock:
   assert downloaded[0]+b-a<6_000_000_000;downloaded[0]+=b-a
  try:
   with local.http.get(url,headers={'Range':f'bytes={a}-{b-1}','Accept-Encoding':'identity'},stream=True,timeout=(10,60)) as response:
    assert response.status_code==206;buf=response.raw.read(b-a+1);assert len(buf)==b-a
   for off,n,dest,take,k in pieces:arrays[k][dest:dest+take]=np.frombuffer(buf,dtype=arrays[k].dtype,count=take,offset=off-a)
   return
  except Exception:
   if attempt==2:raise
with ThreadPoolExecutor(max_workers=3) as executor:
 for i,_ in enumerate(executor.map(fetch,plans)):
  if i%50==0:
   print(i,len(plans),round(downloaded[0]/1e9,3),round(time.time()-start),flush=True);rt.save('artifacts/H1_acquisition_progress.json',dict(bytes=downloaded[0],plans_done=i,plans=len(plans)))
x=sparse.csr_matrix((arrays['data'],arrays['indices'],outptr),shape=(len(rows),len(var)));assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
for name,sel in selections.items():
 ix=np.searchsorted(rows,sel.source_row.to_numpy());assert np.array_equal(rows[ix],sel.source_row.to_numpy());a=ad.AnnData(x[ix].astype(np.int32),obs=sel[['target_gene','batch','evaluation_role','source_row']].copy(),var=var.copy());a.uns['evaluation_only']=True;path=Path('artifacts/tasks')/(name+'.h5ad');a.write_h5ad(path,compression='gzip',compression_opts=1);rt.save(path.with_suffix('.json'),dict(sha256=rt.sha(path),source=url,selected_rows_sha256=rt.sha(Path('artifacts/tasks')/(name+'_selection.parquet')),bytes_total=downloaded[0]));print('Saved',name,a.shape,flush=True)
rt.save('artifacts/H1_acquisition.json',dict(bytes=downloaded[0],planned_bytes=amount,seconds=time.time()-start,source=url));r.http.close()
