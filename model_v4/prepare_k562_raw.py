"""Anonymous bounded extraction from uncompressed dense GWPS HDF5.
200 cells per target + 50 NTC per batch. Exact cell CPM and matched controls.
"""
import runtime,time,threading,requests
from remote_h5 import Remote
import numpy as np,pandas as pd,h5py
from concurrent.futures import ThreadPoolExecutor
rt=runtime.ROOT;out=rt/'artifacts/banks/K562_GWPS_raw.npz'
if out.exists():raise SystemExit('Existing bank')
obs=pd.read_parquet(rt/'artifacts/K562_GWPS_raw.obs.parquet');var=pd.read_parquet(rt/'artifacts/K562_GWPS_raw.var.parquet');targets=np.array(runtime.json.loads((rt/'artifacts/targets.json').read_text()));labels=obs.gene.astype(str).to_numpy();batches=obs.gem_group.astype(str).to_numpy();rng=np.random.default_rng(20261016);selected=[]
for t in targets:
 rows=np.flatnonzero(labels==t)
 if len(rows):selected.extend(rng.choice(rows,min(200,len(rows)),replace=False))
for batch in sorted(set(batches)):
 rows=np.flatnonzero((labels=='non-targeting')&(batches==batch))
 if len(rows):selected.extend(rng.choice(rows,min(50,len(rows)),replace=False))
rows=np.array(sorted(selected),np.int64);chosen=obs.iloc[rows].copy();chosen.to_parquet(rt/'artifacts/K562_GWPS_selected.parquet');genes0=var.gene_name.to_numpy(str);genes=np.array(list(dict.fromkeys(genes0)));gene_map=pd.Index(genes).get_indexer(genes0);raw_ng=len(genes0)
nt=len(targets);batch_names=sorted(set(batches[rows]));bix=pd.Index(batch_names).get_indexer(batches[rows]);tix=pd.Index(targets).get_indexer(labels[rows]);groups=np.where(tix>=0,tix,nt+bix);n=np.bincount(groups,minlength=nt+len(batch_names));shape=(len(n),len(genes));sums=np.zeros(shape);qp=np.zeros(shape);q2=np.zeros(shape);deps=np.zeros(len(n));batch_n=np.zeros((nt,len(batch_names)));batch_lib=np.zeros_like(batch_n)
for i in range(len(rows)):
 if tix[i]>=0:batch_n[tix[i],bix[i]]+=1
r=Remote('https://ndownloader.figshare.com/files/35775507',cap=4_000_000_000,block=65536)
with h5py.File(r,'r') as h:
 ds=h['X'];assert ds.shape==(len(obs),len(var)) and ds.chunks is None and ds.compression is None;offset=ds.id.get_offset();dtype=ds.dtype;assert dtype==np.float32
# Figshare's anonymous redirects expire after ten seconds; refresh in memory.
link_lock=threading.Lock();link={'until':0.,'url':None,'etag':None};control_http=requests.Session();control_http.trust_env=False

def public_link():
 with link_lock:
  if time.monotonic()>=link['until']:
   with control_http.get(r.url+'?fresh='+str(time.time_ns()),headers={'Range':'bytes=0-0'},stream=True,timeout=60) as response:
    assert response.status_code==206;link['url']=response.url;link['etag']=response.headers.get('ETag');assert len(response.raw.read(2))==1
   link['until']=time.monotonic()+3
  return link['url'],link['etag']
rowbytes=raw_ng*4;plans=[]
for i,row in enumerate(rows):
 a=int(offset+row*rowbytes);b=a+rowbytes
 if plans and a<=plans[-1][1]+rowbytes and b-plans[-1][0]<=2*1024**2:plans[-1][1]=b;plans[-1][2].append((i,a))
 else:plans.append([a,b,[(i,a)]])
amount=sum(b-a for a,b,_ in plans);assert amount<3_800_000_000;print('K562 selected',len(rows),'ranges',len(plans),'planned GB',round(amount/1e9,3),flush=True)
threadlocal=threading.local()
lock=threading.Lock();downloaded=[r.bytes];start=time.time()
def fetch(plan):
 a,b,members=plan
 if not hasattr(threadlocal,'http'):
  threadlocal.http=requests.Session();threadlocal.http.trust_env=False
 s=threadlocal.http
 for attempt in range(5):
  with lock:
   assert downloaded[0]+b-a<4_000_000_000;downloaded[0]+=b-a
  try:
   url,etag=public_link()
   with s.get(url,headers={'Range':f'bytes={a}-{b-1}','If-Match':etag,'Accept-Encoding':'identity'},stream=True,timeout=(10,90)) as resp:
    assert resp.status_code==206, f'HTTP {resp.status_code}';buf=resp.raw.read(b-a+1);assert len(buf)==b-a
   return [(i,np.frombuffer(buf,dtype=dtype,count=raw_ng,offset=off-a).astype(float)) for i,off in members]
  except Exception as err:
   print('K562 range retry',attempt,type(err).__name__,str(err)[:80] if isinstance(err,AssertionError) else '',flush=True)
   if attempt==4:raise
   with link_lock:link['until']=0.
   time.sleep(1)
checkpoint=rt/'artifacts/K562_raw_checkpoint.npz';resume=0
if checkpoint.exists():
 with np.load(checkpoint) as cp:
  sums=cp['sums'];qp=cp['qp'];q2=cp['q2'];deps=cp['deps'];batch_lib=cp['batch_lib'];resume=int(cp['next_range']);downloaded[0]+=int(cp['read_bytes'])
with ThreadPoolExecutor(max_workers=16) as pool:
 for k,result in enumerate(pool.map(fetch,plans[resume:]),start=resume):
  runtime.guard()
  for i,x in result:
   x=np.bincount(gene_map,weights=x,minlength=len(genes))
   assert np.isfinite(x).all() and np.all(x>=0) and np.all(x==np.floor(x));d=x.sum();assert d>0;group=groups[i];sums[group]+=x;q=x/d;qp[group]+=q;q2[group]+=q*q;deps[group]+=d
   if tix[i]>=0:batch_lib[tix[i],bix[i]]+=d
  if k%1000==0:
   print('K562',k,len(plans),'GB',round(downloaded[0]/1e9,3),'seconds',round(time.time()-start),flush=True);runtime.save(rt/'artifacts/K562_active_transfer.json',dict(read_bytes=downloaded[0],ranges=k,total_ranges=len(plans)))
  if (k+1)%4096==0:
   temporary=rt/'artifacts/K562_raw_checkpoint.tmp.npz';np.savez(temporary,sums=sums,qp=qp,q2=q2,deps=deps,batch_lib=batch_lib,next_range=k+1,read_bytes=downloaded[0]);temporary.replace(checkpoint)
cr=np.arange(nt,len(n));assert np.all(n[cr]>0);gp=sums[cr].sum(0)/deps[cr].sum();gm=qp[cr].sum(0)/n[cr].sum()*1e6
ng=len(genes);tc=sums[:nt];tn=n[:nt].astype(np.int64);tm=np.divide(qp[:nt]*1e6,tn[:,None],out=np.zeros_like(tc),where=tn[:,None]>0);mp=np.zeros_like(tc);mm=np.zeros_like(tc);dv=np.zeros_like(tc);pv=np.zeros_like(tc);cv=np.zeros_like(tc)
control_p=sums[cr]/deps[cr,None];control_mean=qp[cr]/n[cr,None];control_var=np.maximum(q2[cr]/n[cr,None]-control_mean**2,0)/np.maximum(n[cr,None]-1,1)
for ti,t in enumerate(targets):
 if tn[ti]==0:continue
 lw=batch_lib[ti]/batch_lib[ti].sum();cw=batch_n[ti]/tn[ti];mp[ti]=lw@control_p;mm[ti]=cw@control_mean*1e6
 vp=np.maximum(q2[ti]/tn[ti]-(qp[ti]/tn[ti])**2,0)/max(tn[ti]-1,1)
 # Target CPM variance pooled across batches is conservative under true batch heterogeneity.
 pv[ti]=vp*1e12;cv[ti]=((cw*cw)@control_var)*1e12;dv[ti]=pv[ti]+cv[ti]
np.savez_compressed(out,source=np.array('K562_GWPS_raw'),targets=targets,genes=genes,n_cells=tn,target_count_sums=tc.astype('f4'),target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=gm.astype('f4'),measured_genes=np.ones(ng,bool),delta_cpm_variance=dv.astype('f4'),target_mean_cpm_variance=pv.astype('f4'),matched_control_mean_cpm_variance=cv.astype('f4'),arithmetic_cpm_observed=np.array(True))
runtime.save(rt/'artifacts/K562_raw_transfer.json',dict(read_bytes=downloaded[0],selected_cells=len(rows),targets=int((tn>=20).sum()),seconds=time.time()-start,source='Replogle public GWPS raw dense matrix; 200 cells per target, 50 NTC per batch'))
(rt/'artifacts/K562_active_transfer.json').unlink(missing_ok=True)
print('DONE K562 raw',flush=True)

checkpoint.unlink(missing_ok=True)
