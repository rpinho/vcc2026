"""Anonymous bounded extraction from uncompressed dense GWPS HDF5.
200 cells per target + 50 NTC per batch. Exact cell CPM and matched controls.
"""
import runtime,time,threading,requests,hashlib
from remote_h5 import Remote
from transfer_budget import TransferBudget
import numpy as np,pandas as pd,h5py
from concurrent.futures import ThreadPoolExecutor
rt=runtime.ROOT;out=rt/'artifacts/cd4_banks/K562_GWPS_raw.npz';out.parent.mkdir(exist_ok=True)
if out.exists():raise SystemExit('Existing bank')
obs=pd.read_parquet(rt/'artifacts/GWPS_raw_obs.parquet');var=pd.read_parquet(rt/'artifacts/GWPS_raw_var.parquet');original=dict(np.load(rt.parent/'public_data/model_v4_signature_banks/K562_GWPS_raw.npz'));wanted=runtime.json.loads((rt/'artifacts/CD4_wide_targets.json').read_text());targets=np.array(sorted(set(wanted)-set(original['targets'][original['n_cells']>=200])));labels=obs.gene.astype(str).to_numpy();batches=obs.gem_group.astype(str).to_numpy();rng=np.random.default_rng(20261016);selected=[]
for t in targets:
 rows=np.flatnonzero(labels==t)
 if len(rows):selected.extend(rng.choice(rows,min(200,len(rows)),replace=False))
# Reproduce the ORIGINAL control draw exactly, without re-downloading old target cells.
oldrng=np.random.default_rng(20261016)
for t in runtime.json.loads((rt.parent/'model_v4/artifacts/targets.json').read_text()):
 rows=np.flatnonzero(labels==t)
 if len(rows):oldrng.choice(rows,min(200,len(rows)),replace=False)
for batch in sorted(set(batches)):
 rows=np.flatnonzero((labels=='non-targeting')&(batches==batch))
 if len(rows):selected.extend(oldrng.choice(rows,min(50,len(rows)),replace=False))
rows=np.array(sorted(selected),np.int64);chosen=obs.iloc[rows].copy();selection_path=rt/'artifacts/CD4_source_GWPS_selected.parquet'
if selection_path.exists():assert np.array_equal(pd.read_parquet(selection_path).source_row.to_numpy(),chosen.source_row.to_numpy()),'Resume changed selected source cells'
chosen.to_parquet(selection_path);genes0=var.gene_name.to_numpy(str);genes=np.array(list(dict.fromkeys(genes0)));gene_map=pd.Index(genes).get_indexer(genes0);raw_ng=len(genes0)
nt=len(targets);batch_names=sorted(set(batches[rows]));bix=pd.Index(batch_names).get_indexer(batches[rows]);tix=pd.Index(targets).get_indexer(labels[rows]);groups=np.where(tix>=0,tix,nt+bix);n=np.bincount(groups,minlength=nt+len(batch_names));shape=(len(n),len(genes));sums=np.zeros(shape);qp=np.zeros(shape);q2=np.zeros(shape);deps=np.zeros(len(n));batch_n=np.zeros((nt,len(batch_names)));batch_lib=np.zeros_like(batch_n)
for i in range(len(rows)):
 if tix[i]>=0:batch_n[tix[i],bix[i]]+=1
budget=TransferBudget(rt/'artifacts/CD4_source_GWPS_charges.jsonl',2_900_000_000)
# Precharge the metadata attempt, even if killed before a checkpoint. Remote's
# 20-MB cap plus four possible last-block retries fits below this 25-MB reserve.
budget.charge(25_000_000,'metadata attempt upper bound (20 MB cap plus retries)')
r=Remote('https://ndownloader.figshare.com/files/35775507?fresh='+str(time.time_ns()),cap=20_000_000,block=65536)
budget.bind('source',dict(file_id=35775507,etag=r.etag,size=r.size))
with h5py.File(r,'r') as h:
 ds=h['X'];assert ds.shape==(len(obs),len(var)) and ds.chunks is None and ds.compression is None;offset=ds.id.get_offset();dtype=ds.dtype;assert dtype==np.float32
# Figshare's anonymous redirects expire after ten seconds; refresh in memory.
link_lock=threading.Lock();link={'until':0.,'url':None,'etag':None};control_http=requests.Session();control_http.trust_env=False

def public_link():
 with link_lock:
  if time.monotonic()>=link['until']:
   budget.charge(1,'redirect refresh probe')
   with control_http.get('https://ndownloader.figshare.com/files/35775507?fresh='+str(time.time_ns()),headers={'Range':'bytes=0-0'},stream=True,timeout=60) as response:
    assert response.status_code==206;link['url']=response.url;link['etag']=response.headers.get('ETag');assert link['etag']==r.etag;assert len(response.raw.read(2))==1
   link['until']=time.monotonic()+3
  return link['url'],link['etag']
rowbytes=raw_ng*4;plans=[]
for i,row in enumerate(rows):
 a=int(offset+row*rowbytes);b=a+rowbytes
 if plans and a<=plans[-1][1]+rowbytes and b-plans[-1][0]<=2*1024**2:plans[-1][1]=b;plans[-1][2].append((i,a))
 else:plans.append([a,b,[(i,a)]])
assert r.bytes<25_000_000
amount=sum(b-a for a,b,_ in plans);assert amount+25_000_000<2_900_000_000;print('CD4 source K562 selected',len(rows),'ranges',len(plans),'planned GB',round(amount/1e9,3),flush=True)
budget.bind('selection_and_plan',dict(rows_sha256=hashlib.sha256(rows.tobytes()).hexdigest(),ranges_sha256=hashlib.sha256(np.array([(a,b) for a,b,_ in plans],dtype=np.int64).tobytes()).hexdigest(),original_bank_sha256=runtime.sha(rt.parent/'public_data/model_v4_signature_banks/K562_GWPS_raw.npz')))
restart_receipt=rt/'artifacts/CD4_source_GWPS_restart.json';restart_reserve=runtime.json.loads(restart_receipt.read_text())['uncheckpointed_read_reserve_bytes'] if restart_receipt.exists() else 0
if restart_reserve:
 prefix=np.r_[0,np.cumsum([b-a for a,b,_ in plans])];assert int(np.max(prefix[4096:]-prefix[:-4096]))+2*2*1024**2<=restart_reserve
threadlocal=threading.local()
assert not restart_reserve,'Archive prior-job reserve; live charges come only from durable journal'
start=time.time()
def fetch(plan):
 a,b,members=plan
 if not hasattr(threadlocal,'http'):
  threadlocal.http=requests.Session();threadlocal.http.trust_env=False
 s=threadlocal.http
 for attempt in range(5):
  budget.charge(b-a,'payload range attempt')
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
checkpoint=rt/'artifacts/CD4_source_GWPS_checkpoint.npz';resume=0
if checkpoint.exists():
 with np.load(checkpoint) as cp:
  sums=cp['sums'];qp=cp['qp'];q2=cp['q2'];deps=cp['deps'];batch_lib=cp['batch_lib'];resume=int(cp['next_range']);assert int(cp['read_bytes'])<=budget.total
remaining=sum(b-a for a,b,_ in plans[resume:]);assert budget.total+remaining<=budget.limit,'Insufficient budget to finish from checkpoint'
def consume(k,result):
  runtime.guard()
  for i,x in result:
   x=np.bincount(gene_map,weights=x,minlength=len(genes))
   assert np.isfinite(x).all() and np.all(x>=0) and np.all(x==np.floor(x));d=x.sum();assert d>0;group=groups[i];sums[group]+=x;q=x/d;qp[group]+=q;q2[group]+=q*q;deps[group]+=d
   if tix[i]>=0:batch_lib[tix[i],bix[i]]+=d
  if k%1000==0:
   print('CD4 source K562',k,len(plans),'GB',round(budget.total/1e9,3),'seconds',round(time.time()-start),flush=True);runtime.save(rt/'artifacts/CD4_source_GWPS_active_transfer.json',dict(read_bytes=budget.total,ranges=k,total_ranges=len(plans)))
  if (k+1)%4096==0:
   temporary=rt/'artifacts/CD4_source_GWPS_checkpoint.tmp.npz';np.savez(temporary,sums=sums,qp=qp,q2=q2,deps=deps,batch_lib=batch_lib,next_range=k+1,read_bytes=budget.total);temporary.replace(checkpoint)
with ThreadPoolExecutor(max_workers=16) as pool:
 for batch_start in range(resume,len(plans),64):
  for k,result in enumerate(pool.map(fetch,plans[batch_start:batch_start+64]),start=batch_start):consume(k,result)
cr=np.arange(nt,len(n));assert np.all(n[cr]>0);gp=sums[cr].sum(0)/deps[cr].sum();gm=qp[cr].sum(0)/n[cr].sum()*1e6
ng=len(genes);tc=sums[:nt];tn=n[:nt].astype(np.int64);tm=np.divide(qp[:nt]*1e6,tn[:,None],out=np.zeros_like(tc),where=tn[:,None]>0);mp=np.zeros_like(tc);mm=np.zeros_like(tc);dv=np.zeros_like(tc);pv=np.zeros_like(tc);cv=np.zeros_like(tc)
control_p=sums[cr]/deps[cr,None];control_mean=qp[cr]/n[cr,None];control_var=np.maximum(q2[cr]/n[cr,None]-control_mean**2,0)/np.maximum(n[cr,None]-1,1)
for ti,t in enumerate(targets):
 if tn[ti]==0:continue
 lw=batch_lib[ti]/batch_lib[ti].sum();cw=batch_n[ti]/tn[ti];mp[ti]=lw@control_p;mm[ti]=cw@control_mean*1e6
 vp=np.maximum(q2[ti]/tn[ti]-(qp[ti]/tn[ti])**2,0)/max(tn[ti]-1,1)
 # Target CPM variance pooled across batches is conservative under true batch heterogeneity.
 pv[ti]=vp*1e12;cv[ti]=((cw*cw)@control_var)*1e12;dv[ti]=pv[ti]+cv[ti]
assert np.array_equal(genes,original['genes']) and np.array_equal(gp.astype('f4'),original['global_control_probability']) and np.array_equal(gm.astype('f4'),original['global_control_mean_cpm']),'Original control reproduction failed'
new=dict(source=np.array('K562_GWPS_raw'),targets=targets,genes=genes,n_cells=tn,target_count_sums=tc.astype('f4'),target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=gm.astype('f4'),measured_genes=np.ones(ng,bool),delta_cpm_variance=dv.astype('f4'),target_mean_cpm_variance=pv.astype('f4'),matched_control_mean_cpm_variance=cv.astype('f4'),arithmetic_cpm_observed=np.array(True));alltargets=np.array(sorted(set(original['targets'])|set(targets)));oi=pd.Index(alltargets).get_indexer(original['targets']);ni=pd.Index(alltargets).get_indexer(targets);merged=dict(original)
for key,value in original.items():
 if value.ndim and value.shape[0]==len(original['targets']) and key!='targets':
  assert key in new;z=np.zeros((len(alltargets),)+value.shape[1:],value.dtype);z[oi]=value;z[ni]=new[key];merged[key]=z
merged['targets']=alltargets;temp=out.with_suffix('.partial.npz');np.savez_compressed(temp,**merged);temp.replace(out)
for name in ['H1','HCT116','HEK293T']:
 path=out.parent/(name+'.npz')
 if not path.exists():path.symlink_to('../../../public_data/model_v4_signature_banks/'+name+'.npz')
runtime.save(rt/'artifacts/CD4_source_GWPS_transfer.json',dict(read_bytes=budget.total,accounting='Durable pre-I/O charges, including failed attempts and 25 MB metadata reserve per invocation',charge_journal='artifacts/CD4_source_GWPS_charges.jsonl',source_etag=r.etag,source_size=r.size,selected_cells=len(rows),targets_new=len(targets),targets_total=len(alltargets),seconds=time.time()-start,source='Replogle public GWPS raw dense matrix; 200 cells per new target, original 50 NTC per batch reproduced',sha256=runtime.sha(out),original_controls_bitwise_equal=True));(rt/'artifacts/CD4_source_GWPS_active_transfer.json').unlink(missing_ok=True);print('DONE CD4 source K562 expansion',flush=True)

checkpoint.unlink(missing_ok=True)
