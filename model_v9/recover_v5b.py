"""Reproduce the previously submitted baseline; require archived exact CSR hashes."""
import runtime,sys,gc,types,time
R=runtime.ROOT
runtime.shutil=__import__('shutil')
sys.path.insert(0,str(R.parent/'model_v5'))
import build_v5 as old
import model,v5
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from recovery_worker import initialize,emit
original_control_axes=v5.control_axes
def recorded_control_axes(x,k=50,seed=0):
 V,s=original_control_axes(x,k=k,seed=seed)
 np.savez_compressed(R/f'artifacts/control_axes_{seed}.npz',V=V,s=s)
 return V,s
v5.control_axes=recorded_control_axes
def bounded_template(a,genes,cells=400,pool_k=4,seed=20260910):
 # Original template's multiply(...).mean(0) expands/sorts all 110M entries.
 # This arithmetic mean only feeds the discarded CPM head in v5b; the bulk
 # head, selected real cells, depths, pooling and integer emitter are unchanged.
 assert np.array_equal(a.var_names,genes)
 X=sparse.csr_matrix(a.X);depth=np.asarray(X.sum(1,dtype=np.float64)).ravel()
 p=np.asarray(X.T@(1/depth)).ravel()/X.shape[0]
 b=np.asarray(X.sum(0,dtype=np.float64)).ravel();b/=b.sum()
 selected=np.random.default_rng(seed).choice(len(depth),cells*pool_k,replace=False);selected=selected[np.argsort(depth[selected],kind='stable')]
 q=X[selected].toarray().astype(np.float64)/depth[selected,None];q=q.reshape(cells,pool_k,len(genes)).mean(1)
 d=np.rint(depth[selected].reshape(cells,pool_k).mean(1)).astype(np.int64)
 return q,d,p,b
old.template=bounded_template
genes=pd.read_csv(R.parent/'data/gene_names.csv').gene_name.to_numpy(str)
targets=pd.read_csv(R.parent/'data/pert_counts.csv').target_gene.to_numpy(str)
old.ARGS=types.SimpleNamespace(k=18533)
saved=runtime.json.loads((R.parent/'model_v5/artifacts/v5b_build.json').read_text())
assert runtime.sha(R.parent/'model_v5/model.py')==saved['model_sha256']
assert runtime.sha(R.parent/'model_v5/v5.py')==saved['v5_sha256']
path=R/'out/recovered_v5b.h5ad';writer=old.Writer(path,genes)
tot=[];conf=[];start=time.time()
try:
 for ci,c in enumerate('ABC'):
  runtime.guard();print('recovery loading',c,flush=True);a=ad.read_h5ad(R.parent/f'data/context_{c}.h5ad');print('recovery profile',c,flush=True)
  q,d,L,Lc,mask,info=old.profile('b',targets,genes,a,c,ci,print)
  bf=v5.from_logbulk(v5.rescale(L,saved['contexts'][ci]['alpha'],Lc))
  names=model.load_banks(None,True);src=[model.aligned_bank(n,targets,genes,'log2fc') for n in names];cd=model.cd4_family(targets,genes,'log2fc')
  if cd is not None:src.append(cd)
  sign=np.sign(L-Lc);net=np.zeros(L.shape);zs=np.zeros(L.shape);agree=np.zeros(L.shape,np.int8);oppose=np.zeros(L.shape,np.int8);positive=np.zeros(L.shape,np.int8);negative=np.zeros(L.shape,np.int8)
  for e,m,v in src:
   z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=np.abs(z)>1.96
   aa=sig*(np.sign(z)==sign);oo=sig*(np.sign(z)==-sign)
   positive+=(sig&(z>0)).astype(np.int8);negative+=(sig&(z<0)).astype(np.int8)
   agree+=aa.astype(np.int8);oppose+=oo.astype(np.int8)
   net+=aa.astype(float)-oo.astype(float);zs+=np.minimum(np.abs(z),100)*(np.sign(z)==sign)
  conf.append(net*1e6+zs);del a,src,cd;gc.collect()
  np.savez_compressed(R/f'artifacts/recovery_{c}.npz',q=q,d=d,bf=bf,confidence=conf[-1],agree=agree,oppose=oppose,positive_votes=positive,negative_votes=negative,genes=genes,targets=targets)
  totals=[]
  with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('fork'),initializer=initialize,initargs=(q,d)) as pool:
   pending={};nxt=0
   for i,t in enumerate(targets):
    while nxt<min(i+6,len(targets)):
     pending[nxt]=pool.submit(emit,(bf[nxt],model.atlas.seed_for(c+':'+targets[nxt],20260910)));nxt+=1
    x=pending.pop(i).result();totals.append(np.asarray(x.sum(0,dtype=np.int64)).ravel());writer.append(x,c,t)
    if i%30==0:print('recovery',c,i,'seconds',round(time.time()-start),flush=True)
  tot.append(totals)
finally:writer.close()
cert=old.audit_h5(path,genes,targets)
before=runtime.json.loads((R.parent/'model_v5/artifacts/v5b_structure.json').read_text())
assert cert==before,[(k,before[k],cert[k]) for k in before if before[k]!=cert[k]]
np.savez_compressed(R/'artifacts/stage1_v5b.npz',totals=np.array(tot),confidence=np.array(conf),genes=genes,targets=targets)
runtime.save(R/'artifacts/recovered_v5b.json',dict(status='PASS: exact archived canonical CSR and axes reproduced',structure=cert,source='model_v5/artifacts/v5b_structure.json',seconds=time.time()-start))
print('PASS exact baseline recovery',flush=True)
