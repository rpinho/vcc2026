"""Matched public-context cells, with all responses held out of source confidence."""
import runtime,sys,gc,types
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from de import prepare,Reference
from layout import emit,restrict
from scorer_check import config
from call_budget import budget,supported
from cell_eval2.run import _compute_de_side
R=runtime.ROOT
sys.path.insert(0,str(R.parent/'model_v5'))
import model,v5,build_v5 as old
def template(a,genes,cells=400,pool_k=4,seed=20260910):
 X=sparse.csr_matrix(a.X);d=np.asarray(X.sum(1,dtype=float)).ravel();p=np.asarray(X.T@(1/d)).ravel()/len(d);b=np.asarray(X.sum(0,dtype=float)).ravel();b/=b.sum()
 rows=np.random.default_rng(seed).choice(len(d),cells*pool_k,replace=False);rows=rows[np.argsort(d[rows],kind='stable')]
 q=X[rows].toarray().astype(float)/d[rows,None];q=q.reshape(cells,pool_k,len(genes)).mean(1)
 depths=np.rint(d[rows].reshape(cells,pool_k).mean(1)).astype(np.int64)
 return q,depths,p,b
def main(name):
 runtime.guard();ci={'Jurkat':0,'HepG2':1,'RPE1':2}[name];c='ABC'[ci]
 a=ad.read_h5ad(R/f'out/public_{name}.h5ad');control=a[a.obs.target.astype(str)=='non-targeting'].copy();genes=a.var_names.to_numpy(str)
 targets=np.array(list(runtime.json.loads((R/f'artifacts/public_{name}.json').read_text())['targets']),str)
 assert len(targets)>0 and pd.Index(genes).is_unique
 refdir=R/f'artifacts/ref_public_{name}'
 if not (refdir/'cdf_scale.npy').exists():prepare(sparse.csr_matrix(control.X),refdir,np.random.default_rng(20261021))
 ref=Reference('public_'+name);old.template=template;old.ARGS=types.SimpleNamespace(k=len(genes))
 q,d,L,Lc,mask,info=old.profile('b',targets,genes,control,c,ci,print)
 alpha=runtime.json.loads((R.parent/'model_v5/artifacts/v5b_build.json').read_text())['contexts'][ci]['alpha'];bf=v5.from_logbulk(v5.rescale(L,alpha,Lc))
 src=[model.aligned_bank(n,targets,genes,'log2fc') for n in model.load_banks(None,True)];cd=model.cd4_family(targets,genes,'log2fc')
 if cd is not None:src.append(cd)
 sign=np.sign(L-Lc);ag=np.zeros(L.shape,np.int8);op=np.zeros_like(ag);zs=np.zeros(L.shape);positive=np.zeros_like(ag);negative=np.zeros_like(ag)
 for e,m,v in src:
  z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=abs(z)>1.96;ag+=(sig&(np.sign(z)==sign)).astype(np.int8);op+=(sig&(np.sign(z)==-sign)).astype(np.int8);zs+=np.minimum(abs(z),100)*(np.sign(z)==sign);positive+=(sig&(z>0)).astype(np.int8);negative+=(sig&(z<0)).astype(np.int8)
 cf=(ag-op).astype(float)*1e6+zs;del src,cd;gc.collect()
 matrices={arm:[] for arm in ['v5b','v9a','v9b']};records=[]
 for i,t in enumerate(targets):
  runtime.guard();original=model.emit(q,None,bf[i],d,model.atlas.seed_for(c+':'+t,20260910),'anchor_single_mean').astype(np.int32);own=int(np.flatnonzero(genes==t)[0])
  xa,called,ia=emit(original,ref,cf[i],budget(c),own)
  _,_,aq=ref.test(xa);actual=np.flatnonzero((aq<.05)&(np.arange(len(genes))!=own))
  amean=(xa/(xa.sum(1,dtype=np.int64)[:,None]/1e6)).mean(0);alfc=np.log2((amean+1e-9)/(ref.mean+1e-9))
  xb,bcalled,ib=restrict(original,xa,ref,cf[i],actual,supported(alfc,positive[i],negative[i]),own)
  for arm,x in [('v5b',original),('v9a',xa),('v9b',xb)]:matrices[arm].append(sparse.csr_matrix(x))
  records.append(dict(target=t,selected_a=len(called),selected_b=len(bcalled)));print('PUBLIC BUILD',name,i,t,len(called),len(bcalled),flush=True)
 del q,bf,L,Lc,ref;gc.collect()
 for arm,parts in matrices.items():
  runtime.guard();X=sparse.vstack([sparse.csr_matrix(control.X),*parts],format='csr');labels=['non-targeting']*control.n_obs+np.repeat(targets,400).tolist()
  mat=ad.AnnData(X,obs=pd.DataFrame({'target':labels},index=[str(i) for i in range(len(labels))]),var=pd.DataFrame(index=genes))
  frame=_compute_de_side(mat,cfg=config(),fp=None,store=None,side='pred');frame.write_parquet(R/f'artifacts/public_pred_{arm}_{name}.parquet')
  print('PUBLIC OFFICIAL PRED DE',name,arm,frame.height,flush=True);del mat,X;gc.collect()
 runtime.save(R/f'artifacts/public_build_{name}.json',dict(status='complete matched public-context cells and official prediction-DE producer',targets=targets.tolist(),held_source=name,alpha_from_v5b=alpha,records=records,qualification='Public native-depth and partial-gene proxy, 16 targets selected using source coverage/counts only; no official targets or hidden labels.'))
if __name__=='__main__':main(sys.argv[1])
