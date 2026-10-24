"""A bounded v5b -> v9a -> stricter-call v9b job."""
import runtime
import numpy as np,h5py
from scipy import sparse
from de import Reference
from layout import emit,restrict
from audit_worker import audit
from call_budget import budget,NCONF,supported
R=runtime.ROOT
CACHE={}
def group(job):
 c,ci,i=job;runtime.guard()
 if c not in CACHE:
  z=np.load(R/f'artifacts/recovery_{c}.npz');CACHE[c]=(Reference(c),z['confidence'],z['positive_votes'],z['negative_votes'],z['genes'].astype(str),z['targets'].astype(str))
 ref,cf,positive,negative,genes,targets=CACHE[c]
 with h5py.File(R/'out/recovered_v5b.h5ad') as f:
  X=f['X'];lo=(ci*300+i)*400;ptr=X['indptr'][lo:lo+401];a,b=int(ptr[0]),int(ptr[-1])
  original=sparse.csr_matrix((X['data'][a:b],X['indices'][a:b],ptr-a),shape=(400,len(genes))).toarray()
 own=int(np.flatnonzero(genes==targets[i])[0]);k=budget(c,targets[i])
 a,called,infoa=emit(original,ref,cf[i],k,own)
 _,_,aq=ref.test(a);actual=np.flatnonzero((aq<.05)&(np.arange(len(genes))!=own))
 amean=(a/(a.sum(1,dtype=np.int64)[:,None]/1e6)).mean(0);alfc=np.log2((amean+1e-9)/(ref.mean+1e-9))
 b,bcalled,infob=restrict(original,a,ref,cf[i],actual,supported(alfc,positive[i],negative[i]),own)
 ra,da=audit(a,original,c,cf[i],own,called);rb,db=audit(b,original,c,cf[i],own,bcalled,baseline=da)
 bmask=(db['q']<.05)&(np.arange(len(genes))!=own)
 rb['actual_calls_without_two_agreeing_unopposed_sources']=int(np.sum(bmask&~supported(db['lfc'],positive[i],negative[i])))
 rb['actual_calls_outside_v9a']=int(np.sum(bmask&(da['q']>=.05)))
 rb['common_call_direction_changes_from_v9a']=int(np.sum(bmask&(da['q']<.05)&(np.sign(db['lfc'])!=np.sign(da['lfc']))))
 ra['MSE_sufficient_condition']=bool(ra['jk']>=ra['v5b_jk']);rb['MSE_sufficient_condition']=bool(rb['jk']>=rb['v5b_jk'])
 for arm,rec,values,info in [('v9a',ra,da,infoa),('v9b',rb,db,infob)]:
  rec.update(context=c,target=targets[i],index=i,n_conf_est=NCONF[c],requested_call_budget=k,**info)
  np.savez_compressed(R/f'artifacts/de/{arm}_{c}_{i}.npz',**values)
 return a,b,ra,rb
