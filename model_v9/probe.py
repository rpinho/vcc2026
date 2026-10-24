"""Bounded deterministic full-gene/full-control feasibility probes."""
import runtime,sys,time,gc,zipfile,os
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr
from cell_eval2.moments import jackknife_correction
from de import Reference
from layout import emit,restrict
from recovery_worker import initialize,emit as baseline_emit
from call_budget import budget,supported
R=runtime.ROOT
def wait_inputs(c):
 while True:
  try:
   if not (R/f'artifacts/ref_{c}/receipt.json').exists():time.sleep(5);continue
   z=np.load(R/f'artifacts/recovery_{c}.npz');q=z['q'];bf=z['bf']
   ref=Reference(c)
   return z,q,bf,ref
  except (OSError,ValueError,EOFError,zipfile.BadZipFile):time.sleep(5)
def main():
 sys.path.append(str(R.parent/'model_v5'));import model
 # Recovery uses three group workers. Prepare finishes its one-worker pass
 # before this fourth numerical worker starts.
 while not all((R/f'artifacts/ref_{c}/receipt.json').exists() for c in 'ABC'):time.sleep(5)
 records=[]
 for ci,c in enumerate('ABC'):
  z,q,bf,ref=wait_inputs(c);initialize(q,z['d']);targets=z['targets'].astype(str);genes=z['genes'].astype(str)
  for i in [0,4]:
   runtime.guard();original=baseline_emit((bf[i],model.atlas.seed_for(c+':'+targets[i],20260910))).toarray();own=int(np.flatnonzero(genes==targets[i])[0])
   jk0=float(jackknife_correction(sparse.csr_matrix(original),np.zeros(400,np.intp),1,50000.,chunk=512)[0]);d=original.sum(1);mean0=(original/(d[:,None]/1e6)).mean(0)
   for arm in ['a','b']:
    start=time.time()
    if arm=='a':
     x,called,info=emit(original,ref,z['confidence'][i],budget(c,targets[i]),own,block=int(os.environ.get('V9_BLOCK','10')))
     ax=x.copy();acalled=called.copy()
    else:
     amean=(ax/(ax.sum(1,dtype=np.int64)[:,None]/1e6)).mean(0)
     strict=supported(np.log2((amean+1e-9)/(ref.mean+1e-9)),z['positive_votes'][i],z['negative_votes'][i])
     _,_,aq=ref.test(ax);actual=np.flatnonzero((aq<.05)&(np.arange(len(genes))!=own))
     x,called,info=restrict(original,ax,ref,z['confidence'][i],actual,strict,own,block=int(os.environ.get('V9_BLOCK','10')))
    _,p,adj=ref.test(x);sig=adj<.05;sig[own]=False
    jk=float(jackknife_correction(sparse.csr_matrix(x),np.zeros(400,np.intp),1,50000.,chunk=512)[0]);mean=(x/(d[:,None]/1e6)).mean(0)
    rho=float(spearmanr(z['confidence'][i,sig],-np.log10(np.maximum(p[sig],1e-300))).statistic)
    row=dict(context=c,index=i,arm=arm,calls=int(sig.sum()),selected=len(called),unintended_calls=int(np.sum(sig&~np.isin(np.arange(len(genes)),called))),confidence_rho=rho,jk_v5b=jk0,jk=jk,jk_ratio=jk/jk0,mean_CPM_L1=float(abs(mean-mean0).sum()/1e6),seconds=time.time()-start,**info)
    records.append(row);print(row,flush=True)
    np.savez_compressed(R/f'artifacts/probe_{arm}_{c}_{i}.npz',x=x,original=original,called=called,p=p,q=adj,mean=mean,mean0=mean0)
    runtime.save(R/'artifacts/probes.json',records)
  del z,q,bf,ref;gc.collect()
if __name__=='__main__':main()
