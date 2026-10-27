"""Strongest form of the DM mechanism: continuous non-negative per-cell weights solved by NNLS so the pooled profile
matches control + s * (B0 delta) in variance-stabilised space. If this cannot reach s=1.38 with small orthogonal
error, no integer depth reweighting can."""
import sys,os,json,numpy as np
from scipy.optimize import nnls
H='h2h2_1026';sys.path.insert(0,H);os.chdir(H)
import score_run as sr
j=sr.j;TS=5e4
a,ctrl,real,base,anchors,targets,genes,response,bulk,conf,q,depth,ref,model,layout=sr.setup('H1')
with np.load(j.ROOT/'scratch/confidences.npz') as z:confidence=z['B0']
with np.load(j.ROOT/'results/predictions_B0.npz') as z:r=z['response']
Lc=np.log1p(TS*bulk/bulk.sum());res=[]
for i in [int(x) for x in sys.argv[1:]]:
  t=targets[i];own=np.flatnonzero(genes==t);own=int(own[0]) if len(own) else np.empty(0,dtype=np.intp)
  prob=np.maximum(np.expm1(np.log1p(50000*bulk)+r[i]),0);prob/=prob.sum()
  orig=model.emit(q,None,prob,depth,model.atlas.seed_for('H1:'+t,20260910),'anchor_single_mean').astype(np.int32)
  x,_,_=layout.emit(orig,ref,confidence[i],int(base['renderer']['budget']),own)
  mask=np.ones(len(genes),bool);mask[own]=False
  f0=x.sum(0)/x.sum();L0=np.log1p(TS*f0);d0=(L0-Lc)[mask];n0=np.linalg.norm(d0);u=d0/n0
  keep=(f0>0)&mask
  for s in [1.2,1.38]:
    Ls=Lc+s*(L0-Lc);fs=np.expm1(Ls);fs=fs/fs.sum()
    # rows = genes (variance-stabilising weight 1/sqrt(fs)), cols = cells' count vectors; sum constraint via a heavy row
    D=x.sum(1).astype(np.float64);sc=1/np.sqrt(np.maximum(fs[keep],1e-9))
    A=(x[:,keep].T.astype(np.float64)-np.outer(fs[keep],D))*sc[:,None]
    A=np.vstack([A,1e3*(D/D.mean())[None,:]]);bvec=np.r_[np.zeros(keep.sum()),1e3*len(x)]
    w,_=nnls(A,bvec,maxiter=5000)
    f=(w[:,None]*x).sum(0);f=f/f.sum();d=(np.log1p(TS*f)-Lc)[mask]
    amp=float(d @ u/n0);orth=float(np.linalg.norm(d-(d @ u)*u)/n0)
    wi=np.maximum(1,np.round(w/w[w>0].min()) if (w>0).any() else w)
    res.append(dict(target=t,s=s,amp=round(amp,3),orth=round(orth,3),cos=round(float(d @ d0/np.linalg.norm(d)/n0),4),
      nonzero=int((w>0).sum()),wmax_over_wmin=float(w.max()/w[w>0].min()) if (w>0).any() else None))
    print(json.dumps(res[-1]),flush=True)
json.dump(res,open(j.ROOT/'scratch/dm_probe_nnls.json','w'),indent=1)
