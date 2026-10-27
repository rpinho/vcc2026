"""DM feasibility probe (consult_1026_amp R10, Astra's mechanism): re-emit B0's cells for a few H1 targets
(hash-checked against render_B0.json), then measure how far integer per-cell depth multipliers x'_c = k_c x_c
can amplify the pooled log1p(5e4) delta along B0's own predicted delta. Label-free: uses only our emitted cells and
the control pseudobulk; no real perturbed cells are read here."""
import sys,os,json,time,numpy as np
H='h2h2_1026';sys.path.insert(0,H);os.chdir(H)
import score_run as sr
j=sr.j;TS=5e4
t0=time.time();a,ctrl,real,base,anchors,targets,genes,response,bulk,conf,q,depth,ref,model,layout=sr.setup('H1');print('setup',round(time.time()-t0),'s',flush=True)
with np.load(j.ROOT/'scratch/confidences.npz') as z:confidence=z['B0']
with np.load(j.ROOT/'results/predictions_B0.npz') as z:r=z['response']
rec={x['target']:x for x in json.load(open(j.ROOT/'results/render_B0.json'))['records']}
Lc=np.log1p(TS*bulk/bulk.sum())
idx=[int(x) for x in sys.argv[1:]] or [0,10,20,40,60]
out=[]
for i in idx:
  t=targets[i];own=np.flatnonzero(genes==t);own=int(own[0]) if len(own) else np.empty(0,dtype=np.intp)
  prob=np.maximum(np.expm1(np.log1p(50000*bulk)+r[i]),0);prob/=prob.sum()
  orig=model.emit(q,None,prob,depth,model.atlas.seed_for('H1:'+t,20260910),'anchor_single_mean').astype(np.int32)
  x,called,detail=layout.emit(orig,ref,confidence[i],int(base['renderer']['budget']),own)
  ok=sr.arrsha(x)==rec[t]['counts_sha256']
  mask=np.ones(len(genes),bool);mask[own]=False
  def L(w):s=(w[:,None]*x).sum(0).astype(np.float64);return np.log1p(TS*s/s.sum())
  d0=(L(np.ones(len(x)))-Lc)[mask];u=d0/np.linalg.norm(d0)
  # per-cell projection of its own log-profile deviation onto u
  cl=np.log1p(TS*x/x.sum(1,keepdims=True));proj=((cl-Lc)[:,mask]@u)
  zc=(proj-proj.mean())/proj.std()
  res=[]
  for K in [2,3,4,6]:
    for frac in [0.5,0.3,0.15]:
      w=np.where(zc>=np.quantile(zc,1-frac),K,1).astype(np.int64)
      d=(L(w)-Lc)[mask];amp=float(d @ u/np.linalg.norm(d0));orth=float(np.linalg.norm(d-(d @ u)*u)/np.linalg.norm(d0))
      res.append(dict(K=K,frac=frac,amp=round(amp,3),orth=round(orth,3),cos=round(float(d @ d0/np.linalg.norm(d)/np.linalg.norm(d0)),4),maxdepth=int((w*x.sum(1)).max())))
  best=max(res,key=lambda z:z['amp']-z['orth'])
  out.append(dict(target=t,hash_ok=ok,cells=int(len(x)),depth_med=float(np.median(x.sum(1))),best=best,grid=res));print(json.dumps(out[-1]['best']),t,'hash',ok,'cells',len(x),flush=True)
json.dump(out,open('h2h2_1026/scratch/dm_probe.json','w'),indent=1)
