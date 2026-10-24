"""Public-control predicted-call diagnostic; not fitness, no treated labels.
Column-streamed tie-corrected MW normal approximation, no continuity correction.
Only matched probe indices 4 and12 in each arm; all 18,533 genes, all18,400 controls.
"""
import runtime
import numpy as np,pandas as pd,anndata as ad,json,gc
from scipy import sparse
from scipy.special import ndtr
from scipy.stats import false_discovery_control,mannwhitneyu
R=runtime.ROOT
fs=[p for p in sorted((R/'artifacts').glob('fixture_*.npz')) if int(p.stem.split('_')[-1]) in (4,12)]
ys=[];masks=[]
for f in fs:
 with np.load(f) as z:ys.append(z['x'].astype(float)/z['depths'][:,None]*1e6)
 with np.load(R/'artifacts'/f"fixture_b_{int(f.stem.split('_')[-1])}.npz") as z:masks.append(z['mask'])
assert len(ys)==4
# Keep sparse control values in original float32, normalize one gene at a time.
a=ad.read_h5ad(R.parent/'data/context_A.h5ad');X=sparse.csr_matrix(a.X);dep=np.asarray(X.sum(1)).ravel().astype(float);C=X.tocsc();genes=a.var_names.to_numpy();n=C.shape[0];del a,X;gc.collect()
pvals=np.ones((len(ys),C.shape[1]));cmean=np.zeros(C.shape[1]);checks=0
for j in range(C.shape[1]):
 lo,hi=C.indptr[j:j+2];v=C.data[lo:hi].astype(float)*(1e6/dep[C.indices[lo:hi]]);cmean[j]=v.sum()/n
 u,ct=np.unique(np.r_[0.,v],return_counts=True);ct[0]=n-len(v)
 if ct[0]==0:u=u[1:];ct=ct[1:]
 cum=np.r_[0,np.cumsum(ct)];tie0=float(np.sum(ct.astype(float)**3-ct))
 for f,Y in enumerate(ys):
  y=Y[:,j];m=len(y);uy,cy=np.unique(y,return_counts=True);pos=np.searchsorted(u,uy);present=pos<len(u);present[present]&=u[pos[present]]==uy[present]
  cc=np.zeros(len(uy),int);cc[present]=ct[pos[present]];less=cum[np.minimum(pos,len(u))];U=float(np.sum(cy*(less+.5*cc)));tie=tie0+np.sum((cc+cy).astype(float)**3-(cc+cy))-np.sum(cc.astype(float)**3-cc)
  var=m*n/12*((m+n+1)-tie/((m+n)*(m+n-1)))
  if var>0:pvals[f,j]=2*ndtr(-abs((U-.5*m*n)/np.sqrt(var)))
  if f==0 and checks<8 and var>0:
   exact=mannwhitneyu(y,np.r_[np.zeros(n-len(v)),v],alternative='two-sided',method='asymptotic',use_continuity=False).pvalue
   assert np.isclose(pvals[f,j],exact,rtol=1e-10,atol=1e-14),(j,pvals[f,j],exact)
   checks+=1
 if j%3000==0:print('control-column',j,flush=True)
T=pd.read_csv(R.parent/'data/pert_counts.csv').target_gene.to_numpy(str);out=[]
for f,Y,p,mask in zip(fs,ys,pvals,masks):
 q=false_discovery_control(p);sig=q<.05;target=T[int(f.stem.split('_')[-1])];sig[genes==target]=False
 sgn=np.sign(Y.mean(0)-cmean);v=np.var(Y,axis=0)
 out.append(dict(file=f.name,target=target,predicted_significant=int(sig.sum()),significant_supported=int((sig&mask).sum()),significant_other=int((sig&~mask).sum()),up=int((sig&(sgn>0)).sum()),down=int((sig&(sgn<0)).sum()),median_variance_supported=float(np.median(v[mask])),median_variance_other=float(np.median(v[~mask])) if (~mask).any() else None));print(out[-1],flush=True)
(R/'artifacts/predicted_de.json').write_text(json.dumps(dict(status='diagnostic only; no treated labels; no candidate ranking',control_cells=n,scipy_crosschecks=checks,records=out),indent=2))
