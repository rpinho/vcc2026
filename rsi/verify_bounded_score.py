"""Bitwise compare bounded official reductions with unmodified official scoring."""
import runtime as rt
import numpy as np,pandas as pd,anndata as ad,polars as pl
from scipy import sparse
from dataclasses import replace
import official_scoring as s
from cell_eval2 import compute_metrics
import cell_eval2.anchor as anchor
rng=np.random.default_rng(int(__import__('os').environ.get('RSI_FIXTURE_SEED','10')));genes=np.array(['G'+str(i) for i in range(80)]);targets=genes[:12];base=rng.uniform(1,30,80);x=[];labels=[]
for t in ['non-targeting']+list(targets):
 rate=base.copy()
 if t!='non-targeting':rate[int(t[1:])]=.1;rate[15:35]*=rng.choice([.25,3],20)
 x.append(rng.poisson(rate,(40,80)));labels.extend([t]*40)
real=ad.AnnData(sparse.csr_matrix(np.vstack(x)),obs=pd.DataFrame({'target_gene':labels},index=['r'+str(i) for i in range(len(labels))]),var=pd.DataFrame(index=genes));pred=real.copy();pred.X.data=np.maximum(0,pred.X.data+rng.integers(-1,2,len(pred.X.data)));pred.X.eliminate_zeros();cfg=s.config();a=compute_metrics(pred,real,config=cfg).sort(['metric','perturbation']);one_thread=compute_metrics(pred,real,config=replace(cfg,num_threads=1)).sort(['metric','perturbation']);assert a.equals(one_thread);sp0,an0=anchor.compute_replicate_anchor(real,config=cfg,n_splits=2)
import bounded_de as de,bounded_score as b
de.install();aa=b.resident(pred,real,config=cfg).sort(['metric','perturbation']);assert a.equals(aa),a.join(aa,on=['metric','perturbation'],suffix='_new').filter(pl.col('value')!=pl.col('value_new'))
anchor.compute_metrics=b.resident;sp1,an1=anchor.compute_replicate_anchor(real,config=cfg,n_splits=2);assert sp0.equals(sp1) and an0.equals(an1)
# Omitted pred controls; unequal cell groups; and exact official tiled baseline.
for case in ['no_controls','unequal','baseline']:
 if case=='no_controls':pp=pred[pred.obs.target_gene.astype(str)!='non-targeting'].copy();cc=cfg
 elif case=='unequal':pp=pred[np.arange(pred.n_obs)%3!=0].copy();cc=cfg
 else:
  from cell_eval2.baseline import GenericProfile,build_baseline_prediction
  mean=np.asarray(real.X[np.array(labels)=='non-targeting'].mean(0)).ravel();p=GenericProfile(mean,genes,12,False,0);pp=build_baseline_prediction(p,real,pert_col='target_gene',control='non-targeting',emit='tile');cc=replace(cfg,allow_fractional_counts=True)
 import cell_eval2.run as run
 run._compute_de_side=de.ORIGINAL;ordinary=ad.concat([pp,real[real.obs.target_gene.astype(str)=='non-targeting'].copy()],join='inner',index_unique=':ref') if case=='no_controls' else pp;original=compute_metrics(ordinary,real,config=cc).sort(['metric','perturbation']);de.install()
 adapted=b.resident(pp,real,config=cc).sort(['metric','perturbation']);assert original.equals(adapted),(case,original.join(adapted,on=['metric','perturbation'],suffix='_new').filter(pl.col('value')!=pl.col('value_new')))
 if case=='baseline':
  virtual=b.score_iter(real,b.baseline_iter(real,mean),cc)[0].sort(['metric','perturbation']);assert original.equals(virtual),(case,original.join(virtual,on=['metric','perturbation'],suffix='_new').filter(pl.col('value')!=pl.col('value_new')))

# Production baseline must use 400 prediction cells even where reference groups are smaller.
base400=[]
for t in sorted(targets):base400.append(b.make(np.tile(mean.astype(np.float32),(400,1)),t,genes))
base400.append(real[real.obs.target_gene.astype(str)=='non-targeting'].copy());ordinary400=ad.concat(base400,index_unique=':b');import cell_eval2.run as run
run._compute_de_side=de.ORIGINAL;expected400=compute_metrics(ordinary400,real,config=cc).sort(['metric','perturbation']);de.install();actual400=b.score_iter(real,b.baseline_iter(real,mean,cells=400),cc)[0].sort(['metric','perturbation']);assert expected400.equals(actual400)

rt.save('artifacts/bounded_score_fixture.json',dict(status='PASS bitwise equality',driver_sha256=rt.sha('bounded_score.py'),metric_rows=a.height,anchor_rows=sp0.height,cases=['integer counts','no submitted controls','unequal cell counts','official dense tiled baseline versus virtual groups','official two-split anchors'],official=s.provenance()));print('PASS bounded score bitwise equality',flush=True)
