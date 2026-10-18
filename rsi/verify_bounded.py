import runtime as rt
import numpy as np,pandas as pd,anndata as ad,polars as pl
from scipy import sparse
from dataclasses import replace
import official_scoring as s
from cell_eval2 import compute_metrics
from cell_eval2.anchor import compute_replicate_anchor
rng=np.random.default_rng(int(__import__('os').environ.get('RSI_FIXTURE_SEED','10')));genes=np.array(['G'+str(i) for i in range(80)]);targets=genes[:12];base=rng.uniform(1,30,80);x=[];labels=[]
for t in ['non-targeting']+list(targets):
 rate=base.copy()
 if t!='non-targeting':rate[int(t[1:])]=.1;rate[15:35]*=rng.choice([.25,3],20)
 x.append(rng.poisson(rate,(40,80)));labels.extend([t]*40)
real=ad.AnnData(sparse.csr_matrix(np.vstack(x)),obs=pd.DataFrame({'target_gene':labels},index=['r'+str(i) for i in range(len(labels))]),var=pd.DataFrame(index=genes));pred=real.copy();pred.X.data=np.maximum(0,pred.X.data+rng.integers(-1,2,len(pred.X.data)));pred.X.eliminate_zeros();cfg=s.config();a=compute_metrics(pred,real,config=cfg).sort(['metric','perturbation']);sp0,an0=compute_replicate_anchor(real,config=cfg,n_splits=2)
import bounded_de as b
b.install_unchecked();aa=compute_metrics(pred,real,config=cfg).sort(['metric','perturbation']);sp1,an1=compute_replicate_anchor(real,config=cfg,n_splits=2)
if not a.equals(aa):
 mismatch=a.join(aa,on=['metric','perturbation'],suffix='_new').filter(pl.col('value')!=pl.col('value_new'));mismatch.write_csv('artifacts/bounded_mismatch.csv');raise AssertionError(mismatch)
assert sp0.equals(sp1) and an0.equals(an1)

# Unequal group sizes and fractional tiled baseline exercise production's two other paths.
for fractional in [False,True]:
 import cell_eval2.run as run,cell_eval2.anchor as anchor,cell_eval2.lfc_nmae_ref as ref
 run._compute_de_side=b.ORIGINAL;anchor._compute_de_side=b.ORIGINAL;ref._compute_de_side=b.ORIGINAL
 pp=pred.copy()
 if fractional:
  pp.X=pp.X.astype(float)*.975;cc=replace(cfg,allow_fractional_counts=True)
 else:
  pp=pp[np.arange(pp.n_obs)%3!=0].copy();cc=cfg
 original=compute_metrics(pp,real,config=cc).sort(['metric','perturbation'])
 b.install_unchecked();adapted=compute_metrics(pp,real,config=cc).sort(['metric','perturbation'])
 assert original.equals(adapted),(fractional,original.join(adapted,on=['metric','perturbation'],suffix='_new').filter(pl.col('value')!=pl.col('value_new')))

rt.save('artifacts/bounded_fixture.json',dict(status='PASS bitwise equality',metric_rows=a.height,anchor_rows=sp0.height,fixture='12 targets x 80 genes; equal and unequal group sizes; integer and fractional baseline; two-split anchors',adapter_sha256=rt.sha('bounded_de.py'),cell_eval2=s.provenance()));print('PASS bitwise equality',flush=True)
