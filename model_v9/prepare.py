import runtime,gc,time
import numpy as np,anndata as ad,pandas as pd
from scipy import sparse
from de import prepare,Reference
from cell_eval2.moments import jackknife_correction
R=runtime.ROOT
for ci,c in enumerate('ABC'):
 runtime.guard();out=R/'artifacts'/f'ref_{c}'
 if (out/'receipt.json').exists():continue
 a=ad.read_h5ad(R.parent/f'data/context_{c}.h5ad');x=sparse.csr_matrix(a.X);del a;gc.collect()
 prepare(x,out,np.random.default_rng(20261018+ci))
 jks=[];rng=np.random.default_rng(73591+ci)
 for _ in range(8):
  runtime.guard();z=x[rng.choice(x.shape[0],400,replace=False)];jks.append(float(jackknife_correction(z,np.zeros(400,np.intp),1,50000.,chunk=32)[0]))
 runtime.save(out/'receipt.json',dict(control_cells=x.shape[0],genes=x.shape[1],control400_jk=jks,cap_jk=float(np.mean(jks)),threshold='reference arithmetic mean CPM > 5, before per-target BH',status='controls only'))
 print('REFERENCE',c,'JK',np.mean(jks),flush=True);del x;gc.collect()
