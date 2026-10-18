"""Predeclared depth sensitivity, not another independent lineage vote."""
import runtime as rt
import sys
from pathlib import Path
import numpy as np,anndata as ad,shutil
path=Path('artifacts/tasks/manifest.json');manifest=rt.json.loads(path.read_text())

# Restrict each destination to its measured intersection with the actual official readout list.
# Missing readouts are excluded, never padded with invented zero observations.
import pandas as pd
from scipy import sparse
official=pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:,0].astype(str).to_numpy()
for task,m in list(manifest.items()):
 if m.get('view'):continue
 file=Path('artifacts/tasks')/(task+'.h5ad')
 if not file.exists():continue
 receipt=rt.json.loads(file.with_suffix('.json').read_text())
 if receipt.get('readout_policy')=='measured official intersection':continue
 a=ad.read_h5ad(file,backed='r');cols=pd.Index(a.var_names).get_indexer(official);genes=official[cols>=0];cols=cols[cols>=0];blocks=[]
 for lo in range(0,a.n_obs,512):blocks.append(sparse.csr_matrix(a.X[lo:lo+512])[:,cols])
 b=ad.AnnData(sparse.vstack(blocks,format='csr'),obs=a.obs.copy(),var=pd.DataFrame(index=genes));b.uns['evaluation_only']=True;raw_genes=a.n_vars;a.file.close();rt.guard(2*1024**3);temp=file.with_suffix('.projected.h5ad');b.write_h5ad(temp,compression='gzip',compression_opts=1);temp.replace(file);receipt.update(unprojected_sha256=receipt['sha256'],sha256=rt.sha(file),readout_policy='measured official intersection',original_readouts=raw_genes,scored_readouts=len(genes));rt.save(file.with_suffix('.json'),receipt);m['readouts']=len(genes);m['readout_policy']='measured official intersection';print('Official readout intersection',task,raw_genes,len(genes),flush=True);del a,b,blocks

if '--projection-only' in sys.argv:
 rt.save(path,manifest);raise SystemExit(0)

for split,seed in [('search',815),('holdout',1815)]:
 parent='H1_'+split;task='H1_depth20k_'+split;file=Path('artifacts/tasks')/(task+'.h5ad');m=dict(manifest[parent])
 if not file.exists():
  original=file.with_name(parent+'.h5ad');a=ad.read_h5ad(original);x=a.X.tocsr();depth=np.asarray(x.sum(1)).ravel();rng=np.random.default_rng(seed);x.data=rng.binomial(x.data.astype(np.int64),np.repeat(np.minimum(1,20000/depth),np.diff(x.indptr))).astype(np.int32);x.eliminate_zeros();a.X=x;a.write_h5ad(file,compression='gzip',compression_opts=1);rt.save(file.with_suffix('.json'),dict(sha256=rt.sha(file),parent_sha256=rt.sha(original),thinning='Per-cell binomial min(1,20000/depth)',seed=seed));shutil.copyfile(file.with_name(parent+'_selection.parquet'),file.with_name(task+'_selection.parquet'))
 m.update(source='H1_depth20k',parent_task=parent,view='depth20k',sensitivity_only=True);manifest[task]=m
rt.save(path,manifest)
