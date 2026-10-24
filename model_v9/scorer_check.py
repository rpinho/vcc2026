"""Official DE producer, catalog dispatch, raw metric sensitivity and schema proof."""
import runtime,gc,sys,h5py
import numpy as np,pandas as pd,anndata as ad,polars as pl
from scipy import sparse
from dataclasses import replace
from cell_eval2 import EvalConfig
from cell_eval2.run import _compute_de_side,dispatch_de_metrics
from cell_eval2.de import prepare_de
from de import Reference
R=runtime.ROOT
NAMES=['de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard','de_wilcoxon_lfc_nmae']
TABLE_AXES={}
def config():
 cfg=EvalConfig.from_preset('vcc2026')
 return replace(cfg,device='cpu',num_threads=1,gather_threads=1,de=replace(cfg.de,backend='pdex'))
def table(c,t,p,q,lfc):
 if c not in TABLE_AXES:TABLE_AXES[c]=(pd.read_csv(R.parent/'data/gene_names.csv').gene_name.to_numpy(str),np.load(R/f'artifacts/ref_{c}/mean.npy')>5)
 genes,m=TABLE_AXES[c]
 return pl.DataFrame(dict(target=[t]*int(m.sum()),feature=genes[m],p_value=p[m],p_adj=q[m],log2_fold_change=lfc[m]))
def score(pred,real,return_prepared=False):
 cfg=config();prepared=prepare_de(pred,real,control=cfg.control,sort_by=cfg.de.sort_by,p_adj_threshold=cfg.de.p_adj_threshold,nan_lfc_policy=cfg.de.nan_lfc_policy,min_abs_log2fc=cfg.de.min_abs_log2fc)
 rows=dispatch_de_metrics(NAMES,prepared,cfg)
 return (rows,prepared) if return_prepared else rows
def main():
 records=[];genes=pd.read_csv(R.parent/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(R.parent/'data/pert_counts.csv').target_gene.to_numpy(str)
 for ci,c in enumerate('ABC'):
  runtime.guard();control=ad.read_h5ad(R.parent/f'data/context_{c}.h5ad');ref=Reference(c)
  for arm in ['a','b']:
   if '--emitted' in sys.argv:
    d=np.load(R/f'artifacts/de/v9{arm}_{c}_0.npz')
    with h5py.File(R/f'out/v9{arm}.h5ad') as h:
     X=h['X'];ptr=X['indptr'][ci*120000:ci*120000+401];lo,hi=int(ptr[0]),int(ptr[-1]);x=sparse.csr_matrix((X['data'][lo:hi],X['indices'][lo:hi],ptr-lo),shape=(400,len(genes))).toarray()
   else:d=np.load(R/f'artifacts/probe_{arm}_{c}_0.npz');x=d['x']
   label=targets[0]
   mat=ad.AnnData(sparse.vstack([sparse.csr_matrix(control.X),sparse.csr_matrix(x)],format='csr'),obs=pd.DataFrame({'target':['non-targeting']*control.n_obs+[label]*400},index=[str(i) for i in range(control.n_obs+400)]),var=pd.DataFrame(index=genes))
   frame=_compute_de_side(mat,cfg=config(),fp=None,store=None,side='pred');f=frame.to_pandas().set_index('feature');ix=pd.Index(genes).get_indexer(f.index)
   assert set(f.index)==set(genes[ref.keep])
   assert np.allclose(f.p_value.to_numpy(),d['p'][ix],atol=1e-13,rtol=1e-9)
   assert np.array_equal(f.p_adj.to_numpy()<.05,d['q'][ix]<.05)
   expected=np.log2((d['mean'][ix]+1e-9)/(ref.mean[ix]+1e-9))
   assert np.allclose(f.log2_fold_change.to_numpy(),expected,rtol=1e-8,atol=1e-8)
   frame.write_parquet(R/f'artifacts/scorer_{arm}_{c}.parquet')
   records.append(dict(context=c,arm=arm,genes=len(ix),input='final emitted H5AD' if '--emitted' in sys.argv else 'layout probe',p_max_error=float(abs(f.p_value.to_numpy()-d['p'][ix]).max()),mean_LFC_max_error=float(abs(f.log2_fold_change.to_numpy()-expected).max()),status='PASS official full-control producer, BH masks and arithmetic LFC'))
   print(records[-1],flush=True);runtime.save(R/'artifacts/scorer_check.json',records);del mat;gc.collect()
  del control,ref;gc.collect()
if __name__=='__main__':main()
