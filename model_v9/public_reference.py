"""Held-source DE tables on released real cells, using the official producer.

These sources never enter v5b/v9 confidence. Native public sample sizes and
gene coverage are retained; transported metrics are not hidden score forecasts.
"""
import runtime,gc,sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from dataclasses import replace
from cell_eval2 import EvalConfig
from cell_eval2.run import _compute_de_side
R=runtime.ROOT
def main(name):
 runtime.guard()
 a=ad.read_h5ad(R.parent/f'public_data/{name}.h5ad',backed='r')
 label=a.obs['gene' if 'gene' in a.obs else 'target_gene'].astype(str).to_numpy()
 genes=a.var['gene_name'].astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
 bank=R.parent/'public_data/model_v4_signature_banks'
 shared=set.intersection(*[set(np.load(bank/f'{b}.npz')['targets'].astype(str)) for b in ['K562_GWPS_raw','HCT116','HEK293T']])
 counts=pd.Series(label).value_counts();selected=sorted(t for t in shared&set(genes) if counts.get(t,0)>=40)
 if len(selected)>16:
  chosen=np.random.default_rng(20261020).choice(len(selected),16,replace=False)
  selected=[selected[i] for i in sorted(chosen)]
 assert selected,('no independent public targets with multi-source coverage',name)
 controls=np.flatnonzero(np.char.startswith(label.astype(str),'non-targeting')|np.isin(np.char.lower(label.astype(str)),['ntc','control','ctrl','unperturbed']))
 rng=np.random.default_rng(20261019);controls=np.sort(rng.choice(controls,min(len(controls),18400),replace=False))
 sets=[controls];labs=['non-targeting']*len(controls);sizes={}
 for t in selected:
  rows=np.flatnonzero(label==t);rows=np.sort(rng.choice(rows,min(len(rows),400),replace=False));sets.append(rows);labs.extend([t]*len(rows));sizes[t]=len(rows)
 requested=np.concatenate(sets);order=np.argsort(requested);pieces=[]
 for start in range(0,len(order),512):
  runtime.guard();pieces.append(sparse.csr_matrix(a.X[requested[order[start:start+512]]]))
 X=sparse.vstack(pieces,format='csr');del pieces
 lab=np.asarray(labs)[order];a.file.close();del a;gc.collect()
 # Public feature IDs can map to the same gene symbol. Sum those raw counts
 # before normalization, rather than duplicate or relabel biological genes.
 unique=pd.Index(pd.unique(genes));codes=unique.get_indexer(genes)
 if len(unique)!=len(genes):
  mapping=sparse.csr_matrix((np.ones(len(genes),np.int32),(np.arange(len(genes)),codes)),shape=(len(genes),len(unique)))
  X=(X@mapping).tocsr();genes=unique.to_numpy(str)
 mat=ad.AnnData(X,obs=pd.DataFrame({'target':lab},index=[str(i) for i in range(len(lab))]),var=pd.DataFrame(index=genes))
 mat.write_h5ad(R/f'out/public_{name}.h5ad',compression='gzip',compression_opts=1)
 cfg=EvalConfig.from_preset('vcc2026');cfg=replace(cfg,device='cpu',num_threads=1,gather_threads=1,de=replace(cfg.de,backend='pdex'))
 frame=_compute_de_side(mat,cfg=cfg,fp=None,store=None,side='real')
 frame.write_parquet(R/f'artifacts/public_{name}.parquet')
 runtime.save(R/f'artifacts/public_{name}.json',dict(status='exact official _compute_de_side on held public source',targets=sizes,reference_cells=len(controls),native_genes=len(genes),rows=frame.height,scorer=runtime.scorer_hashes(),limitations='Native public power and gene coverage; context transport. These are released cells, not challenge treated labels.'))
 print('PUBLIC REFERENCE READY',name,len(selected),frame.height,flush=True)
if __name__=='__main__':main(sys.argv[1])
