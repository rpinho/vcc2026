import runtime,sys
sys.path.insert(0,'vendor/deps')
import lance,pandas as pd,numpy as np,h5py,anndata as ad
from remote_h5 import Remote
from pathlib import Path
rt=runtime.ROOT
url='https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/pert_counts_Training.csv'
import requests
s=requests.Session();s.trust_env=False;r=s.get(url);r.raise_for_status();(rt/'artifacts/h1_targets.csv').write_text(r.text)
targets=list(dict.fromkeys(pd.read_csv('../data/pert_counts.csv').target_gene.tolist()+pd.read_csv(rt/'artifacts/h1_targets.csv').target_gene.tolist()));targets=[t for t in targets if t!='non-targeting']
runtime.save(rt/'artifacts/targets.json',targets)
for name in ['HCT116','HEK293T']:
 base=f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/'
 d=lance.dataset(base+'cells.lance',storage_options={'hf_token':''})
 labels=','.join("'"+t+"'" for t in targets+['Non-Targeting'])
 cells=d.to_table(columns=['cell_integer_id','cell_start_index','n_genes_by_counts','total_counts','sample','gene_target','pass_guide_filter'],filter=f'gene_target IN ({labels})').to_pandas()
 cells.to_parquet(rt/f'artifacts/{name}_cells.parquet')
 v=lance.dataset(base+'genes.lance',storage_options={'hf_token':''}).to_table().to_pandas();v.to_parquet(rt/f'artifacts/{name}_genes.parquet')
 print(name,'cells',len(cells),'targets',cells.gene_target.nunique(),'entries',cells.n_genes_by_counts.sum(),'bytes',lance.bytes_read_counter(),flush=True)
for name in ['GWCD4i.pseudobulk_merged.h5ad','GWCD4i.DE_stats.h5ad']:
 remote=Remote('https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/'+name,cap=250_000_000)
 with h5py.File(remote,'r') as h:
  obs=ad.io.read_elem(h['obs']);var=ad.io.read_elem(h['var']);obs.to_parquet(rt/f'artifacts/{name}.obs.parquet');var.to_parquet(rt/f'artifacts/{name}.var.parquet')
  print(name,'shape',len(obs),len(var),'obs',obs.columns.tolist(),'X',h['X'],'metadata bytes',remote.bytes,flush=True)
  print(obs.head(2).to_string(),flush=True)
  if isinstance(h['X'],h5py.Group):print([(k,h['X'][k].shape,h['X'][k].compression,h['X'][k].chunks) for k in h['X']],flush=True)
  if 'layers' in h:print('layers',[(k,h['layers'][k].shape,h['layers'][k].chunks,h['layers'][k].compression) for k in h['layers']],flush=True)
 runtime.save(rt/f'artifacts/{name}.probe.json',dict(size=remote.size,bytes=remote.bytes))
