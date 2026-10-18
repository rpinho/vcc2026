"""Record source identities and coverage; anonymous public metadata only."""
import runtime,requests,sys,importlib.metadata
import numpy as np,pandas as pd
rt=runtime.ROOT;records=[];targets=pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.to_numpy(str);genes=pd.read_csv(rt.parent/'data/gene_names.csv').gene_name.to_numpy(str)
for path in sorted((rt/'artifacts/banks').glob('*.npz')):
 if path.stem.endswith('_blocks'):continue
 with np.load(path) as d:
  if 'n_cells' not in d or d['n_cells'].ndim!=1:continue
  ri=pd.Index(d['targets']).get_indexer(targets);ci=pd.Index(d['genes']).get_indexer(genes);available=ri>=0
  if d['n_cells'].ndim!=1:continue
  available[available]&=d['n_cells'][ri[available]]>=20;measured=ci>=0;measured[measured]&=d['measured_genes'][ci[measured]]
  record=dict(source=path.stem,sha256=runtime.sha(path),source_targets=len(d['targets']),source_genes=len(d['genes']),official_targets=int(available.sum()),official_readouts=int(measured.sum()),target_readout_coordinates=int(available.sum()*measured.sum()),retained_cells=int(d['n_cells'].sum()),arithmetic_CPM_observed=bool(d['arithmetic_cpm_observed']) if 'arithmetic_cpm_observed' in d else True)
  records.append(record)
pd.DataFrame(records).to_csv(rt/'artifacts/source_coverage.csv',index=False)
s=requests.Session();s.trust_env=False;public={}
for name,url in [('cd4_pseudobulk','https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.pseudobulk_merged.h5ad'),('cd4_de','https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad'),('orion_dataset','https://huggingface.co/api/datasets/slaf-project/X-Atlas-Orion')]:
 r=s.get(url,stream=True,timeout=60) if name=='orion_dataset' else s.head(url,timeout=60);r.raise_for_status()
 public[name]=dict(url=url,etag=r.headers.get('ETag'),last_modified=r.headers.get('Last-Modified'),bytes=r.headers.get('Content-Length'))
 if name=='orion_dataset':public[name]['revision']=r.json()['sha']
 r.close()
packages={}
for name in ['numpy','scipy','h5py','anndata','cell-eval2','vcc-cli','numba']:
 try:packages[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:packages[name]='not installed under this distribution name'
runtime.save(rt/'artifacts/source_identity.json',dict(public=public,banks=records,packages=packages,python=sys.version,code={p.name:runtime.sha(p) for p in [rt/'model.py',rt/'vendor/atlas_model.py',rt/'vendor/adaptive_shrink.py']},verification='Hashes identify actual source estimates, not biological correctness'))
print(pd.DataFrame(records).to_string(index=False),flush=True)
