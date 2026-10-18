import runtime
from remote_h5 import Remote
import h5py,anndata as ad
from pathlib import Path
r=Remote('https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad',cap=100_000_000,block=65536)
with h5py.File(r,'r') as h:
 obs=ad.io.read_elem(h['obs']);var=ad.io.read_elem(h['var']);obs.to_parquet('artifacts/H1_full_obs.parquet');var.to_parquet('artifacts/H1_full_var.parquet');print(obs.shape,obs.columns,obs.batch.value_counts());print('X',dict(h['X'].attrs))
 for k in ['data','indices','indptr']:print(k,h['X'][k].shape,h['X'][k].dtype,h['X'][k].chunks,h['X'][k].compression)
 counts=obs.groupby(['target_gene','batch'],observed=True).size().unstack(fill_value=0)
 print('balanced capacity',counts.min(axis=1).describe());print(counts.head())
runtime.save('artifacts/H1_metadata_transfer.json',{'bytes':r.bytes,'url':r.url})
r.http.close()
