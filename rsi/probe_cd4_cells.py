"""Bounded anonymous metadata-only feasibility probe; no model or holdout scoring."""
import runtime as rt
from remote_h5 import Remote
import h5py
url='https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/D4_Stim8hr.assigned_guide.h5ad';r=Remote(url,cap=10_000_000,block=65536)
try:
 with h5py.File(r,'r') as h:
  x=h['X'];layout={k:dict(shape=list(x[k].shape),dtype=str(x[k].dtype),chunks=x[k].chunks,compression=x[k].compression) for k in x} if isinstance(x,h5py.Group) else dict(shape=list(x.shape),dtype=str(x.dtype),chunks=x.chunks,compression=x.compression,offset=x.id.get_offset());obs={k:dict(kind='group' if isinstance(v,h5py.Group) else 'dataset',shape=list(v.shape) if isinstance(v,h5py.Dataset) else {a:list(v[a].shape) for a in v}) for k,v in h['obs'].items()};result=dict(source=url,etag=r.etag,source_size=r.size,x_attributes=dict(x.attrs),layout=layout,obs=obs,var_columns=list(h['var']))
 rt.save('artifacts/CD4_cell_feasibility.json',result);print(result,flush=True)
finally:rt.save('artifacts/CD4_probe_transfer.json',dict(bytes=r.bytes,source=url));r.http.close()
