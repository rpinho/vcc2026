import runtime
from remote_h5 import Remote
import h5py,anndata as ad
r=Remote('https://ndownloader.figshare.com/files/35775507',cap=1_000_000_000,block=1024**2)
print('remote bytes',r.size,flush=True)
with h5py.File(r,'r') as h:
 o=ad.io.read_elem(h['obs']);v=ad.io.read_elem(h['var']);o.to_parquet('artifacts/K562_GWPS_raw.obs.parquet');v.to_parquet('artifacts/K562_GWPS_raw.var.parquet');print('metadata',len(o),len(v),r.bytes,flush=True);print(o.columns.tolist(),flush=True);print([(k,h['X'][k].shape,h['X'][k].compression,h['X'][k].chunks,str(h['X'][k].dtype)) for k in h['X']],flush=True)
runtime.save('artifacts/K562_probe_transfer.json',dict(read_bytes=r.bytes,size=r.size))
