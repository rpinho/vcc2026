"""Check the already-used Replogle GWPS source for a matched destination.
Availability only, no response-based target selection or model scoring.
"""
import runtime as rt
from remote_h5 import Remote
import h5py,anndata as ad,numpy as np,pandas as pd,time,hashlib
base='https://ndownloader.figshare.com/files/35775507'
class Fresh(Remote):
 def fetch(self,start,n):
  self.url=base+'?fresh='+str(time.time_ns());return super().fetch(start,n)
r=Fresh(base+'?fresh='+str(time.time_ns()),cap=250_000_000,block=262144)
try:
 with h5py.File(r,'r') as h:
  # This older H5AD stores categorical codes in columns and categories elsewhere.
  # Decode the DATAFRAME through AnnData; standalone column reads yield integer codes.
  whole=ad.io.read_elem(h['obs']);obs=whole[['gene','gem_group','sgID_AB','UMI_count']].copy();del whole;assert 'non-targeting' in set(obs.gene.astype(str));var=ad.io.read_elem(h['var']);ds=h['X'];layout=dict(offset=ds.id.get_offset(),shape=list(ds.shape),dtype=str(ds.dtype),chunks=ds.chunks,compression=ds.compression);assert ds.chunks is None and ds.compression is None and str(ds.dtype)=='float32'
 obs['source_row']=np.arange(len(obs));obs.to_parquet('artifacts/GWPS_raw_obs.parquet',index=False);var.to_parquet('artifacts/GWPS_raw_var.parquet')
 table=obs.groupby(['gene','gem_group'],observed=True).size().unstack(fill_value=0);table.columns=table.columns.astype(str);batches=sorted(table.columns);available=set(np.load(rt.ROOT.parent/'public_data/model_v4_signature_banks/HCT116.npz')['targets'].astype(str));available&=set(table.index);available.discard('non-targeting');trials=[];plan=None
 # Largest technical groups, chosen solely by sample count. Whole-dataset first.
 largest=obs.gem_group.astype(str).value_counts().index.tolist();subsets=[batches]+[sorted(largest[:n]) for n in [len(batches)//2,4,2,1] if 0<n<len(batches)]
 for use in subsets:
  for n in [200,160,96,64]:
   q=n//len(use)
   if q<1:continue
   eligible=[t for t in available if table.loc[t,use].min()>=q];trials.append(dict(batches=use,quota=q,n_real=q*len(use),eligible=len(eligible)))
   if plan is None and len(eligible)>=300 and table.loc['non-targeting',use].min()>=16*q:
    targets=sorted(sorted(eligible,key=lambda t:hashlib.sha256(('GWPS_panel:'+t).encode()).hexdigest())[:300]);plan=dict(context='K562',targets=targets,batches=use,quota=q,n_real=q*len(use),control_ratio=8)
 rt.save('artifacts/GWPS_matched_feasibility.json',dict(status='FEASIBLE' if plan else 'INSUFFICIENT comparable matched target panel',plan=plan,trials=trials,n_candidates=len(available),batch_sizes=obs.gem_group.astype(str).value_counts().to_dict(),layout=layout,source=base,etag=r.etag,bytes=r.bytes));print('GWPS plan',plan,'trials',trials,flush=True)
finally:rt.save('artifacts/GWPS_census_transfer.json',dict(bytes=r.bytes,source=base));r.http.close()
