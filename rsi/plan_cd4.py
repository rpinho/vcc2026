"""Response-blind census for a genuinely batch-matched CD4 destination.
Raw guide assignments are mapped through the authors' curated pseudobulk guide table.
No response matrix or HOLDOUT score is inspected.
"""
import runtime as rt
from remote_h5 import Remote
from pathlib import Path
import anndata as ad,h5py,numpy as np,pandas as pd,hashlib
base='https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/'
receipts=[]
def columns(filename,fields,cap):
 r=Remote(base+filename,cap=cap,block=65536)
 try:
  with h5py.File(r,'r') as h:
   out={k:ad.io.read_elem(h['obs/'+k]) for k in fields}
  return pd.DataFrame(out)
 finally:
  receipts.append(dict(source=base+filename,etag=r.etag,bytes=r.bytes));rt.save('artifacts/CD4_census_transfer.json',receipts);r.http.close()
cur=columns('GWCD4i.pseudobulk_merged.h5ad',['guide_id','perturbed_gene_name','guide_type'],80_000_000)
mapping=cur[['guide_id','perturbed_gene_name']].astype(str).drop_duplicates();ambiguous=mapping.guide_id.value_counts();mapping=mapping[~mapping.guide_id.isin(ambiguous[ambiguous>1].index)];mapping=mapping.set_index('guide_id').perturbed_gene_name.to_dict()
x=columns('D4_Stim8hr.assigned_guide.h5ad',['guide_id','lane_id','low_quality','guide_type','total_counts'],100_000_000)
x['source_row']=np.arange(len(x));x['target_gene']=x.guide_id.astype(str).map(mapping);x.loc[x.guide_type.astype(str)=='non-targeting','target_gene']='non-targeting';x['batch']=x.lane_id.astype(str)
x=x[(~x.low_quality.astype(bool))&x.target_gene.notna()&(x.total_counts>0)].copy();x.to_parquet('artifacts/CD4_cell_census.parquet',index=False)
table=x.groupby(['target_gene','batch'],observed=True).size().unstack(fill_value=0);batches=sorted(table.columns)
shared=rt.ROOT.parent/'public_data/model_v4_signature_banks';available=set()
for bank in ['K562_GWPS_raw','HCT116']:
 p=shared/(bank+'.npz')
 if p.exists():available.update(np.load(p)['targets'].astype(str))
available&=set(np.load('artifacts/sources/K562_GWPS_bulk.npz')['targets'].astype(str));available&=set(table.index);available.discard('non-targeting')
tests=[];chosen=None
for use in [batches,batches[:12],batches[:8]]:
 for q in [8,4,2,1]:
  eligible=[t for t in available if table.loc[t,use].min()>=q];tests.append(dict(lanes=len(use),quota=q,n_real=q*len(use),eligible=len(eligible)))
  if chosen is None and len(eligible)>=300 and q*len(use)>=64 and table.loc['non-targeting',use].min()>=16*q:
   targets=sorted(sorted(eligible,key=lambda t:hashlib.sha256(('CD4_panel:'+t).encode()).hexdigest())[:300]);chosen=dict(targets=targets,batches=use,quota=q,n_real=q*len(use))
rt.save('artifacts/CD4_matched_feasibility.json',dict(status='FEASIBLE' if chosen else 'INSUFFICIENT matched cells in comparably covered universe',plan=chosen,trials=tests,source_rows=len(x),source_targets=len(table)-1,comparably_covered_targets=len(available),ambiguous_guides_excluded=int((ambiguous>1).sum()),transfer=receipts,policy='Fixed lanes and abundance eligibility only; no response-strength selection. Entire CD4 source family must be excluded before destination scoring.'))
print('CD4 census',len(x),'rows; plan',chosen,'trials',tests,flush=True)
