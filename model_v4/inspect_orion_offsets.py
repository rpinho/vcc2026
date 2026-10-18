import runtime,sys,os
os.environ['RAYON_NUM_THREADS']='1';sys.path.insert(0,'vendor/deps')
import lance,pandas as pd,numpy as np
from public_store import PublicStore
c=pd.read_parquet('artifacts/HCT116_selected.parquet');lo=52608;rows=c.iloc[lo:lo+256];p=PublicStore();e=lance.dataset('s3://public/data/HCT116/expression.lance',version=682,storage_options=p.options(),block_size=4096)
try:
 for k in range(0,len(rows),32):
  rr=rows.iloc[k:k+32];ix=np.ravel(np.column_stack([rr.cell_start_index.to_numpy(),rr.cell_start_index.to_numpy()+rr.n_genes_by_counts.to_numpy()-1]));tab=e.take(ix.astype('uint64'),columns=['cell_integer_id','gene_integer_id','value']);actual=tab['cell_integer_id'].to_numpy().reshape(-1,2);df=rr[['cell_integer_id','cell_start_index','n_genes_by_counts']].copy();df['actual_first']=actual[:,0];df['actual_last']=actual[:,1];bad=df[(df.cell_integer_id!=df.actual_first)|(df.cell_integer_id!=df.actual_last)];print('batch',lo+k,'bad',len(bad),flush=True)
  if len(bad):print(bad.head(20).to_string(),flush=True)
finally:p.close()
