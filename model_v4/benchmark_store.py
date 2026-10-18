import runtime,sys,time
sys.path.insert(0,'vendor/deps')
import lance,numpy as np,pandas as pd
from public_store import PublicStore
c=pd.read_parquet('artifacts/HCT116_selected.parquet');p=PublicStore();e=lance.dataset('s3://public/data/HCT116/expression.lance',version=682,storage_options=p.options(),block_size=4096);start=time.time()
for lo in range(0,1280,128):
 rows=c.iloc[lo:lo+128];ix=np.concatenate([np.arange(s,s+k,dtype=np.uint64) for s,k in zip(rows.cell_start_index,rows.n_genes_by_counts)]);tab=e.take(ix,columns=['cell_integer_id','gene_integer_id','value']);print(lo+128,round(time.time()-start,1),round(p.bytes/1e6),flush=True)
p.close()
