import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())
for name,m in manifest.items():
 if m['source']!='local_raw':continue
 path=Path('artifacts/tasks')/(name+'.h5ad')
 if path.exists():continue
 rt.guard(2*1024**3);sel=pd.read_parquet(path.with_name(name+'_selection.parquet'));a=ad.read_h5ad(rt.ROOT.parent/'public_data'/f'{m["context"]}.h5ad',backed='r');order=np.argsort(sel.source_row.to_numpy());sel=sel.iloc[order].copy();blocks=[]
 for lo in range(0,len(sel),512):blocks.append(sparse.csr_matrix(a.X[sel.source_row.to_numpy()[lo:lo+512]]))
 genes=a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str);x=sparse.vstack(blocks,format='csr');unique=np.array(sorted(set(genes)));proj=sparse.csr_matrix((np.ones(len(genes),np.int32),(np.arange(len(genes)),pd.Index(unique).get_indexer(genes))),shape=(len(genes),len(unique)));x=(x@proj).astype(np.int32);assert np.all(x.data>=0);b=ad.AnnData(x,obs=sel[['target_gene','batch','evaluation_role','source_row']].copy(),var=pd.DataFrame(index=unique));b.uns['evaluation_only']=True;b.write_h5ad(path,compression='gzip',compression_opts=1);rt.save(path.with_suffix('.json'),dict(sha256=rt.sha(path),source='../public_data/'+m['context']+'.h5ad',selected_rows_sha256=m['selection_sha256']));a.file.close();print('Saved',name,b.shape,flush=True)
