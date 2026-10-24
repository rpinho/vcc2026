import runtime
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
class Writer:
 def __init__(self,path,genes):
  self.n=0;self.nnz=0;self.obs=[];self.genes=genes;self.f=h5py.File(path,'w');self.f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
  self.g=self.f.create_group('X');self.g.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0'})
  for k,dt in [('data','int32'),('indices','int32'),('indptr','int64')]:self.g.create_dataset(k,shape=(1 if k=='indptr' else 0,),maxshape=(None,),dtype=dt,chunks=(4096 if k=='indptr' else 262144,),compression='gzip',compression_opts=1,shuffle=True)
  self.g['indptr'][0]=0
 def append(self,x,c,t):
  runtime.guard();x=sparse.csr_matrix(x)
  for k,a in [('data',x.data.astype(np.int32)),('indices',x.indices.astype(np.int32)),('indptr',x.indptr[1:].astype(np.int64)+self.nnz)]:
   z=self.g[k];j=len(z);z.resize((j+len(a),));z[j:]=a
  self.n+=x.shape[0];self.nnz+=x.nnz;self.obs.extend([(c,t)]*x.shape[0])
 def close(self):
  self.g.attrs['shape']=(self.n,len(self.genes));o=pd.DataFrame(self.obs,columns=['context','target_gene'],index=[f'pred_{i}' for i in range(self.n)])
  for k in o:o[k]=pd.Categorical(o[k])
  ad.io.write_elem(self.f,'obs',o);ad.io.write_elem(self.f,'var',pd.DataFrame(index=pd.Index(self.genes)))
  for k in ['uns','obsm','varm','obsp','varp','layers']:ad.io.write_elem(self.f,k,{})
  self.f.close()
