"""Split input controls from a loaded CSR without a second full reference copy."""
import numpy as np,anndata as ad
from scipy import sparse

def split(path):
 a=ad.read_h5ad(path);roles=a.obs.evaluation_role.astype(str).to_numpy();input_mask=roles=='input_control';ctrl=a[input_mask].copy();keep=np.flatnonzero(~input_mask);x=a.X
 assert sparse.isspmatrix_csr(x)
 old=x.indptr;length=np.diff(old)[keep];ptr=np.r_[0,np.cumsum(length,dtype=np.int64)];dest=0
 # Compact forward. A bounded overlapping copy cannot overwrite any future source.
 boundaries=np.r_[0,np.flatnonzero(np.diff(keep)!=1)+1,len(keep)]
 for left,right in zip(boundaries[:-1],boundaries[1:]):
  start=int(old[keep[left]]);stop=int(old[keep[right-1]+1]);size=stop-start
  if dest!=start:
   for offset in range(0,size,1_000_000):
    n=min(1_000_000,size-offset);x.data[dest+offset:dest+offset+n]=x.data[start+offset:start+offset+n];x.indices[dest+offset:dest+offset+n]=x.indices[start+offset:start+offset+n]
  dest+=size
 assert dest==ptr[-1]
 compact=sparse.csr_matrix((x.data[:dest],x.indices[:dest],ptr),shape=(len(keep),a.n_vars),copy=False);real=ad.AnnData(compact,obs=a.obs.iloc[keep].copy(),var=a.var.copy());return ctrl,real
