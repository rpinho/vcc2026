"""Exact tied two-sample normal MWU; controls are sorted once, no labels."""
import runtime
import numpy as np,numba as nb,math
from scipy import sparse
from scipy.special import ndtr
from scipy.stats import false_discovery_control

@nb.njit(cache=True)
def search(a,value,right=False):
 lo=0;hi=len(a)
 while lo<hi:
  mid=(lo+hi)//2
  if a[mid]<value or (right and a[mid]==value):lo=mid+1
  else:hi=mid
 return lo

@nb.njit(cache=True)
def uz(y,ref,n,tie0):
 m=len(y);U=0.;tie=tie0;nz=n-len(ref)
 order=np.sort(y);i=0
 while i<m:
  j=i+1
  while j<m and order[j]==order[i]:j+=1
  v=order[i];k=j-i
  if v==0:less=0;equal=nz
  else:
   l=search(ref,v);r=search(ref,v,True);less=nz+l;equal=r-l
  U+=k*(less+.5*equal)
  tie+=(equal+k)**3-(equal+k)-(equal**3-equal)
  i=j
 var=m*n/12*((m+n+1)-tie/((m+n)*(m+n-1)))
 z=0. if var<=0 else (U-.5*m*n)/math.sqrt(var)
 return U,z,0. if var<=0 else math.sqrt(var)

@nb.njit(cache=True)
def _stats_log(logx,ref,ptr,n,ties,keep):
 g=logx.shape[1];z=np.zeros(g);pv=np.ones(g)
 for j in range(g):
  if not keep[j]:continue
  y=logx[:,j]
  U,zz,s=uz(y,ref[ptr[j]:ptr[j+1]],n,ties[j]);z[j]=zz
  if s>0:pv[j]=math.erfc(max(0.,abs(zz)-.5/s)/math.sqrt(2.))
 return z,pv

def stats(x,depth,ref,ptr,n,ties,keep=None):
 # NumPy log1p, exactly as the scorer's preprocessing, rather than a JIT libm
 # variant which could break a count/depth tie by one ulp.
 logx=np.asfortranarray(np.log1p(np.asarray(x,dtype=np.float64)/(np.asarray(depth,dtype=float)[:,None]/1e6)))
 if keep is None:keep=np.ones(logx.shape[1],bool)
 return _stats_log(logx,ref,ptr,n,ties,keep)

def adjusted(p,keep):
 q=np.ones(len(p));q[keep]=false_discovery_control(p[keep]);return q

class Reference:
 def __init__(self,c):
  d=runtime.ROOT/'artifacts'/f'ref_{c}'
  self.values=np.load(d/'values.npy',mmap_mode='r');self.ptr=np.load(d/'ptr.npy');self.ties=np.load(d/'ties.npy');self.mean=np.load(d/'mean.npy');self.base=np.load(d/'base.npy');self.depth=np.load(d/'depth.npy');self.n=int(np.load(d/'n.npy'));self.keep=self.mean>5
  self.bulk=np.load(d/'bulk.npy');self.base=np.asfortranarray(self.base)
  self.cdf=np.load(d/'cdf.npy',mmap_mode='r');self.cdf_scale=np.load(d/'cdf_scale.npy')
 def test(self,x):
  d=x.sum(1,dtype=np.int64);z,p=stats(np.asarray(x),d,self.values,self.ptr,self.n,self.ties,self.keep)
  return z,p,adjusted(p,self.keep)

def prepare(x,out,rng,base_cells=400):
 """Match scorer counts -> float64 scanpy normalize_total/log1p; retain real-cell ranks."""
 out.mkdir(parents=True,exist_ok=True);x=sparse.csr_matrix(x);x.eliminate_zeros();dep=np.asarray(x.sum(1,dtype=np.float64)).ravel();assert (dep>0).all()
 mean=np.zeros(x.shape[1])
 bulk=np.asarray(x.sum(0,dtype=np.float64)).ravel();bulk=bulk/bulk.sum()*1e6
 rows=np.sort(rng.choice(x.shape[0],base_cells,replace=False));base=x[rows].toarray().astype(float)/dep[rows,None]*1e6
 x=x.tocsc();x.sort_indices()
 ties=np.zeros(x.shape[1]);vals=np.lib.format.open_memmap(out/'values.npy',mode='w+',dtype='float64',shape=(x.nnz,));cdf=np.lib.format.open_memmap(out/'cdf.npy',mode='w+',dtype='float32',shape=(x.shape[1],2049));scales=np.ones(x.shape[1])
 for j in range(x.shape[1]):
  lo,hi=x.indptr[j:j+2];linear=x.data[lo:hi].astype(np.float64)/(dep[x.indices[lo:hi]]/1e6);mean[j]=linear.sum()/x.shape[0];v=np.sort(np.log1p(linear));vals[lo:hi]=v;nz=x.shape[0]-len(v);_,ct=np.unique(v,return_counts=True);ties[j]=nz**3-nz+np.sum(ct.astype(float)**3-ct)
  scales[j]=max(float(v[-1]) if len(v) else 0.,1.);grid=np.linspace(0,scales[j],2049);left=np.searchsorted(v,grid);right=np.searchsorted(v,grid,side='right');cdf[j]=nz+left+.5*(right-left);cdf[j,0]=nz
 for name,value in [('ptr',x.indptr),('ties',ties),('mean',mean),('bulk',bulk),('base',base),('depth',dep[rows]),('n',np.array(x.shape[0]))]:np.save(out/f'{name}.npy',value)
 np.save(out/'cdf_scale.npy',scales);vals.flush();cdf.flush();del vals,cdf
