"""Independent exact integer and rational-interval certificates.
This module does not import the predictor, emitter, or evaluator.
"""
from fractions import Fraction as F
import numpy as np,hashlib,h5py,pandas as pd,collections
S=2**48

def target_floor(values,D):
 ans=np.empty(len(values),np.int64)
 for i,v in enumerate(values):
  n,d=float(v).as_integer_ratio();ans[i]=(n*D)//d
 return ans

def moment_certificate(counts,depths,p,b):
 assert counts.dtype.kind in 'iu' and np.all(counts>=0)
 a=counts.astype(np.int64);depths=np.asarray(depths,np.int64);n,g=a.shape
 assert np.array_equal(a.sum(1),depths) and np.all((depths>0)&(depths<=1_000_000))
 columns=a.sum(0);total=int(columns.sum());reciprocal=S//depths
 # floor(S/depth) <= S/depth <= floor(S/depth)+1, entirely integer arithmetic.
 lo=(a*reciprocal[:,None]).sum(0);hi=lo+columns;D=n*S
 fp=target_floor(p,D)
 upper=int(np.maximum(abs(lo-(fp+1)),abs(hi-fp)).sum())
 lower=int(np.maximum(np.maximum(lo-(fp+1),fp-hi),0).sum())
 fb=target_floor(b,S);blo=columns*(S//total);bhi=blo+columns
 bu=int(np.maximum(abs(blo-(fb+1)),abs(bhi-fb)).sum())
 bl=int(np.maximum(np.maximum(blo-(fb+1),fb-bhi),0).sum())
 return dict(mean_L1_lower=float(F(lower,D)),mean_L1_upper=float(F(upper,D)),bulk_L1_lower=float(F(bl,S)),bulk_L1_upper=float(F(bu,S)),row_depths_exact=True,certificate='Integer reciprocal bounds against exact binary64 requested moments; bulk includes feasibility projection/rounding')

def column(obj):
 if isinstance(obj,h5py.Dataset):return obj.asstr()[:]
 if 'values' in obj:
  assert not np.any(obj['mask'][:]);return column(obj['values'])
 cats=column(obj['categories']);codes=obj['codes'][:];assert np.all(codes>=0);return cats[codes]

def audit_h5(path,genes,targets,contexts=tuple('ABC')):
 hashes={};result={}
 with h5py.File(path,'r') as h:
  x=h['X'];shape=tuple(map(int,x.attrs['shape']));assert shape==(400*len(targets)*len(contexts),len(genes));ptr=x['indptr'][:].astype(np.int64)
  assert ptr[0]==0 and np.all(np.diff(ptr)>=0);nnz=int(ptr[-1]);assert nnz==len(x['data'])==len(x['indices']) and nnz<=4750000000
  vg=h['var'];axis=column(vg[vg.attrs['_index']]);assert list(axis)==list(genes)
  cc=column(h['obs/context']);tt=column(h['obs/target_gene']);groups=collections.Counter(zip(cc,tt));assert set(groups)=={(c,t) for c in contexts for t in targets};assert set(groups.values())=={400}
  maxdepth=0;maxcount=0;total=0
  for lo in range(0,shape[0],128):
   hi=min(lo+128,shape[0]);a,b=int(ptr[lo]),int(ptr[hi]);v=x['data'][a:b];ix=x['indices'][a:b];assert np.isfinite(v).all() and np.all(v>0) and np.all(v==np.floor(v)) and np.all((ix>=0)&(ix<shape[1]))
   ps=np.r_[0,np.cumsum(v,dtype=np.int64)];ls=ps[(ptr[lo+1:hi+1]-a)]-ps[(ptr[lo:hi]-a)];assert np.all((ls>0)&(ls<=1000000));maxdepth=max(maxdepth,int(ls.max()));maxcount=max(maxcount,int(v.max()));total+=int(v.sum(dtype=np.int64))
   for row in range(lo,hi):assert np.all(np.diff(ix[int(ptr[row])-a:int(ptr[row+1])-a])>0)
  for name,dtype in [('data','<u4'),('indices','<u4'),('indptr','<u8')]:
   dig=hashlib.sha256();ds=x[name]
   for lo in range(0,len(ds),2**22):dig.update(ds[lo:lo+2**22].astype(dtype).tobytes())
   hashes[name]=dig.hexdigest()
  for name,values in [('genes',axis),('context',cc),('target_gene',tt)]:hashes[name]=hashlib.sha256('\0'.join(values).encode()).hexdigest()
  result=dict(status='VERIFIED exact integer/structure invariants',shape=list(shape),nnz=nnz,max_cell_depth=maxdepth,max_count=maxcount,total_counts=total,canonical_hashes=hashes)
 return result

def scale_record(values,baseline,replicate,metric_names):
 """Exact rational arithmetic on stored floating-point metric aggregates only."""
 result={}
 for m in metric_names:
  x,b,r=map(lambda x:F(float(x)),[values[m],baseline[m],replicate[m]]);score=(x-b)/(r-b)
  if m=='expr_mse_unbiased_capped_norm':score=min(F(1),max(F(0),score))
  if m=='de_wilcoxon_lfc_nmae':score=max(F(-6),score)
  result[m]=score
 overall=sum(result.values(),F(0))/len(result)
 return result,overall
