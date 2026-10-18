"""Re-run the pinned released 400x40 numerical baseline before improvement."""
import runtime as rt
from fractions import Fraction as F
import numpy as np
from atlas_model import dual_moment_counts
assert rt.sha('vendor/atlas_model.py')=='8ad68112dbfa9c5c8c8b21a717dbe475fe23613f84c25e844fbfbf8bb47cff3b'
rng=np.random.default_rng(7331);n,g=400,40;depths=np.rint(np.linspace(11000,31000,n)).astype('int64');p0=rng.dirichlet(np.full(g,3.));z=depths/depths.mean();comp=p0[None,:]*np.exp((z[:,None]-1)*np.linspace(-.5,.5,g)[None,:]);comp/=comp.sum(1,keepdims=True);p=comp.mean(0);bulk=depths@comp/depths.sum();template=rng.gamma(5,1,(n,g))*p0[None,:];a=dual_moment_counts(template,p,bulk,depths=depths,seed=20260910);b=dual_moment_counts(template,p,bulk,depths=depths,seed=20260910)
assert a.dtype.kind in 'iu' and np.all(a>=0) and np.array_equal(a,b)
assert all(sum(map(int,row))==int(d) for row,d in zip(a,depths))
pe=sum((abs(sum((F(int(a[i,j]),int(depths[i])) for i in range(n)),F(0))/n-F(float(p[j]))) for j in range(g)),F(0));be=sum((abs(F(sum(map(int,a[:,j])),sum(map(int,depths)))-F(float(bulk[j]))) for j in range(g)),F(0));assert pe<F(1,1000) and be<F(1,1000)
rt.save('artifacts/atlas_reproduction.json',dict(status='PASS released numerical fixture; not full source/hidden-context reproduction',source_sha256=rt.sha('vendor/atlas_model.py'),shape=[n,g],exact_integer_depths=True,deterministic=True,mean_L1_exact=str(pe),bulk_L1_exact=str(be)))
print('PASS released AtlasShift numerical fixture; no official score inherited')
