"""Compare pre-calibration projection with released AtlasShift's clipped scalar formula.
No model heads, scorer or hidden labels; separate scalar bisection and feasibility checks.
"""
import runtime,json
import numpy as np
from feasible_bulk import project_bulk
rng=np.random.default_rng(94);P=rng.lognormal(0,1,(12,256));P[:,0]=0;P/=P.sum(1,keepdims=True);B=rng.lognormal(0,1.5,P.shape);B/=B.sum(1,keepdims=True);zmin=.35;zmax=2.8
new,errors=project_bulk(B,P,zmin,zmax);reference=[]
for p,b in zip(P,B):
 lower=(zmin+.01*(1-zmin))*p;upper=(zmax-.01*(zmax-1))*p;lo=0.;hi=1.
 while np.clip(b*hi,lower,upper).sum()<1:hi*=2
 for _ in range(90):
  mid=(lo+hi)/2
  if np.clip(b*mid,lower,upper).sum()<1:lo=mid
  else:hi=mid
 reference.append(np.clip(b*((lo+hi)/2),lower,upper))
reference=np.array(reference);maxdiff=float(np.max(np.abs(new-reference)));assert maxdiff<1e-11
assert np.allclose(new.sum(1),1,rtol=0,atol=1e-12);assert np.all(new>=(zmin+.01*(1-zmin))*P-1e-13);assert np.all(new<=(zmax-.01*(zmax-1))*P+1e-13)
assert np.allclose(errors,np.abs(new-B).sum(1),rtol=0,atol=1e-12)
runtime.save(runtime.ROOT/'artifacts/projection_verification.json',dict(status='PASS released scalar projection equality and depth-feasibility constraints',max_abs_difference=maxdiff,max_mass_error=float(np.max(np.abs(new.sum(1)-1))),shape=list(P.shape)));print('PASS',maxdiff)
