"""Memory-bounded fused version of AtlasShift's alternating dual-moment fit.
Same update/constraints/stopping tolerances; IEEE fp64, no fastmath or parallel workers.
The original MIT code remains in vendor/atlas_model.py and competitors/code/AtlasShift.
"""
import numpy as np
from numba import njit
@njit(cache=True)
def fit_dual(x,p,bulk,z,iterations,tolerance):
 n,g=x.shape;desired=n*p;ratio=np.ones(g)
 for j in range(g):
  if p[j]>0:ratio[j]=bulk[j]/p[j]
 zmin=np.min(z);zmax=np.max(z);cpm_error=0.;bulk_error=0.
 col=np.empty(g);first=np.empty(g);second=np.empty(g);fac=np.empty(g);tilt=np.empty(g);means=np.empty(g)
 for step in range(iterations):
  col[:]=0.;first[:]=0.;second[:]=0.
  for i in range(n):
   zi=z[i]
   for j in range(g):
    xx=x[i,j];col[j]+=xx;first[j]+=zi*xx;second[j]+=zi*zi*xx
  for j in range(g):
   fac[j]=desired[j]/col[j] if col[j]>0 else 0.
   mean=first[j]/col[j] if desired[j]>0 and col[j]>0 else 1.
   var=max(second[j]/col[j]-mean*mean,0.) if desired[j]>0 and col[j]>0 else 0.
   v=(ratio[j]-mean)/var if var>1e-12 else 0.
   v=max(-.95/max(zmax-mean,1e-12),min(.95/max(mean-zmin,1e-12),v));means[j]=mean;tilt[j]=v
  for i in range(n):
   s=0.;zi=z[i]
   for j in range(g):
    xx=x[i,j]*fac[j]*(1.+tilt[j]*(zi-means[j]));x[i,j]=xx;s+=xx
   for j in range(g):x[i,j]/=max(s,1e-30)
  if step%5==4 or step==iterations-1:
   col[:]=0.;first[:]=0.
   for i in range(n):
    for j in range(g):col[j]+=x[i,j];first[j]+=z[i]*x[i,j]
   cpm_error=0.;bulk_error=0.
   for j in range(g):cpm_error+=abs(col[j]/n-p[j]);bulk_error+=abs(first[j]/n-bulk[j])
   if max(cpm_error,bulk_error)<tolerance:break
 return x,step+1,cpm_error,bulk_error
