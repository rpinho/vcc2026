"""AtlasShift depth-feasible bulk projection, before its budget is calibrated.
Same clipped scalar projection as MIT dual_moment_counts; safeguarded Newton root.
"""
import numpy as np
from numba import njit
@njit(cache=True)
def project_bulk(B,P,zmin,zmax):
 t,g=B.shape;out=np.empty_like(B);errors=np.zeros(t)
 a=zmin+.01*(1-zmin);b=zmax-.01*(zmax-1)
 for i in range(t):
  lo=0.;hi=1.
  for k in range(1024):
   mass=0.
   for j in range(g):mass+=min(max(B[i,j]*hi,a*P[i,j]),b*P[i,j])
   if mass>=1.-1e-14:break
   hi*=2.
  assert np.isfinite(hi) and k<1023
  s=min(1.,hi)
  for it in range(90):
   mass=0.;derivative=0.
   for j in range(g):
    v=B[i,j]*s;l=a*P[i,j];u=b*P[i,j];mass+=min(max(v,l),u)
    if v>l and v<u:derivative+=B[i,j]
   if abs(mass-1.)<1e-13:break
   if mass<1.:lo=s
   else:hi=s
   step=s+(1.-mass)/derivative if derivative>0 else .5*(lo+hi)
   s=step if lo<step<hi else .5*(lo+hi)
  err=0.
  for j in range(g):out[i,j]=min(max(B[i,j]*s,a*P[i,j]),b*P[i,j]);err+=abs(out[i,j]-B[i,j])
  errors[i]=err
 return out,errors
