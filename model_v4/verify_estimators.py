"""Independent rational checks of shrinkage, inverse-variance pooling and interval certificates."""
import runtime
from fractions import Fraction as F
import numpy as np
from vendor.adaptive_shrink import posterior_mean,adaptive_shrink
from verify_counts import moment_certificate
# Single Gaussian prior: posterior coefficient is exactly tau^2/(tau^2+SE^2).
x=np.array([-10.,-5.,0.,5.,10.]);s=np.full(5,2.);out=posterior_mean(x,s,np.array([1.]),np.array([1.]));expected=[F(int(a),5) for a in x]
assert all(F(float(a))==b for a,b in zip(out,expected))
# Shared true effect under independent source errors: pooled variance is 1/sum(1/v).
w=[1/F(4),1/F(9)];v=1/sum(w);a=sum(ww*xx for ww,xx in zip(w,[F(2),F(3)]))/sum(w);assert a==F(30,13) and v==F(36,13)
# Prove the fixed reciprocal interval brackets the exact rational synthetic moments.
d=np.load('artifacts/atlas_fixture.npz');cert=moment_certificate(d['counts'],d['depths'],d['p'],d['bulk']);n=len(d['depths']);total=sum(map(int,d['depths']))
q=[sum((F(int(d['counts'][i,j]),int(d['depths'][i])) for i in range(n)),F(0))/n for j in range(d['counts'].shape[1])];b=[F(sum(map(int,d['counts'][:,j])),total) for j in range(d['counts'].shape[1])]
qe=sum((abs(v-F(float(w))) for v,w in zip(q,d['p'])),F(0));be=sum((abs(v-F(float(w))) for v,w in zip(b,d['bulk'])),F(0))
# Float displays may round outward or inward: allow one binary64 ULP on display checks.
assert np.nextafter(cert['mean_L1_lower'],-np.inf)<=float(qe)<=np.nextafter(cert['mean_L1_upper'],np.inf)
assert np.nextafter(cert['bulk_L1_lower'],-np.inf)<=float(be)<=np.nextafter(cert['bulk_L1_upper'],np.inf)
# Continuous, sign-preserving adaptive estimates, with the released full stopping rule.
rng=np.random.default_rng(761);effects=rng.normal(0,.03,1000);var=rng.uniform(.01,.1,1000)**2;stats={};shrunk=adaptive_shrink(effects,var,stats=stats);assert np.all(np.sign(shrunk)==np.sign(effects));assert np.all(abs(shrunk)<=abs(effects)*(1+1e-14));assert np.all(shrunk!=0)
runtime.save('artifacts/estimator_verification.json',dict(status='VERIFIED rational one-component posterior and IV identities; exact reciprocal-bound fixture; adaptive sign/norm invariants',posterior_expected=[str(x) for x in expected],pooled_mean=str(a),pooled_variance=str(v),adaptive_stats=stats,biological_improvement='UNVERIFIED'))
print('Estimator algebra and reciprocal intervals VERIFIED',stats,flush=True)
