"""Independent exact rational verifier for the saved synthetic count artifact."""
from fractions import Fraction as F
import numpy as np,json,pathlib
p=pathlib.Path('artifacts');d=np.load(p/'atlas_fixture.npz');a=d['counts'];depth=d['depths'];n,g=a.shape
assert a.dtype.kind in 'iu' and a.min()>=0
assert all(sum(map(int,row))==int(dep) for row,dep in zip(a,depth))
assert np.array_equal(a,d['repeat'])
# Floats are interpreted as their exact stored binary rationals, never rounded decimals.
cp=[sum((F(int(a[i,j]),int(depth[i])) for i in range(n)),F(0))/n for j in range(g)]
total=sum(map(int,depth));bp=[F(sum(map(int,a[:,j])),total) for j in range(g)]
cperr=sum((abs(x-F(float(y))) for x,y in zip(cp,d['p'])),F(0))
berr=sum((abs(x-F(float(y))) for x,y in zip(bp,d['bulk'])),F(0))
assert cperr<F(1,1000) and berr<F(1,1000)
# Score arithmetic using exact decimal examples; no inaccessible biological truth.
b,r=F(99,100),F(4,100)
u=b+F(1,2)*(r-b);assert u==F(103,200)
# In the noiseless one-direction model min_a ||a*x-y||^2/||y||^2 = 1-cos^2.
# Baseline-like raw .99 to scaled .5 requires raw .515 => squared cosine >= .485.
required_cos2=1-u
result={'verified':'exact integer rows and deterministic equality; rational L1 moment bounds < 1/1000','shape':[n,g],'mean_cpm_L1_float_display':float(cperr),'bulk_L1_float_display':float(berr),'example_raw_for_scaled_mse_half':str(u),'example_required_cos_squared':str(required_cos2),'biological_improvement':'UNVERIFIED; synthetic fixture only'}
(p/'atlas_fixture_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
