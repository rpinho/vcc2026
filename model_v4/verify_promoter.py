"""Independent rational functional check of a real source-guide-supported short-distance pair."""
import runtime
import numpy as np,pandas as pd
from fractions import Fraction as F
from model import promoters
pairs=pd.read_csv(runtime.ROOT/'artifacts/promoter_pairs.csv');row=pairs[pairs.guide_supported&(pairs.guide_distance<=500)&(pairs.guide_count>=2)].iloc[0]
genes=np.array([row.target,row.neighbor,'SYNTHETIC_UNLINKED_MARKER']);targets=np.array([row.target]);q=np.array([[.25,.25,.5]]);control=q[0];result=promoters(q,control,targets,genes,np.zeros(q.shape),True)
# At this measured distance the neighbour retains 15% before final renormalization.
weights=[F(1,4),F(1,4)*F(.15),F(1,2)];z=sum(weights);expected=[w/z for w in weights]
errors=[abs(F(float(x))-e) for x,e in zip(result[0],expected)];assert max(errors)<F(1,2**50)
measured=promoters(q,control,targets,genes,np.full(q.shape,2),True);assert np.array_equal(measured,q)
runtime.save(runtime.ROOT/'artifacts/promoter_verification.json',dict(status='VERIFIED rational functional invariant with real public guide-supported pair; biological destination effect UNVERIFIED',target=row.target,neighbor=row.neighbor,source_guides=int(row.guide_count),guide_distance=float(row.guide_distance),renormalized_expected=[str(e) for e in expected],exact_binary64_error_upper=str(max(errors)),two_family_measured_coordinates_unchanged=True))
print('Guide-prior action and abstention VERIFIED',flush=True)
