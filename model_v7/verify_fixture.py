"""Independent numerical certificate: no model/emitter/scorer imports.
Integer reciprocal interval bounds verify CPM and bulk moments against exact stored floats.
"""
import pathlib,json,hashlib
import numpy as np
from verify_counts import moment_certificate
R=pathlib.Path(__file__).resolve().parent
records=[]
for fn in sorted((R/'artifacts').glob('fixture_*.npz')):
 with np.load(fn) as z:
  x=z['x'];d=z['depths'];p=z['p'];b=z['b'];mask=z['mask']
  assert x.shape==(400,18533);assert x.dtype.kind in 'iu';assert np.all(x>=0);assert np.array_equal(x.sum(1),d)
  assert len(mask)==18533
  c=moment_certificate(x,d,p,b);assert c['mean_L1_upper']<.015 and c['bulk_L1_upper']<.0015
  # Emission is nondegenerate on both sides of the variance mask; no biology claimed.
  v=np.var(x/d[:,None],axis=0);c.update(file=fn.name,sha256=hashlib.sha256(fn.read_bytes()).hexdigest(),variable_genes=int((v>0).sum()),low_variance_mask=int(mask.sum()))
  records.append(c)
assert len(records)==8,f'need 4 fixtures per arm, found {len(records)}'
out=dict(status='PASS exact integer and rational-interval moments, not a biological validation',records=records)
(R/'artifacts/fixture_verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=1))
