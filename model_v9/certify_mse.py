"""Separate rational checker for the sufficient aggregate MSE condition."""
import runtime
from fractions import Fraction as F
R=runtime.ROOT
output={}
for name in ['v9a','v9b']:
 gates=runtime.json.loads((R/f'artifacts/{name}_gates.json').read_text());rows=gates['records'];assert len(rows)==900
 contexts=[]
 for c in 'ABC':
  part=[r for r in rows if r['context']==c];assert len(part)==300
  differences=[F(float(r['jk']))-F(float(r['v5b_jk'])) for r in part]
  failed=[r['target'] for r,diff in zip(part,differences) if diff<0]
  assert all(r['totals_exact'] and r['depths_exact'] and r['nonDE_exact'] for r in part)
  smallest=min(differences)
  contexts.append(dict(context=c,componentwise_JK_dominance=not failed,failed_targets=failed,minimum_JK_gap_rational=[smallest.numerator,smallest.denominator],minimum_JK_gap_float=float(smallest)))
 output[name]=dict(status='PASS sufficient MSE preservation condition on exact recorded official-kernel floats' if all(c['componentwise_JK_dominance'] for c in contexts) else 'FAIL sufficient condition; do not assert guaranteed MSE preservation',contexts=contexts,proof='For every real cap R>=0, C_new>=C_old implies min(C_new,R)>=min(C_old,R). Positive weighted sum and min(B,.) preserve the inequality. Fixed profiles fix distances and B. The enrolled ratio-of-sums denominator is submission-independent. This proves mathematical aggregate monotonicity; floating reduction limits remain explicit.',numeric_qualification='Actual default 512-row vendored jackknife, independently checked by gates_core; exact Fraction comparisons of stored binary64 values. No claim of interval verification of platform log1p or GPU arithmetic.')
runtime.save(R/'artifacts/mse_certificate.json',output)
print({k:v['status'] for k,v in output.items()},flush=True)
