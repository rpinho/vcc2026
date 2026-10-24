"""Independent integer/interval certificate for unavoidable low-count DE calls.

With T total integer transcripts, at most min(T,400) cells can be nonzero.
Each positive cell contributes at most n_control to MWU U; every zero cell
contributes n_control_zero/2. This upper bound applies to EVERY within-group
layout, even if row depths change. No layout/search implementation is imported.
"""
import runtime,sys
from fractions import Fraction as F
import numpy as np,mpmath as mp
R=runtime.ROOT
mp.iv.dps=60
z=np.load(R/'artifacts/stage1_v5b.npz');totals=z['totals'];genes=z['genes'].astype(str)
progress='--progress' in sys.argv
rows=runtime.json.loads((R/'artifacts/build_progress.json').read_text())['records']['v9b'] if progress else runtime.json.loads((R/'artifacts/v9b_gates.json').read_text())['records']
if not progress:assert len(rows)==900
records=[]
for ci,c in enumerate('ABC'):
 n=int(np.load(R/f'artifacts/ref_{c}/n.npy'));ptr=np.load(R/f'artifacts/ref_{c}/ptr.npy');mtest=int((np.load(R/f'artifacts/ref_{c}/mean.npy')>5).sum())
 votes=np.load(R/f'artifacts/recovery_{c}.npz');positive=votes['positive_votes'];negative=votes['negative_votes']
 for row in [v for v in rows if v['context']==c and v['actual_calls_without_two_agreeing_unopposed_sources']]:
  runtime.guard();i=row['index'];values=np.load(R/f'artifacts/de/v9b_{c}_{i}.npz');lfc=values['lfc'];agree=np.where(lfc>0,positive[i],negative[i]);oppose=np.where(lfc>0,negative[i],positive[i]);supported=np.isfinite(lfc)&(lfc!=0)&(agree>=2)&(oppose==0)
  bad=(values['q']<.05)&~supported;bad[genes==row['target']]=False
  assert int(bad.sum())==row['actual_calls_without_two_agreeing_unopposed_sources']
  for j in np.flatnonzero(bad):
   T=int(totals[ci,i,j]);nonzero=min(T,400);nz=n-int(ptr[j+1]-ptr[j]);umax=F(2*nonzero*n+(400-nonzero)*nz,2);delta=F(400*n,2)-umax-F(1,2)
   rec=dict(context=c,index=i,target=row['target'],gene=genes[j],integer_gene_total=T,control_zero_cells=nz,control_cells=n,maximum_positive_cells=nonzero,U_upper_rational=[umax.numerator,umax.denominator],actual_q=float(values['q'][j]),proven_unavoidable=False)
   if delta>0:
    z2=delta*delta/F(400*n*(400+n+1),12)
    # Tie correction can only reduce variance. For x>=0, erfc(x)<=exp(-x*x).
    # BH adjusted p is at most the number of tested genes times its raw p.
    upper=mtest*mp.iv.exp(-(mp.iv.mpf(z2.numerator)/z2.denominator)/2)
    rec.update(minimum_z_squared_rational=[z2.numerator,z2.denominator],BH_q_upper_interval=str(upper),proven_unavoidable=bool(upper.b<mp.iv.mpf(1)/20))
    if rec['proven_unavoidable']:assert mp.iv.mpf(float(values['q'][j]))<=upper.b
   records.append(rec)
runtime.save(R/('artifacts/forced_call_probe.json' if progress else 'artifacts/forced_call_certificate.json'),dict(status='PARTIAL probe' if progress else 'Independent exact integer support bound with interval exponential tail enclosure',groups_checked=len(rows),unsupported_calls=len(records),proven_unavoidable=sum(r['proven_unavoidable'] for r in records),records=records,scope='A proven call cannot be removed by ANY integer within-group layout with this fixed pooled gene total and 400 cells, under the enrolled asymptotic MWU/BH producer. Changing other genes or row depths does not evade the support bound. Unproven rows only mean this conservative bound is inconclusive.'))
print('FORCED CALL BOUND',len(records),sum(r['proven_unavoidable'] for r in records),flush=True)
