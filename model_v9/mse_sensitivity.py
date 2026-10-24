"""Official capped numerator sensitivity to unknown real per-target caps.

Pooled predictions are fixed, so expression distances cancel. The hypothetical
real means below are identical across candidates; only the real JK cap varies.
These are mathematical sensitivity checks, never hidden score estimates.
"""
import runtime
import numpy as np
from cell_eval2.moments import GroupMoments
from cell_eval2.metrics.delta import mse_unbiased_capped,PRED_TRACE_CAP_K
R=runtime.ROOT
def main():
 z=np.load(R/'artifacts/stage1_v5b.npz');genes=z['genes'].astype(str);targets=z['targets'].astype(str);totals=z['totals']
 reports={arm:runtime.json.loads((R/f'artifacts/{arm}_gates.json').read_text()) for arm in ['v9a','v9b']};output=[]
 for ci,c in enumerate('ABC'):
  L=np.log1p(50000*totals[ci]/totals[ci].sum(1,keepdims=True));ctrl=np.log1p(.05*np.load(R/f'artifacts/ref_{c}/bulk.npy'))
  baseline=np.array([r['v5b_jk'] for r in reports['v9a']['records'] if r['context']==c]);pred_bulk=(targets,L);rlabels=np.r_[targets,'non-targeting'];real_bulk=(rlabels,np.vstack([np.broadcast_to(ctrl,L.shape),ctrl]))
  def metric(jk,cap):
   pm=GroupMoments(targets,np.full(300,400.),np.zeros(300),jk=jk)
   rm=GroupMoments(rlabels,np.r_[np.full(300,400.),18400.],np.zeros(301),jk=np.r_[cap/PRED_TRACE_CAP_K,0.])
   return np.array(list(mse_unbiased_capped(pred_bulk=pred_bulk,real_bulk=real_bulk,pred_moments=pm,real_moments=rm,control='non-targeting',comparator='bulk_lognorm',genes=genes).values()))
  for family in ['constant','proportional_to_baseline','reverse_baseline']:
   for factor in [.25,.5,.75,1.,1.25,1.5,2.,4.]:
    if family=='constant':cap=np.full(300,baseline.mean()*factor)
    elif family=='proportional_to_baseline':cap=baseline*factor
    else:cap=baseline[::-1]*factor
    old=metric(baseline,cap)
    for arm in reports:
     jk=np.array([r['jk'] for r in reports[arm]['records'] if r['context']==c]);new=metric(jk,cap);diff=new-old
     output.append(dict(context=c,arm=arm,cap_family=family,cap_factor=factor,numerator_delta_sum=float(diff.sum()),numerator_delta_max_abs=float(abs(diff).max()),status='exact official capped numerator on hypothetical real-cap panel; distance terms cancel'))
 runtime.save(R/'artifacts/mse_sensitivity.json',output)
 print('MSE SENSITIVITY COMPLETE',len(output),flush=True)
if __name__=='__main__':main()
