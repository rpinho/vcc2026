"""Independent integer/rational raw fidelity and reach checker on prepared tables.

No direction metric implementation or cached component frame is imported here.
Preparation is shared with the enrolled path; ranking and counting are separate.
"""
import runtime,math
from fractions import Fraction as F
import polars as pl
from scorer_check import score
R=runtime.ROOT
FID='de_wilcoxon_direction_fidelity_yield_raw'
REACH='de_wilcoxon_direction_reach_raw'
def defined(x):return x is not None and not math.isnan(x) and x!=0
def rank(x):return float('inf') if x is None or math.isnan(x) else x
def independent(prepared):
 real={};pred={};alpha=float(prepared.p_adj_threshold)
 for row in prepared.real_df.to_dicts():
  if row['feature']!=row['target']:real.setdefault(row['target'],{})[row['feature']]=row
 for row in prepared.pred_df.to_dicts():
  if row['feature']!=row['target']:pred.setdefault(row['target'],{})[row['feature']]=row
 out={}
 for target in prepared.perturbations:
  a=pred.get(target,{});b=real.get(target,{})
  nconf=sum(rank(v['p_adj'])<alpha for v in b.values());n=0;k=0;pool=[]
  for gene,v in a.items():
   if gene not in b:continue
   w=b[gene];rl=w['log2_fold_change'];vl=v['log2_fold_change']
   can=defined(rl);match=can and defined(vl) and ((vl>0)==(rl>0))
   sig=rank(v['p_adj'])<alpha
   if sig and can:n+=1;k+=int(match)
   if rank(w['p_adj'])<alpha and can:
    mag=abs(vl) if vl is not None and not math.isnan(vl) else -float('inf')
    pool.append(((not sig,rank(v['p_adj']),rank(v.get('p_value')),-mag,gene),match))
  pool.sort(key=lambda v:v[0]);hit=0;deepest=0
  for depth,(_,match) in enumerate(pool,1):
   hit+=int(match)
   if 10*hit>=9*depth:deepest=depth
  # These raw catalog members retain undefined targets as NaN. Do not replace
  # an absent reference decision with a fabricated zero.
  out[target]={FID:F(k,max(n,nconf)) if max(n,nconf) else None,REACH:F(deepest,nconf) if nconf else None}
 return out
receipts=[]
for source in ['Jurkat','HepG2','RPE1']:
 for condition in ['native','projected400']:
  suffix='' if condition=='native' else '_projected400';real=pl.read_parquet(R/f'artifacts/public_{source}{suffix}.parquet')
  for arm in ['v5b','v9a','v9b']:
   runtime.guard();pred=pl.read_parquet(R/f'artifacts/public_pred_{arm}_{source}.parquet');rows,prepared=score(pred,real,return_prepared=True);exact=independent(prepared);checked=0
   for row in rows:
    if row['metric'] not in [FID,REACH]:continue
    target=row.get('target',row.get('perturbation'));expected=exact[target][row['metric']];value=row['value']
    if expected is None:assert value is None or math.isnan(value)
    else:assert abs(float(expected)-value)<1e-14,(source,arm,target,row,expected)
    checked+=1
   rational_means={}
   for metric in [FID,REACH]:
    finite=[v[metric] for v in exact.values() if v[metric] is not None];mean=sum(finite,F(0))/len(finite) if finite else None
    rational_means[metric]=None if mean is None else dict(numerator=str(mean.numerator),denominator=str(mean.denominator),float=float(mean),finite_targets=len(finite))
   targets={t:{m:None if v is None else [str(v.numerator),str(v.denominator)] for m,v in metrics.items()} for t,metrics in exact.items()}
   receipts.append(dict(source=source,condition=condition,candidate=arm,checked=checked,means=rational_means,per_target_rationals=targets,status='PASS independent integer ranking/counts and rational metrics agree with official per-target values'))
runtime.save(R/'artifacts/public_rational_verification.json',receipts)
print('PASS independent public raw fidelity/reach verification',flush=True)
