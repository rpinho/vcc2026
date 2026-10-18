"""Independent integer verifier for the 300-target availability ceilings.
Uses row counters, not the search's groupby/minimum/swap implementation.
"""
import runtime as rt
from pathlib import Path
from collections import Counter
import numpy as np,pandas as pd
universe=set(rt.json.loads(Path('artifacts/sources/target_universe.json').read_text()));available=set(np.load(rt.ROOT.parent/'public_data/model_v4_signature_banks/HCT116.npz')['targets'].astype(str))&universe;results={}
for context,file,target,batch in [('GWPS','GWPS_raw_obs.parquet','gene','gem_group'),('CD4','CD4_cell_census.parquet','target_gene','batch')]:
 frame=pd.read_parquet('artifacts/'+file,columns=[target,batch]);counter=Counter(zip(frame[target].astype(str),frame[batch].astype(str)));targets=sorted(available&set(frame[target].astype(str))-{'non-targeting'});batches=sorted(set(frame[batch].astype(str)));limits=[];raw_limits=[]
 for b in batches:
  raw_limits.append(sorted((counter[(t,b)] for t in targets),reverse=True)[299]);limits.append(sorted((min(counter[(t,b)],counter[('non-targeting',b)]//16) for t in targets),reverse=True)[299])
 bound=sum(limits);search=rt.json.loads(Path('artifacts/'+context+'_nonuniform_capacity.json').read_text());assert bound==search['integer_upper_bound'];cases=0;largest=0;deficit_budget=bound-64
 if 0<=deficit_budget<=3:
  # Exhaust every quota vector whose sum could reach 64. Each quota is <= its
  # independently counted coordinate ceiling, so total deficits are <= bound-64.
  def visit(position,budget,quotas):
   global cases,largest
   if position==len(batches):
    cases+=1;eligible=sum(all(counter[(t,b)]>=q for b,q in zip(batches,quotas)) for t in targets);largest=max(largest,eligible);assert eligible<300,'A feasible panel exists; cannot certify impossibility';return
   for d in range(min(budget,limits[position])+1):visit(position+1,budget-d,quotas+[limits[position]-d])
  visit(0,deficit_budget,[])
 plan=search['plan']
 if plan:
  assert len(set(plan['targets']))==300 and sum(plan['target_batch_quotas'].values())==plan['n_real']
  for b,q in plan['target_batch_quotas'].items():assert counter[('non-targeting',b)]>=16*q and all(counter[(t,b)]>=q for t in plan['targets'])
 impossible=bound<64 or cases>0
 results[context]=dict(status='PASS exact integer impossibility at >=64 shared real cells/target' if impossible else 'PASS upper bound; feasibility unresolved',census_sha256=rt.sha('artifacts/'+file),candidates=len(targets),coordinate_upper_bound_sum=bound,unlimited_control_coordinate_bound=sum(raw_limits),control_supply_was_not_binding=raw_limits==limits,batches=len(batches),deficit_vectors_exhausted=cases,max_eligible_in_exhaustive_test=largest if cases else None,certified_upper_bound=63 if cases else bound,scope='300 targets in the frozen prepared-source universe; common technical-batch quotas; two disjoint control arms, each 8x quota. If control supply was not binding, the same bound holds with unlimited controls. Not a claim about all public genes or all future source preparations.');print(context,results[context],flush=True)
rt.save('artifacts/matched_capacity_certificate.json',dict(verifier_sha256=rt.sha(__file__),results=results))
