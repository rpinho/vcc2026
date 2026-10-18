"""Integer availability search and certified upper bounds for 300-target matching.
Only label/batch counts are used. This is benchmark design, never model search.
"""
import runtime as rt
import numpy as np,pandas as pd,hashlib
from pathlib import Path
shared=rt.ROOT.parent/'public_data/model_v4_signature_banks';universe=set(rt.json.loads(Path('artifacts/sources/target_universe.json').read_text()));available=set(np.load(shared/'HCT116.npz')['targets'].astype(str))&universe
for context,file,target,batch in [('GWPS','GWPS_raw_obs.parquet','gene','gem_group'),('CD4','CD4_cell_census.parquet','target_gene','batch')]:
 obs=pd.read_parquet('artifacts/'+file);table=obs.groupby([target,batch],observed=True).size().unstack(fill_value=0);table.columns=table.columns.astype(str);targets=sorted(available&set(table.index)-{'non-targeting'});batches=sorted(table.columns);controls=table.loc['non-targeting',batches].to_numpy(np.int64);counts=table.loc[targets,batches].to_numpy(np.int64);counts=np.minimum(counts,controls[None,:]//16);upper=np.sort(counts,axis=0)[-300];bound=int(upper.sum());use=upper>0;c=counts[:,use];history=[];best=None
 if bound>=64:
  # Several deterministic starts; retain only availability counts, never responses.
  for mode in ['total','log','hash']:
   values=c.sum(1) if mode=='total' else np.log1p(c).sum(1) if mode=='log' else np.array([int(hashlib.sha256(('matching:'+t).encode()).hexdigest()[:12],16) for t in targets]);ix=np.argsort(-values,kind='stable')[:300];score=int(c[ix].min(0).sum())
   for step in range(80):
    sub=c[ix];order=np.argsort(sub,axis=0,kind='stable');low=sub[order[0],np.arange(c.shape[1])];second=sub[order[1],np.arange(c.shape[1])];outside=np.setdiff1d(np.arange(len(targets)),ix);candidate=None
    for j in range(len(ix)):
     without=np.where(order[0]==j,second,low);scores=np.minimum(c[outside],without).sum(1);k=int(np.argmax(scores)) if len(scores) else 0
     if len(scores) and int(scores[k])>score and (candidate is None or int(scores[k])>candidate[0]):candidate=(int(scores[k]),j,int(outside[k]))
    if candidate is None:break
    score,j,k=candidate;ix[j]=k
   history.append(dict(start=mode,capacity=score))
   if best is None or score>best[0]:best=(score,ix.copy())
 plan=None
 if best and best[0]>=64:
  ix=best[1];q=counts[ix].min(0);n=min(400,int(q.sum()))
  if q.sum()>n:
   ideal=n*q/q.sum();quota=np.floor(ideal).astype(np.int64)
   for j in np.argsort(-(ideal-quota),kind='stable')[:n-int(quota.sum())]:quota[j]+=1
   q=quota
  assert q.sum()==n and np.all(counts[ix]>=q);plan=dict(context='K562' if context=='GWPS' else 'CD4',targets=sorted(np.array(targets)[ix].tolist()),target_batch_quotas={b:int(v) for b,v in zip(batches,q) if v},n_real=n,control_ratio=8,availability_search_only=True)
 rt.save('artifacts/'+context+'_nonuniform_capacity.json',dict(status='FEASIBLE' if plan else 'IMPOSSIBLE at 64 cells: certified upper bound' if bound<64 else 'No feasible panel found by availability-only search; not an impossibility proof',panel_size=300,comparable_candidates=len(targets),integer_upper_bound=bound,bound_explanation='For any 300-target subset, each shared batch quota is at most the 300th-largest target count, capped by half the control supply divided by 8. Sum is an integer upper bound independent of subset search.',history=history,plan=plan,census_sha256=rt.sha('artifacts/'+file)));print(context,'upper bound',bound,'search',history,'plan',plan,flush=True)
