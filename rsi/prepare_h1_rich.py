"""Availability-only exact 400-cell H1 panels across all 48 batches.
Preserve original 160 perturbation cells and add 240; independent four 3200-cell control arms.
"""
import runtime as rt
from pathlib import Path
import hashlib,shutil
import pandas as pd,numpy as np
folder=Path('artifacts/tasks');manifest=rt.json.loads((folder/'manifest.json').read_text());obs=pd.read_parquet('artifacts/H1_full_obs.parquet');labels=obs.target_gene.astype(str).to_numpy();batch=obs.batch.astype(str).to_numpy();batches=sorted(set(batch));groups=obs.groupby(['target_gene','batch'],observed=True).indices;both=manifest['H1_search']['targets']+manifest['H1_holdout']['targets'];table=obs.groupby(['target_gene','batch'],observed=True).size().unstack(fill_value=0).reindex(columns=batches,fill_value=0);capacity=np.minimum(table.loc[both].min(axis=0).to_numpy(),table.loc['non-targeting'].to_numpy()//32);baseq=np.array([10 if b.startswith('Flex_2_') else 0 for b in batches]);assert np.all(capacity>=baseq) and capacity.sum()>=400;space=capacity-baseq;ideal=240*space/space.sum();added=np.floor(ideal).astype(int)
for i in sorted(range(len(batches)),key=lambda i:(-(ideal[i]-added[i]),batches[i]))[:240-int(added.sum())]:added[i]+=1
quota=baseq+added;assert quota.sum()==400 and np.all(quota<=capacity);rng=np.random.default_rng(20261017);controls={k:[] for k in ['search_input','search_reference','holdout_input','holdout_reference']}
for b,q in zip(batches,quota):
 ix=groups[('non-targeting',b)];chosen=rng.choice(ix,32*int(q),replace=False)
 for j,arm in enumerate(controls):controls[arm].extend(chosen[j*8*q:(j+1)*8*q].tolist())
blueprint={}
for split in ['search','holdout']:
 parent='H1_'+split;task='H1_rich_'+split;old=pd.read_parquet(folder/(parent+'_selection.parquet'));old=old[old.evaluation_role=='perturbation'];selected=[];roles=[]
 for t in manifest[parent]['targets']:
  previous=old.loc[old.target_gene.astype(str)==t,'source_row'].to_numpy(np.int64);pick=list(previous)
  for b,q in zip(batches,quota):
   used=previous[batch[previous]==b];available=groups[(t,b)];available=available[~np.isin(available,used)];need=int(q)-len(used);assert need>=0;pick.extend(rng.choice(available,need,replace=False).tolist())
  assert len(set(pick))==400;selected.extend(pick);roles.extend(['perturbation']*400)
 for role,arm in [('input_control',split+'_input'),('reference_control',split+'_reference')]:selected.extend(controls[arm]);roles.extend([role]*len(controls[arm]))
 sel=obs.iloc[selected].copy();sel['source_row']=selected;sel['evaluation_role']=roles;assert sel.source_row.is_unique;sel.to_parquet(folder/(task+'_selection.parquet'));m=dict(manifest[parent]);m.update(source='H1_full_rich',parent_task=parent,n_real=400,selection_sha256=rt.sha(folder/(task+'_selection.parquet')),reference_precision='Exactly 400 real cells; same frozen targets; all 48 batches; original 160 perturbation cells retained');m['audit']=dict(exact_target_batch_matching=True,target_batch_quotas=dict(zip(batches,map(int,quota))),control_ratio=8,control_batch_quotas=dict(zip(batches,map(lambda q:int(8*q),quota))),per_target_control_mixture_L1={t:0.0 for t in m['targets']},input_reference_disjoint=True,all_four_control_arms_disjoint=True);blueprint[task]=m
rt.save(folder/'h1_rich_blueprint.json',blueprint);print('Exact quota capacity',int(capacity.sum()),'selected',dict(zip(batches,map(int,quota))),'controls per arm',len(controls['search_input']),flush=True)
