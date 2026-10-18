"""Freeze a count-selected, fully lane-matched wider CD4 panel before scoring."""
import runtime as rt
import pandas as pd,numpy as np,hashlib
from pathlib import Path
obs=pd.read_parquet('artifacts/CD4_cell_census.parquet');raw=pd.read_parquet('artifacts/GWPS_raw_obs.parquet');source_counts=raw.gene.astype(str).value_counts();table=obs.groupby(['target_gene','batch'],observed=True).size().unstack(fill_value=0);batches=sorted(table.columns);eligible=set(table.index[table.min(axis=1)>=5])&set(source_counts.index[source_counts>=200]);eligible.discard('non-targeting');assert len(eligible)>=300;targets=sorted(sorted(eligible,key=lambda t:hashlib.sha256(('CD4_wide_20261016:'+t).encode()).hexdigest())[:300]);rng=np.random.default_rng(20261017);groups=obs.groupby(['target_gene','batch'],observed=True).indices;chosen=[];roles=[]
for t in targets:
 for b in batches:chosen.extend(rng.choice(groups[(t,b)],5,replace=False).tolist());roles.extend(['perturbation']*5)
arms={k:[] for k in ['input_control','reference_control']}
for b in batches:
 rows=rng.choice(groups[('non-targeting',b)],80,replace=False)
 for j,role in enumerate(arms):arms[role].extend(rows[40*j:40*(j+1)].tolist())
for role,rows in arms.items():chosen.extend(rows);roles.extend([role]*len(rows))
sel=obs.iloc[chosen].copy();sel['evaluation_role']=roles;extra=[];extra_rng=np.random.default_rng(20261018)
for b,g in obs[obs.target_gene.astype(str)=='non-targeting'].groupby('batch',observed=True):
 used=set(sel.loc[sel.batch==b,'source_row']);pool=g[~g.source_row.isin(used)];assert len(pool)>=60;extra.append(pool.iloc[extra_rng.choice(len(pool),60,replace=False)])
additional=pd.concat(extra,ignore_index=True);additional['evaluation_role']='input_control';sel=pd.concat([sel,additional],ignore_index=True);assert sel.source_row.is_unique;folder=Path('artifacts/tasks');sel.to_parquet(folder/'CD4_matched_search_selection.parquet',index=False)
m=dict(context='CD4',split='search',targets=targets,n_real=120,n_pred=400,exact_matched=True,source='CD4_raw_matched',selection_sha256=rt.sha(folder/'CD4_matched_search_selection.parquet'),eligible_for_primary_calibration=True,seal='SEARCH',source_profile='cd4_wide',audit=dict(target_batch_quotas={b:5 for b in batches},control_ratio=8,control_batch_quotas={b:40 for b in batches},per_target_control_mixture_L1={t:0.0 for t in targets},input_reference_disjoint=True),selection_policy='300 fixed SHA targets among 335 eligible, >=5 cells in each lane and >=200 GWPS source cells; no response-strength or model-score selection',budget_fallback='If exact planned ranges exceed 2.6 GB, reduce each target quota to 4 and both control quotas to 32 before response extraction; target identities unchanged')
m['audit'].update(input_control_ratio=20,reference_control_ratio=8,input_control_batch_quotas={b:100 for b in batches});m['control_revision']='Preserve reference controls; expand input to 2400 real matched cells so the frozen 400x4 template is executable. Budget fallback keeps input >=1920.';rt.save(folder/'cd4_wide_blueprint.json',m);rt.save('artifacts/CD4_wide_targets.json',targets);rt.save('artifacts/CD4_wide_design_event.json',dict(utc=__import__('time').strftime('%Y-%m-%d %H:%M:%S UTC',__import__('time').gmtime()),eligible=len(eligible),targets=300,reference_cells=120,source_cells_min=200,HOLDOUT='Unchanged, sealed',new_source_banks='Isolated; original banks never modified'))
print('Frozen wider CD4 panel',len(targets),'targets',len(sel),'cells',flush=True)
