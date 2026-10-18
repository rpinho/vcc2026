"""Read-only metadata audit, kept separate from model selection."""
import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd
manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());results={}
for task,m in manifest.items():
 sel=pd.read_parquet(Path('artifacts/tasks')/(task+'_selection.parquet'));assert sel.source_row.is_unique;labs=sel.target_gene.astype(str);roles=sel.evaluation_role.astype(str);targets=set(labs[roles=='perturbation']);assert targets==set(m['targets']);assert set(labs[roles!='perturbation'])=={'non-targeting'};assert len(targets)==len(m['targets']);assert set(sel.loc[roles=='input_control','source_row']).isdisjoint(sel.loc[roles=='reference_control','source_row']);counts=labs[roles=='perturbation'].value_counts();expected=pd.Series(m['n_real']).reindex(counts.index) if isinstance(m['n_real'],dict) else m['n_real'];assert (counts==expected).all()
 if 'pass_guide_filter' in sel:assert sel.pass_guide_filter.all()
 if m['exact_matched']:
  b=sorted(sel.batch.astype(str).unique());q=m['audit'].get('target_batch_quotas');reference=None
  for target in m['targets']:
   observed=sel.loc[(roles=='perturbation')&(labs==target),'batch'].astype(str).value_counts().reindex(b,fill_value=0).to_numpy();expected=np.array([q[k] for k in b]) if q else np.repeat(m['n_real']//len(b),len(b));np.testing.assert_array_equal(observed,expected);reference=observed/observed.sum()
  for role in ['input_control','reference_control']:
   observed=sel.loc[roles==role,'batch'].astype(str).value_counts().reindex(b,fill_value=0).to_numpy();np.testing.assert_allclose(observed/observed.sum(),reference,atol=0,rtol=0)
 mismatch=np.array(list(m['audit']['per_target_control_mixture_L1'].values()));results[task]=dict(targets=len(targets),actual_cells_per_target=sorted(counts.unique().tolist()),input_controls=int((roles=='input_control').sum()),reference_controls=int((roles=='reference_control').sum()),exact_batch_match=m['exact_matched'],mixture_L1_median=float(np.median(mismatch)),mixture_L1_max=float(mismatch.max()),all_guide_calls_passed=bool(sel.pass_guide_filter.all()) if 'pass_guide_filter' in sel else None,selection_sha256=rt.sha(Path('artifacts/tasks')/(task+'_selection.parquet')))
assert set(manifest['H1_search']['targets']).isdisjoint(manifest['H1_holdout']['targets'])
if 'H1_rich_search' in manifest:
 arms=[]
 for task in ['H1_rich_search','H1_rich_holdout']:
  selected=pd.read_parquet(Path('artifacts/tasks')/(task+'_selection.parquet'))
  for role in ['input_control','reference_control']:arms.append(set(selected.loc[selected.evaluation_role==role,'source_row']))
 assert all(a.isdisjoint(b) for i,a in enumerate(arms) for b in arms[i+1:])
rt.save('artifacts/protocol_audit.json',dict(status='PASS raw selection/control disjointness and source guide QC',tasks=results,limitations='H1 and CD4 primary panels use exact batch quotas; other panels use pooled controls. Scorer validity and calibration are separate checks.'))
print('PASS protocol metadata audit',results)
