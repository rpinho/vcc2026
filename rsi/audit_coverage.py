"""Availability-only source audit, separate from response metrics and selection."""
import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd
manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());rows=[];banks=rt.ROOT.parent/'public_data/model_v4_signature_banks'
for task,m in manifest.items():
 targets=m['targets'];legacy=np.zeros(len(targets),int);atlas=np.zeros(len(targets),int);parts={}
 for name in ['H1','HepG2','Jurkat','K562','RPE1','K562_GWPS_bulk']:
  if name in ['RPE1','HEK293T'] or name==m['context'] or name.startswith(m['context']+'_'):continue
  with np.load(Path('artifacts/sources')/(name+'.npz')) as d:
   ix=pd.Index(d['targets']).get_indexer(targets);ok=ix>=0;mask=np.zeros(len(targets),bool);mask[ok]=d['n_cells'][ix[ok]]>=30
   if name=='H1':mask&=~np.isin(targets,manifest['H1_holdout']['targets'])
  legacy+=mask;parts[name]=int(mask.sum())
 for name in ['K562_GWPS_raw','HCT116']:
  if name==m['context'] or name.startswith(m['context']+'_'):continue
  with np.load(banks/(name+'.npz')) as d:
   ix=pd.Index(d['targets']).get_indexer(targets);ok=ix>=0;mask=np.zeros(len(targets),bool);mask[ok]=d['n_cells'][ix[ok]]>=20
  atlas+=mask;parts['v4_'+name]=int(mask.sum())
 with np.load(banks/'CD4_DE_statistics.npz') as d:
  ix=pd.Index(d['targets']).get_indexer(targets);ok=ix>=0;mask=np.zeros(len(targets),bool);available=(d['n_cells']>=20)&d['available']&d['quality_pass'];mask[ok]=available[:,ix[ok]].any(0)
 atlas+=mask;parts['v4_CD4']=int(mask.sum());rows.append(dict(task=task,context=m['context'],targets=len(targets),legacy_any_same_target=int((legacy>0).sum()),v4_any_same_target=int((atlas>0).sum()),**parts))
frame=pd.DataFrame(rows);frame.to_csv('artifacts/coverage_audit.csv',index=False);rt.save('artifacts/coverage_audit.json',dict(rows=rows,status='Availability only; does not measure biological transfer reliability',warning='Widened legacy banks reproduce old rows exactly. Adapted v4 banks retain the prepared target universe; low coverage can make nonofficial-target essential panels unsuitable as official-score calibrators.'))
print(frame.to_string(index=False))
