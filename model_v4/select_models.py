"""Predeclared equal-lineage selection from completed native-view six-metric scores.
The H1 20k view is a feasibility/depth sensitivity diagnostic, not a second lineage vote.
"""
import runtime
from fractions import Fraction as F
import pandas as pd
rt=runtime.ROOT
folders={'H1':rt/'artifacts/evaluation_native','HCT116':rt/'artifacts/evaluation_HCT116_native','HEK293T':rt/'artifacts/evaluation_HEK293T_native'}
from evaluate import VARIANTS
usable={}
for fold,folder in folders.items():
 if not (folder/'anchors.json').exists():continue
 records={p.stem:runtime.json.loads(p.read_text()) for p in folder.glob('*.json') if p.stem in VARIANTS}
 if records:usable[fold]=records
assert usable and 'H1' in usable
common=set.intersection(*(set(r) for r in usable.values()));assert common,'No variant completed every native lineage'
rows=[]
for variant in common:
 records=[r[variant] for r in usable.values()];exact=sum((F(r['exact_overall']) for r in records),F(0))/len(records);row=dict(variant=variant,mean_native_overall=float(exact),exact_mean_native=str(exact),lineages=len(records),fingerprint='|'.join(r['prediction_sha256'] for r in records))
 for fold,r in usable.items():row[fold]=r[variant]['overall']
 rows.append(row)
# Prefer requested guide/shrinkage integration only on exact arithmetic ties, never over a measured score difference.
rows.sort(key=lambda r:(F(r['exact_mean_native']),int(r['variant'].endswith('_guide')),int('ash' in r['variant']),r['variant']),reverse=True)
a=rows[0];distinct=[r for r in rows[1:] if r['fingerprint']!=a['fingerprint']];assert distinct
b=next(r for r in rows if r['variant']=='iv_ash_dual_guide');assert b['fingerprint']!=a['fingerprint'];record=dict(status='UNVERIFIED biological quality; public cross-lineage selection only',policy='Highest equal-lineage mean of six scaled native-view metrics; H1 20k retained as depth/feasibility sensitivity only; B is the requested full-integration ablation iv_ash_dual_guide, rather than the numerical runner-up; guide and adaptive shrinkage are retained for an explicit method comparison',a=a,b=b,numerical_runner_up=distinct[0],b_role='full-integration ablation',ranking=rows,native_lineages=list(usable),incomplete_variants=sorted(set(VARIANTS)-common))
runtime.save(rt/'artifacts/selection.json',record);pd.DataFrame(rows).drop(columns=['fingerprint']).to_csv(rt/'artifacts/selection_ranking.csv',index=False);print('v4a',a['variant'],a['mean_native_overall'],'v4b',b['variant'],b['mean_native_overall'],flush=True)
