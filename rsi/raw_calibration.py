"""Checkpoint fixed historical recipe metrics independently of expensive anchors.
This is calibration scheduling, not a model search. HOLDOUT is refused.
"""
import runtime as rt
import argparse,time
from pathlib import Path
import benchmark as b,bounded_de,bounded_score,stat_cache,official_scoring as s

def run(task,recipes,seed):
 contract=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())[task]
 assert contract['split']=='search','Raw calibration cannot access sealed HOLDOUT'
 m,ctrl,real=b.setup(task);reference=s.matrix_id(real);folder=Path('artifacts/raw_scores')/task;folder.mkdir(parents=True,exist_ok=True)
 bounded_de.install();stat_cache.install();signature=__import__('provenance').identity(task)
 assert s.provenance()['official_files']==rt.json.loads(Path('artifacts/bounded_score_fixture.json').read_text())['official']['official_files']
 if m['context']=='CD4':
  assert rt.json.loads(Path('artifacts/CD4_source_GWPS_transfer.json').read_text())['original_controls_bitwise_equal']
  assert all(x['quantized_bitwise_equal'] for x in rt.json.loads(Path('artifacts/cd4_sources/legacy_reproduction.json').read_text()).values())
 for recipe in recipes:
  assert recipe in ['v1','v2b','v3b','v4a','v4b','atlas_adapted']
  path=folder/f'{recipe}_s{seed}.json'
  if path.exists():
   record=rt.json.loads(path.read_text());assert record['implementation_sha256']==signature and record['raw_driver_sha256']==rt.sha(__file__);print('Cached raw',task,recipe,flush=True);continue
  invalid=Path('artifacts/scores')/task/'invalid_scale.json'
  if invalid.exists() and rt.json.loads(invalid.read_text())['reference_sha256']==reference:raise ValueError('Whole task has invalid scale; no further calibration work')
  start=time.time();seen=[]
  def predictions():
   for target,x in b.predict(ctrl,m['targets'],m['context'],recipe,seed):
    assert x.shape==(400,real.n_vars);seen.append(target);rt.guard(512*1024**2)
    if len(seen)%25==1:print(task,recipe,'raw emission',len(seen),'/',len(m['targets']),flush=True)
    yield target,x
  tidy,digest=bounded_score.score_iter(real,predictions(),s.config(cache_real=str(Path('artifacts/scores')/task/'real_cache')),release=b.release_models)
  assert set(seen)==set(m['targets']);record=dict(task=task,recipe=recipe,seed=seed,raw=b.aggregate(tidy),seconds=time.time()-start,implementation_sha256=signature,raw_driver_sha256=rt.sha(__file__),prediction_sha256=digest,reference_sha256=reference,selection_sha256=m['selection_sha256'],n_targets=len(m['targets']),cells_pred=400,cells_real=m['n_real'],exact_target_matched=m['exact_matched'],primary_calibration_eligible=m['eligible_for_primary_calibration'],status='RAW ONLY: no normalized overall until all six official scales validate')
  rt.save(path,record);print('RAW COMPLETE',task,recipe,record['raw'],flush=True);materialize(task,recipe,seed)

def materialize(task,recipe,seed):
 folder=Path('artifacts/scores')/task;anchor=folder/'anchors.json';raw=Path('artifacts/raw_scores')/task/f'{recipe}_s{seed}.json'
 if not anchor.exists() or not raw.exists():return False
 a=rt.json.loads(anchor.read_text());r=rt.json.loads(raw.read_text());assert r['implementation_sha256']==__import__('provenance').identity(task);assert r['raw_driver_sha256']==rt.sha(__file__)
 assert a['real_sha256']==r['reference_sha256'] and a['driver_sha256']==rt.sha('bounded_score.py') and a['scorer_sha256']==rt.sha('official_scoring.py');assert isinstance(a['baseline_construction'],dict)
 r['scaled'],r['overall']=s.scale(r['raw'],a['baseline'],a['replicate']);r['raw_checkpoint_sha256']=rt.sha(raw);r['status']='Public-context calibration estimate; not an official score or biological verification'
 from verify_results import check
 check(r,a);out=folder/f'{recipe}_s{seed}.json'
 if out.exists():
  old=rt.json.loads(out.read_text());assert old['prediction_sha256']==r['prediction_sha256'] and old['raw']==r['raw'],'Refuse conflicting completed recipe'
 else:rt.save(out,r)
 print('SCORE',task,recipe,r['overall'],flush=True);return True

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--recipes',default='v1,v2b,v3b,v4a,v4b');p.add_argument('--seed',type=int,default=42);p.add_argument('--materialize-only',action='store_true');a=p.parse_args()
 if a.materialize_only:
  for recipe in a.recipes.split(','):materialize(a.task,recipe,a.seed)
 else:run(a.task,a.recipes.split(','),a.seed)
 __import__('verify_results').calibrate()
