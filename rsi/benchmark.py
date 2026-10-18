"""Resumable frozen-recipe calibration using bitwise-verified bounded official scoring."""
import runtime as rt
import argparse,time,gc
from pathlib import Path
from dataclasses import replace
import numpy as np,pandas as pd,anndata as ad
from cell_eval2.anchor import compute_replicate_anchor
from cell_eval2 import aggregate_metrics
import official_scoring as scoring
import bounded_score
from recipes import predict
NT='non-targeting'

def setup(task):
 io_receipt=rt.json.loads(Path('artifacts/task_io_fixture.json').read_text());assert io_receipt['status']=='PASS bitwise CSR and label identities' and io_receipt['implementation_sha256']==rt.sha('task_io.py')
 m=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())[task];path=Path('artifacts/tasks')/(task+'.h5ad');receipt=rt.json.loads(path.with_suffix('.json').read_text());assert rt.sha(path)==receipt['sha256'];ctrl,real=__import__('task_io').split(path)
 sel=pd.read_parquet(path.with_name(task+'_selection.parquet'))
 if 'guide_id' in sel:
  guides=dict(zip(sel.source_row,sel.guide_id.astype(str)));real.obs['guide_id']=[guides[i] for i in real.obs.source_row];ctrl.obs['guide_id']=[guides[i] for i in ctrl.obs.source_row]
 assert set(ctrl.obs.source_row).isdisjoint(set(real.obs.source_row));bounded_score.validate(real.X)
 if m['exact_matched']:
  bs=sorted(real.obs.batch.astype(str).unique())
  quotas=m.get('audit',{}).get('target_batch_quotas')
  expected=np.array([quotas[b] for b in bs]) if quotas else None
  for t in m['targets']:
   nn=real.obs[real.obs.target_gene.astype(str)==t].batch.astype(str).value_counts().reindex(bs,fill_value=0).to_numpy();assert np.all(nn==expected) if expected is not None else (nn==nn[0]).all() and nn[0]>0
  for role,part in [('input',ctrl),('reference',real[real.obs.target_gene.astype(str)==NT])]:
   ratio=m['audit'].get(role+'_control_ratio',m['audit'].get('control_ratio'));nn=part.obs.batch.astype(str).value_counts().reindex(bs,fill_value=0).to_numpy();assert np.all(nn==expected*ratio) if expected is not None else (nn==nn[0]).all()
 return m,ctrl,real

def baseline_mean(real):
 mask=real.obs.target_gene.astype(str).to_numpy()==NT;x=real.X[mask].astype(float)
 if 'guide_id' in real.obs:
  groups=real.obs.loc[mask,'guide_id'].astype(str).to_numpy();return np.mean([np.asarray(x[groups==g].mean(0)).ravel() for g in np.unique(groups)],axis=0)
 return np.asarray(x.mean(0)).ravel()

def aggregate(tidy):
 frame=aggregate_metrics(tidy);return dict(zip(frame['metric'].to_list(),frame['mean'].to_list()))

def anchors(task,real,folder):
 policy=rt.json.loads(Path('artifacts/anchor_policy.json').read_text());n_splits=policy.get(task,policy['default_splits']);file=folder/'anchors.json';identity=scoring.matrix_id(real)
 if file.exists():
  r=rt.json.loads(file.read_text());assert r['real_sha256']==identity and r['scorer_sha256']==rt.sha('official_scoring.py') and r['driver_sha256']==rt.sha('bounded_score.py')
  if isinstance(r.get('baseline_construction'),dict):return r
  rt.save(folder/'anchors_control_mean_sensitivity.json',r)
  b,description=response_baseline(task,real,folder);scoring.scale(b,b,r['replicate']);r.update(baseline=b,baseline_construction=description);rt.save(file,r);rescale_existing(folder,r);return r
 cfg=scoring.config(cache_real=str(folder/'real_cache'));b,description=response_baseline(task,real,folder)
 print(task,n_splits,'official split replicate anchors',flush=True);splits,frame=compute_replicate_anchor(real,config=replace(cfg,cache_real=None),n_splits=n_splits);splits.write_csv(folder/'anchor_splits.csv');frame.write_csv(folder/'anchor_summary.csv');r=dict(zip(frame['metric'].to_list(),frame['replicate'].to_list()));scoring.scale(b,b,r);out=dict(baseline=b,replicate=r,real_sha256=identity,scorer_sha256=rt.sha('official_scoring.py'),driver_sha256=rt.sha('bounded_score.py'),version=__import__('cell_eval2').__version__,baseline_construction=description,n_splits=n_splits);rt.save(file,out);return out

def response_baseline(task,real,folder):
 from baseline_protocol import mean_response
 file=folder/'baseline_mean_response.json';description_file=folder/'baseline_mean_response_protocol.json'
 if file.exists():
  description=rt.json.loads(description_file.read_text());current=rt.sha('baseline_protocol.py')
  if description['builder_sha256']!=current:
   fixture=rt.json.loads(Path('artifacts/baseline_profile_fixture.json').read_text());assert fixture['status']=='PASS bitwise profile equality and unchanged reference' and fixture['builder_sha256']==current and description['builder_sha256'] in fixture['compatible_legacy_builder_sha256'];description['equivalence_certificate_sha256']=rt.sha('artifacts/baseline_profile_fixture.json')
  assert description['reference_sha256']==scoring.matrix_id(real);return rt.json.loads(file.read_text()),description
 mean,description=mean_response(real);description['reference_sha256']=scoring.matrix_id(real);print(task,'Arc-documented mean-response tiled baseline',flush=True);tidy,_=bounded_score.score_iter(real,bounded_score.baseline_iter(real,mean,cells=400),replace(scoring.config(cache_real=str(folder/'real_cache')),allow_fractional_counts=True));b=aggregate(tidy);rt.save(description_file,description);rt.save(file,b);return b,description

def rescale_existing(folder,anchors):
 for path in folder.glob('*_s[0-9]*.json'):
  record=rt.json.loads(path.read_text());assert record['reference_sha256']==anchors['real_sha256'];archive=Path('artifacts/control_scale_history')/folder.name/path.name
  if not archive.exists():rt.save(archive,record)
  record['scaled'],record['overall']=scoring.scale(record['raw'],anchors['baseline'],anchors['replicate']);record['scaling_revision']=dict(protocol=anchors['baseline_construction'],reason='Arc primary documentation supersedes public H1 harness control-mean description; original raw metrics unchanged',original_record_sha256=rt.sha(archive));rt.save(path,record)

def release_models():
 import recipes
 recipes.source.cache_clear();recipes.robust_bulk.cache_clear();recipes.v4._POOL_CACHE.clear();recipes.v2._FINE=None;recipes.v2._POST_STATE=None;gc.collect()

def run(task,recipe,seed,holdout_final=False):
 start=time.time();contract=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())[task]
 if contract['split']=='holdout' and not holdout_final:raise RuntimeError('HOLDOUT sealed: final/checkpoint event required')
 if contract['split']=='holdout':
  with open('artifacts/holdout_events.jsonl','a') as event:event.write(rt.json.dumps(dict(utc=time.strftime('%Y-%m-%d %H:%M:%S UTC',time.gmtime()),task=task,recipe=recipe,seed=seed,event='final report audit'))+'\n')
 m,ctrl,real=setup(task);folder=Path('artifacts/scores')/task;folder.mkdir(exist_ok=True,parents=True);path=folder/f'{recipe}_s{seed}.json'
 if path.exists():
  old=rt.json.loads(path.read_text());assert old['implementation_sha256']==__import__('provenance').identity(task);print('Cached',task,recipe,flush=True);return old
 rt.status('1 calibration',f'Task {task}; recipe {recipe}; seed {seed}. Loop gated pending calibration. HOLDOUT '+('FINAL AUDIT' if m['split']=='holdout' else 'sealed.'))
 import bounded_de
 bounded_de.install();assert scoring.provenance()['official_files']==rt.json.loads(Path('artifacts/bounded_score_fixture.json').read_text())['official']['official_files'];__import__('stat_cache').install();bounded_score.install_anchor();invalid=folder/'invalid_scale.json'
 if invalid.exists() and rt.json.loads(invalid.read_text()).get('reference_sha256')==scoring.matrix_id(real):raise ValueError('Cached invalid six-metric scale; redesign task before retry')
 try:a=anchors(task,real,folder)
 except (ValueError,MemoryError) as error:
  rt.save(invalid,dict(task=task,reference_sha256=scoring.matrix_id(real),error=str(error),status='INVALID full six-metric scale; whole task excluded'));raise
 if m['context']=='CD4':
  source_receipt=rt.json.loads(Path('artifacts/CD4_source_GWPS_transfer.json').read_text());assert source_receipt['original_controls_bitwise_equal'];checks=rt.json.loads(Path('artifacts/cd4_sources/legacy_reproduction.json').read_text());assert all(v['quantized_bitwise_equal'] for v in checks.values())
 signature=__import__('provenance').identity(task)
 if (Path('artifacts/raw_scores')/task/f'{recipe}_s{seed}.json').exists():
  assert __import__('raw_calibration').materialize(task,recipe,seed)
  return rt.json.loads(path.read_text())
 seen=[]
 def iterator():
  for t,x in predict(ctrl,m['targets'],m['context'],recipe,seed):
   assert x.shape==(400,real.n_vars);seen.append(t);rt.guard(512*1024**2)
   if len(seen)%25==1:print(task,recipe,'emission',len(seen),'/',len(m['targets']),'seconds',round(time.time()-start),flush=True)
   yield t,x
 tidy,prediction_sha=bounded_score.score_iter(real,iterator(),scoring.config(cache_real=str(folder/'real_cache')),release=release_models);raw=aggregate(tidy);scaled,overall=scoring.scale(raw,a['baseline'],a['replicate']);assert set(seen)==set(m['targets']);record=dict(task=task,recipe=recipe,seed=seed,raw=raw,scaled=scaled,overall=overall,seconds=time.time()-start,implementation_sha256=signature,prediction_sha256=prediction_sha,reference_sha256=a['real_sha256'],selection_sha256=m['selection_sha256'],n_targets=len(m['targets']),cells_pred=400,cells_real=m['n_real'],exact_target_matched=m['exact_matched'],primary_calibration_eligible=m['eligible_for_primary_calibration'],status='Public-context estimate; not an official score or biological verification');rt.save(path,record);print('SCORE',task,recipe,overall,scaled,flush=True);return record

if __name__=='__main__':
 if not Path('artifacts/atlas_reproduction.json').exists():import reproduce_atlas
 p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--recipes',default='v1,v2b,v3b,v4a,atlas_adapted,v4b');p.add_argument('--seed',type=int,default=42);p.add_argument('--final-holdout',action='store_true');args=p.parse_args()
 for recipe in args.recipes.split(','):run(args.task,recipe,args.seed,args.final_holdout)
 __import__('verify_results').calibrate()
