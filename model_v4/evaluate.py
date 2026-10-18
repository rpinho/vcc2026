"""Frozen attribution experiment on real batch-balanced public H1 holdout.
Published hidden A/B/C scores cannot be reproduced without server labels.
"""
import runtime,argparse,time,hashlib
from dataclasses import replace
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from cell_eval2 import EvalConfig,compute_metrics,aggregate_metrics
from cell_eval2.config import DEParams
from model import heads,template,emit,atlas,BANKS
from independent_metrics import verify
from verify_counts import moment_certificate,scale_record
rt=runtime.ROOT;NT='non-targeting'
SIX=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
VARIANTS=['anchor','anchor_restricted','anchor_single_mean','anchor_single_dual','anchor_dual_mean','anchor_nopromoter','old_single_mean','old_single_dual','new_single_mean','new_single_dual','new_dual_mean','new_dual_dual','new_smooth_dual','new_smooth_dual_guide','iv_none_dual','iv_js_dual','iv_ash_dual','iv_ash_dual_guide','new_dual_dual_guide']

def thin(a,seed):
 a=a.copy();x=a.X.tocsr();d=np.asarray(x.sum(1)).ravel();rng=np.random.default_rng(seed);p=np.minimum(1,20000/d);x.data=rng.binomial(x.data.astype(np.int64),np.repeat(p,np.diff(x.indptr))).astype(np.int32);x.eliminate_zeros();a.X=x;return a

def identity(a):
 h=hashlib.sha256()
 for v in [a.X.data,a.X.indices,a.X.indptr]:h.update(v.tobytes())
 h.update('\0'.join(a.var_names).encode());h.update('\0'.join(a.obs_names).encode());return h.hexdigest()

def eval_folder(view,fold='H1'):
 return rt/('artifacts/evaluation_native' if view=='native' else 'artifacts/evaluation') if fold=='H1' else rt/f'artifacts/evaluation_{fold}_{view}'

def setup(view='depth20k',fold='H1'):
 a=ad.read_h5ad(rt.parent/'public_data/v2_diagnostic/H1_balanced.h5ad' if fold=='H1' else rt/f'artifacts/validation_{fold}.h5ad');assert a.uns['evaluation_only'];ctrl=a[a.obs.evaluation_role.astype(str)=='input_control'].copy();real=a[a.obs.evaluation_role.astype(str)!='input_control'].copy();del a
 assert set(ctrl.obs_names).isdisjoint(real.obs_names)
 if view=='depth20k' and fold=='H1':
  ctrl=thin(ctrl,816);real=thin(real,815);previous=rt.parent/'model_v2/artifacts/balanced_h1';anchors=runtime.json.loads((previous/'anchors.json').read_text());assert anchors['input_arrays_sha256']==identity(ctrl);assert anchors['reference_arrays_sha256']==identity(real)
  bases=runtime.json.loads((previous/'baseline_audit.json').read_text());base=bases['tile_oracle']
 else:
  folder=eval_folder(view,fold);folder.mkdir(exist_ok=True);path=folder/'anchors.json'
  if path.exists():
   anchors=runtime.json.loads(path.read_text());assert anchors['input_arrays_sha256']==identity(ctrl) and anchors['reference_arrays_sha256']==identity(real);base=anchors['tiled_baseline']
  else:
   from cell_eval2.anchor import compute_replicate_anchor
   from cell_eval2.baseline import generic_response_profile,build_baseline_prediction
   cfg=replace(EvalConfig.from_preset('vcc2026'),pert_col='target_gene',device='cpu',num_threads=2,gather_threads=2,de=DEParams(backend='pdex'),cache_real=None,allow_fractional_counts=True)
   def aggregate(pred):
    try:return aggregate_metrics(compute_metrics(pred,real,config=cfg)).to_pandas().set_index('metric')['mean'].to_dict()
    except ValueError as error:
     if 'expr_mse_unbiased_capped_norm' in str(error) and 'not positive' in str(error):
      runtime.save(folder/'reference_invalid.json',dict(status='UNDEFINED official MSE reference normalization; excluded whole fold from composite ranking',fold=fold,view=view,error=str(error),targets=sorted(set(real.obs.target_gene.astype(str))-{NT})))
     raise

   profile=generic_response_profile(real,pert_col='target_gene',control=NT)
   base=aggregate(build_baseline_prediction(profile,real,pert_col='target_gene',control=NT,emit='tile'));base['pds_cosine']=.5
   dispersed=aggregate(build_baseline_prediction(profile,real,pert_col='target_gene',control=NT,seed=0))
   splits,summary=compute_replicate_anchor(real,config=replace(cfg,allow_fractional_counts=False),n_splits=5);splits.write_csv(folder/'anchor_splits.csv');rep=dict(zip(summary['metric'].to_list(),summary['replicate'].to_list()))
   anchors=dict(baseline=dispersed,replicate=rep,input_arrays_sha256=identity(ctrl),reference_arrays_sha256=identity(real),targets=sorted(set(real.obs.target_gene.astype(str))-{NT}),tiled_baseline=base)
 folder=eval_folder(view,fold);runtime.save(folder/'anchors.json',dict(**{k:v for k,v in anchors.items() if k not in ('tiled_baseline','view')},tiled_baseline=base,view=view))
 return real,ctrl,anchors,base

def main(args):
 required=['HCT116','HEK293T','CD4_DE_statistics']
 for n in required:assert (BANKS/f'{n}.npz').exists(),f'Missing required bank {n}'
 folder=eval_folder(args.view,args.fold);folder.mkdir(exist_ok=True);real,ctrl,anchor,base=setup(args.view,args.fold);rep=anchor['replicate'];targets=np.array(anchor['targets']);genes=real.var_names.to_numpy();labels=real.obs.target_gene.astype(str).to_numpy();q,d,mean,bulk=template(ctrl,genes,seed=20260910)
 config=replace(EvalConfig.from_preset('vcc2026'),pert_col='target_gene',device='cpu',num_threads=2,gather_threads=2,de=DEParams(backend='pdex'),cache_real=str(folder/'real_cache'))
 runtime.save(folder/'protocol.json',dict(heldout='Entire '+args.fold+' source family excluded',targets=targets.tolist(),cells_per_target=400,input_controls=ctrl.n_obs,reference_controls=int((labels==NT).sum()),variants=VARIANTS,source_policy='Exclude all '+args.fold+' banks; fixed remaining public source families',status='UNVERIFIED public selection estimates; same targets reused across ablations; no hidden A/B/C labels',selection='Highest mean six-metric scaled score under conservative tiled oracle; examine raw MSE and dispersed-baseline sensitivity; runner-up highest distinct completed candidate',seed=20260910,metric_version=__import__('cell_eval2').__version__))
 selected=VARIANTS if not args.variants else args.variants.split(',');records=[];model_sha=runtime.sha(rt/'model.py');source_sha=runtime.sha(rt/'artifacts/source_identity.json')
 for variant in selected:
  path=folder/f'{variant}.json'
  if path.exists():
   previous=runtime.json.loads(path.read_text());assert previous['model_sha256']==model_sha and previous['source_identity_sha256']==source_sha,'Stale attribution cache';records.append(previous);continue
  if variant!='anchor':assert (folder/'anchor.json').exists() or (rt/'artifacts/evaluation_native/anchor.json').exists(),'Anchor must execute first'
  start=time.time();stats={};(p,b),coverage=heads(targets,genes,mean,bulk,variant,exclude=args.fold,stats=stats);chunks=[];obs=[];mom=[]
  try:
   for i,t in enumerate(targets):
    x=emit(q,p[i],b[i],d,atlas.seed_for(args.fold+':'+t,20260910),variant);mom.append(dict(target=t,**moment_certificate(x,d,p[i],b[i])));chunks.append(sparse.csr_matrix(x));obs.extend([t]*400)
  except ValueError as e:
   runtime.save(folder/f'{variant}_failure.json',dict(error=str(e),target=t,variant=variant));raise
  chunks.append(real.X[labels==NT]);obs.extend([NT]*int((labels==NT).sum()));pred=ad.AnnData(sparse.vstack(chunks,format='csr'),obs=pd.DataFrame({'target_gene':obs},index=[f'p{i}' for i in range(len(obs))]),var=real.var.copy());del chunks
  raw=aggregate_metrics(compute_metrics(pred,real,config=config)).to_pandas().set_index('metric')['mean'].to_dict();independent=verify(pred,real,raw,exclusion_panel=pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.to_numpy(str))
  scaled,overall=scale_record(raw,base,rep,SIX);s2,o2=scale_record(raw,anchor['baseline'],rep,SIX)
  record=dict(variant=variant,model_sha256=model_sha,source_identity_sha256=source_sha,prediction_sha256=identity(pred),view=args.view,fold=args.fold,raw=raw,scaled={m:float(v) for m,v in scaled.items()},overall=float(overall),exact_scaled={m:str(v) for m,v in scaled.items()},exact_overall=str(overall),dispersed_baseline_overall=float(o2),independent=independent,moments=mom,adaptive_stats=stats,seconds=time.time()-start,status='UNVERIFIED biological estimate; exact rational scaling and count-constraint intervals verified')
  runtime.save(path,record);records.append(record);del pred
  print(variant,'overall',round(float(overall),6),'raw MSE',round(raw[SIX[1]],6),'PDS',round(raw[SIX[0]],6),'seconds',round(time.time()-start),flush=True)
  rows=[dict(variant=r['variant'],overall=r['overall'],raw_mse=r['raw'][SIX[1]],raw_pds=r['raw'][SIX[0]],**r['scaled']) for r in records];pd.DataFrame(rows).to_csv(folder/'scores_progress.csv',index=False)
 allrecords=[runtime.json.loads(p.read_text()) for p in folder.glob('*.json') if p.stem in VARIANTS]
 df=pd.DataFrame([dict(variant=r['variant'],overall=r['overall'],raw_mse=r['raw'][SIX[1]],raw_pds=r['raw'][SIX[0]],dispersed_baseline_overall=r['dispersed_baseline_overall'],**r['scaled']) for r in allrecords]);df=df.sort_values('overall',ascending=False);df.to_csv(folder/'attribution.csv',index=False)
 print(df.to_string(index=False),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--variants');p.add_argument('--view',choices=['native','depth20k'],default='native');p.add_argument('--fold',choices=['H1','HCT116','HEK293T'],default='H1');main(p.parse_args())
