"""Official scorer validation. Local SEARCH anchors are never hidden-panel anchors."""
import local_runtime as j
j.enforce_io();j.guard()
import argparse,gc,sys
from dataclasses import replace
import numpy as np,pandas as pd,anndata as ad,polars as pl
from scipy import sparse
import cell_eval2
from cell_eval2 import EvalConfig,compute_metrics,aggregate_metrics
from cell_eval2.score import _replicate_entries,_reference_column
SIX=['expr_mse_unbiased_capped_norm','pds_cosine','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
NT='non-targeting'

def config(fractional=False):
 c=EvalConfig.from_preset('vcc2026')
 return replace(c,pert_col='target_gene',device='cpu',num_threads=1,gather_threads=1,
                allow_fractional_counts=fractional,de=replace(c.de,backend='pdex'))

def matrix_id(a):
 import hashlib
 h=hashlib.sha256();x=a.X.tocsr()
 for v in [x.data,x.indices,x.indptr]:h.update(v.tobytes())
 h.update('\0'.join(a.obs.target_gene.astype(str)).encode());h.update('\0'.join(a.var_names).encode());return h.hexdigest()

def aggregate(frame):
 a=aggregate_metrics(frame);return dict(zip(a['metric'].to_list(),a['mean'].to_list()))

def scale(raw,anchors):
 entries=_replicate_entries(anchors['baseline'],pl.DataFrame({'metric':list(anchors['replicate']),'replicate':list(anchors['replicate'].values())}))
 col=_reference_column(SIX+['avg_score'],[raw[m] for m in SIX],SIX,entries,column='from_replicate',label='validated local SEARCH anchors, NOT official leaderboard panel')
 return dict(zip(SIX,col.to_list()[:-1])),float(col[-1])

def make(x,labels,genes):
 return ad.AnnData(sparse.csr_matrix(x),obs=pd.DataFrame({'target_gene':labels},index=['p'+str(i) for i in range(len(labels))]),var=pd.DataFrame(index=genes))

def read_task(task):
 p=j.REPO/'rsi/artifacts/tasks'/f'{task}.h5ad'
 a=ad.read_h5ad(p);mask=a.obs.evaluation_role.astype(str).to_numpy()=='input_control'
 ctrl=a[mask].copy();real=a[~mask].copy();del a;gc.collect()
 return ctrl,real

def bounded():
 j.modules();sys.path.insert(0,str(j.REPO/'rsi'))
 import bounded_score
 return bounded_score

def versions():
 from importlib.metadata import version
 names=['cell-eval2','pdex','numba-mwu','numpy','scipy','anndata','polars','pandas','h5py','numba','torch']
 return {n:version(n) for n in names}

def validate():
 start=j.utc();task='CD4_depth192_search';ctrl,real=read_task(task)
 anchor_path=j.REPO/'rsi/artifacts/scores'/task/'anchors.json';anchors=j.json.loads(anchor_path.read_text())
 identity=matrix_id(real);assert identity==anchors['real_sha256'],'Reference anchor mismatch'
 print('ANCHOR IDENTITY PASS',identity,flush=True)
 labels=real.obs.target_gene.astype(str).to_numpy();targets=sorted(set(labels)-{NT})
 # Small same-gene-axis fixture, selected without using outcomes.
 ix=np.concatenate([np.flatnonzero(labels==t)[:64] for t in targets[:3]]+[np.flatnonzero(labels==NT)[:128]])
 fixture=real[ix].copy()
 from cell_eval2.baseline import generic_response_profile
 def profile(a):
  a=a.copy();labs=a.obs.target_gene.astype(str).to_numpy();col='target_gene'
  if 'guide_id' in a.obs:
   col='baseline_construct';a.obs[col]=np.where(labs==NT,NT,a.obs.guide_id.astype(str).to_numpy())
  return generic_response_profile(a,pert_col=col,control=NT,exclude_target_gene=False).values
 mean=profile(fixture);ft=targets[:3]
 pred=make(sparse.vstack([sparse.csr_matrix(np.tile(mean.astype(np.float32),(400*len(ft),1))),fixture.X[fixture.obs.target_gene.astype(str).to_numpy()==NT]],format='csr'),sum(([t]*400 for t in ft),[])+[NT]*128,fixture.var_names)
 print('UNMODIFIED compute_metrics fixture',flush=True)
 native=compute_metrics(pred,fixture,config=config(True));native_raw=aggregate(native)
 b=bounded();print('BOUNDED fixture equality',flush=True)
 reduced,_=b.score_iter(fixture,((t,sparse.csr_matrix(np.tile(mean.astype(np.float32),(400,1)))) for t in ft),config(True))
 n=native.sort(['perturbation','metric']);r=reduced.sort(['perturbation','metric'])
 assert n.select('perturbation','metric').equals(r.select('perturbation','metric'))
 nv=n['value'].to_numpy();rv=r['value'].to_numpy()
 assert np.array_equal(np.isnan(nv),np.isnan(rv))
 assert np.array_equal(np.isfinite(nv),np.isfinite(rv))
 err=float(np.max(np.abs(nv[np.isfinite(nv)]-rv[np.isfinite(rv)])))
 assert err<1e-12,(err,native_raw,aggregate(reduced))
 assert native_raw['pds_cosine']==0.5
 print('KNOWN PDS 0.5 + REDUCTION PASS',err,flush=True)
 # Recompute the historical full-panel baseline, all metrics, without altering cell_eval2.
 del pred,fixture,native,reduced;gc.collect();full_mean=profile(real)
 print('RECOMPUTE 79-target KNOWN baseline',flush=True)
 tidy,digest=b.score_iter(real,((t,sparse.csr_matrix(np.tile(full_mean.astype(np.float32),(400,1)))) for t in targets),config(True))
 raw=aggregate(tidy);errors={m:raw[m]-anchors['baseline'][m] for m in SIX}
 validated=all(abs(v)<1e-10 for v in errors.values())
 scaled,overall=scale(raw,anchors) if validated else (None,None)
 record=dict(status='PASS' if validated else 'FAIL',utc=j.utc(),started=start,task=task,real_matrix_id=identity,
  fixture=dict(targets=ft,real_cells=320,pred_cells=1328,native_raw=native_raw,reduction_max_abs=err,known_pds=.5),
  known_baseline=anchors['baseline'],recomputed_baseline=raw,errors=errors,local_scaled=scaled,local_score_avg=overall,
  normalization='Local SEARCH baseline/replicate anchors only; NOT hidden official panel',
  prediction_digest=digest,official_challenge_score=None)
 j.save('results/scorer_validation.json',record)
 inputs=[anchor_path,j.REPO/'rsi/artifacts/tasks'/f'{task}.h5ad',j.REPO/'rsi/artifacts/tasks'/f'{task}.json',
  j.REPO/'rsi/official_scoring.py',j.REPO/'rsi/bounded_score.py',j.REPO/'rsi/baseline_protocol.py',j.REPO/'rsi/task_io.py',
  j.REPO/'competitors/raw/lb-2026-10-25T164255Z.json',j.ROOT/'local_runtime.py',j.ROOT/'score_official.py']
 package=j.Path(cell_eval2.__file__).parent
 j.save('results/scoring_provenance.json',dict(utc=j.utc(),command='nice -n 10 python -B score_official.py --validate',
  versions=versions(),config=config().to_dict(),fractional_baseline_config=config(True).to_dict(),
  preset_pert_col=EvalConfig.from_preset('vcc2026').pert_col,local_pert_col='target_gene',inputs={str(p):j.sha(p) for p in inputs},
  package_files={str(p.relative_to(package)):j.sha(p) for p in package.rglob('*.py')},
  scorer_validated=validated,official_panel_reference_scaling_validated=False,process_pid=__import__('os').getpid()))
 print(record['status'],errors,flush=True)

def incumbent_fixture_data():
 # This is a schema/scorer fixture, not the frozen incumbent's v9 count renderer.
 # Its raw source-response vector is exact incumbent inference; fractional tiled
 # cells deliberately make NO claim to reproduce per-cell DE of the incumbent.
 task='CD4_depth192_search';ctrl,full=read_task(task)
 labels=full.obs.target_gene.astype(str).to_numpy();targets=sorted(set(labels)-{NT})[:3]
 ix=np.concatenate([np.flatnonzero(labels==t)[:64] for t in targets]+[np.flatnonzero(labels==NT)[:128]])
 real=full[ix].copy();del full
 genes=ctrl.var_names.to_numpy(str);control=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();control/=control.sum()
 depth=float(np.asarray(ctrl.X.sum(1)).mean());del ctrl;gc.collect()
 env=j.frozen();torch=env['torch'];src=env['Sources'](excluded=['CD4'])
 lookup=pd.Index(src.genes).get_indexer(genes);assert (lookup>=0).all()
 canonical=np.zeros(src.g);canonical[lookup]=control
 features=np.stack([src.features(t,canonical) for t in targets]);learned=[]
 for seed in [1019,2027,4093]:
  net=env['ResponseMLP'](src.g).cpu();net.load_state_dict(torch.load(j.REPO/'model_billmB/checkpoints/CD4'/f'seed{seed}_weights.pt',map_location='cpu',weights_only=True));net.eval()
  with torch.no_grad():learned.append(net(torch.from_numpy(features)))
  del net
 learned=torch.stack(learned).mean(0).numpy();response=env['blend'](features[:,:src.g],learned,.3)[:,lookup]
 L=np.log1p(50000*control)[None,:]+response
 probability=np.maximum(np.expm1(L),0);probability/=probability.sum(1,keepdims=True)
 profile=(probability*depth).astype(np.float32)
 blocks=[sparse.csr_matrix(np.tile(profile[i],(400,1))) for i in range(3)]+[real.X[real.obs.target_gene.astype(str).to_numpy()==NT]]
 pred=make(sparse.vstack(blocks,format='csr'),sum(([t]*400 for t in targets),[])+[NT]*128,real.var_names)
 import hashlib
 meta=dict(task=task,targets=targets,real_cells=real.n_obs,pred_cells=pred.n_obs,genes=real.n_vars,
  real_matrix_id=matrix_id(real),pred_matrix_id=matrix_id(pred),response_sha256=hashlib.sha256(response.tobytes()).hexdigest(),
  scalar=1.,checkpoint_fold='CD4',count_construction='fractional tiled pseudobulk profile from exact incumbent post-blend response, plus real controls',
  qualification='Scorer/layout fixture ONLY. Does not reproduce frozen v9 renderer; not a calibration arm or full held-out-family score')
 del src,features,learned;gc.collect();return pred,real,meta

def score_incumbent_fixture():
 pred,real,meta=incumbent_fixture_data()
 raw=aggregate(compute_metrics(pred,real,config=config(True)))
 j.save('results/incumbent_fixture.json',dict(utc=j.utc(),fixture=meta,raw_unscaled=raw,score_avg=None,
  reason_no_composite='The fixture is a different real matrix from the full-panel anchors; applying those anchors would invent a scale',
  official_challenge_score=None))
 print('INCUMBENT LAYOUT FIXTURE, UNSCALED',raw,flush=True)

def score_count_fixture():
 # s=1 only. No calibrated arm, full-family result or cross-matrix scale.
 from renderer_adapter import Renderer
 task='CD4_depth192_search';ctrl,full=read_task(task)
 labels=full.obs.target_gene.astype(str).to_numpy();targets=sorted(set(labels)-{NT})[:3]
 ix=np.concatenate([np.flatnonzero(labels==t)[:64] for t in targets]+[np.flatnonzero(labels==NT)[:128]])
 real=full[ix].copy();del full;gc.collect()
 env=j.frozen();torch=env['torch'];src=env['Sources'](excluded=['CD4'])
 genes=ctrl.var_names.to_numpy(str);lookup=pd.Index(src.genes).get_indexer(genes);assert (lookup>=0).all()
 control=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();control/=control.sum()
 canonical=np.zeros(src.g);canonical[lookup]=control
 features=np.stack([src.features(t,canonical) for t in targets]);outputs=[]
 for seed in [1019,2027,4093]:
  net=env['ResponseMLP'](src.g).cpu();net.load_state_dict(torch.load(j.REPO/'model_billmB/checkpoints/CD4'/f'seed{seed}_weights.pt',map_location='cpu',weights_only=True));net.eval()
  with torch.no_grad():outputs.append(net(torch.from_numpy(features)))
  del net
 response=env['blend'](features[:,:src.g],torch.stack(outputs).mean(0).numpy(),.3)[:,lookup]
 del src,features,outputs;gc.collect()
 renderer=Renderer(ctrl,targets,response,'CD4');blocks=[];render_info=[]
 try:
  for i,t in enumerate(targets):
   print('INTEGER FIXTURE EMIT',i,t,flush=True)
   block,info=renderer.emit(i,1.);blocks.append(block);render_info.append(info)
  render_inputs=renderer.inputs;budget=renderer.budget_info
 finally:renderer.close()
 blocks.append(real.X[real.obs.target_gene.astype(str).to_numpy()==NT])
 pred=make(sparse.vstack(blocks,format='csr'),sum(([t]*400 for t in targets),[])+[NT]*128,real.var_names)
 assert np.equal(pred.X.data,np.floor(pred.X.data)).all()
 raw=aggregate(compute_metrics(pred,real,config=config(False)))
 paths={}
 for name,a in [('count_fixture_pred',pred),('count_fixture_real',real)]:
  path=j.ROOT/'results'/(name+'.npz');j.guard();sparse.save_npz(path,a.X)
  paths[name]=dict(path=str(path),sha256=j.sha(path),labels=a.obs.target_gene.astype(str).tolist())
 inputs=dict(j.json.loads((j.ROOT/'results/preflight.json').read_text())['inputs'])
 inputs.update(render_inputs)
 for path in [j.REPO/'rsi/artifacts/tasks'/f'{task}.h5ad',j.ROOT/'renderer_adapter.py',j.ROOT/'score_official.py',j.ROOT/'PREREGISTRATION.md']:
  inputs[str(path)]=j.sha(path)
 j.save('results/count_fixture.json',dict(utc=j.utc(),status='PASS integer renderer/scorer fixture ONLY',
  raw_unscaled=raw,score_avg=None,scalar=1.,targets=targets,genes=genes.tolist(),matrices=paths,
  real_matrix_id=matrix_id(real),pred_matrix_id=matrix_id(pred),renderer=render_info,budget=budget,inputs=inputs,
  qualification='Exact requested Sources/MLP/blend response + unchanged v6 integer emitter/v9 layout kernels. Not the production v5b amplitude recipe, a full outer fold, or an arm comparison.',
  command='nice -n 10 python -B score_official.py --count-fixture'))
 print('INTEGER FIXTURE UNSCALED',raw,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--incumbent-fixture',action='store_true');p.add_argument('--count-fixture',action='store_true');args=p.parse_args()
 if args.validate:validate()
 elif args.incumbent_fixture:score_incumbent_fixture()
 elif args.count_fixture:score_count_fixture()
 else:p.error('Choose --validate, --incumbent-fixture or --count-fixture')
