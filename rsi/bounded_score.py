"""Forrest H1 bounded-reduction pattern, using official functions for ALL metrics.
Cell matrices are reduced one target at a time; full-panel PDS and MSE stay global.
Activation requires bitwise equality against unmodified compute_metrics, including anchors.
"""
import runtime as rt
from pathlib import Path
from dataclasses import replace
import hashlib,gc
import numpy as np,pandas as pd,polars as pl,anndata as ad
from scipy import sparse
from cell_eval2.prep import pseudobulk_bulk_lognorm_with_moments
from cell_eval2.moments import GroupMoments
from cell_eval2.catalog import resolve_metrics
import cell_eval2.run as run,cell_eval2.anchor as anchor
NT='non-targeting';SCHEMA={'perturbation':pl.Utf8,'metric':pl.Utf8,'value':pl.Float64}

def validate(x,allow_fractional=False):
 values=x.data if sparse.issparse(x) else np.asarray(x)
 assert np.isfinite(values).all() and np.all(values>=0)
 if not allow_fractional:assert np.equal(values,np.floor(values)).all()
 totals=np.asarray(x.sum(1)).ravel();assert np.all((totals>0)&(totals<=1000000))

def make(x,t,genes):
 x=sparse.csr_matrix(x);x.eliminate_zeros();return ad.AnnData(x,obs=pd.DataFrame({'target_gene':[t]*x.shape[0]},index=[f'{t}:{i}' for i in range(x.shape[0])]),var=pd.DataFrame(index=genes))

def stat(a,cfg):
 _,means,m=pseudobulk_bulk_lognorm_with_moments(a,cfg.pert_col,bulk_target_sum=cfg.bulk_target_sum);return means[0],m.counts[0],m.sumsq[0],m.jk[0]

def assemble(stats):
 order=np.array(sorted(stats));records=[stats[t] for t in order];means=np.stack([r[0] for r in records]);m=GroupMoments(perts=order,counts=np.array([r[1] for r in records]),sumsq=np.array([r[2] for r in records]),jk=np.array([r[3] for r in records]));return {'bulk_lognorm':(order,means)},{'bulk_lognorm':m}

def score_iter(real,predictions,cfg,pred_control=None,release=None):
 assert cfg.pert_col=='target_gene' and cfg.de.backend=='pdex' and cfg.control==NT and cfg.target_sum==1e6 and cfg.de.fdr_scope=='per_pert'
 genes=real.var_names.to_numpy(str);labels=real.obs.target_gene.astype(str).to_numpy();targets=sorted(set(labels)-{NT});control=real[labels==NT].copy();validate(control.X);rstats={NT:stat(control,cfg)};pctrl=control if pred_control is None else pred_control;pstats={NT:stat(pctrl,cfg)};rframes=[];pframes=[];pending=[];digest=hashlib.sha256();seen=[]
 def flush():
  if not pending:return
  ts=[t for t,x in pending];blocks=[x for t,x in pending];labs=[t for t,x in pending for _ in range(x.shape[0])]
  if pred_control is not None:
   blocks.append(pred_control.X);labs.extend([NT]*pred_control.n_obs)
  pred=ad.AnnData(sparse.vstack(blocks,format='csr'),obs=pd.DataFrame({'target_gene':labs},index=['p'+str(i) for i in range(len(labs))]),var=real.var.copy());rr=real[np.isin(labels,ts+[NT])].copy()
  rd,pd_=run._materialize_de_sides(pred,rr,cfg=cfg,de_real=None,de_pred=None,real_fp=None,real_store=None,pred_store=None);rframes.append(rd);pframes.append(pd_);pending.clear();del pred,rr
 for i,(t,x) in enumerate(predictions):
  assert t in targets and t not in seen;seen.append(t);x=sparse.csr_matrix(x);validate(x,cfg.allow_fractional_counts);part=make(x,t,genes);pstats[t]=stat(part,cfg);rpart=real[labels==t].copy();validate(rpart.X);rstats[t]=stat(rpart,cfg);del rpart,part
  for value in [x.data,x.indices,x.indptr]:digest.update(value.tobytes())
  digest.update(t.encode());pending.append((t,x))
  if len(pending)==8:flush();print('Bounded reductions',len(seen),'/',len(targets),flush=True)
 flush();assert set(seen)==set(targets);del control,pctrl
 if release is not None:release()
 gc.collect();rb,rm=assemble(rstats);pb,pm=assemble(pstats);names=list(resolve_metrics(cfg.metrics,version=cfg.version)[0]);rd=pl.concat(rframes);pd_=pl.concat(pframes);prepared=run._prepare_de_cached(pd_,rd,cfg=cfg,real_store=None,pred_store=None,de_real_supplied=False,de_pred_supplied=False);rows=run.dispatch_de_metrics(names,prepared,cfg);rows.extend(run.dispatch_anndata_metrics(names,pb,rb,genes,cfg,comparator='bulk_lognorm',pred_moments=pm,real_moments=rm,driver='RSI bounded official reductions'));return pl.DataFrame(rows,schema=SCHEMA),digest.hexdigest()

def resident(pred,real,*,config,**kwargs):
 assert not kwargs
 labels=pred.obs.target_gene.astype(str).to_numpy();targets=sorted(set(labels)-{NT});pctrl=pred[labels==NT].copy() if NT in labels else None
 if config.control_source=='pred':assert pctrl is not None
 # Pseudobulk keeps submitted controls even under control_source=real; DE injection is official.
 it=((t,pred.X[labels==t]) for t in targets);return score_iter(real,it,config,pctrl)[0]

def baseline_iter(real,mean,cells=None):
 labels=real.obs.target_gene.astype(str).to_numpy()
 for t in sorted(set(labels)-{NT}):yield t,np.tile(np.asarray(mean,dtype=np.float32),(int((labels==t).sum()) if cells is None else int(cells),1))

def install_anchor():
 receipt=rt.json.loads(Path('artifacts/bounded_score_fixture.json').read_text());assert receipt['status']=='PASS bitwise equality' and receipt['driver_sha256']==rt.sha(__file__);anchor.compute_metrics=resident
