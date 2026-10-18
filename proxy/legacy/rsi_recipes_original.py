"""Port archived model recipes to LOCO without reading destination responses.
All reused numerical model kernels remain in vendor, with source hashes in provenance.
External AtlasShift/Sidechain leaderboard labels are NOT assigned to our adapted recipes.
"""
import runtime as rt
import sys,json,functools
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import legacy_pipeline as pipe
sys.modules['pipeline']=pipe;pipe.ROOT=rt.ROOT
import legacy_v2 as v2
import legacy_v3 as v3
import legacy_profiles as profiles
profiles.ROOT=rt.ROOT
import v4_model as v4
from atlas_model import seed_for
NT='non-targeting';SEALED={'HEK293T','RPE1'}
MANIFEST=rt.json.loads((rt.ROOT/'artifacts/tasks/manifest.json').read_text());H1_SEALED=set(MANIFEST['H1_holdout']['targets'])

def allowed(name,destination):return name.split('_')[0] not in SEALED and not (name==destination or name.startswith(destination+'_'))

@functools.lru_cache(maxsize=6)
def source(name,profile='base'):
 folder='cd4_sources' if profile=='CD4' else 'sources';d=dict(np.load(rt.ROOT/f'artifacts/{folder}/{name}.npz'));d['name']=name;d['family']='K562' if name.startswith('K562') else name;d['gene_ids']=np.array([pipe.annotation_aliases().get(g,'') for g in d['genes']]);d['target_ids']=np.array([pipe.annotation_aliases().get(t,'') for t in d['targets']]);return d

def sources(destination):
 out=[]
 for name in ['H1','HepG2','Jurkat','K562','RPE1','K562_GWPS_bulk']:
  if not allowed(name,destination):continue
  d=source(name,'CD4' if destination=='CD4' else 'base').copy();keep=~np.isin(d['targets'],list(H1_SEALED)) if name=='H1' else np.ones(len(d['targets']),bool)
  for k in ['targets','target_ids','lfc','fine','n_cells']:d[k]=d[k][keep]
  d['lfc']=d['lfc'].astype(float)*float(d['lfc_scale']);out.append(d)
 return out

def control_summary(ctrl):
 x=sparse.csr_matrix(ctrl.X,dtype=float);dep=np.asarray(x.sum(1)).ravel();cm=np.asarray(x.mean(0)).ravel();norm=x.multiply((1e4/dep)[:,None]);avg=np.asarray(norm.mean(0)).ravel();var=np.maximum(np.asarray(norm.multiply(norm).mean(0)).ravel()-avg**2,0);theta=np.clip(avg**2/np.maximum(var-avg*np.mean(1e4/dep),1e-4),.1,1000)
 genes=ctrl.var_names.to_numpy(str);return dict(genes=genes,gene_ids=np.array([pipe.annotation_aliases().get(g,'') for g in genes]),mean=cm,libs=dep,theta=theta)

def generic_sources(destination,targets):
 out={}
 for name in ['H1','HepG2','Jurkat','K562','RPE1','K562_GWPS_bulk']:
  if not allowed(name,destination):continue
  d=source(name,'CD4' if destination=='CD4' else 'base');labs=d['raw_labels'];counts=np.r_[d['n_controls'],d['n_cells']];use=np.ones(len(labs),bool)
  if name=='H1':use&=~np.isin(labs,list(H1_SEALED))
  data=dict(genes=d['genes'],targets=labs[use],means=d['raw_means'][use],n=counts[use]);out[name]=profiles.summarize(data,set(targets))
 return out

@functools.lru_cache(maxsize=2)
def robust_bulk(targets,genes,generic_genes,generic_effect):
 # Exact archived robust-residual rule, widened target panel; constructed only from K562.
 a=ad.read_h5ad(rt.ROOT.parent/'public_data/K562_GWPS_bulk.h5ad');X=np.asarray(a.X,dtype=float);n=a.obs.num_cells_filtered.fillna(0).to_numpy(float);import re
 labels=np.array([NT if NT in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$',s).group(1) for s in a.obs_names]);valid=n>=30;cc=(labels==NT)&valid&a.obs.core_control.to_numpy(bool);rawgenes=a.var.gene_name.to_numpy(str);unique=np.array(sorted(set(rawgenes)));proj=sparse.csr_matrix((np.ones(len(rawgenes)),(np.arange(len(rawgenes)),pd.Index(unique).get_indexer(rawgenes))),shape=(len(rawgenes),len(unique)));X=np.asarray(X@proj);rawgenes=unique;assert tuple(rawgenes)==generic_genes
 cpm=X/X.sum(1,keepdims=True)*1e6;rawbase=np.average(X[cc],weights=n[cc],axis=0);base=rawbase/rawbase.sum()*1e6;sd=np.maximum(np.sqrt(np.mean((np.log2((cpm[cc]+5)/(base+5))*np.sqrt(n[cc]/170)[:,None])**2,axis=0)),1e-3);res=np.zeros((len(targets),len(genes)),np.float32);ix=profiles.align(rawgenes,np.array(genes));ok=ix>=0;generic=np.array(generic_effect)
 for i,t in enumerate(targets):
  sel=(labels==t)&valid
  if not sel.any():continue
  mu=np.average(X[sel],weights=n[sel],axis=0);mu=mu/mu.sum()*1e6;lfc=np.log2((mu+5)/(base+5))-generic;z=lfc/(sd*np.sqrt(170/n[sel].sum()));constructs=np.log2((cpm[sel]+5)/(base+5))-generic;signok=np.mean(np.sign(constructs)==np.sign(lfc)[None,:],axis=0)>=.75;good=(np.abs(z)>=4)&(np.abs(lfc)>=.25)&(base>=5)&signok
  if sel.sum()==1:good&=(n[sel].sum()>=100)&(np.abs(z)>=5)
  mapped=np.zeros(len(genes),bool);mapped[ok]=good[ix[ok]];mapped&=np.array(genes)!=t;res[i,mapped]=.3*np.clip(lfc[ix[mapped]],-1.5,1.5)
 return res

def predict(ctrl,targets,destination,recipe,seed=42):
 """Return one target's count matrix at a time; caller streams to disk."""
 s=control_summary(ctrl);s['panel']=np.array(targets);rng=np.random.default_rng(seed);X=sparse.csr_matrix(ctrl.X);genes=s['genes'];n=400
 if recipe in ['v1','v2b']:
  src=sources(destination)
  # v2b's selected recipe is postfine100, not generic t100.
  v2._FINE={d['name']:dict(lfc=d['fine'],genes=d['genes'],targets=d['targets'],baseline=d['baseline'],n_cells=d['n_cells'],n_controls=d['n_controls']) for d in src};v2._POST_STATE=None
  for t in targets:
   if recipe=='v1':mu=pipe.predict_mean(s,t,src,'transfer');x=pipe.sample_counts(mu,s['libs'],s['theta'],n,rng)
   else:mu=v2.predict(s,t,src,'postfine100');x=v2.posterior_transport(mu,X,s['mean'],s['theta'],n,rng)
   yield t,x
 elif recipe=='v3b':
  common=generic_sources(destination,targets);shift,_=profiles.common_effect(common,genes,s['mean'],exclude=destination);res=np.zeros((len(targets),len(genes)))
  if 'K562_GWPS_bulk' in common:
   z=common['K562_GWPS_bulk'];res=robust_bulk(tuple(targets),tuple(genes),tuple(z['genes']),tuple(z['effect']))
  for i,t in enumerate(targets):yield t,v3.emit(X,s['mean'],v3.mean_vector(s['mean'],shift,genes,t,res[i]),0 if seed==42 else seed,i)
 elif recipe in ['v4a','v4b','atlas_adapted']:
  # Original model kernels; sealed destination families are unavailable to SEARCH.
  old=v4.load_banks;oldpath=v4.BANKS;oldcd4=v4.cd4_family;oldalign=v4.atlas.aligned_cd4
  if destination=='CD4':
   v4.BANKS=rt.ROOT/'artifacts/cd4_banks';v4.cd4_family=lambda *args,**kwargs:None
   def no_cd4(path,targets,genes,*args,**kwargs):return np.zeros((len(targets),len(genes)),np.float32),np.zeros((len(targets),len(genes)),bool)
   v4.atlas.aligned_cd4=no_cd4
  def banks(exclude=None,expanded=True,validation_extra=False):return [name for name in old(exclude,expanded,validation_extra) if allowed(name,destination) and name!='H1']
  v4.load_banks=banks
  try:
   effective_seed=20260910 if seed==42 else seed;q,depth,mean,bulk=v4.template(ctrl,genes,seed=effective_seed);variant={'v4a':'anchor_single_mean','v4b':'iv_ash_dual_guide','atlas_adapted':'anchor'}[recipe];(p,b),_=v4.heads(np.array(targets),genes,mean,bulk,variant,exclude=destination,stats={})
   for i,t in enumerate(targets):yield t,sparse.csr_matrix(v4.emit(q,p[i],b[i],depth,seed_for(destination+':'+t,effective_seed),variant))
  finally:v4.load_banks=old;v4.BANKS=oldpath;v4.cd4_family=oldcd4;v4.atlas.aligned_cd4=oldalign
 else:raise ValueError(recipe)
