"""Multi-family, two-summary transfer with post-pooling adaptive shrinkage.
AtlasShift and Sidechain kernels are MIT; see vendor/*_LICENSE.
"""
import runtime,os
import numpy as np,pandas as pd
from pathlib import Path
from scipy import sparse
from vendor import atlas_model as atlas
from vendor.adaptive_shrink import adaptive_shrink
LOCAL_BANKS=runtime.ROOT/'artifacts/banks'
SHARED_BANKS=runtime.ROOT.parent/'public_data/model_v4_signature_banks'
BANKS=LOCAL_BANKS if (LOCAL_BANKS/'K562_GWPS_raw.npz').exists() or not SHARED_BANKS.exists() else SHARED_BANKS
WEIGHTS={'K562_GWPS_raw':2,'K562_GWPS_bulk':2,'H1':2,'HCT116':1,'HEK293T':1,'HepG2':1,'Jurkat':1,'K562':1,'RPE1':1,'CD4':.5}

def load_banks(exclude=None,expanded=True,validation_extra=False):
 raw=(BANKS/'K562_GWPS_raw.npz').exists()
 names=['K562_GWPS_raw' if raw else 'K562_GWPS_bulk','H1']
 if expanded:names+=['HCT116','HEK293T']
 if validation_extra:names+=['HepG2','Jurkat','RPE1','K562']
 names=[n for n in names if (BANKS/f'{n}.npz').exists() and not (exclude and (n==exclude or n.startswith(exclude+'_')))]
 return names

def aligned_bank(name,targets,genes,space,center=True):
 path=BANKS/f'{name}.npz';src=atlas.load_xatlas(path,prior_counts=100000)
 full,mask=atlas.aligned_effect(src,src.targets,genes,space=space,common_subtract=int(center))
 # Delta-method sampling errors are propagated before inverse-variance pooling.
 # CPM-contrast variances are source-specific estimates, not cross-line uncertainty.
 with np.load(path) as d:
  dv=d['delta_cpm_variance'].astype(float)
  ctrl=d['matched_control_mean_cpm'].astype(float);pert=d['target_mean_cpm'].astype(float)
  pseudocount=1 if space=='log2fc' else 20
  # Use the conservative smaller expression denominator when only contrast variance is retained.
  if 'target_mean_cpm_variance' in d:
   var=d['target_mean_cpm_variance'].astype(float)/(pert+pseudocount)**2+d['matched_control_mean_cpm_variance'].astype(float)/(ctrl+pseudocount)**2
  else:var=dv/np.maximum(np.minimum(ctrl,pert)+pseudocount,1e-12)**2
  if space=='log2fc':var/=np.log(2)**2
  # A small numerical clamp matches Sidechain's historically released IV path.
  var=np.maximum(var,1e-6)
 row=pd.Index(src.targets).get_indexer(targets);col=pd.Index(src.genes).get_indexer(genes);r=np.flatnonzero(row>=0);c=np.flatnonzero(col>=0)
 eff=np.zeros((len(targets),len(genes)),np.float32);m=np.zeros_like(eff,bool);v=np.full(eff.shape,np.inf)
 eff[r]=full[row[r]];m[r]=mask[row[r]];v[np.ix_(r,c)]=var[np.ix_(row[r],col[c])];v[~m]=np.inf
 return eff,m,v

def cd4_family(targets,genes,space,center=True):
 """One lineage vote: preserve condition banks, average correlated conditions.
 Published DE LFC/SE overrides pseudobulk approximations on its measured axis.
 """
 names=[p.stem for p in BANKS.glob('CD4_*.npz') if p.stem!='CD4_DE_statistics'];values=[];masks=[];variances=[]
 for name in sorted(names):
  a,m,v=aligned_bank(name,targets,genes,space,center);values.append(a);masks.append(m);variances.append(v)
 if not values:return None
 a=np.stack(values);m=np.stack(masks);v=np.stack(variances);n=m.sum(0)
 effect=np.divide(np.where(m,a,0).sum(0),n,out=np.zeros(a.shape[1:]),where=n>0)
 # Treat conditions as correlated: arithmetic mean of variances, not an unjustified 3x precision bonus.
 var=np.divide(np.where(m,v,0).sum(0),n,out=np.full(a.shape[1:],np.inf),where=n>0)
 path=BANKS/'CD4_DE_statistics.npz'
 if space=='log2fc' and path.exists():
  de,ok=atlas.aligned_cd4(path,targets,genes,common_subtract=int(center),center_scope='source')
  with np.load(path) as d:
   ri=pd.Index(d['targets']).get_indexer(targets);ci=pd.Index(d['genes']).get_indexer(genes);rr=np.flatnonzero(ri>=0);cc=np.flatnonzero(ci>=0)
   vv=d['lfcSE'].astype(float)**2;av=d['available']&d['quality_pass']&(d['n_cells']>=20)
   num=np.where(av[:,:,None],vv,0).sum(0);den=av.sum(0)[:,None];vs=np.divide(num,den,out=np.full(num.shape,np.inf),where=den>0)
   aligned=np.full(effect.shape,np.inf);aligned[np.ix_(rr,cc)]=vs[np.ix_(ri[rr],ci[cc])]
  effect=np.where(ok,de,effect);var=np.where(ok,np.maximum(aligned,1e-6),var);n=np.maximum(n,ok)
 return effect,n>0,var

def _ash_task(item):
 x,var=item;stats={};out=adaptive_shrink(x,var,stats=stats);return out,stats

def _pool(targets,genes,names,space,method='atlas',shrink='none',use_cd4=True,stats=None):
 entries=[(n,aligned_bank(n,targets,genes,space)) for n in names]
 cd=cd4_family(targets,genes,space) if use_cd4 else None
 if cd is not None:entries.append(('CD4',cd))
 num=np.zeros((len(targets),len(genes)));den=np.zeros_like(num);coverage=np.zeros_like(num)
 for name,(effect,mask,var) in entries:
  if method=='iv':w=np.where(mask&np.isfinite(var),1/np.maximum(var,1e-6),0)
  else:w=WEIGHTS.get(name,1)*mask
  num+=w*effect;den+=w;coverage+=mask
 result=np.divide(num,den,out=np.zeros_like(num),where=den>0);variance=np.divide(1,den,out=np.full_like(den,np.inf),where=den>0)
 if shrink!='none':
  assert method=='iv','Pooled SE valid only on unmodified inverse-variance path'
  workers=min(3,int(os.environ.get('VCC_ASH_WORKERS','1')))
  if shrink=='ash' and workers>1 and len(targets)>32:
   import multiprocessing
   from concurrent.futures import ProcessPoolExecutor
   items=[]
   for i,t in enumerate(targets):
    usable=(den[i]>0)&(genes!=t);items.append((result[i,usable],variance[i,usable]))
   with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as executor:
    for i,(values,local_stats) in enumerate(executor.map(_ash_task,items)):
     usable=(den[i]>0)&(genes!=targets[i]);result[i,usable]=values
     if stats is not None:
      for k,value in local_stats.items():stats[k]=stats.get(k,0)+value
  else:
   for i,t in enumerate(targets):
    usable=(den[i]>0)&(genes!=t)
    if shrink=='ash':result[i,usable]=adaptive_shrink(result[i,usable],variance[i,usable],stats=stats)
    elif shrink=='js':
     k=int(usable.sum());energy=float(np.sum(result[i,usable]**2));noise=float(np.sum(variance[i,usable]));factor=max(0,1-(max(k-2,0)/max(k,1))*noise/max(energy,1e-30));result[i,usable]*=factor
    else:raise ValueError(shrink)

 return result,coverage,variance

_POOL_CACHE={}
def pool(targets,genes,names,space,method='atlas',shrink='none',use_cd4=True,stats=None):
 key=(tuple(targets),tuple(genes),tuple(names),space,method,shrink,use_cd4)
 if key not in _POOL_CACHE:
  if len(_POOL_CACHE)>=2:_POOL_CACHE.clear()
  _POOL_CACHE[key]=_pool(targets,genes,names,space,method,shrink,use_cd4,stats)
 return _POOL_CACHE[key]

def desired_cpm(control,effects,amplitude=.6):
 """Invert the same 1-CPM pseudocount used in source LFC estimation."""
 effect=np.clip(amplitude*effects,-3,3);p=np.maximum((control[None,:]*1e6+1)*np.exp2(effect)-1,0);p/=p.sum(1,keepdims=True);return p

def heads(targets,genes,mean,bulk,variant,exclude=None,stats=None,de_amplitude=.6):
 expanded=not variant.startswith('old_') and variant!='anchor_restricted';names=load_banks(exclude,expanded)
 single='_single_' in variant;method='iv' if variant.startswith('iv_') else 'atlas';shrink='ash' if '_ash' in variant else 'js' if '_js' in variant else 'none'
 if variant.startswith('anchor'):
  # Unmodified released fusion with the available, documented source subsets.
  sources=[atlas.load_xatlas(BANKS/f'{n}.npz',prior_counts=100000) for n in names]
  effects={s:atlas.fuse_source_centered(sources,[WEIGHTS[n] for n in names],targets,genes,space=s,common_subtract=1) for s in ['log2fc','bulk_delta']}
  cd,available=atlas.aligned_cd4(BANKS/'CD4_DE_statistics.npz',targets,genes,common_subtract=1,center_scope='source')
  result=[];cov=[]
  for space,control,amp in [('log2fc',mean,.6),('bulk_delta',bulk,.3)]:
   e,m=effects['log2fc' if single else space]
   if variant!='anchor_restricted':e,m=atlas.add_cd4_family(e,m,cd,available,.5,space=('log2fc' if single else space),control_probability=control)
   result.append(atlas.desired_mean(control,e,space=('log2fc' if single else space),amplitude=(de_amplitude if space=='log2fc' else (.6 if single else amp)),clip=3));cov.append(m)
 else:
  fc,fc_cov,_=pool(targets,genes,names,'log2fc',method,shrink,use_cd4=expanded,stats=stats)
  p=desired_cpm(mean,fc) if ('smooth' in variant or method=='iv') else atlas.desired_mean(mean,fc,space='log2fc',amplitude=.6,clip=3)
  if single:b=atlas.desired_mean(bulk,fc,space='log2fc',amplitude=.6,clip=3);bc=fc_cov
  else:
   bd,bc,_=pool(targets,genes,names,'bulk_delta',method,shrink,use_cd4=expanded,stats=stats);b=atlas.desired_mean(bulk,bd,space='bulk_delta',amplitude=.3,clip=3)
  result=[p,b];cov=[fc_cov,bc]
 if (variant.startswith('anchor') and variant!='anchor_nopromoter') or variant.endswith('_guide'):
  result=[promoters(q,c,targets,genes,coverage=cov[i],guide_aware=not variant.startswith('anchor')) for i,(q,c) in enumerate(zip(result,[mean,bulk]))]
 return result,cov

def promoters(prob,control,targets,genes,coverage,guide_aware=True):
 path=runtime.ROOT/'artifacts/promoter_pairs.csv'
 if not guide_aware:return atlas.apply_promoter_prior(prob,control,targets,genes,path,.15)[0]
 pairs=pd.read_csv(path);out=prob.copy();ti={t:i for i,t in enumerate(targets)};gi={g:i for i,g in enumerate(genes)}
 for row in pairs.itertuples():
  if row.target not in ti or row.neighbor not in gi:continue
  i,j=ti[row.target],gi[row.neighbor]
  # Fill unsupported readouts, rather than replacing effects measured in multiple lineages.
  if coverage[i,j]>=2:continue
  # Real guide positions are used only with an annotated matching assembly and promoter support.
  distance=float(row.guide_distance) if row.guide_supported else float(row.distance)
  if not row.guide_supported and not row.divergent:continue
  ramp=np.clip(np.log(max(distance,500)/500)/np.log(10),0,1);remain=.15+.85*ramp
  # Unsupported guide coordinates use a conservative half-strength divergent-TSS prior.
  if not row.guide_supported:remain=np.sqrt(remain)
  out[i,j]=min(out[i,j],control[j]*remain)
 out/=out.sum(1,keepdims=True);return out

def template(a,genes,cells=400,pool_k=4,seed=20260910):
 assert np.array_equal(a.var_names,genes)
 x=sparse.csr_matrix(a.X,dtype=float);depth=np.asarray(x.sum(1)).ravel();assert np.all(depth>0)
 p=np.asarray(x.T@(1/depth)).ravel()/x.shape[0];b=np.asarray(x.sum(0)).ravel();b/=b.sum()
 assert len(depth)>=cells*pool_k
 selected=np.random.default_rng(seed).choice(len(depth),cells*pool_k,replace=False);selected=selected[np.argsort(depth[selected],kind='stable')]
 q=x[selected].toarray()/depth[selected,None];q=q.reshape(cells,pool_k,len(genes)).mean(1);d=np.rint(depth[selected].reshape(cells,pool_k).mean(1)).astype(np.int64)
 return q,d,p,b

def mean_counts(q,p,b,depths,seed):
 """One-moment pooled-profile fitting, same template/depths/rounding as dual arm."""
 x=q.astype(float).copy();x=.999*x+.001*b;desired=depths.sum()*b
 for it in range(100):
  col=depths@x;x*=np.divide(desired,col,out=np.zeros_like(col),where=col>0);x/=np.maximum(x.sum(1,keepdims=True),1e-30)
  if it%5==4 and np.abs(depths@x/depths.sum()-b).sum()<2e-4:break
 exp=x*depths[:,None];integer=np.floor(exp).astype(np.int32);rng=np.random.default_rng(seed)
 for i in range(len(x)):
  rem=int(depths[i]-integer[i].sum())
  if rem:
   cum=np.cumsum(exp[i]-integer[i]);cum*=rem/cum[-1];loc=np.searchsorted(cum,np.arange(rem)+rng.random(),side='right');np.add.at(integer[i],loc,1)
 return atlas.repair_bulk_totals(integer,exp)

def emit(q,p,b,depths,seed,variant):
 if variant.endswith('_mean'):return mean_counts(q,p,b,depths,seed)
 return atlas.dual_moment_counts(q,p,b,depths=depths,seed=seed)
