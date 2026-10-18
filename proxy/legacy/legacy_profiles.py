"""Freeze public common-response priors and sparse, hard-shrunk same-target residuals.
All calculations in this file are estimates, not biological verifications.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,re,pathlib
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad
ROOT=pathlib.Path(__file__).resolve().parent; A=ROOT/'artifacts'; NT='non-targeting'
NAMES=['H1','HepG2','Jurkat','K562','RPE1','K562_GWPS_bulk']

def align(source,dest):
 aliases=json.loads((ROOT/'artifacts/gene_aliases.json').read_text())
 ix=pd.Index(source).get_indexer(dest); ids={aliases[g]:i for i,g in enumerate(source) if g in aliases}
 for j in np.flatnonzero(ix<0):
  if dest[j] in aliases:ix[j]=ids.get(aliases[dest[j]],-1)
 return ix

def summarize(d,panel):
 genes=d['genes'];ts=d['targets'];M=d['means'].astype(float);n=d['n'];ctrl=M[np.flatnonzero(ts==NT)[0]]
 eligible=(ts!=NT)&(n>=30); matched=eligible&np.isin(ts,panel)
 select=matched if matched.sum()>=8 else eligible
 # Equal weight per target; omit each target's own transcript before averaging.
 rows=np.flatnonzero(select); total=M[rows].sum(0);den=np.full(len(genes),len(rows),float)
 gi=align(genes,ts[rows]);ok=gi>=0
 np.add.at(total,gi[ok],-M[rows[ok],gi[ok]]);np.add.at(den,gi[ok],-1)
 avg=total/np.maximum(den,1)
 cpm=ctrl/ctrl.sum()*1e6; pert=avg/avg.sum()*1e6
 # 5 CPM pseudocount shrinks unsupported low-abundance ratios.
 effect=np.clip(np.log2((pert+5)/(cpm+5)),-0.5,0.5)
 return dict(genes=genes,base=cpm,effect=effect,n_targets=int(len(rows)),selection='official-panel' if matched.sum()>=8 else 'all-eligible',selected=ts[rows].tolist(),n_matched=int(matched.sum()))

def common_effect(sources,genes,ctrl,exclude=None,scale=0.5):
 # One K562 vote: prefer broad panel-matched GWPS; exclude its family in K562 folds.
 chosen=[k for k in sources if k!='K562' and (exclude is None or not(k==exclude or k.startswith(exclude+'_')))]
 base=ctrl/ctrl.sum()*1e6; mapped=[];corr=[]
 for name in chosen:
  s=sources[name];ix=align(s['genes'],genes);ok=ix>=0
  mask=ok&(base>=5);src=np.zeros(len(genes));src[ok]=s['base'][ix[ok]]
  mask&=src>=5
  r=float(np.corrcoef(np.log1p(base[mask]),np.log1p(src[mask]))[0,1]) if mask.sum()>10 else 0
  corr.append(r);mapped.append((s,ix,ok))
 weights=np.exp((np.array(corr)-max(corr))/0.15);weights/=weights.sum()
 numerator=np.zeros(len(genes));den=np.zeros(len(genes))
 for w,(s,ix,ok) in zip(weights,mapped):
  numerator[ok]+=w*s['effect'][ix[ok]];den[ok]+=w
 # Missing sources contribute zero: do not amplify one-source genes by dividing out coverage.
 effect=scale*numerator
 return effect,dict(zip(chosen,map(float,weights)))

if __name__=='__main__':
 genes=pd.read_csv(ROOT.parent/'data/gene_names.csv').iloc[:,0].astype(str).to_numpy();targets=pd.read_csv(ROOT.parent/'data/pert_counts.csv').iloc[:,0].astype(str).to_numpy()
 sources={name:summarize(dict(np.load(A/'raw'/f'{name}.npz')),targets) for name in NAMES}
 np.savez_compressed(A/'source_shifts.npz',**{f'{k}_{field}':s[field] for k,s in sources.items() for field in ('genes','base','effect')})
 json.dump({k:{a:b for a,b in s.items() if a not in ('genes','base','effect')} for k,s in sources.items()},open(A/'source_selection.json','w'),indent=2)
 # Diagnostic LOCO shared-profile errors, heldout response never used in its prediction.
 diagnostics=[]
 for name in ('H1','HepG2','Jurkat','K562','RPE1'):
  s=sources[name];ix=align(s['genes'],s['genes']);c=s['base'];real=c*np.exp2(s['effect'])
  def prof(x):return np.log1p(5e4*x/x.sum())
  # The actual target average is saved too; pseudocount/clipping above regularizes prior only.
  d=dict(np.load(A/'raw'/f'{name}.npz'));M=d['means'].astype(float);ts=d['targets'];sel=np.isin(ts,s['selected']);avg=M[sel].mean(0)
  # Per-gene own-target omission, exactly as in summarize.
  cnt=np.full(len(c),sel.sum(),float);tot=M[sel].sum(0);gi=align(d['genes'],ts[sel]);ok=gi>=0
  np.add.at(tot,gi[ok],-M[sel][np.flatnonzero(ok),gi[ok]]);np.add.at(cnt,gi[ok],-1);avg=tot/np.maximum(cnt,1)
  control=prof(c);truth=prof(avg)
  for scale in (0,0.25,0.5,1):
   ef,w=common_effect(sources,s['genes'],c,exclude=name,scale=scale);pred=prof(c*np.exp2(ef));err=float(np.mean((pred-truth)**2))
   diagnostics.append(dict(context=name,scale=scale,shared_profile_squared_error=err,control_error=float(np.mean((control-truth)**2)),weights=json.dumps(w)))
 pd.DataFrame(diagnostics).to_csv(A/'common_loco.csv',index=False)
 # Robust residuals from raw GWPS constructs. Noise calibrated on 503 core-control constructs.
 a=ad.read_h5ad(ROOT.parent/'public_data/K562_GWPS_bulk.h5ad');X=np.asarray(a.X,dtype=float);n=a.obs.num_cells_filtered.fillna(0).to_numpy(float)
 labels=np.array([NT if NT in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$',s).group(1) for s in a.obs_names]);valid=n>=30;cc=(labels==NT)&valid&a.obs.core_control.to_numpy(bool)
 rawgenes=a.var.gene_name.astype(str).to_numpy()
 unique=np.array(sorted(set(rawgenes)));proj=sp.csr_matrix((np.ones(len(rawgenes)),(np.arange(len(rawgenes)),pd.Index(unique).get_indexer(rawgenes))),shape=(len(rawgenes),len(unique)))
 X=np.asarray(X@proj);rawgenes=unique
 cpm=X/X.sum(1,keepdims=True)*1e6;rawbase=np.average(X[cc],weights=n[cc],axis=0);base=rawbase/rawbase.sum()*1e6
 # Use pooled-count centering consistently for per-target means and controls.
 sd=np.maximum(np.sqrt(np.mean((np.log2((cpm[cc]+5)/(base+5))*np.sqrt(n[cc]/170)[:,None])**2,axis=0)),1e-3)
 s=sources['K562_GWPS_bulk'];assert len(set(rawgenes))==len(rawgenes)
 generic=s['effect'];assert np.array_equal(rawgenes,s['genes'])
 row=[];col=[];value=[];stats=[]
 for ti,t in enumerate(targets):
  sel=(labels==t)&valid
  if not sel.any():stats.append(dict(target=t,n_cells=0,constructs=0,robust_genes=0));continue
  mu=np.average(X[sel],weights=n[sel],axis=0);mu=mu/mu.sum()*1e6
  # Remove the same shared source shift used by A to avoid adding it twice.
  lfc=np.log2((mu+5)/(base+5))-generic;z=lfc/(sd*np.sqrt(170/n[sel].sum()))
  constructs=np.log2((cpm[sel]+5)/(base+5))-generic
  sign_ok=np.mean(np.sign(constructs)==np.sign(lfc)[None,:],axis=0)>=0.75
  robust=(np.abs(z)>=4)&(np.abs(lfc)>=0.25)&(base>=5)&sign_ok
  if sel.sum()==1:robust&=(n[sel].sum()>=100)&(np.abs(z)>=5)
  # Align onto official genes, aliases from annotation only.
  ix=align(rawgenes,genes);ok=ix>=0;good=np.zeros(len(genes),bool);good[ok]=robust[ix[ok]];good&=genes!=t
  jj=np.flatnonzero(good);vv=0.3*np.clip(lfc[ix[jj]],-1.5,1.5)
  row.extend([ti]*len(jj));col.extend(jj.tolist());value.extend(vv.tolist())
  stats.append(dict(target=t,n_cells=float(n[sel].sum()),constructs=int(sel.sum()),robust_genes=len(jj)))
 R=sp.csr_matrix((value,(row,col)),shape=(len(targets),len(genes)),dtype=np.float32);sp.save_npz(A/'robust_signatures.npz',R)
 pd.DataFrame(stats).to_csv(A/'signature_coverage.csv',index=False)
 np.savez_compressed(A/'axes.npz',genes=genes.astype(str),targets=targets.astype(str))
 controls={};summary={}
 for c in 'ABC':
  a=ad.read_h5ad(ROOT.parent/'data'/f'context_{c}.h5ad');assert np.array_equal(a.var_names,genes)
  ctrl=np.asarray(a.X.astype(float).mean(0)).ravel();effect,w=common_effect(sources,genes,ctrl)
  controls[f'{c}_ctrl']=ctrl;controls[f'{c}_shift']=effect
  summary[c]=dict(weights=w,shift_lfc_q=list(map(float,np.quantile(effect,[0,.1,.5,.9,1]))),mean_logprofile_displacement=float(np.mean((np.log1p(5e4*ctrl*np.exp2(effect)/(ctrl*np.exp2(effect)).sum())-np.log1p(5e4*ctrl/ctrl.sum()))**2)),zero_control_genes=int((ctrl==0).sum()))
 np.savez_compressed(A/'context_priors.npz',**controls)
 json.dump(dict(common_scale=0.5,signature_scale=0.3,signature_z=4,single_construct_z=5,signature_min_lfc=0.25,signature_min_source_cpm=5,signature_min_dest_cpm=5,knockdown=0.1,contexts=summary),open(A/'model.json','w'),indent=2)
 print(json.dumps(summary,indent=2));print(pd.DataFrame(stats).robust_genes.describe().to_string());print(pd.DataFrame(diagnostics).to_string(index=False))
