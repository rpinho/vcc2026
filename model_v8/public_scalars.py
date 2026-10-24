"""Estimate reference call budgets from public real cells, without challenge labels.
Actual-depth calls use exact MWU. Projection to 400 cells uses normal local power
and normal-mixture shrinkage; it is an uncertain scalar estimate, not a score.
"""
import runtime,gc,sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.special import ndtr
from de import prepare,Reference,stats,adjusted
R=runtime.ROOT
def posterior_effect(z):
 # Empirical-Bayes normal-means deconvolution. Directly scaling observed z
 # would turn sampling noise in small public groups into thousands of false
 # projected calls, including under a true null.
 variance=np.array([0.,.25,1.,4.,16.,64.,256.]);pi=np.ones(len(variance))/len(variance)
 density=np.exp(-z[:,None]**2/(2*(1+variance[None,:])))/np.sqrt(1+variance[None,:])
 for _ in range(60):
  weight=density*pi[None,:];weight/=np.maximum(weight.sum(1,keepdims=True),1e-300);pi=weight.mean(0)
 return z*np.sum(weight*(variance/(1+variance))[None,:],axis=1),pi.tolist()
genes=pd.read_csv(R.parent/'data/gene_names.csv').gene_name.to_numpy(str)
context_keep={c:np.load(R/f'artifacts/ref_{c}/mean.npy')>5 for c in 'ABC'}
out=[]
for ci,name in enumerate(['Jurkat','HepG2','RPE1']):
 runtime.guard();a=ad.read_h5ad(R.parent/f'public_data/{name}.h5ad',backed='r');print('PUBLIC',name,a.shape,list(a.obs.columns),flush=True)
 labels=a.obs['gene' if 'gene' in a.obs else 'target_gene'].astype(str).to_numpy()
 g=a.var['gene_name'].astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
 dest=pd.Index(genes).get_indexer(g);valid=dest>=0
 mapping=sparse.csr_matrix((np.ones(valid.sum(),np.int8),(np.flatnonzero(valid),dest[valid])),shape=(len(g),len(genes)))
 control=np.flatnonzero(np.char.startswith(labels.astype(str),'non-targeting')|np.isin(np.char.lower(labels.astype(str)),['ntc','control','ctrl','unperturbed']))
 assert len(control)>=400,(name,len(control))
 rng=np.random.default_rng(7344+ci);control=np.sort(rng.choice(control,min(len(control),18400),replace=False))
 reference_rng=np.random.default_rng(0);reference_rng.bit_generator.state=rng.bit_generator.state
 # Consume exactly prepare's real-cell selection, even on cache hits, so the
 # target/cell cohort is deterministic across fresh and resumed runs.
 rng.choice(len(control),400,replace=False)
 counts=pd.Series(labels).value_counts();eligible=sorted(t for t,n in counts.items() if n>=40 and not t.lower().startswith('non-targeting') and t.lower() not in ['ntc','control','ctrl','unperturbed'])
 chosen=np.sort(rng.choice(eligible,min(64,len(eligible)),replace=False));rowsets={}
 for t in chosen:
  rows=np.flatnonzero(labels==t);rowsets[t]=np.sort(rng.choice(rows,min(len(rows),400),replace=False))
 target_rows=np.unique(np.concatenate(list(rowsets.values())));requested=np.union1d(control,target_rows)
 print('LOAD SORTED COHORT',name,'cells',len(requested),'targets',len(chosen),flush=True)
 cohort=(sparse.csr_matrix(a.X[requested])@mapping).tocsr();a.file.close()
 refdir=R/'artifacts'/f'ref_public_{name}'
 if not (refdir/'cdf_scale.npy').exists():prepare(cohort[np.searchsorted(requested,control)],refdir,reference_rng)
 target_cohort=cohort[np.searchsorted(requested,target_rows)];del cohort;gc.collect()
 ref=Reference('public_'+name);records=[]
 for k,t in enumerate(chosen):
  runtime.guard();rows=rowsets[t];x=target_cohort[np.searchsorted(target_rows,rows)].toarray().astype(np.int32);d=x.sum(1,dtype=np.int64)
  z,p=stats(x,d,ref.values,ref.ptr,ref.n,ref.ties,ref.keep);q=adjusted(p,ref.keep)
  factor=np.sqrt((1/len(rows)+1/ref.n)/(1/400+1/18400));zlatent,mixture=posterior_effect(z[ref.keep]);z400=np.zeros(len(z));z400[ref.keep]=zlatent*factor;p400=2*ndtr(-np.abs(z400));q400=adjusted(p400,ref.keep)
  if t in genes:q[genes==t]=1;q400[genes==t]=1
  by_context={}
  for c in 'ABC':
   qc=adjusted(p400,ref.keep&context_keep[c]);qc[genes==t]=1.;by_context[c]=int((qc<.05).sum())
  records.append(dict(target=t,n=len(rows),calls_observed=int((q<.05).sum()),calls400_estimate=int((q400<.05).sum()),calls400_by_context=by_context,normal_mixture_weights=mixture))
  if k%8==0:print(name,k,records[-1],flush=True)
 vals=np.array([v['calls400_estimate'] for v in records]);result=dict(source=name,reference_cells=ref.n,n_targets=len(records),cohort_cells=len(requested),selection_seed=7344+ci,call_budget_median=float(np.median(vals)),call_budget_quantiles=np.quantile(vals,[.1,.25,.5,.75,.9]).tolist(),records=records,method='Local-normal n=400 power projection of exact MWU after empirical-Bayes normal-mixture shrinkage. Variance grid [0,.25,1,4,16,64,256], 60 EM steps; posterior mean effect, BH on reference-mean>5 CPM genes. No challenge treated labels. Uncertainty includes shrinkage/power projection and cross-context transport.')
 # The official panel has nonzero cross-source response support for every
 # target. Use a responder-conditional prior rather than treating all 2,394
 # TF-screen targets (including failed/no-response knockdowns) as that panel.
 # This rule is fixed before any of these estimates or challenge predictions
 # are inspected. Retain the unconditional budget for a visible sensitivity check.
 responders=[v for v in records if v['calls_observed']>=5]
 if len(responders)<8:responders=[v for v in records if v['calls_observed']>0]
 assert len(responders)>=4,('insufficient public responders for N_conf estimate',name)
 result['responder_targets']=len(responders);result['conditioning']='At least 5 observed source-depth BH calls; if fewer than 8 targets, at least 1 call. Official baseline has source-sign support on every target; this is a responder prior, not known hidden truth.'
 result['context_budgets']={c:dict(median=float(np.median([v['calls400_by_context'][c] for v in responders])),unconditional_median=float(np.median([v['calls400_by_context'][c] for v in records])),quantiles=np.quantile([v['calls400_by_context'][c] for v in responders],[.1,.25,.5,.75,.9]).tolist()) for c in 'ABC'}
 out.append(result);runtime.save(R/'artifacts/public_scalars.json',dict(status='ESTIMATES, not hidden-truth measurements',sources=out));print('SCALAR',name,result['call_budget_median'],result['context_budgets'],flush=True)
 del a,ref,target_cohort;gc.collect()
