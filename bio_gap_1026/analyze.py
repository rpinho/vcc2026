import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,pandas as pd,time
from scipy.stats import norm,spearmanr
start=time.monotonic()
z=np.load(rt.ROOT/'results/baseline.npz');targets=z['targets'];genes=z['genes'];Y=z['truth'].astype(float);P=z['realized'].astype(float);raw=z['raw'];own=z['own_idx'];n,g=Y.shape
rendered_path=rt.ROOT/'scratch/rendered_bulk.npz'
if rendered_path.exists():
    with np.load(rendered_path) as rz:
        assert np.array_equal(rz['targets'],targets) and np.array_equal(rz['genes'],genes);P=rz['response'].astype(float)
meta=rt.json.loads((rt.ROOT/'results/training_metadata.json').read_text())
ess=z['essential'];kd=z['source_kd'];bulk=z['bulk']
def masked(a):
    a=a.copy();a[np.arange(n)[own>=0],own[own>=0]]=0;return a
Y=masked(Y);P=masked(P)
def moments(a):
    a=masked(a);dot=(a*Y).sum(1);pn=(a*a).sum(1);yn=(Y*Y).sum(1);err=((a-Y)**2).sum(1)
    alpha=dot/np.maximum(pn,1e-30)
    return dict(dot=dot,pred_norm2=pn,truth_norm2=yn,sse=err,relative_error=err/yn,cosine=dot/np.maximum(np.sqrt(pn*yn),1e-30),best_scalar=alpha,norm_ratio=np.sqrt(pn/yn),direction_error=(yn-2*alpha*dot+alpha**2*pn)/yn,magnitude_error=pn*(1-alpha)**2/yn)
base=moments(P);E0=base['sse'].sum();T0=base['truth_norm2'].sum()
def project(a):
    q=np.maximum(np.expm1(np.log1p(50000*bulk)[None,:]+a),0);q/=q.sum(1,keepdims=True)
    return masked(np.log1p(50000*q)-np.log1p(50000*bulk))
def loo_scale(group):
    out=np.empty_like(P);alphas=[]
    for i in range(n):
        tr=(np.arange(n)!=i)&(group==group[i])
        if tr.sum()<8:tr=np.arange(n)!=i
        a=max(base['dot'][tr].sum()/base['pred_norm2'][tr].sum(),0);out[i]=P[i]*a;alphas.append(a)
    return project(out),alphas
preds={'B0':P,'pool':masked(z['pool'].astype(float)),'context_weighted':masked(z['weighted'].astype(float))}
coeff={}
for name,group in [('global_scale',np.zeros(n)),('essential_scale',ess),('source_KD_scale',np.digitize(kd,[.5,.8]))]:preds[name],coeff[name]=loo_scale(group)
# Response-size-only calibration. Weighted scalar regression minimizes training SSE;
# source KD/essentiality/predicted response energy are legal inference features.
X=np.column_stack([np.ones(n),ess.astype(float),np.clip(kd,-1,1),np.log1p(np.sqrt(base['pred_norm2']))])
out=np.empty_like(P);co=[]
for i in range(n):
    tr=np.arange(n)!=i;xm=X[tr,1:].mean(0);sd=X[tr,1:].std(0);sd[sd<1e-8]=1
    A=np.column_stack([np.ones(tr.sum()),(X[tr,1:]-xm)/sd]);q=np.r_[1,(X[i,1:]-xm)/sd]
    w=base['pred_norm2'][tr];ridge=10*w.mean()*np.diag([0,1,1,1]);rhs=A.T @ base['dot'][tr]
    beta=np.linalg.solve(A.T@(w[:,None]*A)+ridge,rhs);alpha=max(float(q @ beta),0);out[i]=alpha*P[i];co.append(alpha)
preds['response_size_scale']=project(out);coeff['response_size_scale']=co
residual=Y-P;preds['LOO_bias_correction']=project(P+(residual.sum(0)[None,:]-residual)/(n-1))
# Mean-only negative control: 167 training targets in every fold.
preds['LOO_mean']=project((Y.sum(0)[None,:]-Y)/(n-1))
family_names=[k[7:] for k in z.files if k.startswith('family_')]
Fs=np.stack([masked(z['family_'+f].astype(float)) for f in family_names])
avail=np.isfinite(Fs).all(2);Fs=np.where(avail[:,:,None],Fs,P[None,:,:]);R=Fs-P
K=len(family_names);Q=np.einsum('kig,lig->ikl',R,R);rhs=np.einsum('kig,ig->ik',R,Y-P)
out=np.empty_like(P);weights=[]
for i in range(n):
    tr=np.arange(n)!=i;gram=Q[tr].sum(0);reg=max(np.trace(gram)/K,1e-12)
    w=np.linalg.solve(gram+reg*np.eye(K),rhs[tr].sum(0));out[i]=P[i]+w @ R[:,i];weights.append(w.tolist())
preds['LOO_family_transfer']=project(out);coeff['family_residual_weights']=weights
# Oracles: diagnostic ceilings, never selectable predictions.
preds['oracle_scalar']=P*np.maximum(base['best_scalar'],0)[:,None]
family_sse=np.sum((Fs-Y[None])**2,2);family_sse[~avail]=np.inf
oracle_family=np.minimum(np.min(family_sse,0),base['sse'])
def pds(a):
    # Fixed panel exclusion and n-1 denominator, vcc2026 preset; equivalent to
    # installed discrimination_score and independently verified in verify.py.
    keep=~np.isin(genes,targets);a=a[:,keep];y=Y[:,keep]
    c=a @ y.T/np.maximum(np.linalg.norm(a,axis=1)[:,None]*np.linalg.norm(y,axis=1)[None,:],1e-30)
    rank=np.sum(c>np.diag(c)[:,None],1)+.5*(np.sum(c==np.diag(c)[:,None],1)-1)
    return 1-rank/(n-1)
summary={};per=[]
for name,a in preds.items():
    st=moments(a);pdsv=pds(a)
    summary[name]=dict(sse=float(st['sse'].sum()),error_ratio_zero=float(st['sse'].sum()/T0),SSE_reduction_vs_B0=float(1-st['sse'].sum()/E0),macro_relative_error=float(st['relative_error'].mean()),mean_cosine=float(st['cosine'].mean()),PDS168=float(pdsv.mean()),centered_energy=float(np.sum((a-a.mean(0))**2)))
    for i,t in enumerate(targets):per.append(dict(target=t,method=name,sse=st['sse'][i],relative_error=st['relative_error'][i],cosine=st['cosine'][i],pds=pdsv[i]))
    print(name,summary[name],flush=True)
pd.DataFrame(per).to_csv(rt.ROOT/'results/fix_loo_per_target.csv',index=False)
# Moment-based n_conf diagnostic: full training raw counts no longer retained.
# This is a BH-adjusted two-sample normal test on arithmetic CPM, NOT Wilcoxon.
src=rt.REPO/'model_sidechain/approx/sources/H1';md=rt.json.loads((src/'meta.json').read_text());gi=pd.Index(md['genes']).get_indexer(genes);ti=pd.Index(md['labels']).get_indexer(targets);ci=md['labels'].index('non-targeting')
assert (ti>=0).all() and (gi>=0).all()
N=np.load(src/'n_cells.npy',mmap_mode='r');S=np.load(src/'cpm_sum.npy',mmap_mode='r');S2=np.load(src/'cpm_sq_sum.npy',mmap_mode='r')
nn=N[ti].astype(float);cn=float(N[ci]);sm=np.asarray(S[ti][:,gi],float);ss=np.asarray(S2[ti][:,gi],float);cm=np.asarray(S[ci][gi],float)/cn;cv=np.maximum((np.asarray(S2[ci][gi],float)-cn*cm**2)/(cn-1),0)
tm=sm/nn[:,None];tv=np.maximum((ss-nn[:,None]*tm**2)/(nn[:,None]-1),0);se=np.sqrt(tv/nn[:,None]+cv/cn)
pv=2*norm.sf(np.abs(tm-cm)/np.maximum(se,1e-30));pv[np.arange(n)[own>=0],own[own>=0]]=1
nconf=[]
for row in pv:
    o=np.sort(row);passed=np.flatnonzero(o<.05*np.arange(1,g+1)/g);nconf.append(int(passed[-1]+1) if len(passed) else 0)
strata=[]
def stratum(kind,label,ind):
    if not np.any(ind):
        strata.append(dict(kind=kind,label=label,n=0));return
    d=base['direction_error'][ind]*base['truth_norm2'][ind];m=base['magnitude_error'][ind]*base['truth_norm2'][ind]
    strata.append(dict(kind=kind,label=label,n=int(ind.sum()),error_mass_pct=float(100*base['sse'][ind].sum()/E0),error_ratio_zero=float(base['sse'][ind].sum()/base['truth_norm2'][ind].sum()),macro_relative_error=float(base['relative_error'][ind].mean()),cosine_mean=float(base['cosine'][ind].mean()),best_scalar_median=float(np.median(base['best_scalar'][ind])),norm_ratio_median=float(np.median(base['norm_ratio'][ind])),direction_pct_of_SSE=float(100*d.sum()/base['sse'][ind].sum()),magnitude_pct_of_SSE=float(100*m.sum()/base['sse'][ind].sum())))
direct=np.array([bool(a) for a in meta['coverage']]);stratum('all','all',np.ones(n,bool))
for value in [True,False]:stratum('coverage','direct' if value else 'none',direct==value);stratum('essential_proxy','screen_member' if value else 'other',ess==value)
for f in family_names:stratum('family_presence',f,np.array([any(a['family']==f for a in row) for row in meta['coverage']]))
for j,label in enumerate(['<50%','50–79%','>=80%']):stratum('source_knockdown',label,np.digitize(kd,[.5,.8])==j)
truthkd=np.array([v if v is not None else np.nan for v in meta['observed_H1_knockdown']])
for j,label in enumerate(['<50%','50–79%','>=80%']):stratum('observed_H1_knockdown',label,np.isfinite(truthkd)&(np.digitize(truthkd,[.5,.8])==j))
stratum('observed_H1_knockdown','unknown',~np.isfinite(truthkd))
for j,label in enumerate(['weak <10','mid 10–999','strong >=1000']):stratum('n_conf_MOMENT_PROXY',label,np.digitize(nconf,[10,1000])==j)
groups=np.array([len({v['family'] for v in row}) for row in meta['coverage']])
for j in np.unique(groups):stratum('number_families',str(j),groups==j)
pd.DataFrame(strata).to_csv(rt.ROOT/'results/strata.csv',index=False)
perbase=pd.DataFrame(dict(target=targets,**{k:v for k,v in base.items() if k!='dot'},essential_proxy=ess,source_knockdown=kd,observed_H1_knockdown=truthkd,n_conf_exact=np.full(n,np.nan),n_conf_MOMENT_PROXY=nconf,n_families=groups,n_cells=z['n_cells']))
perbase['source_families']=[';'.join(sorted({v['family'] for v in a})) for a in meta['coverage']];perbase.to_csv(rt.ROOT/'results/per_target.csv',index=False)
oracle_ceiling=dict(oracle_scalar_gain=float(1-summary['oracle_scalar']['sse']/E0),oracle_existing_family_gain=float(1-oracle_family.sum()/E0),uncovered_error_mass=float(base['sse'][~direct].sum()/E0))
# Fixed paired target bootstrap for exploratory uncertainty; no held-out data.
rng=np.random.default_rng(20261026);boot=rng.integers(n,size=(10000,n));bootstrap={}
for name,a in preds.items():
    if name in ['B0','oracle_scalar']:continue
    e=moments(a)['sse'];gain=(base['sse']-e)[boot].sum(1)/base['sse'][boot].sum(1)
    dp=(pds(a)-pds(P))[boot].mean(1)
    bootstrap[name]=dict(SSE_reduction_CI95=np.quantile(gain,[.025,.975]).tolist(),PDS_delta_CI95=np.quantile(dp,[.025,.975]).tolist())
rt.save('results/analysis.json',dict(summary=summary,strata=strata,oracle_ceilings=oracle_ceiling,coefficients=coeff,bootstrap=bootstrap,limitations=['Exact full-training Wilcoxon n_conf unavailable: raw cell ranks removed; moment-based BH count is explicitly separate.','Essentiality is K562 essential-screen membership, not an H1 fitness measurement.','Baseline uses exact integer emitter bulk when renderer audit exists; proposed fix vectors use continuous bulk projection. No jackknife MSE.','All fits exclude target i; reused single destination H1 is not independent context validation.'],training_rows=168,baseline_integer_emitted=rendered_path.exists(),heldout_outcomes_loaded=False,seconds=time.monotonic()-start))
# Retain only one diagnostic challenger; methods remain reproducible from script.
np.savez_compressed(rt.ROOT/'results/loo_bias.npz',response=preds['LOO_bias_correction'].astype(np.float32),targets=targets,genes=genes)
# Remove large temporary source-family matrices from our own saved baseline.
keep={k:z[k] for k in z.files if not k.startswith('family_') and k not in ['pool','weighted']}
keep['realized']=P.astype(np.float32);z.close();np.savez_compressed(rt.ROOT/'results/baseline_compact.npz',**keep)
(rt.ROOT/'results/baseline.npz').unlink();(rt.ROOT/'results/baseline_compact.npz').rename(rt.ROOT/'results/baseline.npz')
if rendered_path.exists():rendered_path.unlink()
rt.status(f'168-target decomposition and fixed LOO calibration complete in {time.monotonic()-start:.1f}s; exact n_conf unavailable, moment proxy explicitly labeled.')
