"""Independent full-output replay, exact rational SSE, and native PDS check.
Does not import analyze.py or its fitting functions.
"""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,pandas as pd,time,gc
from fractions import Fraction
from threadpoolctl import threadpool_info
start=time.monotonic()
rt.save('results/verify_run.json',dict(pid=rt.os.getpid(),requested_nice=10,effective_nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0)))
audit=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
coef=rt.json.loads((rt.ROOT/'results/coefficients.json').read_text())
for path,h in audit['pins'].items():assert rt.sha(path)==h,path
with np.load(rt.PARENT/'results/baseline.npz') as z:
    y=z['truth'].astype(float);p=z['realized'].astype(float);targets=z['targets'];genes=z['genes'];own=z['own_idx'];bulk=z['bulk']
n,g=y.shape;assert n==168 and g==18077
ix=np.arange(n);valid=own>=0
for i,j in enumerate(own):
    if j>=0:y[i,j]=0;p[i,j]=0

def clean(a):
    a=np.array(a,copy=True)
    for i,j in enumerate(own):
        if j>=0:a[i,j]=0
    return a

def feasible(a):
    b=np.log1p(bulk*50000);q=np.expm1(a+b);q[q<0]=0;q/=q.sum(1)[:,None]
    return clean(np.log1p(q*50000)-b)

# Recompute all nested betas from explicit subcohorts, not global subtraction.
yy=y @ y.T;py=p @ y.T;pp=np.einsum('ij,ij->i',p,p);yp=np.einsum('ij,ij->i',p,y)
yown=np.zeros((n,n));yown[:,valid]=y[:,own[valid]]
pair_cache={};max_beta_err=0

def fit_explicit(ex):
    ids=np.array([i for i in range(n) if i not in ex]);L=len(ids)
    qyy=yy[np.ix_(ids,ids)];qpy=py[np.ix_(ids,ids)]
    mm=(qyy.sum()-2*qyy.sum(1)+np.diag(qyy)-yown[np.ix_(ids,ids)].sum(0)**2)/(L-1)**2
    pm=(qpy.sum(1)-yp[ids])/(L-1)
    ym=(qyy.sum(1)-np.diag(qyy))/(L-1)
    return np.array([np.sum(yp[ids]-pm)/np.sum(pp[ids]),np.sum(yp[ids]-pm*ym/mm)/np.sum(pp[ids]-pm*pm/mm)])

for t in range(n):
    fitted=fit_explicit({t});wanted=np.array([coef['beta_b'][t],coef['beta_c'][t]])
    max_beta_err=max(max_beta_err,float(np.max(abs(fitted-wanted))))
    for j in range(t):
        fitted_pair=fit_explicit({t,j})[1];pair_cache[j,t]=fitted_pair
        max_beta_err=max(max_beta_err,abs(fitted_pair-coef['inner_beta_c'][t][j]),abs(fitted_pair-coef['inner_beta_c'][j][t]))
assert max_beta_err<1e-11,max_beta_err
# Selected folds also use actual gene-space design vectors; excludes t and j
# before computing a mean. This catches a wrong Gram/masking identity.
direct_spots=[]
for ex in [{0},{55},{167},{0,55},{55,167},{0,167}]:
    ids=np.array([i for i in range(n) if i not in ex]);s=y[ids].sum(0);L=len(ids)
    m=(s-y[ids])/(L-1)
    for row,j in enumerate(ids):
        if own[j]>=0:m[row,own[j]]=0
    pred=p[ids];truth=y[ids];norm=np.sum(m*m,1);q=pred-(np.sum(pred*m,1)/norm)[:,None]*m
    b=float(np.sum(pred*(truth-m))/np.sum(pred*pred))
    c=float(np.sum(q*(truth-m))/np.sum(q*q))
    assert np.max(abs(np.array([b,c])-fit_explicit(ex)))<1e-11
    direct_spots.append(sorted(ex))
rt.status(f'Independent nested-fit replay PASS: all outer and 14,028 inner-pair betas; max discrepancy {max_beta_err:.2g}.')

# Independent reconstruction of seven family features from original banks.
manifest=rt.json.loads((rt.REPO/'model_billmB/data/manifest.json').read_text())
families=coef['family_names'];R=np.zeros((len(families),n,g));covered=np.zeros((len(families),n),bool)
for fno,family in enumerate(families):
    accum=np.zeros((n,g),np.float32);den=np.zeros((n,g),np.float32)
    for name,meta in sorted(manifest.items()):
        if meta['family']!=family:continue
        rt.guard(True)
        with np.load(rt.REPO/f'model_billmB/data/{name}.npz') as z:
            rows=pd.Index(z['targets']).get_indexer(targets);columns=pd.Index(z['genes']).get_indexer(genes);measured=z['measured'][columns];response=z['response']
            for i,j in enumerate(rows):
                if j<0:continue
                covered[fno,i]=True;accum[i]+=response[j,columns]*measured;den[i]+=measured
        del response
    realized=feasible(accum/np.maximum(den,1e-8)).astype(np.float32).astype(float)
    realized[~covered[fno]]=p[~covered[fno]]
    R[fno]=realized-p
    del accum,den,realized;gc.collect()

m=np.stack([np.mean(y[ix!=i],axis=0) for i in range(n)]);m=clean(m)
pm=np.sum(p*m,1);mm=np.sum(m*m,1);ym=np.sum(y*m,1);q=p-(pm/mm)[:,None]*m
preds={'B0':p,'a_mean':m,'b_mean_plus_B0':clean(m+np.array(coef['beta_b'])[:,None]*p),'c_mean_plus_orthogonal':clean(m+np.array(coef['beta_c'])[:,None]*q)}
# Fit family residual coefficients afresh using actual gene-space residuals,
# with inner c predictions constructed separately for every outer cohort.
d=np.empty_like(p);e=np.empty_like(p);max_w_error=0
for t in range(n):
    tr=ix!=t;ids=ix[tr];rr=R[:,tr];K=len(families)
    gram=np.einsum('kig,lig->kl',rr,rr,optimize=True);pen=np.trace(gram)/K
    rhs=np.einsum('kig,ig->k',rr,y[tr]-p[tr],optimize=True)
    wd=np.linalg.solve(gram+max(pen,1e-12)*np.eye(K),rhs)
    # Sum only the explicitly selected outer cohort, then subtract j.
    # This computes the same 167 inner means without repeatedly gathering
    # all 166 rows for each mean.
    outer_sum=y[tr].sum(0)
    inner=(outer_sum[None,:]-y[tr])/(len(ids)-1)
    for k,j in enumerate(ids):
        if own[j]>=0:inner[k,own[j]]=0
    inner_q=p[tr]-(np.sum(p[tr]*inner,1)/np.sum(inner*inner,1))[:,None]*inner
    ib=np.array([pair_cache[tuple(sorted((t,int(j))))] for j in ids])
    residual=y[tr]-inner-ib[:,None]*inner_q
    rhs=np.einsum('kig,ig->k',rr,residual,optimize=True)
    we=np.linalg.solve(gram+max(pen,1e-12)*np.eye(K),rhs)
    max_w_error=max(max_w_error,float(np.max(abs(wd-coef['family_d'][t]))),float(np.max(abs(we-coef['family_e'][t]))))
    d[t]=p[t]+np.sum(wd[:,None]*R[:,t],0)
    e[t]=preds['c_mean_plus_orthogonal'][t]+np.sum(we[:,None]*R[:,t],0)
assert max_w_error<1e-10,max_w_error
preds['d_family_correction']=clean(d);preds['e_c_plus_family']=clean(e)
rt.status(f'Family fits independently replayed in gene space for all 168 outer folds; max coefficient discrepancy {max_w_error:.2g}. Exact rational SSE certificates next.')

# Exact rational arithmetic over every coordinate of every primary arm.
# Float64 vectors are finite dyadic rationals. An integer grid common to a
# pair of arrays makes both subtraction and the sum of squares exact.
def exact_squared(a,b=None):
    arrays=[np.ascontiguousarray(a)] if b is None else [np.ascontiguousarray(a),np.ascontiguousarray(b)]
    min_exp=2047
    for ar in arrays:
        bits=ar.view(np.uint64);ex=(bits>>52)&2047;nonzero=bits&((1<<63)-1)!=0
        assert np.all((ex[nonzero]>0)&(ex[nonzero]<2047))
        if np.any(nonzero):min_exp=min(min_exp,int(ex[nonzero].min()))
    k=max(0,1075-min_exp);acc=0
    if b is None:
        for av in a.flat:
            num,den=float(av).as_integer_ratio();v=num<<(k-(den.bit_length()-1));acc+=v*v
    else:
        for av,bv in zip(a.flat,b.flat):
            an,ad=float(av).as_integer_ratio();bn,bd=float(bv).as_integer_ratio()
            v=(an<<(k-(ad.bit_length()-1)))-(bn<<(k-(bd.bit_length()-1)));acc+=v*v
    return Fraction(acc,1<<(2*k))

truth_exact=exact_squared(y);cert={};SSE0=None
for name,pred in preds.items():
    value=exact_squared(pred,y);ratio=value/truth_exact;cert[name]=dict(sse_numerator=str(value.numerator),sse_denominator=str(value.denominator),E_numerator=str(ratio.numerator),E_denominator=str(ratio.denominator),E=float(ratio))
    assert abs(float(ratio)-audit['summary']['vector'][name]['E'])<1e-11
    if name=='B0':SSE0=value
    else:assert value<SSE0
    print('exact',name,float(ratio),flush=True)
rt.save('results/exact_SSE.json',dict(status='PASS',meaning='Exact all-coordinate dyadic-rational SSE of independently replayed float64 primary predictions; not a certificate of biological/official score gain.',truth_norm2_numerator=str(truth_exact.numerator),truth_norm2_denominator=str(truth_exact.denominator),arms=cert))

from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2 import EvalConfig
cfg=EvalConfig.from_preset('vcc2026').discrimination
labels=np.concatenate([targets,['non-targeting']]);real=np.vstack([y,np.zeros(g)])
table=pd.read_csv(rt.ROOT/'results/per_target.csv');native={};max_metric_error=0
B0=float(np.sum((p-p.mean(0))**2));rng=np.random.default_rng(20261026);draws=rng.integers(n,size=(2000,n));base_sse=np.sum((p-y)**2,1)
for space in ['vector','projected']:
    for name,primary in preds.items():
        pred=feasible(primary) if space=='projected' and name!='B0' else primary
        values=discrimination_score(pred_bulk=(labels,np.vstack([pred,np.zeros(g)])),real_bulk=(labels,real),genes=genes,distance='cosine',rank_denominator=cfg.rank_denominator,tie_policy=cfg.tie_policy,exclude_target_gene=cfg.exclude_target_gene,exclusion_scope=cfg.exclusion_scope,control_source='real')
        pdsv=np.array([values[t] for t in targets]);want=audit['summary'][space][name]
        subset=table[(table.space==space)&(table.arm==name)]
        assert subset.target.tolist()==targets.tolist()
        sse=np.sum((pred.astype(np.longdouble)-y.astype(np.longdouble))**2,1).astype(float)
        assert np.max(abs(sse-subset.sse))<1e-9
        assert np.max(abs(pdsv-subset.pds))<1e-12
        mean_cos=np.mean(np.sum(pred*y,1)/np.sqrt(np.sum(pred*pred,1)*np.sum(y*y,1)))
        B=float(np.sum((pred-pred.mean(0))**2))
        for val,key in [(sse.sum()/float(truth_exact),'E'),(pdsv.mean(),'PDS'),(mean_cos,'mean_cosine'),(B/B0,'B_relative_B0')]:
            max_metric_error=max(max_metric_error,abs(float(val)-want[key]));assert abs(val-want[key])<1e-10
        native[space+'__'+name]=float(pdsv.mean())
        base_pdsv=table[(table.space=='vector')&(table.arm=='B0')].pds.to_numpy()
        dE=(sse-base_sse)[draws].sum(1)/np.sum(y*y,1)[draws].sum(1)
        dP=(pdsv-base_pdsv)[draws].mean(1)
        assert np.max(abs(np.quantile(dE,[.025,.975])-want['delta_E_CI95']))<1e-10
        assert np.max(abs(np.quantile(dP,[.025,.975])-want['delta_PDS_CI95']))<1e-12

# Orthogonal decomposition, cross-fitted shared fraction, calibrated residual.
yperp=y-(ym/mm)[:,None]*m
shared_exact_spots=[]
for t in [0,55,167]:
    mq=[Fraction(float(v)) for v in m[t]];pq=[Fraction(float(v)) for v in p[t]];yq=[Fraction(float(v)) for v in y[t]]
    mn=sum(v*v for v in mq);mp=sum(a*b for a,b in zip(mq,pq));my=sum(a*b for a,b in zip(mq,yq));pn=sum(v*v for v in pq);yn=sum(v*v for v in yq)
    assert pn-mp*mp/mn>=0 and yn-my*my/mn>=0
    shared_exact_spots.append(dict(target=str(targets[t]),orthogonal_energy_nonnegative_exact=True))
sh=audit['shared'];assert abs(np.sum(ym*ym/mm)/np.sum(y*y)-sh['projection_energy_fraction'])<1e-12
assert abs(np.sum((q-yperp)**2)/np.sum(yperp*yperp)-sh['orthogonal_E'])<1e-12
calibrated=np.array(coef['beta_c'])[:,None]*q
extra=dict(nested_scaled_orthogonal_E=float(np.sum((calibrated-yperp)**2)/np.sum(yperp*yperp)),projection_fraction_per_target_quantiles=np.quantile(ym*ym/mm/np.sum(y*y,1),[0,.25,.5,.75,1]).tolist(),truth_shared_cosine_quantiles=np.quantile(ym/np.sqrt(mm*np.sum(y*y,1)),[0,.25,.5,.75,1]).tolist(),shared_mean_improvement_vs_zero=float(1-audit['summary']['vector']['a_mean']['E']))
assert all(v['num_threads']==1 for v in threadpool_info())
rt.save('results/verification.json',dict(status='PASS',all_input_hashes_unchanged=True,outer_folds_replayed=168,unique_inner_pairs_replayed=14028,direct_gene_space_beta_spots=direct_spots,max_beta_abs_error=max_beta_err,max_family_weight_abs_error=max_w_error,max_summary_abs_error=max_metric_error,native_PDS=native,bootstrap_replicates_replayed=2000,all_primary_SSE_exact_rational=True,exact_orthogonal_energy_spots=shared_exact_spots,extra=extra,nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0),blas=threadpool_info(),free_gib=rt.guard(),seconds=time.monotonic()-start))
rt.status(f'Independent verification PASS in {time.monotonic()-start:.1f}s: nested fits, exact all-coordinate primary SSE, native PDS, paired bootstrap CIs, and spread. Official MSE remains unmeasurable.')
