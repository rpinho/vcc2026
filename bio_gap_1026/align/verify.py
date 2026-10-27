"""Independent fit/metric replay and exact rational alignment certificate.

Does not import analyze.py. RRR is replayed in primal source coordinates using
Cholesky whitening, independent of the search's kernel eigensystem formula.
"""
import runtime as rt
rt.enforce(); rt.guard(True)
import numpy as np, pandas as pd, gc, time
from scipy.linalg import eigh, cholesky, solve_triangular
from fractions import Fraction
from threadpoolctl import threadpool_info

start=time.monotonic()
rt.save('results/verify_run.json',dict(pid=rt.os.getpid(),free_gib=rt.guard(True)))
audit=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
choices=rt.json.loads((rt.ROOT/'results/choices.json').read_text())
for path,wanted in audit['pins'].items():
    assert rt.sha(path)==wanted,path
with np.load(rt.PARENT/'results/baseline.npz') as z:
    Y=z['truth'].astype(float); P=z['realized'].astype(float)
    targets=z['targets']; genes=z['genes']; own=z['own_idx']
n,g=Y.shape; assert (n,g)==(168,18077)
ii=np.arange(n); valid=own>=0
Y[ii[valid],own[valid]]=0; P[ii[valid],own[valid]]=0
YY=Y @ Y.T; PY=P @ Y.T; PP=P @ P.T
yy=np.sum(Y*Y,axis=1); py=np.sum(Y*P,axis=1); pp=np.sum(P*P,axis=1)
Yown=np.zeros((n,n)); Pown=np.zeros((n,n))
Yown[:,valid]=Y[:,own[valid]]; Pown[:,valid]=P[:,own[valid]]
keep=~np.isin(genes,targets)

# Independently rebuild the seven-family Gram; no cached truth-derived input.
manifest=rt.json.loads((rt.REPO/'model_billmB/data/manifest.json').read_text())
sums={}; counts={}
for bank,meta in sorted(manifest.items()):
    if meta['family']=='H1': continue
    rt.guard(True)
    with np.load(rt.REPO/f'model_billmB/data/{bank}.npz') as z:
        source={str(t):i for i,t in enumerate(z['targets'])}
        axis={str(t):i for i,t in enumerate(z['genes'])}
        gi=np.array([axis[str(t)] for t in genes])
        measured=z['measured'][gi]
        response=z['response']; f=meta['family']
        sums.setdefault(f,np.zeros_like(Y)); counts.setdefault(f,np.zeros_like(Y))
        for j,t in enumerate(targets):
            if str(t) in source:
                sums[f][j]+=response[source[str(t)],gi]*measured
                counts[f][j]+=measured
    del response; gc.collect()
raw=np.zeros((n,n))
for f in sums:
    a=sums[f]/np.maximum(counts[f],1)
    a[:,~keep]=0; raw+=a @ a.T
norm=np.sqrt(np.diag(raw))
K=raw/np.maximum(norm[:,None]*norm[None,:],1e-30)
with np.load(rt.ROOT/'results/source_kernel.npz') as z:
    assert np.max(np.abs(K-z['K']))<1e-12
    assert np.allclose(raw,z['raw_gram'],atol=1e-8,rtol=1e-12)
del sums,counts; gc.collect()
rt.status('Independent source kernel and all input hashes PASS. Replaying outer fits in primal coordinates.')
rt.guard(True)

# Any fixed square root of the source kernel represents identical features.
d,U=eigh(K); good=d>1e-12
X=U[:,good]*np.sqrt(d[good])
# Preserve structural zero source vectors exactly; arbitrary eigensolver
# roundoff there would turn exact retrieval ties into arbitrary rankings.
X[np.diag(K)==0]=0
assert np.max(np.abs(X @ X.T-K))<1e-10
KS=(3,5,10,20); RS=(5,10,20,40,80,0); LS=(.01,.1,1.,10.,100.)
HYPERS=[(l,r) for l in LS for r in RS]
WS=[0.,.05,.1,.2,.35,.5,.75,1.]

def primal(train,query,lam,rank):
    train=np.asarray(train); xt=X[train]
    penalty=lam*max(np.sum(xt*xt)/len(train),1e-12)
    if rank==0:
        # Direct linear-system ridge rather than spectral shrinkage.
        coef=np.linalg.solve(K[np.ix_(train,train)]+penalty*np.eye(len(train)),K[np.ix_(train,query)]).T
    else:
        L=cholesky(xt.T @ xt+penalty*np.eye(X.shape[1]),lower=True)
        H=solve_triangular(L,xt.T,lower=True)
        cov=H @ YY[np.ix_(train,train)]@H.T
        _,V=eigh(cov,subset_by_index=[max(0,len(cov)-rank),len(cov)-1])
        B=solve_triangular(L.T,V@(V.T @ H),lower=False)
        coef=X[query]@B
    out=np.zeros((len(query),n)); out[:,train]=coef
    return out

def neighbors(train,query,k):
    b=np.zeros((len(query),n))
    for row,j in enumerate(query):
        pairs=sorted(((float(K[j,i]),int(i)) for i in train),key=lambda p:(-p[0],p[1]))[:k]
        weights=np.array([max(p[0],0) for p in pairs])
        if weights.sum()>1e-15:
            b[row,[p[1] for p in pairs]]=weights/weights.sum()
        else: b[row,train]=1/len(train)
    return b

def metric(a,b,rows):
    """Independently form the complete masked prediction Gram matrix."""
    rows=np.asarray(rows)
    cross=a*py[rows]+np.sum(b*YY[rows],axis=1)
    mix=(a[:,None]*PY[rows])@b.T
    gram=a[:,None]*a[None,:]*PP[np.ix_(rows,rows)]+mix+mix.T+b @ YY @ b.T
    coord=a[:,None]*Pown[np.ix_(rows,rows)]+b @ Yown[:,rows]
    rem=np.diag(coord).copy()
    gram-=coord*rem[None,:]+coord.T*rem[:,None]
    gram[np.diag_indices(len(rows))]+=rem**2
    B=float(np.trace(gram)-gram.sum()/len(rows))
    return float(cross.sum()),B

def b0spread(rows):
    return float(pp[rows].sum()-PP[np.ix_(rows,rows)].sum()/len(rows))

coef=np.load(rt.ROOT/'results/coefficients.npz')
fit_errors={}; preds={'B0':P}
for arm in choices:
    a=coef[arm+'__a']; b=coef[arm+'__b']; expected_a=np.zeros(n); expected_b=np.zeros((n,n))
    for t,ch in enumerate(choices[arm]):
        tr=ii[ii!=t]; hp=ch['hyper']
        if 'knn' in arm:
            base=neighbors(tr,[t],hp['k'])[0]
        else:
            base=primal(tr,[t],hp['lambda'],hp['rank'])[0]
        expected_a[t]=ch['final_scale']*(1-ch['weight'])
        expected_b[t]=base*ch['final_scale']*ch['weight']*ch['base_scale']
    err=float(max(np.max(np.abs(expected_a-a)),np.max(np.abs(expected_b-b))))
    assert err<2e-7,(arm,err)
    assert np.max(np.abs(np.diag(b)))==0
    fit_errors[arm]=err
    q=a[:,None]*P+b @ Y; q[ii[valid],own[valid]]=0
    preds[arm]=q
    rt.status(f'Outer {arm} replay PASS: max coefficient discrepancy {err:.3g}.')
gene=np.empty_like(P)
for t in range(n):
    v=np.var(Y[ii!=t],axis=0,ddof=1)
    row=P[t]*v
    gene[t]=row*np.linalg.norm(P[t])/max(np.linalg.norm(row),1e-15)
assert np.allclose(np.sum(gene*gene,axis=1),pp,atol=1e-8,rtol=1e-12)
preds['gene_variance']=gene

# Check all stored candidate maxima; independently rebuild all candidates for
# three complete nested outer folds (all 167 inner exclusions each).
objectives=pd.read_csv(rt.ROOT/'results/inner_objectives.csv',low_memory=False)
for arm in choices:
    for t,ch in enumerate(choices[arm]):
        sub=objectives[(objectives.arm==arm)&(objectives.outer==t)]
        if '_PDS' in arm:
            sub=sub[sub.feasible.eq(True)]
        assert abs(float(sub.objective.max())-ch['objective'])<1e-12
        match=sub[(sub.h==ch['h'])&(abs(sub.weight-ch['weight'])<1e-12)]
        assert len(match)==1 and abs(float(match.iloc[0].objective)-ch['objective'])<1e-12
nested_error=0.
for t in (0,55,167):
    rt.guard(True)
    rows=ii[ii!=t]; B0=b0spread(rows)
    ik=np.zeros((len(KS),len(rows),n)); im=np.zeros((len(HYPERS),len(rows),n))
    for pos,j in enumerate(rows):
        train=ii[(ii!=t)&(ii!=j)]
        for h,k in enumerate(KS): ik[h,pos]=neighbors(train,[j],k)[0]
        for h,(lam,r) in enumerate(HYPERS): im[h,pos]=primal(train,[j],lam,r)[0]
    for arm,inner,blend in (('knn',ik,False),('map',im,False),('B0_knn',ik,True),('B0_map',im,True)):
        sub=objectives[(objectives.arm==arm)&(objectives.outer==t)]
        for h,b in enumerate(inner):
            _,B=metric(np.zeros(len(rows)),b,rows); bs=np.sqrt(B0/max(B,1e-20))
            for w in (WS if blend else [1.]):
                c,spread=metric(np.full(len(rows),1-w),b*w*bs,rows)
                fs=np.sqrt(B0/max(spread,1e-20)); val=c*fs/yy[rows].sum()
                wanted=sub[(sub.h==h)&(abs(sub.weight-w)<1e-12)].iloc[0]
                nested_error=max(nested_error,abs(val-wanted.objective))
                assert abs(val-wanted.objective)<2e-7,(t,arm,h,w,val,wanted.objective)
    # Independently replay both added PDS-constrained blend grids. Their
    # retrieval candidate set is rows=S, excluding outer t completely.
    yr=Y[:,keep]; pr=P[:,keep]
    yyretr=yr @ yr.T; pyretr=pr @ yr.T
    ynorm=np.linalg.norm(yr,axis=1)
    def rank_scores(m):
        diag=np.diag(m)[:,None]
        return 1-((m>diag).sum(1)+.5*((m==diag).sum(1)-1))/(len(rows)-1)
    base_scores=pyretr[np.ix_(rows,rows)]/ynorm[rows][None,:]
    base_pds=rank_scores(base_scores)
    for arm,inner in (('B0_knn_PDS',ik),('B0_map_PDS',im),('B0_knn_PDS_CI',ik),('B0_map_PDS_CI',im)):
        if arm not in choices: continue
        sub=objectives[(objectives.arm==arm)&(objectives.outer==t)]
        for h,b in enumerate(inner):
            _,B=metric(np.zeros(len(rows)),b,rows); bs=np.sqrt(B0/max(B,1e-20))
            source_scores=(b @ yyretr)[:,rows]/ynorm[rows][None,:]
            for w in audit['constrained_blends']['weights']:
                c,spread=metric(np.full(len(rows),1-w),b*w*bs,rows)
                fs=np.sqrt(B0/max(spread,1e-20)); val=c*fs/yy[rows].sum()
                differences=rank_scores((1-w)*base_scores+w*bs*source_scores)-base_pds
                dp=float(differences.mean())
                lower=float(dp-1.96*differences.std(ddof=1)/np.sqrt(len(rows)))
                wanted=sub[(sub.h==h)&(abs(sub.weight-w)<1e-12)].iloc[0]
                nested_error=max(nested_error,abs(val-wanted.objective))
                assert abs(val-wanted.objective)<2e-7,(t,arm,h,w,val,wanted.objective)
                assert abs(dp-wanted.inner_delta_PDS)<1e-12,(t,arm,h,w,dp,wanted.inner_delta_PDS)
                if arm.endswith('_CI'):
                    assert abs(lower-wanted.inner_PDS_lower_N95)<1e-12
                assert bool((lower if arm.endswith('_CI') else dp)>=-.005-1e-15)==bool(wanted.feasible)
    # Label poisoning: held-out row/column can change arbitrarily without
    # changing a selected fit, since only the training principal block enters.
    j=int(rows[0]); train=ii[(ii!=t)&(ii!=j)]
    before=primal(train,[j],1.,20)
    oldrow=YY[t].copy(); oldcol=YY[:,t].copy()
    YY[t,:]=1e20; YY[:,t]=1e20
    after=primal(train,[j],1.,20)
    YY[t]=oldrow; YY[:,t]=oldcol
    assert np.array_equal(before,after)
    rt.status(f'Independent full inner grid replay PASS at outer t={t}; outer-label poisoning invariant.')

from cell_eval2.metrics.discrimination import discrimination_score
rt.guard(True)
from cell_eval2 import EvalConfig
cfg=EvalConfig.from_preset('vcc2026').discrimination
labels=np.concatenate([targets,['non-targeting']])
real=np.vstack([Y,np.zeros(g)])
per=pd.read_csv(rt.ROOT/'results/per_target.csv')
native={}; metrics={}; maxerr=0.
for arm,q in preds.items():
    rt.guard()
    out=discrimination_score(pred_bulk=(labels,np.vstack([q,np.zeros(g)])),
        real_bulk=(labels,real),genes=genes,distance='cosine',
        rank_denominator=cfg.rank_denominator,tie_policy=cfg.tie_policy,
        exclude_target_gene=cfg.exclude_target_gene,exclusion_scope=cfg.exclusion_scope,
        control_source='real')
    pds=np.array([out[t] for t in targets])
    G=q @ q.T; dot=np.sum(q.astype(np.longdouble)*Y.astype(np.longdouble),axis=1).astype(float)
    energy=np.sum(q.astype(np.longdouble)**2,axis=1).astype(float)
    B=float(np.sum((q-q.mean(0))**2))
    metrics[arm]=dict(gram=G,dot=dot,energy=energy,B=B,pds=pds)
    scale=np.sqrt(metrics['B0']['B']/B); want=audit['summary'][arm]
    values=dict(cross=float(scale*dot.sum()/yy.sum()),
        aligned_per_spread=float(dot.sum()/np.sqrt(energy.sum()*yy.sum())),
        PDS=float(pds.mean()),B=metrics['B0']['B']/yy.sum())
    for k,v in values.items():
        err=abs(v-want[k]); maxerr=max(maxerr,err); assert err<1e-10,(arm,k,v,want[k])
    rows=per[per.arm==arm]
    assert rows.target.tolist()==targets.tolist()
    assert np.max(np.abs(rows.pds-pds))<1e-12
    assert np.max(np.abs(rows.cross-scale*dot))<1e-7
    native[arm]=float(pds.mean())

draws=np.random.default_rng(20261026).integers(n,size=(2000,n))
counts=np.array([np.bincount(row,minlength=n) for row in draws],dtype=float)
Tboot=counts @ yy; base=metrics['B0']
Bboot0=counts @ base['energy']-np.einsum('bi,ij,bj->b',counts,base['gram'],counts,optimize=True)/n
cboot0=(counts @ base['dot'])/Tboot
pboot0=(counts @ base['pds'])/n
fixed_sensitivity={}; bootstrap_error=0.
for arm,st in metrics.items():
    Bboot=counts @ st['energy']-np.einsum('bi,ij,bj->b',counts,st['gram'],counts,optimize=True)/n
    cboot=np.sqrt(Bboot0/np.maximum(Bboot,1e-30))*(counts @ st['dot'])/Tboot
    pboot=(counts @ st['pds'])/n
    checks=dict(cross_CI95=np.quantile(cboot,[.025,.975]),
        delta_cross_CI95=np.quantile(cboot-cboot0,[.025,.975]),
        PDS_CI95=np.quantile(pboot,[.025,.975]),
        delta_PDS_CI95=np.quantile(pboot-pboot0,[.025,.975]))
    for k,v in checks.items():
        err=float(np.max(np.abs(v-audit['summary'][arm][k])))
        bootstrap_error=max(bootstrap_error,err); assert err<1e-10,(arm,k,err)
    alpha=np.sqrt(base['B']/st['B'])
    fixed=(alpha*st['dot']-base['dot'])[draws].sum(1)/yy[draws].sum(1)
    fixed_sensitivity[arm]=np.quantile(fixed,[.025,.975]).tolist()
rt.status('Independent gene-space metrics, native PDS and all 2000 paired rematched-spread bootstrap intervals PASS. Exact rational certificates next.')

def exact_moments(q,y):
    """Sum products/centered energy exactly on a common integer lattice."""
    minexp=2047
    for arr in (q,y):
        bits=np.ascontiguousarray(arr).view(np.uint64)
        nonzero=(bits&((1<<63)-1))!=0; exp=(bits>>52)&2047
        assert np.all((exp[nonzero]>0)&(exp[nonzero]<2047))
        if nonzero.any(): minexp=min(minexp,int(exp[nonzero].min()))
    power=max(0,1075-minexp)
    qq=0; qy=0; yenergy=0; columns=[0]*q.shape[1]
    for qr,yr in zip(q,y):
        for col,(qv,yv) in enumerate(zip(qr,yr)):
            a,b=float(qv).as_integer_ratio(); c,d=float(yv).as_integer_ratio()
            qi=a<<(power-(b.bit_length()-1)); yi=c<<(power-(d.bit_length()-1))
            qq+=qi*qi; qy+=qi*yi; yenergy+=yi*yi; columns[col]+=qi
    denom=1<<(2*power)
    return dict(cross=Fraction(qy,denom),energy=Fraction(qq,denom),
        truth=Fraction(yenergy,denom),B=Fraction(n*qq-sum(v*v for v in columns),n*denom))

cert={}; exact={}
for arm,q in preds.items():
    rt.guard(True)
    moments=exact_moments(q,Y); exact[arm]=moments
    want=audit['summary'][arm]; m0=exact['B0']
    match=float(moments['cross']/moments['truth'])*np.sqrt(float(m0['B']/moments['B']))
    assert abs(match-want['cross'])<1e-10
    # Exact sign of the ideal matched-B cross difference; both cross terms
    # must be positive before squaring. This avoids a floating sqrt comparison.
    assert moments['cross']>0 and m0['cross']>0
    difference=moments['cross']**2*m0['B']-m0['cross']**2*moments['B']
    exact_sign=1 if difference>0 else (-1 if difference<0 else 0)
    if abs(want['delta_cross'])>1e-10:
        assert exact_sign==int(np.sign(want['delta_cross']))
    cert[arm]=dict(**{k:dict(numerator=str(v.numerator),denominator=str(v.denominator)) for k,v in moments.items()},
                  exact_matched_cross_delta_sign=exact_sign,matched_cross_float=match)
    rt.status(f'Exact rational full-coordinate cross/energy/B PASS: {arm}, matched delta sign={exact_sign}.')
rt.save('results/exact_alignment.json',dict(status='PASS',arms=cert,
    scope='Exact moments and matched-cross point-difference sign of verified float64 outputs; bootstrap intervals numerically replayed, not exact rational intervals; no official MSE claim.'))
assert all(p['num_threads']==1 for p in threadpool_info())
rt.save('results/verification.json',dict(status='PASS',all_input_hashes_unchanged=True,
    source_kernel_independently_rebuilt=True,outer_folds_replayed=168,
    outer_coefficient_max_errors=fit_errors,full_nested_outer_folds_replayed=[0,55,167],
    selected_inner_fits_replayed=3*167*30,nested_objective_max_error=nested_error,
    all_stored_candidate_optima_checked=True,outer_label_poisoning_tests=[0,55,167],
    all_primary_full_coordinate_moments_exact_rational=True,
    native_PDS=native,max_summary_error=maxerr,bootstrap_max_error=bootstrap_error,
    bootstrap_replicates_replayed=2000,fixed_normalization_delta_cross_CI95=fixed_sensitivity,
    seconds=time.monotonic()-start,free_gib=rt.guard(),blas=threadpool_info(),
    effective_nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0)))
rt.status(f'Independent verification PASS in {time.monotonic()-start:.1f}s. Official MSE is not measurable locally.')
