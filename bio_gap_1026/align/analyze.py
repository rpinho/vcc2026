"""Source-only nested target-LOO alignment analysis; no held-out79 access."""
import runtime as rt
rt.enforce(); rt.guard(True)
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from threadpoolctl import threadpool_info
import time, gc

START = time.monotonic()
N = 168
KS = (3, 5, 10, 20)
RANKS = (5, 10, 20, 40, 80, 0)  # 0 = unrestricted ridge
LAMBDAS = (.01, .1, 1., 10., 100.)
WEIGHTS = np.array([0., .05, .1, .2, .35, .5, .75, 1.])
HYPERS = [(lam, rank) for lam in LAMBDAS for rank in RANKS]
PINS = {}

def pin(path):
    path = rt.Path(path); PINS[str(path)] = rt.sha(path); return path

def load():
    global targets, genes, Y, P, own, YY, PY, PP, Yown, Pown, yy, py, pp, T, keep
    with np.load(pin(rt.PARENT/'results/baseline.npz')) as z:
        targets = z['targets']; genes = z['genes']; own = z['own_idx']
        Y = z['truth'].astype(float); P = z['realized'].astype(float)
    assert Y.shape == P.shape == (N, 18077)
    valid = own >= 0
    assert len(set(own[valid])) == valid.sum()
    Y[np.arange(N)[valid], own[valid]] = 0
    P[np.arange(N)[valid], own[valid]] = 0
    with np.load(pin(rt.PARENT/'shared_resp/results/sufficient.npz')) as z:
        assert np.array_equal(targets, z['targets']) and np.array_equal(genes, z['genes'])
        YY = z['YY']; PY = z['PY']; pp = z['PP']; Yown = z['Yown']
        # Baseline comparator replay; these similarities are NEVER features.
        saved_similarity = z['vector__B0']
    PP = P @ P.T
    Pown = np.zeros((N,N)); Pown[:,valid] = P[:,own[valid]]
    assert np.allclose(Y @ Y.T, YY, atol=1e-8, rtol=1e-12)
    assert np.allclose(P @ Y.T, PY, atol=1e-8, rtol=1e-12)
    assert np.allclose(np.diag(PP), pp, atol=1e-8, rtol=1e-12)
    yy = np.diag(YY); py = np.diag(PY); T = float(yy.sum())
    keep = ~np.isin(genes, targets)
    native = P[:,keep]@Y[:,keep].T / (np.linalg.norm(P[:,keep],axis=1)[:,None]*np.linalg.norm(Y[:,keep],axis=1)[None,:])
    assert np.allclose(native, saved_similarity, atol=1e-12, rtol=1e-12)
    rt.status(f'PASS B0/sufficient-statistic replay. Truth energy={T:.6f}; source-only features next.')

def sources():
    """Read only pre-existing non-H1 source rows matching training168 identities."""
    manifest = rt.json.loads(pin(rt.REPO/'model_billmB/data/manifest.json').read_text())
    sums = {}; den = {}; rows = []
    for name, entry in sorted(manifest.items()):
        if entry['family'] == 'H1':
            continue
        rt.guard(True)
        path = pin(rt.REPO/f'model_billmB/data/{name}.npz')
        assert PINS[str(path)] == entry['sha256']
        with np.load(path) as z:
            ti = pd.Index(z['targets']).get_indexer(targets)
            gi = pd.Index(z['genes']).get_indexer(genes)
            assert (gi>=0).all()
            good = np.flatnonzero(ti>=0); measured = z['measured'][gi]
            response = z['response']
            fam = entry['family']
            sums.setdefault(fam, np.zeros_like(Y)); den.setdefault(fam, np.zeros_like(Y))
            sums[fam][good] += response[ti[good]][:,gi]*measured
            den[fam][good] += measured
            rows.append(dict(bank=name, family=fam, matching_rows=len(good)))
        del response; gc.collect()
    gram = np.zeros((N,N)); coverage = np.zeros(N, int)
    for fam in sorted(sums):
        a = sums[fam]/np.maximum(den[fam], 1.)
        # Remove direct-KD coordinates from every source row, using the same
        # fixed training-panel gene exclusion as retrieval. No fitted transform.
        a[:,~keep] = 0
        gram += a @ a.T
        coverage += np.any(den[fam]>0, axis=1)
    norms = np.sqrt(np.maximum(np.diag(gram), 0))
    K = gram/np.maximum(norms[:,None]*norms[None,:], 1e-30)
    assert np.isfinite(K).all() and np.max(np.abs(K-K.T))<1e-12
    np.savez_compressed(rt.ROOT/'results/source_kernel.npz', K=K, raw_gram=gram,
                        targets=targets, coverage=coverage)
    rt.save('results/source_rows.json', rows)
    rt.status(f'Source kernel complete: {len(rows)} banks, {len(sums)} families, {np.sum(norms==0)} targets without a source row. Saved H1 retrieval similarities not used.')
    return K

def cross_energy(a, b, rows):
    """Rows of q = a_i P_i + b_i Y, then zero own output gene.

    Return exact gene-space cross, energy and centered spread using reused
    cross-products. All own indices are distinct.
    """
    rows = np.asarray(rows)
    a = np.asarray(a)
    by = b @ YY
    removed = np.sum(b*Yown[:,rows].T, axis=1)
    cross = a*py[rows] + by[np.arange(len(rows)),rows]
    norm = a*a*pp[rows] + 2*a*np.sum(b*PY[rows],axis=1) + np.sum(by*b,axis=1)-removed**2
    av = np.zeros(N); av[rows] = a
    bv = b.sum(0)
    sum_norm = av @ PP @ av + 2*av @ PY @ bv + bv @ YY @ bv
    sum_own = av @ Pown[:,rows] + bv @ Yown[:,rows]
    sum_norm += -2*removed @ sum_own + removed @ removed
    B = float(norm.sum()-sum_norm/len(rows))
    assert B > -1e-6
    return cross, norm, max(B,1e-20)

def baseline_B(rows):
    rows = np.asarray(rows)
    return float(pp[rows].sum()-PP[np.ix_(rows,rows)].sum()/len(rows))

def knn(K, train, tests):
    sim = K[np.ix_(tests,train)]
    order = np.argsort(-sim, axis=1, kind='stable')
    out = []
    for k in KS:
        inds = order[:,:k]
        weights = np.maximum(np.take_along_axis(sim, inds, axis=1), 0)
        denom = weights.sum(1)
        b = np.zeros((len(tests),N))
        for i in range(len(tests)):
            if denom[i] > 1e-15:
                b[i,np.asarray(train)[inds[i]]] = weights[i]/denom[i]
            else:
                b[i,train] = 1/len(train)
        out.append(b)
    return np.asarray(out)

def maps(K, train, tests):
    """Exact ridge reduced-rank regression, including unrestricted ridge.

    No intercept: zero source-response features map to zero delta. Rank is
    imposed on the coefficient matrix in the ridge-penalized least-squares
    problem, rather than fitting an output PCA on the validation labels.
    """
    train = np.asarray(train)
    k = K[np.ix_(train,train)]
    d, U = eigh(k, check_finite=False, driver='evr')
    good = d > max(d.max()*1e-11, 1e-12)
    d = d[good]; U = U[:,good]
    yyrot = U.T @ YY[np.ix_(train,train)]@U
    query = K[np.ix_(tests,train)]@U
    penalty_scale = max(float(np.trace(k)/len(train)), 1e-12)
    out = []
    for lam in LAMBDAS:
        penalty = lam*penalty_scale
        D = np.sqrt(d/(d+penalty))
        e, V = eigh(D[:,None]*yyrot*D[None,:], check_finite=False, driver='evr')
        L = query/np.sqrt(d*(d+penalty))
        for rank in RANKS:
            if rank == 0:
                coef = (query/(d+penalty))@U.T
            else:
                Vr = V[:,-min(rank,len(d)):]
                coef = ((L @ Vr)@Vr.T*D)@U.T
            b = np.zeros((len(tests),N)); b[:,train] = coef
            out.append(b)
    return np.asarray(out)

def select(inner, rows, blend=False):
    """Choose all hyperparameters jointly on inner held-out cross at matched B.

    The source prediction is first scaled from its inner OOF B to B0's B.
    Blend weights therefore compare components on a common spread scale.
    A final scale for the selected blend is also learned from inner OOF B.
    """
    B0 = baseline_B(rows)
    best = None; scores = []
    for h, b in enumerate(inner):
        _, _, B = cross_energy(np.zeros(len(rows)), b, rows)
        base_scale = np.sqrt(B0/B)
        for w in (WEIGHTS if blend else [1.]):
            a = np.full(len(rows), 1-w)
            bb = b*(w*base_scale)
            c, energy, spread = cross_energy(a,bb,rows)
            scale = np.sqrt(B0/spread)
            value = float(c.sum()*scale/yy[rows].sum())
            scores.append((h,float(w),value,float(base_scale),float(scale)))
            if best is None or value > best['objective']+1e-14:
                best = dict(h=h, weight=float(w), base_scale=float(base_scale),
                            final_scale=float(scale), objective=value)
    return best, scores

def apply(choice, candidate):
    w=choice['weight']; s=choice['final_scale']
    return s*(1-w), s*w*choice['base_scale']*candidate[choice['h'],0]

def nested(K):
    ii = np.arange(N)
    names = ('knn','map','B0_knn','B0_map')
    aa = {name:np.zeros(N) for name in names}
    bb = {name:np.zeros((N,N)) for name in names}
    choices = {name:[] for name in names}
    all_scores=[]
    for t in range(N):
        rt.guard()
        rows = ii[ii!=t]
        ik = np.zeros((len(KS), len(rows), N))
        im = np.zeros((len(HYPERS), len(rows), N))
        for pos,j in enumerate(rows):
            train = ii[(ii!=t)&(ii!=j)]
            ik[:,pos] = knn(K,train,[j])[:,0]
            im[:,pos] = maps(K,train,[j])[:,0]
        outer_k = knn(K,rows,[t]); outer_m = maps(K,rows,[t])
        for name,inner,outer,blend in (
            ('knn',ik,outer_k,False), ('map',im,outer_m,False),
            ('B0_knn',ik,outer_k,True), ('B0_map',im,outer_m,True)):
            choice, scores = select(inner,rows,blend)
            aa[name][t],bb[name][t] = apply(choice,outer)
            choice['target']=str(targets[t])
            choice['hyper'] = dict(k=KS[choice['h']]) if 'knn' in name else dict(zip(('lambda','rank'),HYPERS[choice['h']]))
            choices[name].append(choice)
            for h,w,v,bs,fs in scores:
                all_scores.append(dict(outer=t,arm=name,h=h,weight=w,objective=v,base_scale=bs,final_scale=fs))
        if t%4==0 or t==N-1:
            rt.status(f'Nested LOO {t+1}/{N}: elapsed {time.monotonic()-START:.1f}s; no target t or inner validation j labels enter its fit.')
    pd.DataFrame(all_scores).to_csv(rt.ROOT/'results/inner_objectives.csv',index=False)
    rt.save('results/choices.json',choices)
    np.savez_compressed(rt.ROOT/'results/coefficients.npz',**{name+'__a':aa[name] for name in names},**{name+'__b':bb[name] for name in names})
    return aa,bb

def gene_prior():
    # Unbiased across-target variance, estimated on the other 167 only.
    sy=Y.sum(0); sy2=np.sum(Y*Y,axis=0)
    out=np.empty_like(P)
    for t in range(N):
        variance=np.maximum((sy2-Y[t]**2-(sy-Y[t])**2/(N-1))/(N-2),0)
        q=P[t]*variance
        out[t]=q*np.sqrt(pp[t]/max(float(q @ q),1e-30))
    return out

def stats(q):
    gram=q @ q.T
    dot=np.sum(q*Y,axis=1)
    energy=np.diag(gram).copy()
    B=float(energy.sum()-gram.sum()/N)
    c=q[:,keep]@Y[:,keep].T / np.maximum(np.linalg.norm(Y[:,keep],axis=1)[None,:],1e-30)
    diag=np.diag(c)[:,None]
    ranks=(c>diag).sum(1)+.5*((c==diag).sum(1)-1)
    pds=1-ranks/(N-1)
    centered_cross=float(dot.sum()-q.sum(0)@Y.sum(0)/N)
    return dict(gram=gram,dot=dot,energy=energy,B=B,pds=pds,centered_cross=centered_cross)

def evaluate(aa,bb):
    rt.guard(True)
    preds={'B0':P}
    for name in aa:
        q=aa[name][:,None]*P+bb[name]@Y
        good=own>=0; q[np.arange(N)[good],own[good]]=0
        c,e,B=cross_energy(aa[name],bb[name],np.arange(N))
        assert np.allclose(np.sum(q*Y,axis=1),c,atol=1e-7,rtol=1e-11)
        assert np.allclose(np.sum(q*q,axis=1),e,atol=1e-6,rtol=1e-10)
        assert np.isclose(np.sum((q-q.mean(0))**2),B,rtol=1e-10)
        preds[name]=q
    preds['gene_variance']=gene_prior()
    ss={name:stats(q) for name,q in preds.items()}
    B0=ss['B0']['B']
    assert abs(ss['B0']['pds'].mean()-.8558953521528372)<1e-12
    rng=np.random.default_rng(20261026)
    draws=rng.integers(N,size=(2000,N))
    counts=np.zeros((len(draws),N))
    for i,row in enumerate(draws):
        counts[i]=np.bincount(row,minlength=N)
    truthboot=counts @ yy
    base=ss['B0']
    Bboot0=counts @ base['energy']-np.sum((counts @ base['gram'])*counts,axis=1)/N
    crossboot0=(counts @ base['dot'])/truthboot
    pdsboot0=(counts @ base['pds'])/N
    summary={}; per=[]; archives={}; small={}
    for name,s in ss.items():
        alpha=float(np.sqrt(B0/s['B']))
        C=float(alpha*s['dot'].sum()/T)
        Q=float(alpha**2*s['energy'].sum()/T)
        Bboot=counts @ s['energy']-np.sum((counts @ s['gram'])*counts,axis=1)/N
        alphaboot=np.sqrt(Bboot0/np.maximum(Bboot,1e-30))
        crossboot=alphaboot*(counts @ s['dot'])/truthboot
        pdsboot=(counts @ s['pds'])/N
        dc=crossboot-crossboot0; dp=pdsboot-pdsboot0
        ci=lambda a:np.quantile(a,[.025,.975]).tolist()
        item=dict(cross=C,aligned_per_spread=C/np.sqrt(Q),prediction_energy=Q,
                  B=B0/T,B_relative_B0=1.,raw_B_relative_B0=s['B']/B0,
                  final_match_scale=alpha,PDS=float(s['pds'].mean()),
                  delta_cross=C-summary['B0']['cross'] if name!='B0' else 0.,
                  delta_PDS=float((s['pds']-base['pds']).mean()),
                  cross_CI95=ci(crossboot),PDS_CI95=ci(pdsboot),
                  delta_cross_CI95=ci(dc),delta_PDS_CI95=ci(dp),
                  centered_cross=float(alpha*s['centered_cross']/T),
                  SSE=1+Q-2*C)
        item['qualifies_point_PDS']=bool(name!='B0' and item['delta_cross_CI95'][0]>0 and item['delta_PDS']>=-.005)
        item['qualifies_CI_PDS']=bool(item['qualifies_point_PDS'] and item['delta_PDS_CI95'][0]>=-.005)
        summary[name]=item
        for i in range(N):
            per.append(dict(arm=name,target=targets[i],truth_norm2=yy[i],
                            cross=alpha*s['dot'][i],prediction_norm2=alpha**2*s['energy'][i],
                            pds=s['pds'][i]))
        archives.update({name+'__'+key:s[key] for key in ('gram','dot','energy','pds')})
        small[name]=dict(gram=s['gram'],alpha=alpha)
        rt.status(f'{name}: matched cross={C:.6f}, cross/spread={item["aligned_per_spread"]:.6f}, PDS={item["PDS"]:.6f}, delta cross CI={item["delta_cross_CI95"]}, delta PDS CI={item["delta_PDS_CI95"]}')
    pd.DataFrame(per).to_csv(rt.ROOT/'results/per_target.csv',index=False)
    np.savez_compressed(rt.ROOT/'results/metric_stats.npz',yy=yy,draws=draws,**archives)
    rt.save('results/analysis.json',dict(summary=summary,training_rows=N,
        grids=dict(k=KS,rank=RANKS,lambda_relative=LAMBDAS,blend_weight=WEIGHTS.tolist()),
        pins=PINS,seconds=time.monotonic()-START,free_gib=rt.guard(),
        requested_nice=10,effective_nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0),blas=threadpool_info(),
        heldout79_accessed=False,rendering=False,official_MSE='not measurable locally',
        bootstrap=dict(replicates=2000,seed=20261026,unit='paired targets',
            method='percentile; fixed fitted LOO predictions and retrieval panel; recompute B matching within each resample'),
        normalization='Cross C=dot/T; Q=prediction_energy/T; aligned_per_spread=C/sqrt(Q); B=centered_energy/T. Entire deltas multiplied by sqrt(B0/B). Final normalization descriptive; no labels used to choose scale.',
        selection='Maximize inner OOF cross at matched centered B; jointly tune component hyperparameters and blend weight. Each fit excludes outer t and inner validation j. Inner OOF spread also fixes component/final blend scales.'))
    return preds, summary

def main():
    (rt.ROOT/'results').mkdir(exist_ok=True)
    rt.save('results/run.json',dict(pid=rt.os.getpid(),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),free_gib=rt.guard(True)))
    load()
    K=sources()
    aa,bb=nested(K)
    evaluate(aa,bb)
    rt.status(f'Analysis complete in {time.monotonic()-START:.1f}s; independent verification next.')

if __name__=='__main__':
    main()
