"""Descriptive paired bootstrap for arbitrary arm pairs, using analyze.py's registered machinery verbatim
(aggregate, same 20,000 indices from seed 20261026, official score_one scaling). Usage: paired_boot.py A B [A B ...]
Writes results/paired_boot_<A>_<B>.json. Does not alter analysis.json or the preregistered verdict."""
import sys, numpy as np, polars as pl
import analyze as an
j=an.j; SIX=an.SIX
def sample(arm,indices):
    rec=j.json.loads((j.ROOT/f'results/official_{arm}.json').read_text())
    with np.load(j.ROOT/f'results/target_metrics_{arm}.npz') as z: cols=z['columns'].tolist(); mat=z['values']
    raw=an.aggregate(cols,mat); assert all(abs(raw[m]-rec['raw'][m])<1e-12 for m in SIX)
    an_=rec['anchors']; entries=an._replicate_entries(an_['baseline'],pl.DataFrame({'metric':list(an_['replicate']),'replicate':list(an_['replicate'].values())}))
    vec={m:np.empty(an.N) for m in SIX}
    for lo in range(0,an.N,200):
        v=an.aggregate(cols,mat,indices[lo:lo+200])
        for m in SIX: vec[m][lo:lo+200]=v[m]
    sc=np.column_stack([[an.score_one(float(x),entries[m].base,entries[m].scoring) for x in vec[m]] for m in SIX])
    return rec,{SIX[3]:vec[SIX[3]],SIX[1]:vec[SIX[1]],'score_avg':sc.mean(1)}
idx=np.random.default_rng(an.SEED).integers(79,size=(an.N,79)); cache={}
args=sys.argv[1:]
for a,b in zip(args[::2],args[1::2]):
    for x in (a,b):
        if x not in cache: cache[x]=sample(x,idx)
    (ra,sa),(rb,sb)=cache[a],cache[b]; c={}
    for m in [SIX[3],SIX[1],'score_avg']:
        point=ra['score_avg']-rb['score_avg'] if m=='score_avg' else ra['raw'][m]-rb['raw'][m]
        d=sa[m]-sb[m]; lo,hi=np.quantile(d,[.025,.975],method='linear'); c[m]=dict(delta=float(point),ci95=[float(lo),float(hi)],lower_bound_positive=bool(lo>0),p_delta_le_0=float(np.mean(d<=0)))
    c['raw_mse_ratio']=ra['raw']['expr_mse_unbiased_capped_norm']/rb['raw']['expr_mse_unbiased_capped_norm']
    j.save(f'results/paired_boot_{a}_{b}.json',dict(pair=f'{a} minus {b}',draws=an.N,seed=an.SEED,comparison=c,note=__doc__)); print(a,'-',b,j.json.dumps(c))
