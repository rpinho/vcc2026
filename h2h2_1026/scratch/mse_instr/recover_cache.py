"""Recover retained moments by exact stored-score identities, never fit predictions."""
import runtime_local as j
j.guard(launch=True)
import numpy as np, csv
from pathlib import Path

def tidy(arm):
    p=j.VCC/'final_push/artifacts/de_check/H1_rich_search'/arm/'tidy.csv'
    rows={}
    for r in csv.DictReader(p.open()):
        rows.setdefault(r['metric'],{})[r['perturbation']]=float(r['value'])
    return rows

if __name__=='__main__':
    z=np.load(j.ROOT/'reference.npz'); ts=z['targets'][:-1];r=z['means'][:-1];c=z['jk'][:-1];genes=z['genes']
    own={i:int(np.flatnonzero(genes==t)[0]) for i,t in enumerate(ts) if np.any(genes==t)}
    g=np.full(len(ts),len(genes),float)
    for i in own:g[i]-=1
    arms=['v5b','v8b_like','levers_knob_none']
    wanted={a:np.array([tidy(a)['expr_mse_unbiased'][t] for t in ts]) for a in arms}
    found={a:{} for a in arms};sources={a:{} for a in arms}
    def batch(items):
        if not items:return
        m=np.stack([v[1] for v in items]);p=np.array([v[2] for v in items])
        e=(m*m).sum(1)[:,None]+(r*r).sum(1)[None,:]-2*m @ r.T
        for i,col in own.items():e[:,i]-=(m[:,col]-r[i,col])**2
        values=(e-p[:,None]-c[None,:])/g[None,:]
        for a in arms:
            rr,cc=np.where(abs(values-wanted[a][None,:])<1e-12)
            for k,i in zip(rr,cc):
                exact=float(((m[k]-r[i])**2).sum())
                if i in own:exact-=float(m[k,own[i]]-r[i,own[i]])**2
                if abs((exact-p[k]-c[i])/g[i]-wanted[a][i])<1e-13:
                    found[a][str(ts[i])]=(m[k],p[k]);sources[a][str(ts[i])]=str(items[k][0])
    items=[]
    files=sorted((j.VCC/'rsi/artifacts/stat_cache').glob('*.npz'))
    for k,path in enumerate(files):
        if k%1000==0:j.guard();j.status(f'Scanning retained moment cache {k}/{len(files)}; matches '+str({a:len(v) for a,v in found.items()}))
        with np.load(path) as d:
            if d['mean'].shape!=r.shape[1:] or int(d['count'])!=400:continue
            items.append((path,d['mean'],float(d['jk'])))
        if len(items)>=128:batch(items);items=[]
    batch(items)
    for a,values in found.items():
        order=sorted(values)
        if order:np.savez_compressed(j.ROOT/f'cached_{a}.npz',targets=order,means=np.stack([values[t][0] for t in order]),jk=[values[t][1] for t in order])
    j.save('cache_recovery.json',dict(matches={a:len(v) for a,v in found.items()},sources=sources,
        rule='Stored expr_mse_unbiased identity, residual <1e-13, count=400, same gene dimension; not a model fit.'))
    j.status('Retained-moment scan complete: '+str({a:len(v) for a,v in found.items()}))
