"""Recipient family and outer fold never enter any donor feature or normalizer."""
import common as rt
import numpy as np
from collections import Counter

class Sources:
    def __init__(self, excluded=()):
        self.excluded=set(excluded);self.banks=[]
        for path in sorted((rt.ROOT/'data').glob('*.npz')):
            z=np.load(path); family=str(z['family'])
            if family in self.excluded: continue
            d={key:z[key] for key in z.files};d['name']=path.stem;d['family']=family
            d['index']={t:i for i,t in enumerate(d['targets'])}
            d['mean_response']=d['response'].mean(0);d['mean_control']=d['control'].mean(0)
            self.banks.append(d)
        assert self.banks
        self.genes=self.banks[0]['genes'];self.g=len(self.genes)
        assert all(np.array_equal(d['genes'],self.genes) for d in self.banks)
        self.gi={t:i for i,t in enumerate(self.genes)}
        targets=sorted({t for d in self.banks for t in d['targets']});self.ti={t:i for i,t in enumerate(targets)}
    def features(self, target, control, recipient=None):
        donors=[(d,d['index'][target]) for d in self.banks if d['family']!=recipient and target in d['index']]
        if not donors:
            # Shared source response fallback, computed solely from eligible donors.
            donors=[(d,None) for d in self.banks if d['family']!=recipient]
        if not donors: raise ValueError('No donors after context exclusion')
        counts=Counter(d['family'] for d,_ in donors)
        denom=np.zeros(self.g,np.float32);total=np.zeros_like(denom);sq=np.zeros_like(denom)
        compat=np.zeros_like(denom);cd=np.zeros_like(denom)
        dst=np.log1p(50000*control).astype(np.float32)
        for d,i in donors:
            measured=d['measured']
            r=d['response'][i] if i is not None else d['mean_response']
            ctrl=d['control'][i] if i is not None else d['mean_control']
            similarity=float(np.exp(-np.mean(np.abs(dst[measured]-np.log1p(50000*ctrl[measured])))))
            w=1/counts[d['family']];m=measured*w
            total+=r*m;sq+=r*r*m;denom+=m;compat+=r*m*similarity;cd+=m*similarity
        pooled=total/np.maximum(denom,1e-8)
        sd=np.sqrt(np.maximum(sq/np.maximum(denom,1e-8)-pooled**2,0))
        weighted=compat/np.maximum(cd,1e-8)
        support=np.minimum(denom/len(counts),1)
        own=self.gi.get(target)
        if own is not None:
            # Remove own-target shortcuts from response channels and all losses.
            pooled[own]=0;weighted[own]=0;sd[own]=0
        return np.concatenate([pooled,weighted,sd,dst/5,support]).astype(np.float32)
    def rows(self, families=None):
        return [(b,i) for b,d in enumerate(self.banks) if families is None or d['family'] in families for i in range(len(d['targets']))]
    def batch(self, rows):
        xs=[];ys=[];ms=[];ids=[]
        for b,i in rows:
            d=self.banks[b];t=d['targets'][i]
            xs.append(self.features(t,d['control'][i],d['family']))
            y=d['response'][i];mask=d['measured'].astype(np.float32)
            own=self.gi.get(t)
            if own is not None: mask[own]=0
            ys.append(y);ms.append(mask);ids.append(self.ti[t])
        return np.stack(xs),np.stack(ys),np.stack(ms),np.array(ids,dtype=np.int64)
