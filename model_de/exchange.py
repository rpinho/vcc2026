"""Equal-depth count exchanges: fixed bulk AND arithmetic CPM means.
Standard margin-preserving exchanges (prior art: Wilcoxon-constrained count
exchanges in COMPETITORS.md); no new-method claim. Most useful when the emitter
uses repeated/quantized depths. Does not quantize or change existing row depths.
"""
import numpy as np

def exchange_equal_depth(x,called,strength=.75,preferred_buffer=None):
    if not 0<=strength<=1:raise ValueError('strength must lie in [0,1]')
    x=np.asarray(x)
    if x.ndim!=2 or not np.issubdtype(x.dtype,np.integer) or (x<0).any():
        raise ValueError('nonnegative integer cell matrix required')
    called=np.asarray(called,bool)
    if called.shape!=(x.shape[1],):raise ValueError('gene mask axis differs')
    y=x.copy();depth=x.sum(1,dtype=np.int64);buffers=np.flatnonzero(~called)
    if not len(buffers):return y
    preferred=preferred_buffer if preferred_buffer is not None and preferred_buffer in buffers else None
    for d in np.unique(depth):
        rows=np.flatnonzero(depth==d)
        if len(rows)<2:continue
        for g in np.flatnonzero(called):
            col=y[rows,g].astype(np.int64);total=int(col.sum())
            desired=(1-strength)*col+strength*total/len(rows)
            goal=np.floor(desired).astype(np.int64);remainder=total-int(goal.sum())
            if remainder:goal[np.argsort(-(desired-goal),kind='stable')[:remainder]]+=1
            delta=col-goal;donors=np.flatnonzero(delta>0);receivers=np.flatnonzero(delta<0)
            for j in receivers:
                for i in donors:
                    while delta[j]<0 and delta[i]>0:
                        donor=int(rows[i]);receiver=int(rows[j])
                        h=preferred if preferred is not None and y[receiver,preferred]>0 else int(buffers[np.argmax(y[receiver,buffers])])
                        units=min(int(delta[i]),int(-delta[j]),int(y[receiver,h]))
                        if units==0:break
                        y[donor,g]-=units;y[receiver,g]+=units
                        y[receiver,h]-=units;y[donor,h]+=units
                        delta[i]-=units;delta[j]+=units
                    if delta[j]==0:break
    return y
