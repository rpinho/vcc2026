"""Label-free DE selection and bulk-preserving selective dispersion.

Inputs are public source effects/uncertainties and destination controls/predictions.
No reference treated cells are accepted by this module. Ideas: AtlasShift, SIVC,
Chen Zhang, eye-contact; see PRIOR_ART.md. No unlicensed source code used.
"""
import numpy as np

def rank_sources(sources, predicted_lfc, eligible, genes=None, targets=None):
    """One vote per source family (caller must combine screens of a family first).
    Return confidence, signed agreement, consensus direction, source-derived budgets.
    Budgets are median source significant counts on the destination control gate;
    nominal |z|>1.96 is an uncorrected source proxy, not a BH call or calibrated posterior.
    """
    shape=predicted_lfc.shape
    pos=np.zeros(shape);neg=np.zeros(shape);strength=np.zeros(shape);counts=[]
    for effect, available, variance in sources:
        good=available & np.isfinite(variance) & (variance>0) & np.isfinite(effect)
        z=np.divide(effect,np.sqrt(np.maximum(variance,1e-12)),out=np.zeros(shape),where=good)
        sig=good & (np.abs(z)>1.96) & eligible[None,:]
        pos+=sig & (z>0);neg+=sig & (z<0)
        strength+=np.where(sig,np.minimum(np.abs(z),20),0)
        ct=sig.sum(1).astype(float);ct[~good.any(1)]=np.nan;counts.append(ct)
    consensus=np.sign(pos-neg)
    agree=np.where(predicted_lfc>=0,pos-neg,neg-pos)
    confidence=agree*1e6+strength
    confidence[:,~eligible]=-np.inf
    if genes is not None and targets is not None:
        gi={g:j for j,g in enumerate(genes)}
        for i,t in enumerate(targets):
            if t in gi: confidence[i,gi[t]]=-np.inf
    budgets=np.full(shape[0],300,dtype=int)
    if counts:
        for i in range(shape[0]):
            c=np.array(counts)[:,i];c=c[np.isfinite(c)]
            if c.size: budgets[i]=int(np.clip(np.median(c),50,1000))
    return confidence,agree,consensus,budgets

def select(confidence, agreement, k=300, min_net=1):
    out=np.zeros(confidence.shape,bool)
    for i in range(len(out)):
        ix=np.flatnonzero(np.isfinite(confidence[i]) & (agreement[i]>=min_net))
        n=int(k if np.isscalar(k) else k[i])
        out[i,ix[np.argsort(-confidence[i,ix],kind='stable')[:n]]]=True
    return out

def selective_template(realistic, pooled, called):
    """Use pool4 on supported genes, pool1 elsewhere. Renormalize each cell.
    Pass through an existing exact-total emitter afterwards. 'called' is a desired
    support mask, not a promise of actual Wilcoxon significance.
    """
    q=np.where(called[None,:],pooled,realistic)
    return q/np.maximum(q.sum(1,keepdims=True),1e-30)

def redistribute(x, called, strength=0.75):
    """Compress selected columns toward depth-proportional counts, preserving ALL
    gene totals exactly. Recompute row depths before DE; CPM is compositional.
    A separate verifier checks totals/validity. No bulk-accuracy/PDS change in exact
    arithmetic; jackknife/MSE correction and CPM means can still change.
    """
    if not 0<=strength<=1: raise ValueError('strength must be in [0,1]')
    y=np.asarray(x).copy();depth=y.sum(1).astype(float);w=depth/depth.sum()
    for j in np.flatnonzero(called):
        total=int(y[:,j].sum());a=(1-strength)*y[:,j]+strength*w*total
        col=np.floor(a).astype(y.dtype);rem=total-int(col.sum())
        if rem: col[np.argsort(-(a-col),kind='stable')[:rem]]+=1
        y[:,j]=col
    return y
