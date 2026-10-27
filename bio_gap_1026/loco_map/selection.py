"""Nested selection, with labels supplied only by inner validation recipients."""
import runtime as rt
import numpy as np
from models import HYPERS, WEIGHTS, mask_own

def spread(q):
    q=q.astype(float)
    return max(float(np.sum(q*q)-q.sum(0)@q.sum(0)/len(q)),1e-30)

def ranks(s):
    n=len(s); out=np.empty(n)
    for lo in range(0,n,128):
        hi=min(n,lo+128); block=s[lo:hi]
        diag=block[np.arange(hi-lo),np.arange(lo,hi)][:,None]
        out[lo:hi]=1-((block>diag).sum(1)+.5*((block==diag).sum(1)-1))/(n-1)
    return out

def blend_pds(a,b,w,scale):
    n=len(a); out=np.empty(n)
    for lo in range(0,n,128):
        hi=min(n,lo+128)
        block=(1-w)*a[lo:hi]+w*scale*b[lo:hi]
        diag=block[np.arange(hi-lo),np.arange(lo,hi)][:,None]
        out[lo:hi]=1-((block>diag).sum(1)+.5*((block==diag).sum(1)-1))/(n-1)
    return out

class Validation:
    def __init__(self,data,r,p):
        self.r=r; self.targets=data.recips[r]['targets']; self.genes=data.genes
        self.p=p.astype(float); self.y=data.recips[r]['truth'].astype(float)
        self.n=len(p); self.T=float(np.sum(self.y*self.y)); self.B=spread(self.p)
        self.dot=float(np.sum(self.p*self.y))
        self.keep=~np.isin(self.genes,self.targets)
        self.yn=self.y[:,self.keep]/np.maximum(np.linalg.norm(self.y[:,self.keep],axis=1)[:,None],1e-30)
        self.baseS=self.p[:,self.keep]@self.yn.T
        self.basepds=ranks(self.baseS)
        self.basecross=self.dot/self.T

    def grid(self,candidates):
        out=[]
        for h,(L,W) in enumerate(candidates):
            rt.guard()
            if h in (0,3,6): rt.status(f'Inner {self.r}: scoring hyperparameters {h+1}-{h+3}/9 against plain transfer.')
            q=mask_own(L @ W,self.targets,self.genes)
            Bq=spread(q); scale=np.sqrt(self.B/Bq)
            dot=float(np.sum(q*self.y))
            pq=float(np.sum(self.p*q)-self.p.sum(0)@q.sum(0)/self.n)
            sim=L@(W[:,self.keep]@self.yn.T)
            for w in WEIGHTS:
                Bblend=(1-w)**2*self.B+2*(1-w)*w*scale*pq+w*w*scale*scale*Bq
                C=float(((1-w)*self.dot+w*scale*dot)*np.sqrt(self.B/max(Bblend,1e-30))/self.T)
                pds=self.basepds if w==0 else blend_pds(self.baseS,sim,w,scale)
                dp=pds-self.basepds; mean=float(dp.mean()); se=float(dp.std(ddof=1)/np.sqrt(self.n))
                out.append(dict(recipient=self.r,h=h,weight=w,cross=C,delta_cross=C-self.basecross,
                    PDS=float(pds.mean()),delta_PDS=mean,PDS_lower_normal=mean-1.96*se,component_scale=float(scale)))
            del q,sim
        return out

def choose(records):
    contexts=sorted(set(x['recipient'] for x in records)); pure=[]; bounded=[]
    for h in range(len(HYPERS)):
        for w in WEIGHTS:
            rr=[x for x in records if x['h']==h and x['weight']==w]
            assert len(rr)==len(contexts)
            score=float(np.mean([x['cross'] for x in rr]))
            entry=dict(h=h,lambda_relative=HYPERS[h][0],rank=HYPERS[h][1],weight=w,
                       training_macro_cross=score,training_recipients=contexts,
                       minimum_training_PDS_bound=min(x['PDS_lower_normal'] for x in rr))
            if w==1: pure.append(entry)
            if all(x['PDS_lower_normal']>=-.005-1e-14 for x in rr): bounded.append(entry)
    # Fixed first-in-grid tie resolution (weight zero appears before positive weights).
    best=lambda items:max(enumerate(items),key=lambda pair:(round(pair[1]['training_macro_cross'],14),-pair[0]))[1]
    return dict(map=best(pure),bounded=best(bounded))
