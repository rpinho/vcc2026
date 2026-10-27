"""Label-safe construction and ridge reduced-rank models."""
import runtime as rt
import gc
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from sklearn.utils.extmath import randomized_svd

RECIPIENTS=('H1','K562','RPE1','HepG2','Jurkat')
RECIP_BANK=dict(H1='H1_full',K562='K562_GWPS_full',RPE1='RPE1',HepG2='HepG2',Jurkat='Jurkat')
HYPERS=[(lam,rank) for lam in (.1,1.,10.) for rank in (20,80,0)]
WEIGHTS=(0.,.005,.01,.02,.03,.05,.1,.2,.5,1.)

def mask_own(q,targets,genes):
    q=np.array(q,copy=True)
    ix=pd.Index(genes).get_indexer(targets); ok=ix>=0
    q[np.flatnonzero(ok),ix[ok]]=0
    return q

class Data:
    def __init__(self):
        rt.guard(True)
        manifest_path=rt.REPO/'model_billmB/data/manifest.json'
        manifest=rt.json.loads(manifest_path.read_text()); self.pins={str(manifest_path):rt.sha(manifest_path)}
        metadata=rt.ROOT.parent/'results/provenance.json'
        forbidden=set(rt.json.loads(metadata.read_text())['heldout_targets_metadata_only'])
        assert len(forbidden)==79
        self.pins[str(metadata)]=rt.sha(metadata)
        names=sorted(set(manifest)-{'H1_cached'}); common=None
        for name in names:
            with np.load(rt.REPO/f'model_billmB/data/{name}.npz') as z:
                gs=set(z['genes'][z['measured']]); common=gs if common is None else common&gs
        with np.load(rt.ROOT.parent/'results/baseline.npz') as z:
            self.h1targets=z['targets']; self.fullgenes=z['genes']; self.fullY=z['truth'].astype(np.float32)
            self.fullB0=z['realized'].astype(np.float32)
        self.genes=np.array([g for g in self.fullgenes if g in common])
        assert len(self.genes)==6006
        self.banks={}
        for name in names:
            rt.guard()
            path=rt.REPO/f'model_billmB/data/{name}.npz'; digest=rt.sha(path)
            assert digest==manifest[name]['sha256']; self.pins[str(path)]=digest
            with np.load(path) as z:
                ix=pd.Index(z['genes']).get_indexer(self.genes)
                targets=z['targets']; response=z['response'][:,ix].astype(np.float32)
                controls=z['control'][:,ix].astype(np.float32)
            allow=~np.isin(targets,list(forbidden))
            targets=targets[allow]; response=response[allow]; controls=controls[allow]
            assert not set(targets)&forbidden
            if name=='H1_full': assert np.array_equal(targets,self.h1targets)
            self.banks[name]=dict(name=name,family=manifest[name]['family'],targets=targets,
                response=response,controls=controls,mean_response=response.mean(0,dtype=np.float64).astype(np.float32),
                mean_control=controls.mean(0,dtype=np.float64).astype(np.float32),index=pd.Index(targets))
            rt.status(f'Pinned {name}: {len(targets)} rows, {len(self.genes)} common genes.')
            gc.collect()
        self.recips={}
        for r in RECIPIENTS:
            b=self.banks[RECIP_BANK[r]]
            ctrl=b['mean_control'].astype(float)
            w=np.sqrt(1+np.log1p(50000*ctrl)); w/=np.sqrt(np.mean(w*w))
            self.recips[r]=dict(targets=b['targets'],control=ctrl,weight=w,
                               truth=mask_own(b['response'],b['targets'],self.genes))
        ix=pd.Index(self.fullgenes).get_indexer(self.genes)
        assert np.allclose(self.recips['H1']['truth'],mask_own(self.fullY[:,ix],self.h1targets,self.genes),atol=1e-7)
        self.B0=mask_own(self.fullB0[:,ix],self.h1targets,self.genes)
        self.fullY=mask_own(self.fullY,self.h1targets,self.fullgenes)
        self.fullB0=mask_own(self.fullB0,self.h1targets,self.fullgenes)
        self.fullix=ix
        self.lookups={}; self.similarities={}
        # All of these weights use control probabilities only.
        for r in RECIPIENTS:
            dst=np.log1p(50000*self.recips[r]['control'])
            for name,b in self.banks.items():
                self.lookups[r,name]=b['index'].get_indexer(self.recips[r]['targets'])
                self.similarities[r,name]=np.exp(-np.mean(np.abs(np.log1p(50000*b['controls'])-dst),axis=1))
        # Row control matrices are no longer needed once control-only weights are fixed.
        for b in self.banks.values(): del b['controls']
        gc.collect()
        rt.save('results/input_receipt.json',dict(pins=self.pins,genes=self.genes.tolist(),
            recipients={r:dict(bank=RECIP_BANK[r],n=len(self.recips[r]['targets'])) for r in RECIPIENTS},
            H1_cached_opened=False,H1_heldout79_opened=False,heldout79_identities_excluded_every_bank=True))

    def plain(self,r,excluded=(),return_family=False):
        """All excludes are family-level, including every K562 variant."""
        excluded=set(excluded)|{r}; targets=self.recips[r]['targets']; n=len(targets); g=len(self.genes)
        banks=[b for b in self.banks.values() if b['family'] not in excluded]
        assert banks and not any(b['family'] in excluded for b in banks)
        counts={f:np.zeros(n,int) for f in {b['family'] for b in banks}}
        for b in banks: counts[b['family']]+=self.lookups[r,b['name']]>=0
        sums={f:np.zeros((n,g),np.float32) for f in counts} if return_family else None
        nums=np.zeros((n,g),float); denom=np.zeros(n,float); coverage=np.zeros(n,int)
        for b in banks:
            ix=self.lookups[r,b['name']]; good=np.flatnonzero(ix>=0)
            coverage[good]+=1
            if len(good):
                raw=b['response'][ix[good]]
                w=self.similarities[r,b['name']][ix[good]]/counts[b['family']][good]
                nums[good]+=raw*w[:,None]; denom[good]+=w
                if return_family: sums[b['family']][good]+=raw/counts[b['family']][good,None]
        no=np.flatnonzero(coverage==0)
        if len(no):
            famcounts={f:sum(b['family']==f for b in banks) for f in counts}
            dst=np.log1p(50000*self.recips[r]['control']); tot=np.zeros(g); den=0.
            for b in banks:
                w=np.exp(-np.mean(np.abs(dst-np.log1p(50000*b['mean_control']))))/famcounts[b['family']]
                tot+=w*b['mean_response']; den+=w
            nums[no]=tot; denom[no]=den
        out=mask_own(nums/np.maximum(denom[:,None],1e-30),targets,self.genes).astype(np.float32)
        if return_family:
            # Kernel construction excludes the WHOLE recipient target panel from every source vector.
            keep=~np.isin(self.genes,self.h1targets)
            for f in sums: sums[f][:,~keep]=0
            return out,coverage,sums
        return out,coverage

def features(p,weight,conditional):
    x=p.astype(float)
    if conditional: x*=weight
    x/=np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-30)
    return x

def train_map(data,recipients,excluded,conditional):
    """Only explicitly supplied training recipients may contribute labels."""
    assert not set(recipients)&set(excluded)
    rt.guard(True)
    total=sum(len(data.recips[r]['targets']) for r in recipients)
    x=np.empty((total,len(data.genes)),np.float32); y=np.empty_like(x)
    weights=[]; record=[]; offset=0
    for r in recipients:
        p,coverage=data.plain(r,excluded)
        n=len(p); rows=slice(offset,offset+n)
        x[rows]=features(p,data.recips[r]['weight'],conditional)
        y[rows]=data.recips[r]['truth']
        if conditional: y[rows]/=data.recips[r]['weight']
        x[rows]/=np.sqrt(n); y[rows]/=np.sqrt(n)
        record.append(dict(recipient=r,rows=n,source_exclusions=sorted(set(excluded)|{r}),covered=int((coverage>0).sum())))
        offset+=n
        del p
    u,s,vt=randomized_svd(x,n_components=min(256,min(x.shape)-1),n_oversamples=16,n_iter=2,random_state=20261026,flip_sign=True)
    retained=float((s @ s)/np.sum(x*x))
    # Reproject instead of assuming the approximate randomized U,S are exact singular factors.
    z=(x @ vt.T).astype(float)
    # Float32 storage, double precision sufficient-statistic arithmetic in gene blocks.
    gram=z.T @ z; rhs=np.empty((z.shape[1],y.shape[1]),float)
    for lo in range(0,y.shape[1],256): rhs[:,lo:lo+256]=z.T @ y[:,lo:lo+256].astype(float)
    d,U=eigh(gram,check_finite=False); d=np.maximum(d,0)
    model=dict(basis=vt.astype(float),d=d,U=U,rhs=rhs,penalty_scale=float(np.trace(gram)/len(d)),
               conditional=conditional,retained_energy=retained,training=record)
    from verify import verify_fit
    model['fit_verification']=verify_fit(z,y,model)
    return model

def predict_candidates(model,p,weight,targets,genes):
    x=features(p,weight,model['conditional'])@model['basis'].T
    d=model['d']; U=model['U']; rhs=model['rhs']; out=[]
    for lam in (.1,1.,10.):
        penalty=lam*model['penalty_scale']; invroot=(U/np.sqrt(d+penalty))@U.T
        M=invroot @ rhs; evals,V=eigh(M @ M.T,check_finite=False)
        left=x @ invroot
        for rank in (20,80,0):
            if rank:
                v=V[:,-min(rank,len(d)):]; L=left @ v; W=v.T @ M
            else: L=left; W=M
            if model['conditional']: W=W*weight
            out.append((L,W))
    return out

def direct_h1(data,r):
    """Exact dual ridge/RRR kernel, with outer recipient absent even from source features."""
    assert r!='H1'
    _,_,trainfam=data.plain('H1',{r},True)
    _,_,testfam=data.plain(r,{'H1'},True)
    n=len(data.h1targets); m=len(data.recips[r]['targets'])
    k=np.zeros((n,n)); query=np.zeros((m,n)); tn=np.zeros(n); qn=np.zeros(m)
    assert set(trainfam)==set(testfam)
    for f in trainfam:
        a=trainfam[f].astype(float); b=testfam[f].astype(float)
        k+=a @ a.T; query+=b @ a.T; tn+=np.sum(a*a,1); qn+=np.sum(b*b,1)
    norms=np.sqrt(tn); qnorms=np.sqrt(qn)
    k/=np.maximum(norms[:,None]*norms[None,:],1e-30)
    query/=np.maximum(qnorms[:,None]*norms[None,:],1e-30)
    d,U=eigh(k,check_finite=False); good=d>max(1e-12,d.max()*1e-11); d=d[good];U=U[:,good]
    y=data.recips['H1']['truth'].astype(float); penalty=.1*max(np.trace(k)/n,1e-12)
    D=np.sqrt(d/(d+penalty)); yy=U.T@(y @ y.T)@U
    _,V=eigh(D[:,None]*yy*D[None,:],check_finite=False); v=V[:,-min(40,len(d)):]
    coef=(((query @ U/np.sqrt(d*(d+penalty)))@v)@v.T*D)@U.T
    q=mask_own(coef @ y,data.recips[r]['targets'],data.genes).astype(np.float32)
    return q,dict(lambda_relative=.1,rank=40,training_recipient='H1',source_excluded=['H1',r],
                  training_targets=168,zero_query_features=int((qn==0).sum()))
