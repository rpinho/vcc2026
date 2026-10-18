"""Conservative public-response transfer and expectation-preserving count transport."""
import os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import functools
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pipeline import transfer_effect, annotation_aliases, NT

VARIANTS=['v1','kd','t010','t025','t050','rob025','rob050','shared','mix010','mix025','mix050','nb025','t100','panelshared','panel010','panel025','panel050','nb050','nb100','pois025','pois050','pois100','center025','center050','center100','nbcenter025','nbcenter050','nbcenter100','eb025','eb050','eb100','nbeb025','nbeb050','nbeb100','gate050','gate100','gateeb050','gateeb100','nbgateeb050','nbgateeb100','fine025','fine050','fine100','nbfine025','nbfine050','nbfine100','fineeb025','fineeb050','fineeb100','nbfineeb025','nbfineeb050','nbfineeb100','finegateeb100','nbfinegateeb100','post050','post100','postfine050','postfine100','postfineeb050','postfineeb100','rank025','rank050','rank100','nbrank025','nbrank050','nbrank100']

def align(s,genes):
    aliases=annotation_aliases(); ix=pd.Index(s['genes']).get_indexer(genes)
    ids={gid:j for j,gid in enumerate(s['gene_ids']) if gid}
    for j in np.flatnonzero(ix<0):ix[j]=ids.get(aliases.get(str(genes[j]),''),-1)
    return ix

_SHARED={}
def shared_effect(sources,genes,panel=None):
    key=(tuple(s['name'] for s in sources),tuple(genes),None if panel is None else tuple(sorted(panel)))
    if key in _SHARED:return _SHARED[key]
    # Equal weight per cell line, not per source file; essential source preferred over GWPS.
    effects=[];families=set()
    for s in sources:
        if s['family'] in families:continue
        selected=np.arange(len(s['targets'])) if panel is None else np.flatnonzero(np.isin(s['targets'],list(panel))&(s['n_cells']>=30))
        precomputed=panel is None and 'shared_logratio_full' in s
        if (int(s['full_target_count']) if precomputed else len(selected))<8:continue
        families.add(s['family']);ix=align(s,genes);ok=ix>=0
        # Arithmetic mean of perturbation/control ratios, self target excluded.
        if precomputed:avg=s['shared_logratio_full']
        else:
            ratios=np.exp2(np.clip(s['lfc'][selected].astype(float),-3,3))
            counts=np.full(len(s['genes']),ratios.shape[0],dtype=float);total=ratios.sum(axis=0)
            gi=pd.Index(s['genes']).get_indexer(s['targets'][selected]);valid=gi>=0
            np.add.at(total,gi[valid],-ratios[np.flatnonzero(valid),gi[valid]])
            np.add.at(counts,gi[valid],-1)
            avg=np.log2(np.maximum(total/np.maximum(counts,1),1e-6))
        e=np.zeros(len(genes));e[ok]=avg[ix[ok]]
        # Do not extrapolate a shared shift from unexpressed source genes.
        e[ok]*=s['baseline'][ix[ok]]/(s['baseline'][ix[ok]]+5)
        effects.append(e)
    result=np.clip(np.mean(effects,axis=0),-.5,.5) if effects else np.zeros(len(genes))
    _SHARED[key]=result
    return result



_FINE=None
_FINE_DIR=None
def fine_sources(sources):
    global _FINE
    if _FINE is None:
        from pipeline import ROOT
        _FINE={p.stem:dict(np.load(p)) for p in (_FINE_DIR or ROOT/'artifacts/fine_signatures').glob('*.npz')}
    aliases=annotation_aliases();result=[]
    for source in sources:
        s=dict(source)
        if s['name'] in _FINE:
            s.update(_FINE[s['name']])
            s['target_ids']=np.array([aliases.get(t,'') for t in s['targets']],dtype=str)
        result.append(s)
    return result


_LR=None
_LR_DIR=None
def rank_sources(sources):
    global _LR
    if _LR is None:
        from pipeline import ROOT
        _LR={p.stem:dict(np.load(p)) for p in (_LR_DIR or ROOT/'artifacts/lowrank_signatures').glob('*.npz')}
    missing={s['name'] for s in sources}-set(_LR)
    if missing:raise FileNotFoundError(f'Build lowrank_signatures.py first: {missing}')
    aliases=annotation_aliases();result=[]
    for source in sources:
        s=dict(source)
        if s['name'] in _LR:
            s.update(_LR[s['name']]);s['target_ids']=np.array([aliases.get(t,'') for t in s['targets']],dtype=str)
        result.append(s)
    return result

_NOISE=None
class _NoiseShrunk:
    def __init__(self,source,variance):self.source=source;self.variance=variance
    def __getitem__(self,key):
        row,columns=key;s=self.source;e=s['lfc'][row,columns].astype(float)
        var=self.variance[columns]*(1/max(float(s['n_cells'][row]),1)+1/max(float(s.get('n_controls',1000)),1))
        weight=np.maximum(0,1-np.divide(var,e**2,out=np.full_like(e,np.inf),where=e!=0))
        return e*weight

def noise_sources(sources):
    global _NOISE
    if _NOISE is None:
        from pipeline import ROOT
        _NOISE=dict(np.load(ROOT/'artifacts/source_noise.npz'))
    result=[]
    for source in sources:
        s=dict(source)
        if s['name'] in _NOISE:s['lfc']=_NoiseShrunk(source,_NOISE[s['name']])
        result.append(s)
    return result

def predict(summary,target,sources,variant):
    from pipeline import predict_mean
    if variant=='v1':return predict_mean(summary,target,sources,'transfer')
    cm=np.asarray(summary['mean'],float);genes=summary['genes'];base=cm/cm.sum()*1e6
    if variant=='kd':e=np.zeros(len(genes))
    else:
        active_sources=rank_sources(sources) if 'rank' in variant else fine_sources(sources) if 'fine' in variant else sources
        e=transfer_effect(noise_sources(active_sources) if 'eb' in variant else active_sources,target,genes,base,'transfer',summary.get('gene_ids'))
        strength={'t010':.1,'t025':.25,'t050':.5,'rob025':.25,'rob050':.5,'mix010':.1,'mix025':.25,'mix050':.5,'nb025':.25,'t100':1.,'panel010':.1,'panel025':.25,'panel050':.5,'nb050':.5,'nb100':1.,'pois025':.25,'pois050':.5,'pois100':1.,'center025':.25,'center050':.5,'center100':1.,'nbcenter025':.25,'nbcenter050':.5,'nbcenter100':1.,'eb025':.25,'eb050':.5,'eb100':1.,'nbeb025':.25,'nbeb050':.5,'nbeb100':1.,'gate050':.5,'gate100':1.,'gateeb050':.5,'gateeb100':1.,'nbgateeb050':.5,'nbgateeb100':1.,'post050':.5,'post100':1.}.get(variant,0)
        if 'fine' in variant or 'rank' in variant:strength=int(variant[-3:])/100
        if variant.startswith('rob'):
            # Robust magnitude screen with source and destination expression support.
            # Same-target effects must clear the quantization/noise floor, not just be nonzero.
            support=np.zeros(len(genes));hits=np.zeros(len(genes))
            for s in sources:
                if target not in s['targets']:continue
                row=list(s['targets']).index(target);ix=align(s,genes);ok=ix>=0
                supported=np.zeros(len(genes),bool)
                supported[ok]=(s['baseline'][ix[ok]]>=5)&(np.abs(s['lfc'][row,ix[ok]])>=.24)
                support+=supported;hits+=ok
            e*= (base>=5)&(np.abs(e)>=.24)&(support>=np.maximum(1,np.ceil(hits*.5)))
        if 'center' in variant:
            e=strength*(e-shared_effect(sources,genes))
        elif variant.startswith('panel'):
            shared=shared_effect(sources,genes,summary['panel'])
            e=shared+strength*(e-shared)
        elif variant.startswith('mix'):
            shared=shared_effect(sources,genes)
            e=shared+strength*(e-shared)
        elif variant=='shared':e=shared_effect(sources,genes)
        else:e*=strength
    if 'gate' in variant:
        jdest=np.flatnonzero(genes==target)
        if len(jdest):
            activity=base[jdest[0]]/(base[jdest[0]]+5);source_activities=[];used=set()
            for src in sources:
                if src['family'] in used or target not in src['targets']:continue
                jj=np.flatnonzero(src['genes']==target)
                if len(jj):
                    val=src['baseline'][jj[0]];source_activities.append(val/(val+5));used.add(src['family'])
            if source_activities:e*=min(1,activity/max(float(np.mean(source_activities)),.05))
    mu=cm*np.exp2(np.clip(e,-1.5,1.5))
    idx=np.flatnonzero(genes==target)
    if not len(idx):
        aliases=annotation_aliases();tid=aliases.get(target,'');ids=summary.get('gene_ids',np.array([aliases.get(g,'') for g in genes]));idx=np.flatnonzero(ids==tid) if tid else []
    # Preserve total expected counts while keeping target EXACTLY 10% of control.
    if len(idx):
        j=idx[0];mu[j]=0;mu*= (cm.sum()-.1*cm[j])/mu.sum();mu[j]=.1*cm[j]
    else:mu*=cm.sum()/mu.sum()
    return mu

def transport(mu,controls,cm,n,rng):
    """E[Y_g]=mu_g, no cutoff and no random vector normalization.
    For sampled raw control X, thin decreases and Poisson-augment increases proportional
    to X. For cm_g=0, add independent counts scaled by sampled library. Mixtures preserve
    empirical covariance/zero inflation. Library totals vary naturally after perturbation.
    """
    original=controls[rng.integers(0,controls.shape[0],n)].toarray().astype(np.int64)
    ratio=np.divide(mu,cm,out=np.zeros_like(mu),where=cm>0)
    y=rng.binomial(original,np.minimum(ratio,1)[None,:])
    up=ratio>1
    if up.any():y[:,up]+=rng.poisson(original[:,up]*(ratio[up]-1)[None,:])
    new=(cm==0)&(mu>0)
    if new.any():y[:,new]+=rng.poisson(original.sum(axis=1)[:,None]/cm.sum()*mu[new][None,:])
    assert np.all(y.sum(axis=1)<=1000000)
    return sp.csr_matrix(y.astype(np.int32))

def nb_counts(mu,libs,theta,n,rng):
    # Independent gamma-Poisson has EXACT target expectation; no normalization after gamma.
    depth=rng.choice(libs,n)/np.mean(libs)
    rates=rng.gamma(theta[None,:],mu[None,:]/theta[None,:],size=(n,len(mu)))*depth[:,None]
    return sp.csr_matrix(rng.poisson(rates).astype(np.int32))

_POST_STATE=None
def posterior_calibration(controls,cm,theta):
    """Calibrate the empirical-pool average Gamma posterior mean to cm exactly."""
    global _POST_STATE
    if _POST_STATE is not None:
        previous,oldcm,oldtheta,cal=_POST_STATE
        if previous is controls and np.array_equal(oldcm,cm) and np.array_equal(oldtheta,theta):return cal
    libs=np.asarray(controls.sum(axis=1)).ravel();meanlib=float(libs.mean());sums=np.zeros(len(cm))
    for lo in range(0,controls.shape[0],128):
        x=controls[lo:lo+128].toarray();mean=cm[None,:]*libs[lo:lo+128,None]/meanlib
        sums+=((theta[None,:]+x)*mean/(theta[None,:]+mean)).sum(axis=0)
    average=sums/controls.shape[0];cal=np.divide(cm,average,out=np.ones_like(cm),where=average>0)
    _POST_STATE=(controls,cm.copy(),theta.copy(),cal)
    return cal

def posterior_transport(mu,controls,cm,theta,n,rng):
    """Empirical thinning plus Gamma-posterior augmentation with corrected mean.
    Under a true NB control model this reduces to NB scaling, while the pool-level
    calibration also preserves means when empirical controls deviate from that model.
    """
    calibration=posterior_calibration(controls,cm,theta)
    x=controls[rng.integers(0,controls.shape[0],n)].toarray().astype(np.int64);libs=x.sum(axis=1)
    ratio=np.divide(mu,cm,out=np.zeros_like(mu),where=cm>0);y=rng.binomial(x,np.minimum(ratio,1)[None,:]);up=ratio>1
    if up.any():
        cellmean=cm[up][None,:]*libs[:,None]/cm.sum();shape=theta[up][None,:]+x[:,up]
        latent=rng.gamma(shape,cellmean/(theta[up][None,:]+cellmean)*calibration[up][None,:])
        y[:,up]+=rng.poisson((ratio[up]-1)[None,:]*latent)
    new=(cm==0)&(mu>0)
    if new.any():y[:,new]+=rng.poisson(libs[:,None]/cm.sum()*mu[new][None,:])
    return sp.csr_matrix(y.astype(np.int32))
