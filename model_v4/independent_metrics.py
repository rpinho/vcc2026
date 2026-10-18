"""Independent direct NumPy implementation of official bulk MSE and PDS definitions.
No cell_eval2 code is called here. Finite precision comparison is not an exact improvement proof.
"""
import numpy as np
from scipy import sparse
NT='non-targeting'

def profiles(a):
    labels=a.obs.target_gene.astype(str).to_numpy();out={}
    for t in sorted(set(labels)):
        x=a.X[labels==t].astype(np.float64);tot=np.asarray(x.sum(axis=0)).ravel();S=tot.sum()
        bulk=np.log1p(50000*tot/S);n=x.shape[0];C=0.
        if n>=2:
            # Welford in float64; summing squared deviations avoids catastrophic cancellation.
            avg=np.zeros(x.shape[1]);ss=np.zeros(x.shape[1]);k=0
            for lo in range(0,n,64):
                z=x[lo:lo+64].toarray();den=S-z.sum(axis=1)
                vals=np.log1p(50000*(tot[None,:]-z)/den[:,None])
                for v in vals:
                    k+=1;d=v-avg;avg+=d/k;ss+=d*(v-avg)
            C=float(ss.sum()*(n-1)/n)
        out[t]=(bulk,C)
    return out

def verify(pred,real,official,exclusion_panel=None):
    pp=profiles(pred);rr=profiles(real);genes=real.var_names.to_numpy();targets=sorted(set(rr)-{NT});G=len(genes);P=np.stack([pp[t][0] for t in targets]);R=np.stack([rr[t][0] for t in targets]);ctrl=rr[NT][0]
    keep=np.ones(P.shape,bool)
    for i,t in enumerate(targets):keep[i,genes==t]=False
    weights=keep/keep.sum(axis=1)[:,None]
    avg=(P*weights).sum(axis=0)/weights.sum(axis=0)
    spread=float(((P-avg)**2*weights).sum())
    cap=np.array([min(pp[t][1],rr[t][1])/keep[i].sum() for i,t in enumerate(targets)])
    rho=min(1,spread/cap.sum()) if cap.sum()>0 else 1
    nums=[];dens=[]
    for i,t in enumerate(targets):
        mask=keep[i];d=mask.sum()
        nums.append(np.sum((P[i,mask]-R[i,mask])**2)/d-rho*cap[i]-rr[t][1]/d)
        dens.append((np.sum((R[i,mask]-ctrl[mask])**2)-rr[t][1]-rr[NT][1])/d)
    mse=float(np.sum(nums)/np.sum(dens))
    pkeep=~np.isin(genes,targets);dp=P[:,pkeep]-ctrl[pkeep];dr=R[:,pkeep]-ctrl[pkeep]
    pn=np.linalg.norm(dp,axis=1);rn=np.linalg.norm(dr,axis=1);den=pn[:,None]*rn[None,:]
    cos=np.divide(dp@dr.T,den,out=np.zeros(den.shape),where=den>0);dist=1-cos;rank=[]
    for i,row in enumerate(dist):rank.append(1-((row<row[i]).sum()+.5*((row==row[i]).sum()-1))/(len(targets)-1))
    pds=float(np.mean(rank));got={'pds_cosine':pds}
    if 'expr_mse_unbiased_capped_norm' in official:got['expr_mse_unbiased_capped_norm']=mse
    err={m:abs(v-official[m]) for m,v in got.items()}
    if 'expr_mse_unbiased_capped_norm' in err:assert err['expr_mse_unbiased_capped_norm']<2e-6,err
    assert err['pds_cosine']<1e-7,err
    extra={}
    if exclusion_panel is not None:
        mask=~np.isin(genes,exclusion_panel);a=P[:,mask]-ctrl[mask];b=R[:,mask]-ctrl[mask];den=np.linalg.norm(a,axis=1)[:,None]*np.linalg.norm(b,axis=1)[None,:];similarity=np.divide(a@b.T,den,out=np.zeros(den.shape),where=den>0);distance=1-similarity;rrank=[]
        for i,row in enumerate(distance):rrank.append(1-((row<row[i]).sum()+.5*((row==row[i]).sum()-1))/(len(targets)-1))
        extra={'pds_full_official_target_gene_exclusion_raw':float(np.mean(rrank)),'diagnostic_note':'Same available competitors, all 300 official target genes excluded; not a full 300-competitor server score'}
    return {'status':'VERIFIED independent floating-point definitions agree (not exact biological verification)','independent':got,'official':{m:official[m] for m in got},'absolute_errors':err,'rho':rho,'numerator_sum':float(sum(nums)),'denominator_sum':float(sum(dens)),'prediction_correction_budget_sum':float(cap.sum()),'prediction_spread_budget':spread,'prediction_correction_used':float(min(spread,cap.sum())),'real_jackknife_correction_sum':float(sum(rr[t][1]/keep[i].sum() for i,t in enumerate(targets))),'squared_profile_error_sum':float(sum(np.sum((P[i,keep[i]]-R[i,keep[i]])**2)/keep[i].sum() for i,t in enumerate(targets))),**extra}
