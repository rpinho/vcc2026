"""Independent numerical verifier. Does not import models/selection/analyze."""
import runtime as rt
import numpy as np
from scipy.linalg import solve, cholesky, solve_triangular, svd
from decimal import Decimal, localcontext
import math

def verify_fit(z,y,model):
    G=z.T @ z; H=np.empty((z.shape[1],y.shape[1]),float)
    for lo in range(0,y.shape[1],192): H[:,lo:lo+192]=z.T @ y[:,lo:lo+192].astype(float)
    reconstructed=(model['U']*model['d'])@model['U'].T
    err=float(np.max(np.abs(G-reconstructed)))
    assert np.allclose(G,reconstructed,rtol=1e-9,atol=1e-10)
    assert np.allclose(H,model['rhs'],rtol=1e-12,atol=1e-12)
    # Independent Cholesky whitening, direct solve and SVD instead of eigen RRR.
    lam=model['penalty_scale']; A=G+lam*np.eye(len(G)); L=cholesky(A,lower=True)
    M=solve_triangular(L,H,lower=True)
    V,s,_=svd(M,full_matrices=False,check_finite=False)
    rank=min(80,len(G)); coef=solve_triangular(L.T,(V[:,:rank]@V[:,:rank].T)@M,lower=False)
    normal=solve(A,H,assume_a='pos')
    normal_error=float(np.linalg.norm(A @ normal-H)/max(np.linalg.norm(H),1e-30))
    d=model['d']; U=model['U']; root=(U/np.sqrt(d+lam))@U.T
    em=root @ H; ee,ev=np.linalg.eigh(em @ em.T); rv=ev[:,-rank:]
    other=root @ rv @ rv.T @ em
    rrr_error=float(np.linalg.norm(other-coef)/max(np.linalg.norm(coef),1e-30))
    assert normal_error<1e-9 and rrr_error<1e-7
    return dict(gram_max_error=err,normal_equation_relative=normal_error,independent_cholesky_SVD_RRR_relative=rrr_error)

def parts(a):
    a=np.asarray(a,np.float32); b=a.view(np.uint32)
    exponent=((b>>23)&255).astype(np.int16)
    mantissa=(b&0x7fffff).astype(np.int64)
    mantissa=np.where(exponent>0,mantissa+(1<<23),mantissa)
    mantissa=np.where((b>>31)>0,-mantissa,mantissa)
    shift=np.maximum(exponent-1,0)
    assert (exponent<255).all()
    # Exact value is signed mantissa * 2**shift * 2**-149.
    return mantissa,shift

def exact_product_sum(a,b):
    a=a.ravel(); b=b.ravel()
    accum=[0]*509
    # At most 2**15 products of two <2**24 integers: signed int64 cannot overflow.
    for lo in range(0,len(a),32768):
        hi=min(len(a),lo+32768); bins=np.zeros(509,np.int64)
        ma,ea=parts(a[lo:hi]); mb,eb=parts(b[lo:hi])
        np.add.at(bins,ea+eb,ma*mb)
        for j in np.flatnonzero(bins): accum[j]+=int(bins[j])
    return sum(x<<i for i,x in enumerate(accum))

def exact_column_sums(a):
    g=a.shape[1]; out=[0]*g; bins=np.zeros((255,g),np.int64)
    for lo in range(0,len(a),64):
        m,e=parts(a[lo:lo+64]); np.add.at(bins,(e,np.arange(g)[None,:]),m)
    for power in np.flatnonzero(np.any(bins!=0,axis=1)):
        ss=bins[power]
        for j in np.flatnonzero(ss): out[j]+=int(ss[j])<<int(power)
    return out

def cert(a,y,truth_exact=None):
    n=len(a); cross=exact_product_sum(a,y); energy=exact_product_sum(a,a)
    truth=truth_exact if truth_exact is not None else exact_product_sum(y,y)
    sums=exact_column_sums(a); sumsq=sum(x*x for x in sums)
    B_numerator=n*energy-sumsq
    assert B_numerator>0
    # Physical dot/energy have common denominator 2**298; B adds denominator n.
    return dict(cross_numerator=str(cross),energy_numerator=str(energy),truth_numerator=str(truth),
                B_numerator=str(B_numerator),dyadic_exponent=-298,B_extra_denominator=n)

def pds_native(q,y,targets,genes):
    keep=~np.isin(genes,targets)
    yt=y[:,keep].astype(float); yp=q[:,keep].astype(float)
    yt/=np.maximum(np.sqrt(np.einsum('ij,ij->i',yt,yt))[:,None],1e-30)
    n=len(q); ranks=np.empty(n); sample=[]
    # Streamed, exact float64 native-comparator cosine ranking.
    for lo in range(0,n,128):
        hi=min(n,lo+128); ss=yp[lo:hi]@yt.T
        diag=ss[np.arange(hi-lo),np.arange(lo,hi)]
        rank=np.count_nonzero(ss>diag[:,None],axis=1)+.5*(np.count_nonzero(ss==diag[:,None],axis=1)-1)
        ranks[lo:hi]=1-rank/(n-1)
        if lo==0:
            # Independent scalar math.fsum check of ordering for first row.
            explicit=np.array([math.fsum(float(a)*float(b) for a,b in zip(yp[0],row)) for row in yt])
            v=explicit[0]; scalar=1-(sum(explicit>v)+.5*(sum(explicit==v)-1))/(n-1)
            assert abs(scalar-ranks[0])<1e-12
            sample.append(dict(row=0,PDS=float(scalar)))
    return ranks,sample

def evaluate(q,baseline,y,targets,genes,seed=20261026,baseline_cache=None):
    """All inputs frozen float32; exact certificate independently checks ordinary metrics."""
    rt.guard(True)
    q=np.asarray(q,np.float32); baseline=np.asarray(baseline,np.float32); y=np.asarray(y,np.float32)
    n,g=q.shape
    assert q.shape==baseline.shape==y.shape and np.isfinite(q).all()
    idx={s:i for i,s in enumerate(genes)}
    for i,t in enumerate(targets):
        if t in idx: assert q[i,idx[t]]==baseline[i,idx[t]]==y[i,idx[t]]==0
    yy=np.empty(n)
    for lo in range(0,n,128):
        block=y[lo:lo+128].astype(float)
        yy[lo:lo+128]=np.einsum('ij,ij->i',block,block)
    del block
    T=float(yy.sum())
    rng=np.random.default_rng(seed); draws=rng.integers(n,size=(2000,n))
    counts=np.zeros((2000,n),float)
    for i,inds in enumerate(draws): counts[i]=np.bincount(inds,minlength=n)
    del draws
    truthboot=counts @ yy
    def raw(a,known_cert=None):
        dot=np.empty(n); energy=np.empty(n); sa=np.zeros(g); sy=np.zeros(g)
        for lo in range(0,n,128):
            af=a[lo:lo+128].astype(float); yf=y[lo:lo+128].astype(float)
            dot[lo:lo+128]=np.einsum('ij,ij->i',af,yf)
            energy[lo:lo+128]=np.einsum('ij,ij->i',af,af)
            sa+=af.sum(0); sy+=yf.sum(0)
        B=float(energy.sum()-sa @ sa/n)
        sumsquares=np.zeros(2000)
        for lo in range(0,g,256):
            summed=counts @ a[:,lo:lo+256].astype(float)
            sumsquares+=np.einsum('ij,ij->i',summed,summed)
        bb=counts @ energy-sumsquares/n
        del summed,af,yf
        pds,replay=pds_native(a,y,targets,genes)
        certificate=known_cert if known_cert else cert(a,y)
        cross_exact=int(certificate['cross_numerator'])/2.**298
        energy_exact=int(certificate['energy_numerator'])/2.**298
        B_exact=int(certificate['B_numerator'])/(n*2.**298)
        assert np.isclose(dot.sum(),cross_exact,rtol=1e-12,atol=1e-10)
        assert np.isclose(energy.sum(),energy_exact,rtol=1e-12,atol=1e-10)
        assert np.isclose(B,B_exact,rtol=1e-12,atol=1e-10)
        return dict(dot=dot,energy=energy,B=B,Bboot=bb,pds=pds,certificate=certificate,
                    center_cross=float(dot.sum()-sa @ sy/n),scalar_PDS_replay=replay)
    base=baseline_cache if baseline_cache is not None else raw(baseline)
    ss=base if np.array_equal(q,baseline) else raw(q)
    alpha=float(np.sqrt(base['B']/ss['B'])); C=alpha*float(ss['dot'].sum())/T
    Q=alpha**2*float(ss['energy'].sum())/T
    cb=np.sqrt(base['Bboot']/np.maximum(ss['Bboot'],1e-30))*(counts @ ss['dot'])/truthboot
    basecb=(counts @ base['dot'])/truthboot
    dpboot=(counts@(ss['pds']-base['pds']))/n
    ci=lambda v:np.quantile(v,[.025,.975]).tolist()
    item=dict(n=n,genes=g,cross=C,aligned_per_spread=C/np.sqrt(Q),prediction_energy=Q,B=base['B']/T,
        raw_B=ss['B']/T,match_scale=alpha,PDS=float(ss['pds'].mean()),
        delta_cross=C-float(base['dot'].sum())/T,delta_PDS=float(np.mean(ss['pds']-base['pds'])),
        delta_cross_CI95=ci(cb-basecb),delta_PDS_CI95=ci(dpboot),cross_CI95=ci(cb),
        centered_cross=alpha*ss['center_cross']/T)
    with localcontext() as ctx:
        ctx.prec=70
        qc=ss['certificate']; bc=base['certificate']
        scale=(Decimal(bc['B_numerator'])/Decimal(qc['B_numerator'])).sqrt()
        exactC=scale*Decimal(qc['cross_numerator'])/Decimal(qc['truth_numerator'])
        assert abs(float(exactC)-C)<1e-12
        item['exact_matched_cross_70digit']=str(exactC)
    item['qualifies']=bool(item['delta_cross_CI95'][0]>0 and item['delta_PDS']>=-.005)
    item['qualifies_PDS_CI']=bool(item['qualifies'] and item['delta_PDS_CI95'][0]>=-.005)
    frozen=dict(dot=ss['dot'],energy=ss['energy'],pds=ss['pds'],Bboot=ss['Bboot'],
                crossboot=cb,delta_PDS_boot=dpboot,truth_energy=yy)
    return item,ss['certificate'],frozen,base

if __name__=='__main__':
    rt.enforce(); rt.guard()
    # Independent exact-arithmetic self-test against Fraction on nontrivial floats.
    from fractions import Fraction
    a=np.array([[.1,-1.25,0.],[1e-20,3.75,-.02],[2.,-4.,.125]],np.float32)
    y=np.array([[3.,.13,-1.],[2e-15,-.005,8.],[.01,0.,9.]],np.float32)
    c=cert(a,y)
    assert Fraction(int(c['cross_numerator']),2**298)==sum(Fraction(float(x))*Fraction(float(z)) for x,z in zip(a.flat,y.flat))
    assert Fraction(int(c['energy_numerator']),2**298)==sum(Fraction(float(x))**2 for x in a.flat)
    exactB=sum(Fraction(float(x))**2 for x in a.flat)-sum(sum(Fraction(float(a[i,j])) for i in range(3))**2 for j in range(3))/3
    assert Fraction(int(c['B_numerator']),3*2**298)==exactB
    # Nontrivial rank-80 case (not merely unrestricted regression).
    from scipy.linalg import eigh
    rng=np.random.default_rng(8127); z=rng.normal(size=(256,128)); yytest=rng.normal(size=(256,192))
    G=z.T @ z; H=z.T @ yytest; d,U=eigh(G)
    fitcheck=verify_fit(z,yytest,dict(U=U,d=d,rhs=H,penalty_scale=float(np.trace(G)/len(G))))
    rt.save('results/verifier_selftest.json',dict(PASS=True,test='Exact float32 dyadic dot, norm and centered B versus independent Fraction',
        nontrivial_rank80=fitcheck))
    rt.status('PASS exact dyadic verifier self-test.')
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--replay',action='store_true'); args=parser.parse_args()
    if args.replay:
        receipts=[]
        for path in sorted((rt.ROOT/'results').glob('*_metrics.json')):
            outer=path.name.removesuffix('_metrics.json')
            metrics=rt.json.loads(path.read_text())
            certs=rt.json.loads((rt.ROOT/f'results/{outer}_exact.json').read_text())
            with np.load(rt.ROOT/f'results/{outer}_stats.npz') as z:
                for name,item in metrics.items():
                    base='B0_original' if name.endswith('_original') else 'plain'
                    c=certs[name]; bc=certs[base]; n=item['n']
                    truth=np.asarray(z[name+'__truth_energy']); dot=np.asarray(z[name+'__dot'])
                    energy=np.asarray(z[name+'__energy']); pds=np.asarray(z[name+'__pds'])
                    assert np.isclose(dot.sum(),int(c['cross_numerator'])/2.**298,rtol=1e-12)
                    assert np.isclose(energy.sum(),int(c['energy_numerator'])/2.**298,rtol=1e-12)
                    with localcontext() as ctx:
                        ctx.prec=70
                        value=(Decimal(bc['B_numerator'])/Decimal(c['B_numerator'])).sqrt()*Decimal(c['cross_numerator'])/Decimal(c['truth_numerator'])
                    assert abs(float(value)-item['cross'])<1e-12
                    assert abs(float(pds.mean())-item['PDS'])<1e-12
                    draws=np.random.default_rng(20261026).integers(n,size=(2000,n))
                    # Independent bootstrap replay by indexing rows, avoiding the search's counts @ vector path.
                    Tb=truth[draws].sum(1); dc=dot[draws].sum(1)
                    bdot=np.asarray(z[base+'__dot']); bpds=np.asarray(z[base+'__pds'])
                    Bq=np.asarray(z[name+'__Bboot']); Bb=np.asarray(z[base+'__Bboot'])
                    cb=np.sqrt(Bb/Bq)*dc/Tb
                    diff=cb-bdot[draws].sum(1)/Tb
                    pdiff=(pds[draws]-bpds[draws]).mean(1)
                    assert np.allclose(np.quantile(diff,[.025,.975]),item['delta_cross_CI95'],rtol=1e-9,atol=1e-12)
                    assert np.allclose(np.quantile(pdiff,[.025,.975]),item['delta_PDS_CI95'],rtol=1e-9,atol=1e-12)
                    receipts.append(dict(recipient=outer,arm=name,exact_totals=True,bootstrap_replayed=2000))
        rt.save('results/replay_verification.json',dict(PASS=True,arms=receipts,certificate_domain='frozen float32 predictions',
            coverage='Exact dot/energy/B certified at evaluation; 70-digit matched cross and all 2000 saved-stat bootstrap replicates independently replayed'))
        rt.status(f'PASS standalone replay of {len(receipts)} final arms.')
