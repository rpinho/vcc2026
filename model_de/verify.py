"""Independent exact count invariants + rational metric fixture + official-driver parity.
Does not share selection/search logic. Rational checks certify fixture mathematics;
real-data results remain floating-point official-code measurements, not exact claims.
"""
import os
os.environ['VCC_DE_THREADS']='2'
import runtime as rt
import sys
sys.path.extend([str(rt.VCC/'proxy/scorer')])
import numpy as np,pandas as pd,polars as pl,anndata as ad
from scipy import sparse
from fractions import Fraction
from dataclasses import replace
from strategy import redistribute,selective_template,select,rank_sources
from cell_eval2.de import prepare_de
from cell_eval2.metrics.direction import de_direction_fidelity_yield_raw,de_direction_reach
from cell_eval2.metrics.de import de_sig_jaccard,de_lfc_nmae
from cell_eval2 import compute_metrics
import official_scoring as scoring,bounded_score

def main():
    rng=np.random.default_rng(27);x=rng.poisson(5,(400,37)).astype(np.int32)
    mask=np.arange(37)%3==0;y=redistribute(x,mask)
    # Python integer sums are independent of vectorized implementation arithmetic.
    assert all(sum(map(int,x[:,j]))==sum(map(int,y[:,j])) for j in range(x.shape[1]))
    assert (y>=0).all() and np.issubdtype(y.dtype,np.integer)
    assert np.array_equal(y[:,~mask],x[:,~mask])
    before=x.sum(0)/x.sum();after=y.sum(0)/y.sum()
    assert np.array_equal(before,after)
    from exchange import exchange_equal_depth
    from exchange_verify import certify
    small=rng.integers(0,8,(12,8),dtype=np.int32);ds=np.repeat([60,80,100],4)
    small[:,-1]=ds-small[:,:-1].sum(1)
    call=np.arange(8)<3;balanced=exchange_equal_depth(small,call,preferred_buffer=7)
    exchange_certificate=certify(small,balanced)
    assert exchange_certificate['changed_entries']>0
    for g in range(small.shape[1]):
        exact_before=sum(Fraction(int(small[i,g]),int(ds[i])) for i in range(len(ds)))
        exact_after=sum(Fraction(int(balanced[i,g]),int(ds[i])) for i in range(len(ds)))
        assert exact_before==exact_after
    genes=['T0']+[f'g{i:02}' for i in range(20)]
    def table(pred):
        return pl.DataFrame(dict(target=['T0']*21,feature=genes,
            log2_fold_change=[1.]+[(-1. if pred and i==0 else 1.) for i in range(20)],
            p_adj=[.001]+[(.01 if pred or i<10 else 1.) for i in range(20)],
            p_value=[.001]+[.001+i*.0001 for i in range(20)]))
    prep=prepare_de(table(True),table(False),control='non-targeting')
    f=de_direction_fidelity_yield_raw(prep)['T0']
    r=de_direction_reach(prep,corrected=False)['T0']
    j=de_sig_jaccard(prep)['T0'];n=de_lfc_nmae(prep)['T0']
    expected=dict(fid=Fraction(19,20),reach=Fraction(1),jac=Fraction(1,2),nmae=Fraction(1,5))
    actual=dict(fid=f,reach=r,jac=j,nmae=n)
    for k in expected:assert abs(actual[k]-float(expected[k]))<1e-14,(k,actual,expected)
    # Compare bounded dispatcher against a fresh UNMODIFIED official compute_metrics call.
    gs=['T0','T1']+[f'g{k}' for k in range(24)];parts=[];labs=[]
    for t in ['non-targeting','T0','T1']:
        xx=rng.poisson(12,(40,len(gs))).astype(np.int32)
        if t!='non-targeting':xx[:,2:18]+= 7 if t=='T0' else 3
        parts.append(sparse.csr_matrix(xx));labs.extend([t]*40)
    real=ad.AnnData(sparse.vstack(parts),obs=pd.DataFrame({'target_gene':labs},index=[f'r{i}' for i in range(120)]),var=pd.DataFrame(index=gs))
    pred=real.copy()
    cfg=replace(scoring.config(),num_threads=1,gather_threads=1)
    a=compute_metrics(pred,real,config=cfg).sort(['metric','perturbation'])
    exchanged=pred.copy();xx=pred.X.toarray();lab=pred.obs.target_gene.to_numpy()
    for target in ['T0','T1']:
        rows=np.flatnonzero(lab==target);call=np.zeros(len(gs),bool);call[2:8]=True
        yy=exchange_equal_depth(xx[rows],call,preferred_buffer=gs.index(target));certify(xx[rows],yy);xx[rows]=yy
    exchanged.X=sparse.csr_matrix(xx)
    ev=compute_metrics(exchanged,real,config=cfg).sort(['metric','perturbation'])
    for name in ['pds_cosine','de_wilcoxon_lfc_nmae']:
        av=a.filter(pl.col('metric')==name)['value'].to_numpy();bv=ev.filter(pl.col('metric')==name)['value'].to_numpy()
        np.testing.assert_allclose(av,bv,rtol=0,atol=0 if name=='pds_cosine' else 2e-13,equal_nan=True)
    two=compute_metrics(pred,real,config=replace(cfg,num_threads=2)).sort(['metric','perturbation'])
    np.testing.assert_allclose(a['value'].to_numpy(),two['value'].to_numpy(),rtol=0,atol=0,equal_nan=True)
    b=bounded_score.resident(pred,real,config=cfg).sort(['metric','perturbation'])
    assert a.select('metric','perturbation').equals(b.select('metric','perturbation'))
    np.testing.assert_allclose(a['value'].to_numpy(),b['value'].to_numpy(),rtol=0,atol=0,equal_nan=True)
    from audit import prediction_de
    audit=prediction_de(parts[1],real[real.obs.target_gene.to_numpy()=='non-targeting'].copy(),'T0')
    assert set(audit['target'].to_list())=={'T0'} and 'T0' not in audit['feature'].to_list()
    assert audit.height==len(gs)-1
    singleton=ad.AnnData(sparse.csr_matrix([[5,5]]),var=pd.DataFrame(index=['T0','g']))
    one=prediction_de(sparse.csr_matrix([[9,1]]),singleton,'T0')
    assert one.height==1 and one['feature'].to_list()==['g']
    rt.save(rt.ROOT/'artifacts/verification.json',dict(status='PASS',rational_expected={k:str(v) for k,v in expected.items()},actual=actual,
        count_invariants='exact Python integer column sums; untouched other columns; nonnegative integer values; unchanged bulk proportions',
        scorer_parity='bitwise equality, all six metrics and diagnostics; synthetic integer cells; CPU DE threads1/2 identical',
        tested_strategy_sha256=rt.sha(rt.ROOT/'strategy.py'),tested_harness_sha256=rt.sha(rt.ROOT/'harness.py'),
        tested_audit_sha256=rt.sha(rt.ROOT/'audit.py'),prediction_audit='complete gated table, own target and control rows excluded; singleton groups accepted',
        exchange=dict(**exchange_certificate,official_pds='bitwise unchanged',official_nmae='unchanged within2e-13; rational means exactly equal'),exchange_sha256=rt.sha(rt.ROOT/'exchange.py'),
        exchange_verifier_sha256=rt.sha(rt.ROOT/'exchange_verify.py'),
        fast_de=dict(status='EXCLUDED: upstream license unspecified')))
    print('PASS',actual)
if __name__=='__main__':main()
