"""Independent retained-output verifier; does not import search/analysis code."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,pandas as pd,math
from fractions import Fraction
z=np.load(rt.ROOT/'results/baseline.npz');p=z['realized'].astype(np.longdouble);y=z['truth'].astype(np.longdouble);own=z['own_idx'];targets=z['targets'];genes=z['genes'];n=len(targets)
for i,j in enumerate(own):
    if j>=0:p[i,j]=0;y[i,j]=0
table=pd.read_csv(rt.ROOT/'results/per_target.csv');assert table.target.tolist()==targets.tolist() and table.n_conf_exact.isna().all()
dot=np.sum(p*y,1);pn=np.sum(p*p,1);yn=np.sum(y*y,1);sse=np.sum((p-y)**2,1);alpha=dot/pn
direction=yn-dot*dot/pn;magnitude=pn*(1-alpha)**2
checks={}
for key,value in [('sse',sse),('truth_norm2',yn),('pred_norm2',pn),('best_scalar',alpha),('cosine',dot/np.sqrt(pn*yn)),('relative_error',sse/yn),('direction_error',direction/yn),('magnitude_error',magnitude/yn)]:
    err=float(np.max(np.abs(table[key].to_numpy()-value)));assert err<1e-9,(key,err);checks[key]=err
assert np.max(abs(sse-direction-magnitude))<64*np.finfo(np.longdouble).eps*np.max(sse+direction+magnitude)
exact=[]
for i in [0,55,167]:
    pp=list(map(lambda v:Fraction(float(v)),p[i]));yy=list(map(lambda v:Fraction(float(v)),y[i]))
    ep=sum((a-b)**2 for a,b in zip(pp,yy));en=sum(a*a for a in pp);ey=sum(b*b for b in yy);ed=sum(a*b for a,b in zip(pp,yy));a=ed/en
    assert (ey-ed*ed/en)+en*(1-a)**2==ep
    assert abs(float(ep)-float(sse[i]))<1e-10
    exact.append(dict(target=str(targets[i]),exact_decomposition_identity=True,exact_SSE_numerator=str(ep.numerator),exact_SSE_denominator=str(ep.denominator)))
analysis=rt.json.loads((rt.ROOT/'results/analysis.json').read_text());bulk=z['bulk'];P=np.array(p,dtype=float);Y=np.array(y,dtype=float)
assert analysis['baseline_integer_emitted']
renderer=rt.json.loads((rt.ROOT/'results/renderer_audit.json').read_text());assert renderer['status']=='PASS' and renderer['groups']==168 and renderer['reference_treated_rows_read']==0 and renderer['control_rows']==3200
for path,h in renderer['source_pins'].items():assert rt.sha(path)==h
raw=z['raw'];q=np.maximum(np.expm1(np.log1p(bulk*50000)+raw),0);q/=q.sum(1,keepdims=True);continuous=np.log1p(50000*q)-np.log1p(50000*bulk)
for i,j in enumerate(own):
    if j>=0:continuous[i,j]=0
continuous_sse=float(np.sum((continuous-Y)**2));emitted_sse=float(sse.sum())
def project(a):
    base=np.log1p(bulk*50000);counts=np.expm1(base+a);counts[counts<0]=0;counts/=counts.sum(1)[:,None]
    result=np.log1p(counts*50000)-base
    for i,j in enumerate(own):
        if j>=0:result[i,j]=0
    return result
projections={};loo_checks=[]
for method,group in [('global_scale',np.zeros(n)),('essential_scale',z['essential']),('source_KD_scale',np.digitize(z['source_kd'],[.5,.8]))]:
    coeff=[]
    for i in range(n):
        rows=[j for j in range(n) if j!=i and group[j]==group[i]]
        if len(rows)<8:rows=[j for j in range(n) if j!=i]
        a=max(float(sum(dot[j] for j in rows)/sum(pn[j] for j in rows)),0);coeff.append(a)
        # The omitted outcome can be arbitrary; no summand can include its index.
        assert i not in rows and len(rows)<=167
    pred=project(P*np.array(coeff)[:,None]);e=float(np.sum((pred-Y)**2));want=analysis['summary'][method]['sse'];assert abs(e-want)<1e-7,(method,e,want)
    assert np.max(abs(np.array(coeff)-analysis['coefficients'][method]))<1e-12
    projections[method]=pred;loo_checks.append(method)
pred_bias=project(np.stack([P[i]+np.mean(np.delete(Y-P,i,axis=0),axis=0) for i in range(n)]))
with np.load(rt.ROOT/'results/loo_bias.npz') as f:assert np.max(abs(pred_bias-f['response']))<2e-7
assert abs(float(np.sum((pred_bias-Y)**2))-analysis['summary']['LOO_bias_correction']['sse'])<1e-7
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2 import EvalConfig
cfg=EvalConfig.from_preset('vcc2026');dc=cfg.discrimination
labels=np.concatenate([targets,np.array(['non-targeting'])]);real=np.vstack([Y,np.zeros(len(genes))]);pdss={}
for method,pred in [('B0',P),*projections.items()]:
    value=discrimination_score(pred_bulk=(labels,np.vstack([pred,np.zeros(len(genes))])),real_bulk=(labels,real),genes=genes,distance='cosine',rank_denominator=dc.rank_denominator,tie_policy=dc.tie_policy,exclude_target_gene=dc.exclude_target_gene,exclusion_scope=dc.exclusion_scope,control_source='real')
    mean=float(np.mean(list(value.values())));assert abs(mean-analysis['summary'][method]['PDS168'])<1e-12,(method,mean)
    pdss[method]=mean
strata=pd.read_csv(rt.ROOT/'results/strata.csv')
for kind in ['coverage','essential_proxy','source_knockdown','observed_H1_knockdown','n_conf_MOMENT_PROXY','number_families']:
    sub=strata[strata.kind==kind];assert sub.n.sum()==168 and abs(sub.error_mass_pct.sum()-100)<1e-9
for kind in ['all','coverage','essential_proxy','source_knockdown','observed_H1_knockdown','n_conf_MOMENT_PROXY','number_families']:
    sub=strata[(strata.kind==kind)&(strata.n>0)];assert np.max(abs(sub.direction_pct_of_SSE+sub.magnitude_pct_of_SSE-100))<1e-9
meta=rt.json.loads((rt.ROOT/'results/training_metadata.json').read_text());assert sum(bool(v) for v in meta['coverage'])==161
inventory=rt.json.loads((rt.ROOT/'results/source_inventory.json').read_text())
for path,h in inventory['pins'].items():assert rt.sha(path)==h
for row in inventory['rows']:
    assert abs(row['training168_covered']/168-row['training168_fraction'])<1e-15
    assert set(row['new_targets'])<=set(inventory['uncovered_targets'])
for label in ['one','two_a','two_b']:
    v=rt.json.loads((rt.ROOT/f'results/benchmark_{label}.json').read_text());assert v['nice']>=10 and v['torch_threads']==1 and v['measured_steps']==100
    assert all(t['num_threads']==1 for t in v['blas_threads'])
baseline_budget=float(np.sum((P-P.mean(0))**2));budget={k:float(np.sum((a-a.mean(0))**2)/baseline_budget) for k,a in projections.items()}
rt.save('results/verification.json',dict(status='PASS',independent_per_target_max_abs_errors=checks,exact_rational_spots=exact,LOO_replays=loo_checks,official_PDS_function=pdss,PDS_space='Measured 18077-gene pseudobulk with all 168 panel targets excluded, native installed discrimination comparator. B0 is integer-emitted; proposed fixes use continuous projection. Not a full six-metric evaluation.',baseline_integer_emitted=True,continuous_SSE=continuous_sse,emitted_SSE=emitted_sse,emitter_relative_SSE_change=(emitted_sse-continuous_sse)/continuous_sse,baseline_centered_energy=baseline_budget,centered_energy_ratios=budget,limitations=['Exact full-training Wilcoxon n_conf unavailable, all n_conf_exact are missing; proxy strata explicitly named.','Essentiality annotation is essential-screen membership only.','New-source metadata checks certify overlap, not response quality.','Bootstrap CIs are exploratory and condition on reused H1 labels/fits.'],all_floats_finite=bool(np.isfinite(p).all() and np.isfinite(y).all()),free_gib=rt.guard()))
rt.status('Independent verifier PASS: long-double all-target replay, exact rational decomposition spots, LOO replay, native PDS parity, resource/metadata checks. Official MSE remains unmeasured.')
