"""Separate check through installed native kernels and independent aggregation."""
import runtime_local as j
j.guard(launch=True)
import numpy as np, csv, math
from fractions import Fraction
from cell_eval2.moments import GroupMoments
from cell_eval2.metrics import delta
from cell_eval2 import aggregate_metrics
import polars as pl

def load(path,ts):
    with np.load(path) as z:
        index={str(t):i for i,t in enumerate(z['targets'])}
        return z['means'][[index[str(t)] for t in ts]]

if __name__=='__main__':
    result=j.json.loads((j.ROOT/'results.json').read_text())
    z=np.load(j.ROOT/'reference.npz');ts=z['targets'][:-1];genes=z['genes'];real=z['means'];rts=z['targets']
    rjk=z['jk'];rmo=GroupMoments(perts=rts,counts=z['counts'],sumsq=np.zeros(len(rts)),jk=rjk)
    base=load(j.ROOT/'cached_v5b.npz',ts);amp=load(j.ROOT/'bulk_v8b_like.npz',ts);bill=load(j.ROOT/'bulk_levers_knob_none.npz',ts)
    source=dict(v5b=base,v9a_885=base,v9a_694=base,v8b_like=amp,levers_knob_none=bill,incumbent=bill)
    checks={};pjk=np.load(j.ROOT/'prediction_jk.npz')
    w=np.full((len(ts),len(genes)),1/len(genes))
    for i,t in enumerate(ts):
        where=np.flatnonzero(genes==t)
        if len(where):w[i,:]=1/(len(genes)-1);w[i,where[0]]=0
    for arm,pred in source.items():
        j.guard()
        pm=GroupMoments(perts=ts,counts=np.full(len(ts),400),sumsq=np.zeros(len(ts)),jk=pjk[arm])
        kw=dict(pred_bulk=(ts,pred),real_bulk=(rts,real),pred_moments=pm,real_moments=rmo,
            comparator='bulk_lognorm',genes=genes,driver='mse_instr independent native verification',control='non-targeting')
        native_u=delta.mse_unbiased(**kw);native_c=delta.mse_unbiased_capped(**kw);native_d=delta.distance_unbiased(**kw)
        original={}
        path=j.VCC/'final_push/artifacts/de_check/H1_rich_search'/arm/'tidy.csv'
        for r in csv.DictReader(path.open()):original[(r['perturbation'],r['metric'])]=float(r['value'])
        errors={}
        rows=[]
        for metric,values in [('expr_mse_unbiased',native_u),('expr_mse_unbiased_capped',native_c),('expr_distance_unbiased',native_d)]:
            errors[metric]=max(abs(v-original[(t,metric)]) for t,v in values.items())
            rows.extend([(t,metric,v) for t,v in values.items()])
        assert max(errors.values())<2e-12,(arm,errors)
        frame=pl.DataFrame(rows,schema=['perturbation','metric','value'],orient='row')
        agg=aggregate_metrics(frame);agg=dict(zip(agg['metric'].to_list(),agg['mean'].to_list()))
        ratio=agg['expr_mse_unbiased_capped_norm']
        assert abs(ratio-result['packages'][arm]['raw']['capped_clamp_lifted'])<2e-12
        # Separate per-coordinate weighted-mean implementation of B. Avoid the
        # scorer's shifted two-pass algorithm entirely.
        center=np.sum(w*pred,axis=0)/np.sum(w,axis=0)
        B=float(np.sum(w*(pred-center)**2))
        assert abs(B-result['packages'][arm]['B'])<1e-12
        # Exact rational sum/division of the returned binary float rows.
        num=sum(map(Fraction.from_float,native_c.values()),Fraction())
        den=sum(map(Fraction.from_float,native_d.values()),Fraction())
        exact_ratio=float(num/den)
        assert abs(exact_ratio-ratio)<2e-12
        checks[arm]=dict(native_row_errors=errors,aggregate_raw=ratio,independent_B_error=abs(B-result['packages'][arm]['B']),rational_ratio_error=abs(exact_ratio-ratio))
    packages=result['packages'];a=packages['v8b_like'];b=packages['v9a_694']
    decomposition=dict(plugin_error_increase=a['E']-b['E'],prediction_credit_increase=min(a['B'],a['claim_labels'])-min(b['B'],b['claim_labels']),
        raw_normalized_increase=a['raw']['budget_labels']-b['raw']['budget_labels'])
    # Genuine C6 is absent regardless of whether a surrogate happens to match.
    assert all(not g['qualified'] for g in result['gate'].values())
    assert all(not cs['C4']['sign_ok'] for cs in result['contrasts'].values())
    from importlib.metadata import version
    package=j.Path(delta.__file__).parent.parent
    j.save('verification.json',dict(status='PASS',checks=checks,C4_decomposition=decomposition,
        installed_version=version('cell-eval2'),scorer_files={str(p.relative_to(package)):j.sha(p) for p in package.rglob('*.py')},
        reference_sha256=j.sha(j.ROOT/'reference.npz'),bulk_reconstruction_sha256=j.sha(j.ROOT/'bulk_reconstruction.json'),
        results_sha256=j.sha(j.ROOT/'results.json'),free_gib=j.guard(),nice=j.os.getpriority(j.os.PRIO_PROCESS,0)))
    j.status('Independent installed-kernel, aggregate, rational ratio and per-coordinate budget checks PASS for all six analog packages.')
    print(j.json.dumps(dict(checks=checks,C4_decomposition=decomposition),indent=2))
