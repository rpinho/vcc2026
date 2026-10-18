"""Seven-arm, per-component calibration; frozen gates in ACCEPTANCE_V2.md."""
import os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import warnings
import numpy as np
from scipy import stats

P=Path(__file__).resolve().parent
ARMS=('v1','v2b','v3b','v4a','v4b','v5a','v5b')
METRICS=dict(pds='pds_cosine',mse='expr_mse_unbiased_capped_norm',nmae='de_wilcoxon_lfc_nmae',
             fid='de_wilcoxon_direction_fidelity_yield_raw',reach='de_wilcoxon_direction_reach_raw',jac='de_wilcoxon_sig_jaccard')


def official():
    result={}
    for line in (P.parent/'OFFICIAL_SCORES.md').read_text().splitlines():
        fields=[v.strip().replace('**','') for v in line.strip().strip('|').split('|')]
        if len(fields)<9:
            continue
        arm=fields[0].split()[0] if fields[0] else ''
        if arm in ARMS:
            result[arm]=dict(zip(('overall',*METRICS),map(float,fields[1:8])))
    if set(result)!=set(ARMS):
        raise ValueError('official table must have all seven arms')
    return result


def association(x,y,method='spearman'):
    x=np.round(np.asarray(x,float),12);y=np.round(np.asarray(y,float),12)
    if np.ptp(x)<1e-14 or np.ptp(y)<1e-14:
        return None
    function=stats.spearmanr if method=='spearman' else stats.kendalltau
    return float(function(x,y).statistic)


def fit(x,y):
    dx=x-x.mean();den=dx@dx
    slope=float(dx@(y-y.mean())/den) if den>1e-24 else 0.0
    constrained=max(0.0,slope)
    return float(y.mean()-constrained*x.mean()),constrained,slope


def directions(x,y,arms):
    rows=[]
    for a,b in itertools.combinations(('v4a','v5a','v5b'),2):
        if a not in arms or b not in arms:
            continue
        i,j=arms.index(a),arms.index(b);delta=float(y[j]-y[i])
        if abs(delta)<.001-1e-12:
            continue
        px=float(x[j]-x[i]);rows.append(dict(pair=f'{a}->{b}',official_delta=delta,proxy_delta=px,correct=bool(px*delta>0)))
    return rows


def evaluate(x,y,arms=ARMS):
    x=np.asarray(x,float);y=np.asarray(y,float);held=[];folds=[]
    for i,arm in enumerate(arms):
        mask=np.arange(len(arms))!=i
        intercept,slope,unconstrained=fit(x[mask],y[mask])
        prediction=intercept+slope*x[i];held.append(prediction)
        folds.append(dict(held_out=arm,intercept=intercept,slope=slope,unconstrained_slope=unconstrained,prediction=float(prediction),absolute_error=float(abs(prediction-y[i]))))
    intercept,slope,unconstrained=fit(x,y)
    raw_direction=directions(x,y,list(arms));held_direction=directions(np.array(held),y,list(arms))
    mae=float(np.mean(np.abs(np.array(held)-y)));rho=association(held,y)
    nonconstant=bool(np.ptp(x)>=1e-14)
    direction_pass=all(r['correct'] for r in raw_direction+held_direction)
    passed=bool(nonconstant and rho is not None and rho>=.8 and mae<=.05 and unconstrained>=0 and direction_pass)
    return dict(proxy=dict(zip(arms,x.tolist())),official=dict(zip(arms,y.tolist())),loo_prediction=dict(zip(arms,map(float,held))),
                raw_spearman=association(x,y),raw_kendall=association(x,y,'kendall'),loo_spearman=rho,
                loo_kendall=association(held,y,'kendall'),loo_mae=mae,
                full_fit=dict(intercept=intercept,slope=slope,unconstrained_slope=unconstrained),folds=folds,
                raw_directions=raw_direction,loo_directions=held_direction,direction_pass=direction_pass,passes_component_rule=passed)


def calibrate(folder,weights):
    off=official();rows={};missing=[];provenance_errors=[]
    for dest in weights:
        rows[dest]={}
        for arm in ARMS:
            path=folder/dest/f'{arm}.json'
            if not path.exists():
                missing.append(str(path));continue
            record=json.loads(path.read_text())
            if any(k not in record['scaled'] for k in METRICS.values()):
                raise ValueError(f'incomplete vector: {path}')
            rows[dest][arm]=record
    if missing:
        return dict(status='INCOMPLETE',missing=missing,weights=weights,arms=ARMS)
    for dest in weights:
        path=folder/dest/'anchors.json'
        anchor=json.loads(path.read_text()) if path.exists() else {}
        if anchor.get('status')!='VALID' or anchor.get('n_splits')!=5:
            provenance_errors.append(f'{dest}: not a valid five-split anchor')
        for arm in ARMS:
            if rows[dest][arm].get('task_sha256')!=anchor.get('real_sha256'):
                provenance_errors.append(f'{dest}/{arm}: task identity mismatch')
    registry=P/'artifacts/preregistration_v2.json'
    if registry.exists():
        registration=json.loads(registry.read_text())
        if registration['acceptance_sha256']!=hashlib.sha256((P/'ACCEPTANCE_V2.md').read_bytes()).hexdigest():
            provenance_errors.append('Acceptance document changed after registration')
        if registration['weights']!=weights:
            provenance_errors.append('Destination weights differ from pre-registration')
    fixture_path=P/'artifacts/full_adapter_fixture.json'
    fixture=json.loads(fixture_path.read_text()) if fixture_path.exists() else {}
    if fixture.get('status')!='PASS bitwise official full scorer and two-split anchors':
        provenance_errors.append('Combined official scorer/anchor fixture has not passed')
    else:
        from execution_v2 import numeric_adapter_sha256
        if fixture.get('numeric_adapter_sha256')!=numeric_adapter_sha256():
            provenance_errors.append('Numerical adapter functions changed after bitwise verification')
        if fixture.get('sorted_mwu_fixture',{}).get('status')!='official installed backend; sorted-reference optimization excluded':
            provenance_errors.append('Expected official installed MWU backend')
    result=dict(status='FAIL',weights=weights,arms=ARMS,components={},official_source='OFFICIAL_SCORES.md (reported precision)',
                acceptance_sha256=hashlib.sha256((P/'ACCEPTANCE_V2.md').read_bytes()).hexdigest(),provenance_errors=provenance_errors)
    result['per_destination']={}
    for dest in weights:
        result['per_destination'][dest]={short:evaluate([rows[dest][a]['scaled'][long] for a in ARMS],[off[a][short] for a in ARMS]) for short,long in METRICS.items()}
    for short,long in METRICS.items():
        x=[sum(weight*rows[d][a]['scaled'][long] for d,weight in weights.items()) for a in ARMS]
        result['components'][short]=evaluate(x,[off[a][short] for a in ARMS])
    prediction=np.mean([[r['loo_prediction'][a] for a in ARMS] for r in result['components'].values()],axis=0)
    y=np.array([off[a]['overall'] for a in ARMS]);raw=[sum(w*rows[d][a]['overall'] for d,w in weights.items()) for a in ARMS]
    total=dict(loo_prediction=dict(zip(ARMS,prediction.tolist())),official=dict(zip(ARMS,y.tolist())),loo_spearman=association(prediction,y),
               loo_kendall=association(prediction,y,'kendall'),loo_mae=float(np.mean(np.abs(prediction-y))),
               loo_abs_error=dict(zip(ARMS,np.abs(prediction-y).tolist())),raw_spearman=association(raw,y))
    result['total']=total
    result['overall_fit_diagnostic']=evaluate(raw,y)
    result['trusted_components']=[c for c,r in result['components'].items() if r['passes_component_rule']]
    passed=bool(not provenance_errors and total['loo_spearman'] is not None and total['loo_spearman']>=.8 and total['loo_mae']<=.03 and
                all(r['direction_pass'] and r['loo_spearman'] is not None and r['loo_spearman']>=0 for r in result['components'].values()))
    if provenance_errors:result['trusted_components']=[]
    result['status']='PASS' if passed else ('COMPONENT_ONLY' if result['trusted_components'] else 'FAIL')
    result['full_proxy_passes']=passed
    result['seconds_per_candidate']=float(sum(np.median([rows[d][a]['seconds'] for a in ARMS if not rows[d][a].get('runtime_is_prior_run')]) for d in weights))
    result['seconds_by_arm']={a:sum(rows[d][a]['seconds'] for d in weights) for a in ARMS}
    # Same fixed model on six arms, excluding v3b ONLY as labelled sensitivity.
    subset=tuple(a for a in ARMS if a!='v3b');sens={}
    for c,r in result['components'].items():
        sens[c]=evaluate([r['proxy'][a] for a in subset],[off[a][c] for a in subset],subset)
    sp=np.mean([[r['loo_prediction'][a] for a in subset] for r in sens.values()],axis=0)
    sy=[off[a]['overall'] for a in subset]
    result['sensitivity_excluding_v3b']=dict(components=sens,total_loo_spearman=association(sp,sy),total_loo_mae=float(np.mean(np.abs(sp-sy))))
    return result


def render(result):
    lines=['# Proxy V2 calibration','',f"Verdict: **{result['status']}**.",'',
           'Rule frozen in ACCEPTANCE_V2.md. All seven arms retained; no destination selected by score.']
    if result['status']=='INCOMPLETE':
        return '\n'.join(lines+['','Missing:']+['- '+x for x in result['missing']])+'\n'
    lines+=['',f"Destination weights: {result['weights']}. Official values retain their reported precision.",'',
            '| Component | Raw rho | LOO rho | LOO Kendall | LOO MAE | Direction gate | Trusted |',
            '|---|---:|---:|---:|---:|---|---|']
    def f(x):return 'undefined' if x is None else f'{x:.4f}'
    for c,r in result['components'].items():
        lines.append(f"| {c} | {f(r['raw_spearman'])} | {f(r['loo_spearman'])} | {f(r['loo_kendall'])} | {f(r['loo_mae'])} | {r['direction_pass']} | {r['passes_component_rule']} |")
    r=result['total'];lines += [f"| **total** | {f(r['raw_spearman'])} | **{f(r['loo_spearman'])}** | {f(r['loo_kendall'])} | **{f(r['loo_mae'])}** | all component gates required | {result['full_proxy_passes']} |",'',
        'The total is the equal mean of six held-out component predictions, not a fit to the held-out overall.',
        f"Runtime: {result['seconds_per_candidate']:.1f} seconds per candidate (sum of per-destination median timings); setup excluded.",'',
        '| Arm | Official total | Held-out total | Absolute error |','|---|---:|---:|---:|']
    for a in ARMS:lines.append(f"| {a} | {r['official'][a]:.4f} | {r['loo_prediction'][a]:.4f} | {r['loo_abs_error'][a]:.4f} |")
    lines+=['',f"Sensitivity excluding v3b (not the acceptance result): total LOO rho {f(result['sensitivity_excluding_v3b']['total_loo_spearman'])}, MAE {f(result['sensitivity_excluding_v3b']['total_loo_mae'])}.",'',
            'Full per-component per-arm errors, fold coefficients, direction pairs, source times, and per-line results are in the JSON alongside this report. A retrospective pass still needs prospective confirmation.']
    return '\n'.join(lines)+'\n'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--scores',type=Path,default=P/'artifacts/scores_v2')
    parser.add_argument('--weights',type=Path,default=P/'artifacts/control_audit.json')
    parser.add_argument('--out',type=Path,default=P/'artifacts/calibration_v2.json')
    args=parser.parse_args();weights=json.loads(args.weights.read_text())['destination_weights']
    assert abs(sum(weights.values())-1)<1e-12 and all(v>0 for v in weights.values())
    result=calibrate(args.scores,weights)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    args.out.with_suffix('.md').write_text(render(result))
    print(render(result))


if __name__=='__main__':main()
