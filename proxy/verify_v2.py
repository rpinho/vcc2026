"""Independent rational verification of saved scores and calibration, no model imports."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

P=Path(__file__).resolve().parent
ARMS=['v1','v2b','v3b','v4a','v4b','v5a','v5b']
NAMES={'pds':'pds_cosine','mse':'expr_mse_unbiased_capped_norm','nmae':'de_wilcoxon_lfc_nmae',
       'fid':'de_wilcoxon_direction_fidelity_yield_raw','reach':'de_wilcoxon_direction_reach_raw','jac':'de_wilcoxon_sig_jaccard'}


def close(a,b,tolerance=1e-9):
    assert abs(float(a)-float(b))<=tolerance*max(1,abs(float(a)),abs(float(b))),(a,b)


def mean(x):return sum(x)/len(x)


def rational_fit(x,y):
    xx=[F(v) for v in x];yy=[F(v) for v in y]
    mx=mean(xx);my=mean(yy);den=sum((v-mx)**2 for v in xx)
    raw=sum((u-mx)*(v-my) for u,v in zip(xx,yy))/den if den>F(1e-24) else F(0)
    slope=max(F(0),raw)
    return my-slope*mx,slope,raw


def rational_fit_predict(x,y,i):
    intercept,slope,_=rational_fit([v for j,v in enumerate(x) if j!=i],[v for j,v in enumerate(y) if j!=i])
    return intercept+slope*F(x[i])


def rank(values):
    groups={}
    for i,value in enumerate(values):groups.setdefault(round(float(value),12),[]).append(i)
    out=[F(0)]*len(values);start=1
    for value,indices in sorted(groups.items()):
        mid=F(2*start+len(indices)-1,2)
        for i in indices:out[i]=mid
        start+=len(indices)
    return out


def rho(x,y):
    xx=rank(x);yy=rank(y);mx=mean(xx);my=mean(yy)
    numerator=sum((u-mx)*(v-my) for u,v in zip(xx,yy));dx=sum((u-mx)**2 for u in xx);dy=sum((v-my)**2 for v in yy)
    if not dx or not dy:return None
    return float(numerator)/math.sqrt(float(dx*dy))


def official_table():
    result={}
    for line in (P.parent/'OFFICIAL_SCORES.md').read_text().splitlines():
        if not line.startswith('| v'):continue
        parts=[v.strip().replace('**','') for v in line.split('|')[1:-1]];a=parts[0].split()[0]
        if a not in ARMS:continue
        result[a]=dict(zip(['overall',*NAMES],map(float,parts[1:8])))
    assert set(result)==set(ARMS)
    return result


def verify_scores(folder,weights):
    evidence={}
    for dest in weights:
        anchor=json.loads((folder/dest/'anchors.json').read_text())
        assert anchor['status']=='VALID' and anchor['n_splits']==5
        import polars as pl
        splits=pl.read_parquet(folder/dest/'anchor_splits.parquet')
        prefix=anchor.get('cached_prefix')
        if prefix:
            source=P/prefix['source']
            assert hashlib.sha256(source.read_bytes()).hexdigest()==prefix['source_sha256']
            assert prefix['cached_seed_prefix']+prefix['new_seeds']==anchor['derived_seeds']
            assert prefix['n_cached']==3
        for long in NAMES.values():
            part=splits.filter(pl.col('metric')==long)
            assert sorted(part['seed'].unique().to_list())==sorted(anchor['derived_seeds'][3:] if prefix else anchor['derived_seeds'])
            values=[F(v) for v in part['value'].to_list()]
            if prefix:
                exact=(3*F(prefix['cached_mean'][long])+sum(values))/5
            else:exact=mean(values)
            close(exact,anchor['replicate'][long],1e-12)
        for key,long in NAMES.items():
            b=anchor['baseline'][long];a=anchor['replicate'][long]
            assert math.isfinite(a) and math.isfinite(b)
            assert (b-a>0 if key in ('mse','nmae') else a-b>0),(dest,key,a,b)
        for arm in ARMS:
            row=json.loads((folder/dest/(arm+'.json')).read_text())
            assert row['task_sha256']==anchor['real_sha256']
            values=[]
            for key,long in NAMES.items():
                # Finite affine anchored policy. MSE clamps [0,1]; NMAE floors at -6.
                val=(F(row['raw'][long])-F(anchor['baseline'][long]))/(F(anchor['replicate'][long])-F(anchor['baseline'][long]))
                if key=='mse':val=min(F(1),max(F(0),val))
                if key=='nmae':val=max(F(-6),val)
                close(val,row['scaled'][long],1e-12);values.append(val)
            close(mean(values),row['overall'],1e-12)
        evidence[dest]=dict(real_sha256=anchor['real_sha256'],anchor_sha256=hashlib.sha256((folder/dest/'anchors.json').read_bytes()).hexdigest())
    return evidence


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--calibration',type=Path,default=P/'artifacts/calibration_v2.json')
    parser.add_argument('--scores',type=Path,default=P/'artifacts/scores_v2');args=parser.parse_args()
    result=json.loads(args.calibration.read_text());assert result['status']!='INCOMPLETE'
    off=official_table();weights=result['weights'];evidence=verify_scores(args.scores,weights)
    preds={};trusted=[];all_sign=True;all_nonnegative=True
    for key,long in NAMES.items():
        reported=result['components'][key];x=[];y=[off[a][key] for a in ARMS]
        for a in ARMS:
            value=sum(F(w)*F(json.loads((args.scores/d/(a+'.json')).read_text())['scaled'][long]) for d,w in weights.items())
            close(value,reported['proxy'][a],1e-12);x.append(float(value))
        pred=[rational_fit_predict(x,y,i) for i in range(7)]
        for i,a in enumerate(ARMS):
            close(pred[i],reported['loo_prediction'][a]);close(abs(pred[i]-F(y[i])),reported['folds'][i]['absolute_error'])
            coefficients=rational_fit([v for j,v in enumerate(x) if j!=i],[v for j,v in enumerate(y) if j!=i])
            for name,value in zip(('intercept','slope','unconstrained_slope'),coefficients):close(value,reported['folds'][i][name])
        coefficients=rational_fit(x,y)
        for name,value in zip(('intercept','slope','unconstrained_slope'),coefficients):close(value,reported['full_fit'][name])
        measured=rho([float(v) for v in pred],y)
        if measured is None:assert reported['loo_spearman'] is None
        else:close(measured,reported['loo_spearman'],1e-12)
        mae=mean([abs(v-F(u)) for v,u in zip(pred,y)]);close(mae,reported['loo_mae'])
        sign_ok=True
        for a,b in [('v4a','v5a'),('v4a','v5b'),('v5a','v5b')]:
            i,j=ARMS.index(a),ARMS.index(b);dy=y[j]-y[i]
            if abs(dy)>=.001-1e-12:
                sign_ok=sign_ok and (x[j]-x[i])*dy>0 and float(pred[j]-pred[i])*dy>0
        assert bool(sign_ok)==reported['direction_pass'];all_sign=all_sign and sign_ok
        all_nonnegative=all_nonnegative and measured is not None and measured>=0
        if measured is not None and measured>=.8 and float(mae)<=.05 and sign_ok and max(x)-min(x)>=1e-14 and coefficients[2]>=0:trusted.append(key)
        preds[key]=pred
    total=[mean([preds[c][i] for c in NAMES]) for i in range(7)];y=[off[a]['overall'] for a in ARMS]
    for i,a in enumerate(ARMS):close(total[i],result['total']['loo_prediction'][a]);close(abs(total[i]-F(y[i])),result['total']['loo_abs_error'][a])
    total_rho=rho([float(v) for v in total],y);total_mae=mean([abs(v-F(u)) for v,u in zip(total,y)])
    if total_rho is None:assert result['total']['loo_spearman'] is None
    else:close(total_rho,result['total']['loo_spearman'],1e-12)
    close(total_mae,result['total']['loo_mae'])
    if result['provenance_errors']:trusted=[]
    assert trusted==result['trusted_components']
    passes=bool(total_rho is not None and total_rho>=.8 and float(total_mae)<=.03 and all_sign and all_nonnegative and not result['provenance_errors'])
    assert passes==result['full_proxy_passes']
    certificate=dict(status='PASS independent rational arithmetic and tied-rank verification',score_evidence=evidence,
                     calibration_sha256=hashlib.sha256(args.calibration.read_bytes()).hexdigest(),
                     total_rho=total_rho,total_mae=float(total_mae),trusted_components=trusted,full_proxy_passes=passes)
    output=args.calibration.with_name(args.calibration.stem+'_verification.json');output.write_text(json.dumps(certificate,indent=2))
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
