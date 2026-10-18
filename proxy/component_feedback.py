"""Apply the accepted seven-arm component mapping to a complete fixed-weight macro."""
import argparse
import hashlib
import json
from pathlib import Path
P=Path(__file__).resolve().parent
M={'pds':'pds_cosine','mse':'expr_mse_unbiased_capped_norm'}
FULL={**M,'nmae':'de_wilcoxon_lfc_nmae','fid':'de_wilcoxon_direction_fidelity_yield_raw',
      'reach':'de_wilcoxon_direction_reach_raw','jac':'de_wilcoxon_sig_jaccard'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',nargs='+',required=True,help='DESTINATION=fast_score.json for every fixed destination')
    p.add_argument('--candidate',required=True);p.add_argument('--out',type=Path)
    p.add_argument('--full',action='store_true',help='six-component total; requires full pre-registered proxy pass')
    a=p.parse_args();cal=json.loads((P/'artifacts/calibration_v2.json').read_text())
    verification=json.loads((P/'artifacts/calibration_v2_verification.json').read_text())
    assert verification['status']=='PASS independent rational arithmetic and tied-rank verification'
    assert verification['calibration_sha256']==hashlib.sha256((P/'artifacts/calibration_v2.json').read_bytes()).hexdigest(),'calibration changed after verification'
    assert cal['status'] in ('COMPONENT_ONLY','PASS'),'No proxy component has passed calibration'
    assert hashlib.sha256((P/'ACCEPTANCE_V2.md').read_bytes()).hexdigest()==cal['acceptance_sha256'],'acceptance rule changed after calibration'
    weights=cal['weights'];paths=dict(value.split('=',1) for value in a.scores)
    assert set(paths)==set(weights),'every fixed destination is required'
    records={d:json.loads(Path(path).read_text()) for d,path in paths.items()}
    for d,row in records.items():
        anchor=P/f'artifacts/scores_v2/{d}/anchors.json'
        assert row.get('destination',row.get('dest'))==d
        assert row['anchor_sha256']==hashlib.sha256(anchor.read_bytes()).hexdigest()
        assert row['task_sha256']==json.loads(anchor.read_text())['real_sha256']
    mapping=FULL
    if a.full:assert cal.get('full_proxy_passes'),'full proxy did not pass its pre-registered rule'
    available=[c for c,long in FULL.items() if all(long in row['scaled'] for row in records.values())]
    if a.full:assert set(available)==set(FULL),'full total needs all six components on every destination'
    eligible=list(FULL) if a.full else [c for c in cal['trusted_components'] if c in available]
    assert eligible,'Inputs contain no passing component; use the full scorer if passing DE components are missing'
    components={}
    for short in eligible:
        value=sum(w*records[d]['scaled'][mapping[short]] for d,w in weights.items())
        fit=cal['components'][short]['full_fit']
        components[short]=dict(proxy_macro=value,official_component_estimate=fit['intercept']+fit['slope']*value,
                               retrospective_loo_mae=cal['components'][short]['loo_mae'])
    result=dict(candidate=a.candidate,status='COMPONENT_ONLY',components=components,weights=weights,
                missing_passing_components=[c for c in cal['trusted_components'] if c not in available],
                seconds=sum(r['seconds'] for r in records.values()),
                calibration_sha256=hashlib.sha256((P/'artifacts/calibration_v2.json').read_bytes()).hexdigest(),
                limitation='No total estimate; seven-arm retrospective calibration still requires prospective official confirmation.')
    if a.full:
        result.update(status='FULL_PROXY',official_total_estimate=sum(v['official_component_estimate'] for v in components.values())/6,
                      limitation='Seven-arm retrospective calibration still requires prospective official confirmation.')
    if a.out:a.out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
