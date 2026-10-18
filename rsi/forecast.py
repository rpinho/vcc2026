"""Freeze a prospective affine score forecast before v4 official labels appear.
This is calibration, not model search. With three training entries it is UNVERIFIED.
"""
import runtime as rt
from pathlib import Path
import argparse,time,numpy as np
from verify_results import official
p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--seed',type=int,default=42);args=p.parse_args();dest=Path('artifacts/forecasts')/(args.task+'.json')
if dest.exists():raise SystemExit('Forecast is immutable; use the saved receipt')
numbers=official();train=['v1','v2b','v3b'];test=['v4a','v4b']
assert not any(k in numbers for k in test),'Prospective forecast refused: v4 official label already available'
folder=Path('artifacts/scores')/args.task;anchors=rt.json.loads((folder/'anchors.json').read_text());assert isinstance(anchors['baseline_construction'],dict),'Primary baseline required'
records={k:rt.json.loads((folder/f'{k}_s{args.seed}.json').read_text()) for k in train+test};assert 'holdout' not in args.task.lower()
x=np.array([records[k]['overall'] for k in train]);y=np.array([numbers[k] for k in train]);design=np.column_stack([np.ones(len(x)),x]);beta=np.linalg.lstsq(design,y,rcond=None)[0];loo=[]
for i in range(len(x)):
 keep=np.arange(len(x))!=i;b=np.linalg.lstsq(design[keep],y[keep],rcond=None)[0];loo.append(dict(recipe=train[i],prediction=float(design[i]@b),actual=float(y[i])))
result=dict(status='UNVERIFIED prospective forecast, not a validated leaderboard predictor',frozen_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC',time.gmtime()),task=args.task,method='OLS official_overall = intercept + slope * primary local overall; no tuning',training={k:dict(local=records[k]['overall'],official=numbers[k]) for k in train},intercept=float(beta[0]),slope=float(beta[1]),forecast={k:dict(local=records[k]['overall'],official_forecast=float(beta[0]+beta[1]*records[k]['overall'])) for k in test},leave_one_out=loo,n=len(train),score_sha256={k:rt.sha(folder/f'{k}_s{args.seed}.json') for k in records},official_file_sha256=rt.sha(rt.ROOT.parent/'OFFICIAL_SCORES.md'),limitations='Three dependent recipe families; one public lineage; no reliable uncertainty interval; extreme v3b can dominate Pearson/OLS. Future v4 labels must be treated as prospective tests, not retroactive training.')
rt.save(dest,result);print(result)
