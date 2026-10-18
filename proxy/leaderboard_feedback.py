"""Official-panel empirical envelopes, conditional Pareto sets and candidate sanity checks."""
import os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'competitors'))
from monitor import compatible, entries
from calibrate_v2 import METRICS
C=list(METRICS)


def load(path):
    rows=[]
    for row in entries(json.loads(path.read_text())):
        if compatible(row) and all(isinstance(row.get('score_'+c),(int,float)) for c in C):
            rows.append(dict(id=row['id'],team_id=row.get('team_id'),team=row.get('team_name'),model=row.get('model_name'),
                             submission_date=row.get('submission_date',''),
                             overall=row['score_avg'],description=row.get('description'),**{c:row['score_'+c] for c in C}))
    frame=pd.DataFrame(rows).drop_duplicates('id')
    if not np.isfinite(frame[C+['overall']].to_numpy()).all():
        raise ValueError('nonfinite published components')
    residual=frame[C].mean(axis=1)-frame.overall
    if residual.abs().max()>1e-8:
        raise ValueError('score_avg is not the six scaled component mean')
    return frame


def frontier(frame,columns=('mse','fid','reach','pds')):
    x=frame[list(columns)].to_numpy();keep=[]
    for i,row in enumerate(x):
        dominated=np.all(x>=row,axis=1)&np.any(x>row,axis=1)
        if not dominated.any():keep.append(i)
    return frame.iloc[keep]


def sanity(frame,values):
    x=frame[C].to_numpy();q25,q75=np.quantile(x,[.25,.75],axis=0);scale=np.maximum(q75-q25,.01)
    distance=np.sqrt(np.mean(((x-values)/scale)**2,axis=1));ix=np.argsort(distance)[:8]
    return dict(candidate=dict(zip(C,map(float,values))),nearest=[dict(**frame.iloc[int(i)].drop('description').to_dict(),iqr_rms_distance=float(distance[i])) for i in ix],
                outside_observed_range=[c for j,c in enumerate(C) if values[j]<x[:,j].min() or values[j]>x[:,j].max()],
                note='Empirical support check only; a novel feasible trade-off can lie beyond the observed frontier. Not a predictor of official score.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--snapshot',type=Path)
    parser.add_argument('--candidate',type=Path);parser.add_argument('--out',type=Path,default=P/'artifacts/leaderboard_feedback.json')
    args=parser.parse_args();paths=[args.snapshot] if args.snapshot else sorted((P.parent/'competitors/raw').glob('lb-*.json'))
    frames=[load(path).assign(source_snapshot=path.name) for path in paths]
    f=pd.concat(frames,ignore_index=True).drop_duplicates('id',keep='last')
    latest=f.sort_values('submission_date').drop_duplicates('team_id',keep='last')
    snapshot=paths[-1]
    out=dict(source=snapshot.name,sources=[p.name for p in paths],raw_snapshot_rows=sum(len(v) for v in frames),
             n=len(f),teams=int(f.team_id.nunique()),descriptions=int(f.description.notna().sum()),
                            component_fields=['score_'+c for c in C],six_mean_max_error=float((f[C].mean(axis=1)-f.overall).abs().max()))
    out['quantiles']={c:{str(q):float(f[c].quantile(q)) for q in (0,.01,.1,.25,.5,.75,.9,.99,1)} for c in C}
    out['spearman']=f[C].corr(method='spearman').to_dict()
    out['latest_per_team_spearman']=latest[C].corr(method='spearman').to_dict()
    out['latest_per_team_quantiles']={c:{str(q):float(latest[c].quantile(q)) for q in (.1,.5,.9)} for c in C}
    front=frontier(f);out['pareto_ids']=front.id.tolist();out['pareto_n']=len(front)
    out['conditional']=[]
    for variable,edges in [('nmae',[-10,0,.05,.10,.15,.25,.4,.6,10]),('pds',[-10,.5,.6,.65,.7,.75,.85,10]),('overall',[-10,0,.1,.15,.2,.25,.3,10])]:
        for lo,hi in zip(edges[:-1],edges[1:]):
            subset=f[(f[variable]>=lo)&(f[variable]<hi)]
            if len(subset)==0:continue
            ff=frontier(subset);out['conditional'].append(dict(variable=variable,lo=lo,hi=hi,n=len(subset),pareto_n=len(ff),
                component_median={c:float(subset[c].median()) for c in C},component_q90={c:float(subset[c].quantile(.9)) for c in C},
                component_max={c:float(subset[c].max()) for c in C},pareto_ids=ff.id.tolist()))
    ref=f[(f.pds>=.65)&(f.pds<.7)&(f.nmae>=.05)&(f.nmae<.15)]
    out['near_v5b_regime']=dict(n=len(ref),frontier=frontier(ref).drop(columns='description').to_dict('records') if len(ref) else [])
    # Descriptions provide hypotheses, not treatments; record broad tagged associations transparently.
    tags={'dispersion':r'dispers|variance|noise|sampling|spread|moment|budget',
          'transfer':r'transfer|atlas|same.target|signature',
          'DE_gate':r'gate|differential|p_de|significan|wilcoxon|threshold'}
    out['description_associations']={}
    for tag,pattern in tags.items():
        mask=f.description.fillna('').str.contains(pattern,case=False,regex=True)
        out['description_associations'][tag]=dict(pattern=pattern,n=int(mask.sum()),median=f.loc[mask,C].median().to_dict())
    if args.candidate:
        record=json.loads(args.candidate.read_text());record=record.get('scaled',record)
        values=np.array([record[c] if c in record else record[METRICS[c]] for c in C])
        out['candidate_sanity']=sanity(f,values)
    args.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    front.drop(columns='description').to_csv(args.out.with_name('leaderboard_pareto.csv'),index=False)
    lines=['# Leaderboard as a feedback dataset','',f"Union of {len(paths)} snapshots through competitors/raw/{snapshot.name}: {len(f)} distinct comparable submissions, {out['teams']} unique teams, {out['descriptions']} descriptions.",
           'Filter: validation panel vcc2026-val-1 and exactly valA/B/C-r4. Use score_pds/score_mse/score_nmae/score_fid/score_reach/score_jac; raw fields have different scales and MSE direction. All six scaled fields are higher-is-better.',
           f"Six-component mean verified against score_avg (max discrepancy {out['six_mean_max_error']:.3g}).",'',
           f"The {out['raw_snapshot_rows']} snapshot rows are deduplicated by immutable submission ID. Latest payload/description wins for repeated IDs. Historical distinct submissions remain useful evidence of feasible trade-offs. JSON also reports correlations/quantiles on only the latest submission per team ({len(latest)} rows), since active teams otherwise carry more weight.",
           'These are cross-team observations, not randomized experiments. The API exposes a current entry per team, including latest entries that may be worse than earlier ones; repeated snapshots overlap. Do not count repeated IDs as new experiments or infer causal effects from descriptions.',
           f"Four-component Pareto frontier (mse/fid/reach/pds): {len(front)} entries. Full set in artifacts/leaderboard_pareto.csv. This is observed feasibility, not a biological upper bound.",'',
           '| Component | Median | 90th percentile | Maximum |','|---|---:|---:|---:|']
    for c in C:lines.append(f"| {c} | {f[c].median():.4f} | {f[c].quantile(.9):.4f} | {f[c].max():.4f} |")
    lines+=['','At PDS 0.65–0.70 and NMAE 0.05–0.15 (our current regime):',f"{len(ref)} entries. Joint feasible frontier:",'',
            '| Team/model | Total | PDS | MSE | NMAE | FID | Reach | Jaccard |','|---|---:|---:|---:|---:|---:|---:|---:|']
    if len(ref):
        for _,r in frontier(ref).sort_values('overall',ascending=False).iterrows():
            lines.append('| '+str(r.team).replace('|','/')+'/'+str(r.model).replace('|','/')+' | '+' | '.join(f'{r[c]:.4f}' for c in ['overall']+C)+' |')
    lines+=['','Conditional trade-offs (90th percentile of each component; maxima need not coexist in one entry):','',
            '| NMAE interval | n | MSE q90 | FID q90 | Reach q90 | PDS q90 | Pareto entries |','|---|---:|---:|---:|---:|---:|---:|']
    for row in out['conditional']:
        if row['variable']=='nmae':lines.append(f"| [{row['lo']}, {row['hi']}) | {row['n']} | "+' | '.join(f"{row['component_q90'][c]:.4f}" for c in ['mse','fid','reach','pds'])+f" | {row['pareto_n']} |")
    lines+=['','Use `python proxy/leaderboard_feedback.py --candidate candidate.json` to inspect a predicted scaled vector against nearest official vectors (IQR-normalized), ranges, and conditional frontiers. Distance is a sanity check; it cannot validate a proxy on a new cell context. JSON contains conditional sets at fixed NMAE, PDS and total intervals, description associations, and component rank correlations.']
    (P/'LEADERBOARD_FEEDBACK.md').write_text('\n'.join(lines)+'\n')
    print(f"{len(f)} entries; {len(front)} frontier points; {len(ref)} in current regime",flush=True)


if __name__=='__main__':main()
