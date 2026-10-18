"""Read-only leaderboard trade-off analysis, filtered to the same scoring stamp.
The leaderboard contains metric aggregates, not cells or DE tables: official DE
cannot be rerun on it. We analyze its measured official outputs instead.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json, math, pathlib, argparse, statistics, hashlib
P=pathlib.Path(__file__).resolve().parent
KEYS=['pds','mse','nmae','fid','reach','jac']
def corr(x,y):
    def ranks(a):
        out=[0.]*len(a);order=sorted(range(len(a)),key=a.__getitem__);i=0
        while i<len(order):
            j=i+1
            while j<len(order) and a[order[j]]==a[order[i]]:j+=1
            for k in order[i:j]:out[k]=(i+j-1)/2
            i=j
        return out
    x=ranks(x);y=ranks(y);mx=statistics.mean(x);my=statistics.mean(y)
    den=math.sqrt(sum((v-mx)**2 for v in x)*sum((v-my)**2 for v in y))
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else None
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot',default=str(P.parent/'competitors/raw/lb-2026-10-18.json'));args=ap.parse_args()
    path=pathlib.Path(args.snapshot);data=json.loads(path.read_text())
    rows=data if isinstance(data,list) else data['entries']
    rows=[r for r in rows if all(isinstance(r.get('score_'+k),(int,float)) and math.isfinite(r['score_'+k]) for k in KEYS)]
    stamp=max(set((r.get('panel_id'),r.get('anchor_version')) for r in rows),key=lambda s:sum((r.get('panel_id'),r.get('anchor_version'))==s for r in rows))
    rows=[r for r in rows if (r.get('panel_id'),r.get('anchor_version'))==stamp]
    residual=max(abs(sum(r['score_'+k] for k in KEYS)/6-r['score_avg']) for r in rows)
    assert residual<1e-12,('unexpected aggregate convention',residual)
    rows.sort(key=lambda r:-r['score_avg'])
    def compact(r):return dict(team=r.get('team_name'),model=r.get('model_name'),score=r['score_avg'],id=r['id'],**{k:r['score_'+k] for k in KEYS})
    cohorts={'all':rows,'top50':rows[:50],'top11_50':rows[10:50],
        'pds_065_072':[r for r in rows if .65<=r['score_pds']<=.72],
        'pds_ge065_mse_ge030':[r for r in rows if r['score_pds']>=.65 and r['score_mse']>=.30]}
    out=dict(snapshot=str(path.relative_to(P.parent)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),stamp=stamp,n=len(rows),leader=compact(rows[0]),cohorts={},score_identity_max_error=residual)
    for name,rs in cohorts.items():
        out['cohorts'][name]=dict(n=len(rs),median={k:statistics.median(r['score_'+k] for r in rs) for k in KEYS},
            spearman={a+'__'+b:corr([r['score_'+a] for r in rs],[r['score_'+b] for r in rs]) for a,b in [('mse','fid'),('mse','reach'),('fid','reach'),('fid','jac'),('reach','jac')]},
            best= [compact(r) for r in sorted(rs,key=lambda r:-sum(r['score_'+k] for k in ['mse','fid','reach','jac']))[:5]])
    axes=['mse','fid','reach','jac']
    # Frontiers are measured points, not independently attainable coordinate maxima.
    near=cohorts['pds_065_072'];out['near_pds_frontier']=[]
    for r in near:
        if not any(all(s['score_'+k]>=r['score_'+k] for k in axes) and any(s['score_'+k]>r['score_'+k] for k in axes) for s in near):
            out['near_pds_frontier'].append(compact(r))
    out['limitations']='Observational best-entry snapshots; team selection/confounding. No causal lever estimates. No cells supplied, so DE cannot be recomputed.'
    target=P/'artifacts/leaderboard.json';target.parent.mkdir(exist_ok=True,parents=True);target.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
