"""Final anonymous public source receipts; no novelty claim."""
import runtime,requests,concurrent.futures,time
s=requests.Session();s.trust_env=False
urls={'leaderboard':'https://virtualcellchallenge.org/api/leaderboard?get_final=false','atlas_commit':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/commits?per_page=1','atlas_open_issues_prs':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/issues?state=open&per_page=10','sidechain_commit':'https://api.github.com/repos/saberhq/sidechain/commits?per_page=1'}
def fetch(item):
 name,url=item;r=s.get(url,timeout=30);r.raise_for_status();return name,dict(url=url,body=r.json())
receipts=dict(concurrent.futures.ThreadPoolExecutor(4).map(fetch,urls.items()));stamp=time.strftime('%Y-%m-%d %H:%M UTC',time.gmtime());runtime.save(runtime.ROOT/'artifacts/final_public_recheck.json',dict(checked_utc=stamp,receipts=receipts));entries=receipts['leaderboard']['body']['entries'];best=max(entries,key=lambda r:r.get('score_avg') or -999)
with (runtime.ROOT/'PRIOR_ART.md').open('a') as f:f.write('\nFinal public recheck '+stamp+': official table best '+str(best['team_name'])+' '+str(best['score_avg'])+'. Current AtlasShift/Sidechain commit heads and Atlas open issues/PRs retained in `artifacts/final_public_recheck.json`. No novelty claim.\n')
print(stamp,best['team_name'],best['score_avg'],flush=True)
