"""Late compliance audit, public sources only; no novelty claim."""
import runtime,requests,concurrent.futures
s=requests.Session();s.trust_env=False
urls={'leaderboard':'https://virtualcellchallenge.org/api/leaderboard?get_final=false','code_recent':'https://api.github.com/search/repositories?q=virtual+cell+challenge+2026&sort=updated&per_page=10','atlas_issues':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/issues?state=open&per_page=10','sidechain':'https://api.github.com/repos/saberhq/sidechain/commits?per_page=1','labs':'https://api.github.com/search/code?q=VCC+repo:google-deepmind/alphaevolve_results','atlas':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/commits?per_page=1'}
def get(pair):
 k,u=pair
 try:
  r=s.get(u,timeout=30);return k,dict(url=u,status=r.status_code,body=r.json())
 except Exception as e:return k,dict(url=u,error=str(e))
r=dict(concurrent.futures.ThreadPoolExecutor(4).map(get,urls.items()));runtime.save(runtime.ROOT/'artifacts/prior_art_final.json',r)
rows=r['leaderboard'].get('body',{}).get('entries',[]);best=max(rows,key=lambda x:x.get('score_avg') or -999) if rows else {};print('Best',best.get('team_name'),best.get('score_avg'))
text='''# Prior art — late audit, 2026-10-16

**Process exception:** this required standalone file was not written before construction. Public method/configuration/license and leaderboard receipts were retained during construction, but that does not satisfy RULES.md’s first-file requirement. This late audit does not repair that ordering violation. No novelty or record claim is made.

The method anchor is [AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026), public participant-label overall 0.1545618019 (September 10), scaled MSE about 0.0471. Its released kernel and recipe are reproduced locally with documented resource-adapted H1/K562 statistics, not hidden-context artifacts. [Sidechain](https://github.com/saberhq/sidechain) reports SER-17 0.1657 on October 8; its October 3 SER-14 adaptive method reports 0.1626. Neither number is our local validation. Both code licenses are MIT.

Live [official table](https://virtualcellchallenge.org/api/leaderboard?get_final=false), Atlas open issues (including PRs), recent GitHub repository search and current project commits are captured in artifacts/prior_art_final.json. The best table entry on this read was TEAM, SCORE. The strongest public recipe actually reproduced here is AtlasShift; no claim of reproducing the leaderboard leader is made.

Primary papers checked in search: [SCALE](https://arxiv.org/abs/2603.17380), [in-the-wild benchmark](https://arxiv.org/abs/2604.27646), [VCLMU](https://arxiv.org/abs/2610.04475), [MRFlow](https://arxiv.org/abs/2610.13966). These address broader prediction/transport rather than proving this submission’s performance.

Community checked: [Arc competitor showcase, October 15](https://arcinstitute.org/news/vcc-competitor-showcase-2026), [Sidechain public dated ledger](https://saberhq.com/sidechain/models/), [public participant writeup](https://yraj.in/ai/virtual-cell/crispri/single-cell/2026/09/09/predicting-the-unseen-vcc-2026.html). No account creation or contacts.

AI-lab/agent catalogs required by RULES: BakeLab/Epi-results, google-deepmind alphaevolve_results and repository_of_problems, openai/math, EinsteinArena, plasma-ai/erdos, erdosproblems.ai, kisonecat/openai-math. The GitHub code search audit requires authentication (401); no credentials used. These are not established VCC sources. Coverage of these catalogs and awesome-lists is incomplete; we do not claim an exhaustive prior-art review. The user supplied the ranked VCC research/SOURCES inventory, retained as artifacts/source_inventory_review.json, and implementation follows those identified primary projects.
'''.replace('TEAM',str(best.get('team_name'))).replace('SCORE',str(best.get('score_avg')))
(runtime.ROOT/'PRIOR_ART.md').write_text(text)
