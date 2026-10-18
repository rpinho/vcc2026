import runtime as rt
import requests,time
from pathlib import Path
urls={
 'leaderboard':'https://virtualcellchallenge.org/api/leaderboard?get_final=false',
 'github_recent':'https://api.github.com/search/repositories?q=virtual+cell+challenge+2026&sort=updated&per_page=30',
 'scorer_issues':'https://api.github.com/repos/ArcInstitute/cell-eval2/issues?state=open&per_page=30',
 'atlas_issues':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/issues?state=open&per_page=30',
 'sidechain_issues':'https://api.github.com/repos/saberhq/sidechain/issues?state=open&per_page=30',
 'scorer_head':'https://api.github.com/repos/ArcInstitute/cell-eval2/commits/main',
 'atlas_head':'https://api.github.com/repos/kaipengm2/Virtual-Cell-Challenge-2026/commits/main',
 'sidechain_head':'https://api.github.com/repos/saberhq/sidechain/commits/main',
 'scorer_metrics':'https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics.md',
 'h1_asset':'https://github.com/forrestsheldon/vcc2026-h1-benchmark/releases/download/v0.4.0/vcc2026-h1-benchmark-v2.zip',
}
http=requests.Session();http.trust_env=False;receipts=[]
for name,url in urls.items():
 p=Path('artifacts/prior')/name;p.parent.mkdir(exist_ok=True,parents=True)
 try:
  r=http.get(url,timeout=45);p.write_bytes(r.content);receipts.append(dict(name=name,url=url,status=r.status_code,bytes=len(r.content),sha256=rt.sha(p)))
 except Exception as e:receipts.append(dict(name=name,url=url,error=str(e)))
rt.save('artifacts/prior/receipts.json',receipts)
print(receipts)
