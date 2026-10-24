"""Exact >=30-cell response-feature family coverage, one row per mock target."""
import common as job
import io_guard;io_guard.install()
import numpy as np,pandas as pd,argparse
from pathlib import Path
from collections import Counter

def coverage(folder,targets):
    found={t:set() for t in targets}
    for p in folder.glob('*.npz'):
        with np.load(p) as z:
            f=str(z['family'])
            for t,n in zip(z['targets'].astype(str),z['n_cells']):
                if t in found and n>=30:found[t].add(f)
    return found

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--targets',default=str(job.REPO/'model_v6/mock_final/panel.csv'))
    ap.add_argument('--data',default=str(job.ROOT/'data'));ap.add_argument('--before',default=str(job.REPO/'model_billm/data'))
    ap.add_argument('--out-prefix',default=str(job.ROOT/'artifacts/mock_family_coverage'));a=ap.parse_args()
    prefix=Path(a.out_prefix).resolve();assert prefix.is_relative_to(job.ROOT);prefix.parent.mkdir(parents=True,exist_ok=True)
    frame=pd.read_csv(a.targets);targets=sorted(set(frame['target_gene' if 'target_gene' in frame else frame.columns[0]].astype(str))-{'non-targeting'})
    before=coverage(Path(a.before),targets);after=coverage(Path(a.data),targets)
    rows=[dict(target_gene=t,before=len(before[t]),after=len(after[t]),before_families=','.join(sorted(before[t])),after_families=','.join(sorted(after[t]))) for t in targets]
    pd.DataFrame(rows).to_csv(prefix.with_suffix('.csv'),index=False)
    job.save(prefix.with_suffix('.json'),dict(targets=len(targets),threshold=30,
        before_histogram=dict(Counter(map(len,before.values()))),after_histogram=dict(Counter(map(len,after.values()))),
        before_generic_fallback=sum(not f for f in before.values()),after_generic_fallback=sum(not f for f in after.values()),
        panel_sha256=job.sha(a.targets),csv_sha256=job.sha(prefix.with_suffix('.csv'))))
    print('COVERAGE',job.json.loads(prefix.with_suffix('.json').read_text()),flush=True)
if __name__=='__main__':main()
