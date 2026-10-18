"""Re-scale every local variant record two ways (exact rational arithmetic via verify_counts.scale_record):
  official-def : model_v4/v2 cached 11-target tiled generic baseline (official construction, self-share 1/11)
  loo          : leave-one-perturbation-out tiled baseline (loo_baseline.py; mimics 300-target self-share)
Replicate anchor identical in both."""
import runtime, json, sys
import pandas as pd
from verify_counts import scale_record
import evaluate_v5 as E
SIX = E.SIX
view = sys.argv[1] if len(sys.argv) > 1 else 'depth20k'
_, _, anchors = E.load_fold(view)
base, rep = anchors['tiled_baseline'], anchors['replicate']
loo_p = E.OUT / 'loo_baseline.json'
loo = json.loads(loo_p.read_text()) if loo_p.exists() else None
rows = []
for p in sorted(E.OUT.glob('v5_*.json')):
    r = json.loads(p.read_text()); raw = r['raw']
    s1, o1 = scale_record(raw, base, rep, SIX)
    row = dict(variant=r['variant'], overall=float(o1), **{k.replace('de_wilcoxon_', '').replace('expr_mse_unbiased_capped_norm', 'mse'): round(float(v), 4) for k, v in s1.items()},
               raw_mse_norm=round(raw[SIX[1]], 4), alpha=r.get('alpha'))
    if loo:
        s2, o2 = scale_record(raw, loo, rep, SIX)
        row.update(overall_loo=float(o2), mse_loo=round(float(s2[SIX[1]]), 4))
    rows.append(row)
df = pd.DataFrame(rows).sort_values('overall', ascending=False)
pd.set_option('display.width', 250); print(df.to_string(index=False))
df.to_csv(E.OUT / 'rescored.csv', index=False)
