"""Independent exact rational verification of score scaling and winning aggregate."""
import runtime
from fractions import Fraction as F
rt=runtime.ROOT;s=runtime.json.loads((rt/'artifacts/selection.json').read_text());scores={}
for fold in s['native_lineages']:
 folder=rt/('artifacts/evaluation_native' if fold=='H1' else f'artifacts/evaluation_{fold}_native');anchors=runtime.json.loads((folder/'anchors.json').read_text());base=anchors['tiled_baseline'];rep=anchors['replicate'];scores[fold]={}
 for row in s['ranking']:
  r=runtime.json.loads((folder/(row['variant']+'.json')).read_text());values=[]
  for metric,text in r['exact_scaled'].items():
   raw,b,c=map(lambda x:F(float(x)),[r['raw'][metric],base[metric],rep[metric]]);v=(raw-b)/(c-b)
   if metric=='expr_mse_unbiased_capped_norm':v=max(F(0),min(F(1),v))
   if metric=='de_wilcoxon_lfc_nmae':v=max(F(-6),v)
   assert v==F(text);values.append(v)
  assert len(values)==6;overall=sum(values,F(0))/6;assert overall==F(r['exact_overall']);scores[fold][row['variant']]=overall
for row in s['ranking']:
 v=sum((d[row['variant']] for d in scores.values()),F(0))/len(scores);assert v==F(row['exact_mean_native'])
a=s['a'];assert all(F(a['exact_mean_native'])>=F(r['exact_mean_native']) for r in s['ranking']);eligible=[r for r in s['ranking'] if r['fingerprint']!=a['fingerprint']];assert s['b_role']=='full-integration ablation' and s['b']['variant']=='iv_ash_dual_guide';assert s['b'] in eligible;assert F(s['numerical_runner_up']['exact_mean_native'])==max(F(r['exact_mean_native']) for r in eligible)
runtime.save(rt/'artifacts/selection_verification.json',dict(status='VERIFIED exact binary64-rational score scaling and selection arithmetic',variants=len(s['ranking']),lineages=list(scores),biological_quality='UNVERIFIED; metric aggregates remain floating-point public estimates'))
print('Exact scaling and selection VERIFIED',flush=True)
