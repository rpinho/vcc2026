"""Compare identical finite targets and average the three public contexts equally."""
import runtime
from fractions import Fraction as F
R=runtime.ROOT
proof=runtime.json.loads((R/'artifacts/public_rational_verification.json').read_text())
assert len(proof)==18 and all(r['status'].startswith('PASS') for r in proof)
arms=['v5b','v9a','v9b'];sources=['Jurkat','HepG2','RPE1'];metrics=list(proof[0]['means'])
lookup={(r['source'],r['condition'],r['candidate']):r for r in proof}
details=[];summary=[]
for condition in ['native','projected400']:
 for metric in metrics:
  source_means={a:[] for a in arms};counts={}
  for source in sources:
   sets=[{t for t,v in lookup[source,condition,a]['per_target_rationals'].items() if v[metric] is not None} for a in arms]
   common=sorted(set.intersection(*sets));assert common;counts[source]=len(common);values={}
   for arm in arms:
    table=lookup[source,condition,arm]['per_target_rationals'];mean=sum((F(*map(int,table[t][metric])) for t in common),F(0))/len(common);source_means[arm].append(mean)
    values[arm]=dict(value=float(mean),rational=[str(mean.numerator),str(mean.denominator)])
   details.append(dict(source=source,condition=condition,metric=metric,common_targets=common,values=values))
  means={a:sum(source_means[a],F(0))/3 for a in arms};values={a:float(v) for a,v in means.items()}
  deltas={a:float(means[a]-means['v5b']) for a in arms if a!='v5b'}
  row=dict(condition=condition,metric=metric,matched_targets_per_source=counts,raw_equal_context_means=values,delta_from_v5b=deltas,v9b_minus_v9a=float(means['v9b']-means['v9a']),exact_rational_means={a:[str(v.numerator),str(v.denominator)] for a,v in means.items()})
  summary.append(row);print(row,flush=True)
runtime.save(R/'artifacts/public_matched_summary.json',dict(status='PASS exact rational metrics on identical finite cohorts; equal-context aggregation',summary=summary,per_context=details,qualification='Native values use observed held public cells. Projected400 reference values are estimated empirical-Bayes power, not measured 400-cell labels. Source coverage and transport differ from the official panel; these are RAW proxy metrics, not scaled hidden-score forecasts.'))
start=runtime.json.loads((R/'artifacts/public_final_manifest.json').read_text())['source_sha256']
assert all(start[n]==runtime.sha(R/n) for n in ['layout.py','call_budget.py','public_build.py','public_projection.py','public_compare.py'])
runtime.save(R/'artifacts/public_final_complete.json',dict(status='PASS final layout on three held public panels, official metric dispatch, independent rational raw metrics and matched-cohort comparison',final_layout_sha256=runtime.sha(R/'layout.py'),verifier_sha256=runtime.sha(R/'public_verify.py'),tables_sha256={p.name:runtime.sha(p) for p in R.glob('artifacts/public_*.parquet')}))
