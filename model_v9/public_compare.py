"""Official producer -> prepare_de -> enrolled dispatcher -> native aggregation."""
import runtime
import polars as pl
from cell_eval2.run import aggregate_metrics
from scorer_check import score,NAMES
from cell_eval2.metrics.direction import _components
R=runtime.ROOT
results=[]
for source in ['Jurkat','HepG2','RPE1']:
 for condition in ['native','projected400']:
  runtime.guard();suffix='' if condition=='native' else '_projected400';real=pl.read_parquet(R/f'artifacts/public_{source}{suffix}.parquet')
  for arm in ['v5b','v9a','v9b']:
   pred=pl.read_parquet(R/f'artifacts/public_pred_{arm}_{source}.parquet');rows,prepared=score(pred,real,return_prepared=True);tidy=pl.DataFrame(rows);agg=aggregate_metrics(tidy,metrics=NAMES)
   values={r['metric']:r['mean'] for r in agg.to_dicts()};counts=tidy.group_by('metric').agg(pl.col('value').is_finite().sum().alias('finite_targets')).to_dicts()
   rec=dict(source=source,condition=condition,candidate=arm,public_targets=real['target'].n_unique(),native_metrics=values,finite_cohorts=counts,official_count_components=_components(prepared).select(['target','n_conf','n_pred','k']).to_dicts(),qualification='Exact official prepare_de/dispatch_de_metrics/aggregate_metrics path. Native input: both full official DE producers. Projected400 input: ESTIMATED EB real-side power from public_projection.py, not observed 400-cell labels. Public partial-gene proxy; not an official scaled score prediction.')
   results.append(rec);print(rec,flush=True);runtime.save(R/'artifacts/public_comparison.json',results);tidy.write_csv(R/f'artifacts/public_metrics_{arm}_{source}_{condition}.csv')
