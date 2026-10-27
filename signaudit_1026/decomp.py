import numpy as np,polars as pl,json
from pathlib import Path
H=Path('../h2h2_1026')
p=pl.read_parquet(H/'artifacts/de_pred_B0.parquet');r=pl.read_parquet(H/'artifacts/de_real.parquet')
d=p.join(r,on=['target','feature'],suffix='_r')
with np.load(H/'results/predictions_B0.npz') as z:resp=z['response'];genes=z['genes'].astype(str);targets=list(z['targets'].astype(str))
G=pl.DataFrame({'feature':genes.tolist()}).with_row_index('gi');T=pl.DataFrame({'target':targets}).with_row_index('ti')
d=d.join(G,on='feature',how='left').join(T,on='target',how='left')
d=d.with_columns(pl.Series('intended',resp[d['ti'].to_numpy(),d['gi'].fill_null(0).to_numpy()]))
d=d.with_columns([(pl.col('p_adj')<0.05).alias('psig'),(pl.col('p_adj_r')<0.05).alias('rsig'),
 (pl.col('log2_fold_change').sign()==pl.col('log2_fold_change_r').sign()).alias('agree'),
 (pl.col('log2_fold_change').sign()==pl.col('intended').sign()).alias('kept')])
per=d.group_by('target').agg([pl.col('psig').sum().alias('n_pred'),pl.col('rsig').sum().alias('n_real'),(pl.col('psig')&pl.col('agree')).sum().alias('k')]).with_columns([
 (pl.col('k')/pl.max_horizontal('n_pred','n_real')).fill_nan(0).alias('fid'),(pl.col('k')/pl.col('n_pred')).fill_nan(0).alias('prec')])
per=per.with_columns(pl.when(pl.col('n_real')<10).then(pl.lit('weak')).when(pl.col('n_real')<1000).then(pl.lit('mid')).otherwise(pl.lit('strong')).alias('stratum'))
print(per.group_by('stratum').agg([pl.len(),pl.col('fid').mean(),pl.col('prec').mean(),pl.col('n_pred').median(),pl.col('n_real').median(),(pl.col('n_pred')==0).sum().alias('zero_calls')]).sort('stratum'))
print('targets n_real>n_pred:',(per['n_real']>per['n_pred']).sum(),' zero-call targets:',(per['n_pred']==0).sum())
c=d.filter('psig')
print('called kept-intent agree:',c.filter('kept')['agree'].mean(),' flipped agree:',c.filter(~pl.col('kept'))['agree'].mean(),' flipped frac:',1-c['kept'].mean())
# is the flipped set the selected layout calls or incidental?
import json as J
rec={x['target']:set(x['selected_genes']) for x in J.load(open(H/'results/render_B0.json'))['records']}
c=c.with_columns(pl.struct(['target','feature']).map_elements(lambda s:s['feature'] in rec.get(s['target'],()),return_dtype=pl.Boolean).alias('layout_call'))
print(c.group_by('layout_call').agg([pl.len(),pl.col('agree').mean(),pl.col('kept').mean()]))
# best-case: if every called gene took the real-majority-consistent intended sign... precision of intended among real-sig
print('real-sig genes: intended sign right',d.filter('rsig').select((pl.col('intended').sign()==pl.col('log2_fold_change_r').sign()).mean()).item())
