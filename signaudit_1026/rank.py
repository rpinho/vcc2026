import numpy as np,polars as pl,json
from pathlib import Path
H=Path('../h2h2_1026')
p=pl.read_parquet(H/'artifacts/de_pred_B0.parquet');r=pl.read_parquet(H/'artifacts/de_real.parquet')
d=p.join(r,on=['target','feature'],suffix='_r').with_columns([(pl.col('p_adj')<0.05).alias('psig'),(pl.col('p_adj_r')<0.05).alias('rsig'),
 (pl.col('log2_fold_change').sign()==pl.col('log2_fold_change_r').sign()).alias('agree')])
nreal=d.group_by('target').agg(pl.col('rsig').sum().alias('n_real'))
rows=[]
for x in json.load(open(H/'results/render_B0.json'))['records']:
  for k,g in enumerate(x['selected_genes']):rows.append((x['target'],g,k))
R=pl.DataFrame(rows,schema=['target','feature','rank'],orient='row')
c=R.join(d,on=['target','feature']).join(nreal,on='target')
c=c.with_columns(pl.when(pl.col('n_real')<1000).then(pl.lit('weak+mid')).otherwise(pl.lit('strong')).alias('st'))
edges=[0,50,100,200,300,450,600,750,885]
c=c.with_columns(pl.col('rank').cut(edges[1:-1],left_closed=True).alias('bin'))
print(c.group_by(['st','bin']).agg([pl.len(),pl.col('psig').mean().alias('realised'),pl.col('agree').filter(pl.col('psig')).mean().alias('prec_called'),pl.col('agree').mean().alias('sign_agree_all')]).sort(['st','bin']))
# simulated budget cut on weak+mid: precision of top-N realised calls + incidental calls
