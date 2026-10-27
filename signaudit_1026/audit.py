"""Sign audit of B0's DE calls on the h2h2 H1 harness (analysis only; reads captured official DE tables)."""
import json,sys,numpy as np,polars as pl
from pathlib import Path
H=Path(__file__).resolve().parent.parent/'h2h2_1026'
arm=sys.argv[1] if len(sys.argv)>1 else 'B0'
p=pl.read_parquet(H/f'artifacts/de_pred_{arm}.parquet');r=pl.read_parquet(H/'artifacts/de_real.parquet')
d=p.join(r,on=['target','feature'],suffix='_r')
with np.load(H/'results/predictions_B0.npz') as z:resp=z['response'];genes=z['genes'].astype(str);targets=list(z['targets'].astype(str))
gi={g:i for i,g in enumerate(genes)};ti={t:i for i,t in enumerate(targets)}
d=d.with_columns(pl.struct(['target','feature']).map_elements(lambda s:float(resp[ti[s['target']],gi[s['feature']]]) if s['feature'] in gi else np.nan,return_dtype=pl.Float64).alias('intended'))
d=d.with_columns([(pl.col('p_adj')<0.05).alias('psig'),(pl.col('p_adj_r')<0.05).alias('rsig'),
  (pl.col('log2_fold_change').sign()==pl.col('log2_fold_change_r').sign()).alias('agree_real'),
  (pl.col('log2_fold_change').sign()==pl.col('intended').sign()).alias('render_keeps_intent'),
  (pl.col('intended').sign()==pl.col('log2_fold_change_r').sign()).alias('intent_right')])
# per-target fidelity reproduction: k/max(n_pred,n_real), k = pred-sig genes whose sign agrees with real
per=d.group_by('target').agg([pl.col('psig').sum().alias('n_pred'),pl.col('rsig').sum().alias('n_real'),
  (pl.col('psig')&pl.col('agree_real')).sum().alias('k')]).with_columns((pl.col('k')/pl.max_horizontal('n_pred','n_real')).fill_nan(0).alias('fid'))
c=d.filter('psig')
def frac(x):return f"{x.mean():.4f} (n={len(x)})"
out={'arm':arm,'fid_reproduced_mean':float(per['fid'].mean()),'n_pred_median':float(per['n_pred'].median()),'n_real_median':float(per['n_real'].median()),
 'called_sign_agree_real_all':float(c['agree_real'].mean()),
 'called_and_real_sig':float(c.filter('rsig')['agree_real'].mean()),'called_not_real_sig':float(c.filter(~pl.col('rsig'))['agree_real'].mean()),
 'n_called':len(c),'n_called_real_sig':int(c['rsig'].sum()),
 'render_keeps_intent':float(c['render_keeps_intent'].mean()),'intent_right_on_called':float(c['intent_right'].mean()),
 'called_up_agree':float(c.filter(pl.col('log2_fold_change')>0)['agree_real'].mean()),'called_down_agree':float(c.filter(pl.col('log2_fold_change')<0)['agree_real'].mean()),
 'frac_called_up':float((c['log2_fold_change']>0).mean()),
 'real_nonsig_frac_negative_lfc':float((d.filter(~pl.col('rsig'))['log2_fold_change_r']<0).mean()),
}
# by real-expression proxy: |real lfc| not available as expression; use pred p-value bins instead
for lo,hi in [(0,1e-10),(1e-10,1e-4),(1e-4,0.05)]:
  s=c.filter((pl.col('p_adj')>=lo)&(pl.col('p_adj')<hi));out[f'agree_padj_{lo:g}_{hi:g}']=f"{s['agree_real'].mean():.4f} n={len(s)}"
print(json.dumps(out,indent=1));Path(f'audit_{arm}.json').write_text(json.dumps(out,indent=1))
