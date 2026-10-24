"""Public 400/18400 power sensitivity; estimated reference DE, not measured cells.

The fixed empirical-Bayes normal-mixture procedure is copied from v8's public
budget estimator. Official metric dispatch remains exact on the resulting tables.
"""
import runtime
import numpy as np,polars as pl
from scipy.special import ndtri,ndtr
from scipy.stats import false_discovery_control
R=runtime.ROOT
def posterior_effect(z):
 variance=np.array([0.,.25,1.,4.,16.,64.,256.]);pi=np.ones(len(variance))/len(variance)
 density=np.exp(-z[:,None]**2/(2*(1+variance[None,:])))/np.sqrt(1+variance[None,:])
 for _ in range(60):
  weight=density*pi[None,:];weight/=np.maximum(weight.sum(1,keepdims=True),1e-300);pi=weight.mean(0)
 return z*np.sum(weight*(variance/(1+variance))[None,:],axis=1),pi.tolist()
receipts=[]
for source in ['Jurkat','HepG2','RPE1']:
 runtime.guard();meta=runtime.json.loads((R/f'artifacts/public_{source}.json').read_text());real=pl.read_parquet(R/f'artifacts/public_{source}.parquet');parts=[];records=[]
 for target,native_n in meta['targets'].items():
  f=real.filter(pl.col('target')==target);p=f['p_value'].to_numpy();valid=np.isfinite(p)&(p>=0)&(p<=1)
  z=-ndtri(np.maximum(p[valid],1e-300)/2);p400=np.ones(len(p));q400=np.ones(len(p));weights=None
  if valid.any():
   latent,weights=posterior_effect(z);factor=np.sqrt((1/native_n+1/meta['reference_cells'])/(1/400+1/18400));p400[valid]=2*ndtr(-np.abs(latent*factor));q400[valid]=false_discovery_control(p400[valid])
  projected=f.with_columns(pl.Series('p_value',p400),pl.Series('p_adj',q400));parts.append(projected)
  not_own=f['feature'].to_numpy()!=target
  records.append(dict(target=target,native_cells=native_n,native_conf=int(np.sum((f['p_adj'].to_numpy()<.05)&not_own)),projected_conf=int(np.sum((q400<.05)&not_own)),mixture=weights))
 pl.concat(parts).write_parquet(R/f'artifacts/public_{source}_projected400.parquet')
 receipts.append(dict(source=source,records=records,qualification='ESTIMATED reference power at 400 treated/18400 control cells; EB normal-mixture posterior means from observed MWU magnitudes. Native LFC signs retained. Does not fabricate observed cells or reproduce hidden labels; sign noise, transport and partial coverage remain. Prediction tables retain their actual matched-native-control producer values.'))
runtime.save(R/'artifacts/public_projection.json',receipts)
