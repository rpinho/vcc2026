"""Raw diagnostics on the HCT panel with undefined normalized MSE; no composite score."""
import runtime,time
from dataclasses import replace
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from cell_eval2 import EvalConfig,compute_metrics,aggregate_metrics
from cell_eval2.config import DEParams
from evaluate import VARIANTS,SIX,identity
from model import heads,template,emit,atlas
from independent_metrics import verify
from verify_counts import moment_certificate
rt=runtime.ROOT;fold='HCT116';NT='non-targeting';folder=rt/'artifacts/weak_HCT116';folder.mkdir(exist_ok=True)
a=ad.read_h5ad(rt/'artifacts/validation_HCT116.h5ad');ctrl=a[a.obs.evaluation_role.astype(str)=='input_control'].copy();real=a[a.obs.evaluation_role.astype(str)!='input_control'].copy();del a;genes=real.var_names.to_numpy();labels=real.obs.target_gene.astype(str).to_numpy();targets=np.array(sorted(set(labels)-{NT}));q,d,mean,bulk=template(ctrl,genes);cfg=replace(EvalConfig.from_preset('vcc2026'),pert_col='target_gene',device='cpu',num_threads=2,gather_threads=2,de=DEParams(backend='pdex'),cache_real=str(folder/'real_cache'));metrics=SIX+['expr_mse_unbiased_capped','expr_distance_unbiased'];official_panel=pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.to_numpy(str)
for variant in VARIANTS:
 path=folder/f'{variant}.json'
 if path.exists():continue
 start=time.time();stats={};(p,b),cov=heads(targets,genes,mean,bulk,variant,exclude=fold,stats=stats);chunks=[];obs=[];mom=[]
 try:
  for i,t in enumerate(targets):
   x=emit(q,p[i],b[i],d,atlas.seed_for(fold+':'+t,20260910),variant);mom.append(dict(target=t,**moment_certificate(x,d,p[i],b[i])));chunks.append(sparse.csr_matrix(x));obs.extend([t]*400)
 except ValueError as error:
  runtime.save(folder/f'{variant}_failure.json',dict(variant=variant,error=str(error),target=t));print('Failure',variant,str(error),flush=True);continue
 chunks.append(real.X[labels==NT]);obs.extend([NT]*int((labels==NT).sum()));pred=ad.AnnData(sparse.vstack(chunks,format='csr'),obs=pd.DataFrame({'target_gene':obs},index=[f'p{i}' for i in range(len(obs))]),var=real.var.copy());del chunks
 tidy=compute_metrics(pred,real,config=cfg);raw=aggregate_metrics(tidy,metrics=metrics).to_pandas().set_index('metric')['mean'].to_dict();ind=verify(pred,real,raw,exclusion_panel=official_panel);assert abs(ind['numerator_sum']/len(targets)-raw['expr_mse_unbiased_capped'])<2e-8
 rec=dict(variant=variant,fold=fold,status='UNVERIFIED raw diagnostics only; composite undefined because half-panel replicate-anchor MSE denominator nonpositive; full-panel raw normalization positive but nearly zero',overall=None,raw=raw,independent=ind,moments=mom,adaptive_stats=stats,prediction_sha256=identity(pred),model_sha256=runtime.sha(rt/'model.py'),source_identity_sha256=runtime.sha(rt/'artifacts/source_identity.json'),seconds=time.time()-start);runtime.save(path,rec);print(variant,'raw capped loss',raw['expr_mse_unbiased_capped'],'PDS',raw['pds_cosine'],'seconds',round(time.time()-start),flush=True);del pred
records=[runtime.json.loads(p.read_text()) for p in folder.glob('*.json') if p.stem in VARIANTS];pd.DataFrame([dict(variant=r['variant'],**r['raw']) for r in records]).to_csv(folder/'raw_diagnostics.csv',index=False)
print('Weak-panel diagnostics complete',flush=True)
