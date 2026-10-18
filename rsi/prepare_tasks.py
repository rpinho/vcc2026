"""Freeze availability-only panels, split manifests, and disjoint real controls."""
import runtime as rt
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,anndata as ad
NT='non-targeting';ROOT=rt.ROOT

def order(items,salt):return sorted(items,key=lambda x:hashlib.sha256((salt+str(x)).encode()).digest())
def choose(obs,label,batch,targets,n,exact=False):
 rng=np.random.default_rng(20261016);labs=obs[label].astype(str).replace({'Non-Targeting':NT}).to_numpy();bs=obs[batch].astype(str).to_numpy();rows=[];roles=[];wanted={};batchnames=sorted(set(bs));q=n//len(batchnames)
 for t in targets:
  ix=np.flatnonzero(labs==t)
  if exact:
   sel=np.concatenate([rng.choice(ix[bs[ix]==b],q,replace=False) for b in batchnames])
  else:sel=rng.choice(ix,n,replace=False)
  rows.extend(sel);roles.extend(['perturbation']*len(sel));wanted[t]=pd.Series(bs[sel]).value_counts().to_dict()
 # Controls match the pooled perturbation mixture. Target-specific discrepancy audited.
 nt=np.flatnonzero(labs==NT);mixture=pd.Series(bs[rows]).value_counts().reindex(batchnames,fill_value=0).to_numpy(float);mixture/=mixture.sum();cap=np.array([int((bs[nt]==b).sum())//2 for b in batchnames]);num=min(2000,int(cap.sum()));quota=np.minimum(np.floor(num*mixture).astype(int),cap)
 while quota.sum()<num:
  options=np.flatnonzero(quota<cap);i=max(options,key=lambda i: (num*mixture[i]-quota[i],-i));quota[i]+=1
 if exact:
  # Same quota per batch, not an approximate rounding of a pooled distribution.
  per=min(100,int(cap.min()));quota[:]=per
 for b,k in zip(batchnames,quota):
  if not k:continue
  sel=rng.choice(nt[bs[nt]==b],2*k,replace=False);rows.extend(sel);roles.extend(['input_control']*k+['reference_control']*k)
 selected=obs.iloc[rows].copy();selected['target_gene']=labs[rows];selected['batch']=bs[rows];selected['evaluation_role']=roles;selected['source_row']=rows
 distribution=quota/quota.sum();l1={t:float(np.abs(np.array([wanted[t].get(b,0) for b in batchnames])/n-distribution).sum()) for t in targets}
 return selected,{'exact_target_batch_matching':exact,'per_target_control_mixture_L1':l1,'control_batch_quotas':dict(zip(batchnames,map(int,quota))),'input_reference_disjoint':True}

def main():
 out=ROOT/'artifacts/tasks';out.mkdir(exist_ok=True)
 manifest={};meta=ROOT/'artifacts/H1_full_obs.parquet'
 if meta.exists():
  obs=pd.read_parquet(meta);obs=obs[obs.batch.astype(str).str.startswith('Flex_2_')].copy();obs['source_row_original']=obs.index # Preserve actual physical row via original metadata below.
  allobs=pd.read_parquet(meta);obs['source_row_original']=pd.Index(allobs.index).get_indexer(obs.index)
  counts=obs.groupby(['target_gene','batch'],observed=True).size().unstack(fill_value=0);eligible=counts.index[(counts.min(axis=1)>=10)&(counts.index!=NT)].tolist();eligible=order(eligible,'H1_partition:');hold=set(eligible[:20])
  for split in ['search','holdout']:
   targets=sorted(set(eligible)-hold if split=='search' else hold);sel,audit=choose(obs,'target_gene','batch',targets,160,True);sel['source_row']=sel.source_row_original.to_numpy();name='H1_'+split;sel.to_parquet(out/(name+'_selection.parquet'));manifest[name]=dict(context='H1',split=split,targets=targets,n_real=160,n_pred=400,exact_matched=True,source='H1_full',audit=audit)
 for context,n in [('K562',200),('HepG2',64),('Jurkat',100),('RPE1',100)]:
  path=ROOT.parent/'public_data'/f'{context}.h5ad';a=ad.read_h5ad(path,backed='r');obs=a.obs.copy();col='gene' if 'gene' in obs else 'target_gene';counts=obs[col].astype(str).value_counts();genes=a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str);eligible=set(counts.index[counts>=n])&set(genes)-{NT};targets=order(eligible,context+':panel')[:300];assert len(targets)==300
  split='holdout' if context=='RPE1' else 'search';name=context+'_'+split;sel,audit=choose(obs,col,'gem_group',targets,n);sel.to_parquet(out/(name+'_selection.parquet'));manifest[name]=dict(context=context,split=split,targets=sorted(targets),n_real=n,n_pred=400,exact_matched=False,source='local_raw',audit=audit);a.file.close()
 for context in ['HCT116','HEK293T']:
  obs=pd.read_parquet(ROOT/f'artifacts/{context}_cells.parquet');offset=np.load(ROOT/f'artifacts/{context}_offset_mapping.npz');obs=obs[offset['logical_start'][obs.cell_integer_id.to_numpy()//5000]>=0].copy();counts=obs.gene_target.astype(str).value_counts();available=set(np.load(ROOT.parent/f'public_data/model_v4_signature_banks/{context}.npz')['targets']);targets=order((set(counts.index[counts>=64])&available)-{'Non-Targeting'},context+':panel')[:300];assert len(targets)==300,(context,len(targets));split='holdout' if context=='HEK293T' else 'search';name=context+'_'+split;sel,audit=choose(obs,'gene_target','sample',targets,64);sel.to_parquet(out/(name+'_selection.parquet'));manifest[name]=dict(context=context,split=split,targets=sorted(targets),n_real=64,n_pred=400,exact_matched=False,source='orion_raw',audit=audit)
 for name,m in manifest.items():m['selection_sha256']=rt.sha(out/(name+'_selection.parquet'));m['eligible_for_primary_calibration']=m['exact_matched'];m['seal']='HOLDOUT never scored by calibration or search' if m['split']=='holdout' else 'SEARCH'
 rt.save(out/'manifest.json',manifest);print({k:(len(v['targets']),v['n_real'],v['exact_matched']) for k,v in manifest.items()})
if __name__=='__main__':main()
