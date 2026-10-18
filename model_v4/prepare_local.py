"""Matched-control summaries from existing public single-cell files.
All arithmetic-CPM moments computed before aggregation; no local fold responses enter their own source pool.
"""
import runtime,re,time
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
rt=runtime.ROOT;out=rt/'artifacts/banks';out.mkdir(exist_ok=True)
targets=list(runtime.json.loads((rt/'artifacts/targets.json').read_text()))
# Additional existing public holdout targets permit honest cross-line diagnostics.
p=rt.parent/'model_v2/artifacts/validation/panels.json'
if p.exists():
 for panel in runtime.json.loads(p.read_text()).values():
  if isinstance(panel,list):targets+=panel
targets=np.array(list(dict.fromkeys(targets)));tset=set(targets)
for name in ['H1','HepG2','Jurkat','K562','RPE1','K562_GWPS_bulk']:
 path=out/f'{name}.npz'
 if path.exists():continue
 a=ad.read_h5ad(rt.parent/'public_data'/f'{name}.h5ad',backed='r');start=time.time()
 genes=a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
 if name.endswith('bulk'):
  lab=np.array(['non-targeting' if 'non-targeting' in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$',s).group(1) for s in a.obs_names]);n0=a.obs.num_cells_filtered.fillna(0).to_numpy(float);valid=n0>=20;valid&=(lab!='non-targeting')|a.obs.core_control.to_numpy(bool);batches=np.repeat('unmatched',a.n_obs)
 else:
  col='gene' if 'gene' in a.obs else 'target_gene';lab=a.obs[col].astype(str).to_numpy();lab=np.array(['non-targeting' if s.lower() in ['non-targeting','non-targeting control','ntc','control','ctrl','unperturbed'] or s.startswith('non-targeting') else s for s in lab]);n0=np.ones(a.n_obs);valid=np.ones(a.n_obs,bool);bcol=next((c for c in ['gem_group','batch','sample'] if c in a.obs),None);batches=a.obs[bcol].astype(str).to_numpy() if bcol else np.repeat('global',a.n_obs)
 valid&=np.isin(lab,list(tset|{'non-targeting'}));rows=np.flatnonzero(valid);keys=list(dict.fromkeys(zip(batches[rows],lab[rows])));lookup={k:i for i,k in enumerate(keys)};group=np.array([lookup[k] for k in zip(batches[rows],lab[rows])]);n=np.bincount(group,weights=n0[rows],minlength=len(keys));shape=(len(keys),len(genes));sums=np.zeros(shape);qp=np.zeros(shape);q2=np.zeros(shape);deps=np.zeros(len(keys))
 for lo in range(0,len(rows),512):
  rr=rows[lo:lo+512];x=sparse.csr_matrix(a.X[rr],dtype=float);gg=group[lo:lo+512]
  if name.endswith('bulk'):x=x.multiply(n0[rr,None]).tocsr()
  else:assert np.all(x.data==np.floor(x.data))
  dep=np.asarray(x.sum(1)).ravel();assert np.all(dep>0);member=sparse.csr_matrix((np.ones(len(rr)),(gg,np.arange(len(rr)))),shape=(len(keys),len(rr)));sums+=(member@x).toarray();proportion=x.multiply((1/dep)[:,None]).tocsr()
  if name.endswith('bulk'):proportion=proportion.multiply(n0[rr,None]).tocsr()
  qp+=(member@proportion).toarray()
  if not name.endswith('bulk'):q2+=(member @ proportion.multiply(proportion)).toarray()
  np.add.at(deps,gg,dep)
  if lo%51200==0:print(name,lo,len(rows),flush=True)
 if len(set(genes))!=len(genes):
  unique=np.array(list(dict.fromkeys(genes)));proj=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),pd.Index(unique).get_indexer(genes))),shape=(len(genes),len(unique)));sums=sums@proj;qp=qp@proj;q2=q2@proj;genes=unique
 controls={b:i for i,(b,t) in enumerate(keys) if t=='non-targeting'};cr=np.array(list(controls.values()),dtype=int);assert len(cr)>0;gp=sums[cr].sum(0)/deps[cr].sum();gm=qp[cr].sum(0)/n[cr].sum()*1e6
 nt=len(targets);ng=len(genes);tc=np.zeros((nt,ng));tm=np.zeros_like(tc);mp=np.zeros_like(tc);mm=np.zeros_like(tc);dv=np.zeros_like(tc);tn=np.zeros(nt,np.int64)
 for ti,t in enumerate(targets):
  rr=[i for i,(b,l) in enumerate(keys) if l==t]
  if not rr:continue
  tn[ti]=round(n[rr].sum());tc[ti]=sums[rr].sum(0);tm[ti]=qp[rr].sum(0)/n[rr].sum()*1e6
  for i in rr:
   ci=controls.get(keys[i][0]);lw=deps[i]/deps[rr].sum();cw=n[i]/n[rr].sum()
   cp=gp if ci is None else sums[ci]/deps[ci];cm=gm if ci is None else qp[ci]/n[ci]*1e6;mp[ti]+=lw*cp;mm[ti]+=cw*cm
   if ci is not None and not name.endswith('bulk'):
    vp=np.maximum(q2[i]/n[i]-(qp[i]/n[i])**2,0)/max(n[i]-1,1);vc=np.maximum(q2[ci]/n[ci]-(qp[ci]/n[ci])**2,0)/max(n[ci]-1,1);dv[ti]+=cw*cw*(vp+vc)*1e12
   else:
    # Conservative Poisson floor, explicitly a heuristic when cell moments are unavailable.
    dv[ti]+=cw*cw*(tm[ti]/max(tn[ti],1)+cm/max(n[cr].sum(),1))
 np.savez_compressed(path,source=np.array(name),targets=targets,genes=genes,n_cells=tn,target_count_sums=tc.astype('f4'),target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),matched_control_mean_cpm=mm.astype('f4'),global_control_probability=gp.astype('f4'),global_control_mean_cpm=gm.astype('f4'),measured_genes=np.ones(ng,bool),delta_cpm_variance=dv.astype('f4'),arithmetic_cpm_observed=np.array(not name.endswith('bulk')))
 runtime.save(rt/f'artifacts/{name}_local_summary.json',dict(cells=int(n.sum()),targets=int((tn>=20).sum()),groups=len(keys),genes=len(genes),seconds=time.time()-start,arithmetic_CPM_observed=not name.endswith('bulk'),matched_controls=not name.endswith('bulk')))
 a.file.close();print('DONE',name,round(time.time()-start),flush=True)
