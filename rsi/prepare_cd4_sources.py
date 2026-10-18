"""Rebuild released transfer recipes on wider targets; never use own-context banks."""
import runtime as rt
import numpy as np,pandas as pd,anndata as ad,re,time
from scipy import sparse
from pathlib import Path
NT='non-targeting';OUT=rt.ROOT/'artifacts/cd4_sources';OUT.mkdir(exist_ok=True)
manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());needed=set(rt.json.loads(Path('artifacts/sources/target_universe.json').read_text()));needed.update(rt.json.loads(Path('artifacts/CD4_wide_targets.json').read_text()));needed=np.array(sorted(needed));rt.save(OUT/'target_universe.json',needed.tolist())
for name in ['H1','HepG2','Jurkat','K562','K562_GWPS_bulk']:
 path=OUT/f'{name}.npz'
 if path.exists():continue
 start=time.time();a=ad.read_h5ad(rt.ROOT.parent/'public_data'/f'{name}.h5ad',backed='r');genes=a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str);bulk=name.endswith('_bulk')
 if bulk:
  labels=np.array([NT if NT in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$',s).group(1) for s in a.obs_names]);n=a.obs.num_cells_filtered.fillna(0).to_numpy(float);valid=(n>=30)&((labels!=NT)|a.obs.core_control.to_numpy(bool));batches=np.repeat('unmatched',len(labels))
 else:
  labels=a.obs['gene' if 'gene' in a.obs else 'target_gene'].astype(str).to_numpy();n=np.ones(len(labels));valid=np.ones(len(labels),bool);batches=a.obs['gem_group' if 'gem_group' in a.obs else 'batch'].astype(str).to_numpy()
 targets=np.array(sorted(set(labels[valid])&set(needed)-{NT}));alllabels=np.concatenate(([NT],targets));codes=pd.Index(alllabels).get_indexer(labels);select=valid&(codes>=0);rows=np.flatnonzero(select);sums=np.zeros((len(alllabels),len(genes)));raws=np.zeros_like(sums);counts=np.bincount(codes[rows],weights=n[rows],minlength=len(alllabels));bybatch={};batchn={}
 # Store raw pooled summaries for exact v3 prior construction as well as mean-cell CPM.
 for lo in range(0,len(rows),512):
  rr=rows[lo:lo+512];x=sparse.csr_matrix(a.X[rr],dtype=float);assert np.isfinite(x.data).all() and np.all(x.data>=0)
  if not bulk:assert np.equal(x.data,np.floor(x.data)).all()
  depths=np.asarray(x.sum(1)).ravel();norm=x.multiply((1e6/np.maximum(depths,1))[:,None]).tocsr();cc=codes[rr];grp=sparse.csr_matrix((n[rr],(cc,np.arange(len(rr)))),shape=(len(alllabels),len(rr)));sums+=(grp@norm).toarray();raws+=(grp@x).toarray()
  if lo%51200==0:print(name,lo,len(rows),flush=True)
 # Bulk legacy code pools raw means THEN normalizes, not average of normalized constructs.
 if bulk:sums=raws/np.maximum(raws.sum(1,keepdims=True),1e-30)*1e6*counts[:,None]
 if len(set(genes))!=len(genes):
  unique=np.array(sorted(set(genes)));proj=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),pd.Index(unique).get_indexer(genes))),shape=(len(genes),len(unique)));sums=sums@proj;raws=raws@proj;genes=unique
 means=sums/np.maximum(counts[:,None],1);fine=np.clip(np.log2((means[1:]+.5)/(means[0]+.5)),-5,5).astype(np.float32);q04=np.clip(np.rint(fine/.04),-125,125).astype(np.int8);q08=np.clip(np.rint(q04.astype(np.float32)*float(np.array(.04,dtype=np.float32))/.08),-125,125).astype(np.int8)
 np.savez_compressed(path,genes=genes.astype(str),targets=targets.astype(str),fine=fine,lfc=q08,lfc_scale=np.array(.08,np.float32),baseline=means[0].astype(np.float32),n_cells=counts[1:].astype(np.int64),n_controls=np.array(int(counts[0])),cell_line=np.array('K562' if bulk else name),is_pseudobulk=np.array(bulk),raw_means=(raws/np.maximum(counts[:,None],1)).astype(np.float32),raw_labels=alllabels)
 rt.save(path.with_suffix('.json'),dict(source='../public_data/'+name+'.h5ad',targets=len(targets),genes=len(genes),sha256=rt.sha(path),seconds=time.time()-start));a.file.close();print('DONE',name,len(targets),round(time.time()-start),flush=True)
# Retained official-target rows must match; no silent historical-score reassignment.
results={}
for p in OUT.glob('*.npz'):
 d=np.load(p);old=np.load(rt.ROOT.parent/'model/artifacts/signatures'/p.name);ri=pd.Index(d['targets']).get_indexer(old['targets']);ci=pd.Index(d['genes']).get_indexer(old['genes']);assert np.all(ri>=0) and np.all(ci>=0)
 got=d['lfc'][np.ix_(ri,ci)];results[p.stem]=dict(n_old_targets=len(ri),quantized_bitwise_equal=np.array_equal(got,old['lfc']),max_q_difference=int(np.abs(got.astype(int)-old['lfc'].astype(int)).max()) if got.size else 0,baseline_bitwise_equal=np.array_equal(d['baseline'][ci],old['baseline']))
rt.save(OUT/'legacy_reproduction.json',results);print(results,flush=True)
