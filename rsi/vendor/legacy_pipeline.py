"""Zero-shot perturbation transfer and streaming submission construction."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[k]='1'
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('NUMBA_NUM_THREADS','4')
import argparse, hashlib, json, pathlib, subprocess, sys, time, shutil, functools
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import h5py
ROOT=pathlib.Path(__file__).resolve().parent
NT='non-targeting'
BASE_VARIANTS=['mean','knockdown','transfer','shrink','adapt','adapt_strong','nearest','bulk_transfer','bulk_shrink']

def csv_labels(path):
 a=pd.read_csv(path,header=None,dtype=str).iloc[:,0].str.strip().tolist()
 if a and a[0].lower() in ('gene','gene_name','gene_names','target_gene','perturbation','gene_symbol'):a=a[1:]
 if not a or len(set(a))!=len(a):raise ValueError(f'Empty or duplicate labels: {path}')
 return np.array(a)

def gene_symbols(a):
 for key in ('gene_name','gene_symbol','gene_symbols','feature_name','gene'):
  if key in a.var:return a.var[key].astype(str).to_numpy()
 return a.var_names.astype(str).to_numpy()

def raw_matrix(a):
 # Do not infer counts from a filename; choose only verified count space.
 for key in ('counts','count','raw_counts'):
  if key in a.layers:return a.layers[key]
 return a.X

def check_counts(x):
 v=x.data if sp.issparse(x) else np.asarray(x)
 if np.any(~np.isfinite(v)) or np.any(v<0) or np.any(v!=np.floor(v)):raise ValueError('Input is not raw integer counts')

def control_summary(path,genes,expected=None):
 a=ad.read_h5ad(path,backed='r')
 if 'context' not in a.obs:raise ValueError(f'{path}: context column required')
 if 'target_gene' in a.obs and set(a.obs.target_gene.astype(str))!={NT}:raise ValueError(f'{path}: non-control cells in a control file')
 labels=a.obs.context.astype(str).unique()
 if len(labels)!=1:raise ValueError(f'{path}: ambiguous context {labels}')
 context=labels[0]
 if expected and context!=expected:raise ValueError(f'Filename/obs context disagreement: {path} != {context}')
 syms=a.var_names.astype(str).to_numpy()
 if not set(genes).issubset(set(syms)):syms=gene_symbols(a)
 ix=pd.Index(syms).get_indexer(genes)
 if np.any(ix<0):raise ValueError('Control missing official genes')
 if len(syms)!=len(genes):raise ValueError('Control gene set differs from official list')
 sums=np.zeros(len(genes)); sq=np.zeros(len(genes)); normsum=np.zeros(len(genes)); libs=[]
 for i in range(0,a.n_obs,256):
  x=a.X[i:i+256,:][:,ix];check_counts(x)
  sums+=np.asarray(x.sum(axis=0)).ravel();l=np.asarray(x.sum(axis=1)).ravel();libs.extend(l)
  y=x.multiply((1e4/np.maximum(l,1))[:,None]) if sp.issparse(x) else x*(1e4/np.maximum(l,1))[:,None]
  normsum+=np.asarray(y.sum(axis=0)).ravel()
  sq+=np.asarray(y.multiply(y).sum(axis=0)).ravel() if sp.issparse(y) else (y*y).sum(axis=0)
 if 'gene_id' in a.var:gene_ids=a.var['gene_id'].astype(str).to_numpy()[ix]
 elif all(str(g).startswith('ENSG') for g in a.var_names):gene_ids=a.var_names.astype(str).to_numpy()[ix]
 else:gene_ids=np.array([annotation_aliases().get(g,'') for g in genes])
 gene_ids=np.array([str(g).split('.')[0] for g in gene_ids])
 a.file.close();mean=sums/a.n_obs;libs=np.array(libs)
 if np.any(libs>1e6) or np.any(libs<=0):raise ValueError('Invalid control library sizes')
 normmean=normsum/a.n_obs
 var=np.maximum(sq/a.n_obs-normmean**2,0)
 # Estimate overdispersion in normalized count space; gamma mixing plus multinomial sampling.
 theta=np.clip(normmean**2/np.maximum(var-normmean*np.mean(1e4/libs),1e-4),.1,1000)
 return dict(context=context,genes=genes,gene_ids=gene_ids,n_controls=int(a.n_obs),mean=mean,libs=libs,theta=theta,path=str(path),sha256=sha256(path))

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

def sample_counts(mean,libs,theta,n,rng):
 mean=np.asarray(mean,dtype=float)
 if np.any(~np.isfinite(mean)) or np.any(mean<0) or mean.sum()<=0:raise ValueError('Invalid predicted mean')
 # Suppress sub-1-CPM mass only if needed for the global density ceiling.
 p=mean/mean.sum();p=np.where(p>=1e-6,p,0);p/=p.sum()
 rows=[]
 for _ in range(n):
  total=int(rng.choice(libs));w=rng.gamma(theta,1/theta)*p
  if w.sum()<=0:w=p.copy()
  rows.append(sp.csr_matrix(rng.multinomial(total,w/w.sum()).astype(np.int32)[None,:]))
 return sp.vstack(rows,format='csr')

class CSRWriter:
 def __init__(self,path,genes):
  self.path=path;self.genes=genes;self.n=0;self.nnz=0;self.obs=[]
  self.f=h5py.File(path,'w');self.f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
  self.g=self.f.create_group('X');self.g.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0'})
  for key,dtype in [('data','float32'),('indices','int32'),('indptr','int64')]:
   self.g.create_dataset(key,shape=(1 if key=='indptr' else 0,),maxshape=(None,),dtype=dtype,chunks=(4096 if key=='indptr' else 262144,),compression='gzip',compression_opts=1)
  self.g['indptr'][0]=0
 def append(self,x,context,target):
  x=x.tocsr();x.eliminate_zeros();check_counts(x)
  if shutil.disk_usage('/').free-2*(x.data.nbytes+x.indices.nbytes)<60*1024**3:raise RuntimeError('60 GiB disk reserve')
  if np.max(np.asarray(x.sum(axis=1)))>1e6:raise ValueError('Per-cell count cap')
  if self.nnz+x.nnz>4_750_000_000:raise ValueError('Density cap')
  for k,v in [('data',x.data),('indices',x.indices),('indptr',x.indptr[1:]+self.nnz)]:
   d=self.g[k];start=len(d);d.resize((start+len(v),));d[start:]=v
  self.nnz+=x.nnz;self.n+=x.shape[0]
  self.obs.extend([(context,target)]*x.shape[0])
 def close(self):
  self.g.attrs['shape']=(self.n,len(self.genes))
  obs=pd.DataFrame(self.obs,columns=['context','target_gene'],index=[f'pred_{i}' for i in range(self.n)])
  for k in obs:obs[k]=pd.Categorical(obs[k])
  ad.io.write_elem(self.f,'obs',obs);ad.io.write_elem(self.f,'var',pd.DataFrame(index=self.genes))
  for k in ('uns','obsm','varm','obsp','varp','layers'):ad.io.write_elem(self.f,k,{})
  self.f.close()

@functools.lru_cache(maxsize=1)
def annotation_aliases():
 path=ROOT/'artifacts/gene_aliases.json'
 return json.loads(path.read_text()) if path.exists() else {}

def transfer_effect(sources,target,genes,baseline,variant,gene_ids=None):
 # Each held-out cell line is excluded by the caller; expression adaptation uses controls only.
 avail=[];weights=[];used_families=set()
 aliases=annotation_aliases()
 if gene_ids is None:gene_ids=np.array([aliases.get(g,'') for g in genes])
 target_ids=np.flatnonzero(genes==target)
 target_id=gene_ids[target_ids[0]] if len(target_ids) else aliases.get(target,'')
 ordered_sources=sorted(sources,key=lambda s: not s.get('is_pseudobulk',False)) if variant.startswith('bulk_') else sources
 for s in ordered_sources:
  if target in s['targets']:target_row=list(s['targets']).index(target)
  elif target_id and 'target_ids' in s and np.any(s['target_ids']==target_id):target_row=int(np.flatnonzero(s['target_ids']==target_id)[0])
  else:continue
  if 'n_cells' in s and s['n_cells'][target_row]<30:continue
  family=s.get('family',s.get('name',''))
  if family and family in used_families:continue
  gi=pd.Index(s['genes']).get_indexer(genes)
  if 'gene_ids' in s:
   id_to_column={gid:j for j,gid in enumerate(s['gene_ids']) if gid}
   for j in np.flatnonzero(gi<0):
    if gene_ids[j]:gi[j]=id_to_column.get(gene_ids[j],-1)
  ok=gi>=0
  if not np.any(ok):continue
  if family:used_families.add(family)
  effect=np.zeros(len(genes));effect[ok]=s['lfc'][target_row,gi[ok]]
  sourcebase=np.zeros(len(genes));sourcebase[ok]=s['baseline'][gi[ok]]
  if variant.startswith('adapt'):
   response=np.clip(np.sqrt((baseline+.5)/(sourcebase+.5)),.25,2)
   target_index=np.flatnonzero(genes==target)
   if not len(target_index) and target_id:target_index=np.flatnonzero(gene_ids==target_id)
   if len(target_index):
    j=target_index[0];dest_activity=baseline[j]/(baseline[j]+10);source_activity=sourcebase[j]/(sourcebase[j]+10)
    target_scale=np.clip(np.sqrt(dest_activity/max(source_activity,1e-6)),0,1.5)
   else:target_scale=1.
   reliability=np.sqrt(sourcebase/(sourcebase+5))
   effect*=response*target_scale*reliability
  valid=ok & (baseline>1) & (sourcebase>1)
  sim=np.corrcoef(np.log1p(baseline[valid]),np.log1p(sourcebase[valid]))[0,1] if valid.sum()>10 else 0
  avail.append(effect);weights.append(np.exp(4*np.nan_to_num(sim)))
 if not avail:return np.zeros(len(genes))
 if variant=='nearest':return avail[int(np.argmax(weights))]
 e=np.average(avail,axis=0,weights=weights if variant.startswith('adapt') else None)
 shrink={'transfer':1.,'shrink':.5,'adapt':.5,'adapt_strong':1.,'bulk_shrink':.5}.get(variant,1.)
 return np.clip(e*shrink,-3,3)

def predict_mean(summary,target,sources,variant):
 variant=variant.removesuffix('+empirical')
 mu=summary['mean'].copy();genes=summary['genes'];base=mu/mu.sum()*1e6
 if variant not in ('mean','knockdown'):
  e=transfer_effect(sources,target,genes,base,variant,summary.get('gene_ids'));mu*=np.exp2(e)
 target_index=np.flatnonzero(genes==target)
 if not len(target_index):
  ids=summary.get('gene_ids',np.array([annotation_aliases().get(g,'') for g in genes]));tid=annotation_aliases().get(target,'')
  if tid:target_index=np.flatnonzero(ids==tid)
 if variant!='mean' and len(target_index):mu[target_index[0]]=summary['mean'][target_index[0]]*.1
 return mu

def load_sources(directory):
 result=[]
 for p in sorted(pathlib.Path(directory).glob('*.npz')):
  if p.stem.endswith(('.tmp','.repack','.compact')):continue
  s=dict(np.load(p,allow_pickle=False))
  if 'lfc_scale' in s:s['lfc']=s['lfc'].astype(np.float32)*float(s['lfc_scale'])
  if s['lfc'].shape!=(len(s['targets']),len(s['genes'])):raise ValueError(f'{p}: signature axes disagree')
  if len(set(s['genes']))!=len(s['genes']):raise ValueError(f'{p}: duplicate source genes')
  aliases=annotation_aliases()
  s['gene_ids']=np.array([aliases.get(g,'') for g in s['genes']],dtype=str)
  s['target_ids']=np.array([aliases.get(g,'') for g in s['targets']],dtype=str)
  s['name']=p.stem;s['family']=str(s.get('cell_line',p.stem));result.append(s)
 return result

def build(args):
 d=pathlib.Path(args.data);out=pathlib.Path(args.output);out.parent.mkdir(exist_ok=True,parents=True)
 genes=csv_labels(d/'gene_names.csv');targets=csv_labels(d/'pert_counts.csv')
 if not args.synthetic and (len(genes)!=18533 or len(targets)!=300):raise ValueError('Unexpected official dimensions')
 phase=args.contexts.split(',')
 if set(phase) not in ({'A','B','C'},{'D','E','F'}):raise ValueError('Specify one exact phase: A,B,C or D,E,F')
 paths=[d/f'context_{c}.h5ad' for c in phase]
 if not all(p.exists() for p in paths):raise ValueError('All three requested control files are required')
 summaries=[control_summary(p,genes,p.stem.split('_')[-1]) for p in paths]
 contexts=[s['context'] for s in summaries]
 if not args.synthetic and any(s['n_controls']!=18400 for s in summaries):raise ValueError('Expected 18,400 official controls per context')
 if set(contexts) not in ({'A','B','C'},{'D','E','F'}):raise ValueError(f'Invalid phase {contexts}')
 sources=load_sources(args.signatures);rng=np.random.default_rng(args.seed)
 means=None
 if args.means:
  means=dict(np.load(args.means,allow_pickle=False))
  if not np.array_equal(means['genes'],genes):raise ValueError('Mean archive gene order mismatch')
  if len(set(means['contexts']))!=3 or set(means['contexts'])!=set(contexts):raise ValueError('Mean archive contexts disagree with phase')
  if len(set(means['targets']))!=len(targets) or set(means['targets'])!=set(targets):raise ValueError('Mean archive targets disagree with official list')
  if means['means'].shape!=(3,len(targets),len(genes)):raise ValueError('Mean archive dimensions disagree')
 writer=CSRWriter(out,genes)
 try:
  for c,s in enumerate(summaries):
   control_x=None
   if args.variant.endswith('+empirical'):
    control=ad.read_h5ad(s['path']);names=control.var_names.astype(str).to_numpy()
    if not set(genes).issubset(set(names)):names=gene_symbols(control)
    col=pd.Index(names).get_indexer(genes);control_x=control.X if np.array_equal(col,np.arange(len(genes))) else control.X[:,col]
    del control
   for i,t in enumerate(targets):
    if means is None:mu=predict_mean(s,t,sources,args.variant)
    else:
     ci=list(means['contexts']).index(s['context']);ti=list(means['targets']).index(t);mu=means['means'][ci,ti]
    if control_x is None:x=sample_counts(mu,s['libs'],s['theta'],400,rng)
    else:
     from samplers import empirical_counts
     x=empirical_counts(mu,control_x,s['mean'],s['theta'],400,rng)
    writer.append(x,s['context'],t)
    if i%20==0:print(s['context'],i,'/',len(targets),'nnz',writer.nnz,flush=True)
   del control_x
 finally:writer.close()
 manifest={'variant':args.variant,'seed':args.seed,'shape':[writer.n,len(genes)],'nnz':writer.nnz,'sources':[s['name'] for s in sources], 'pipeline_sha256':sha256(__file__), 'mean_archive_sha256':sha256(args.means) if args.means else None, 'controls':[{k:s[k] for k in ('context','path','sha256')} for s in summaries]}
 out.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2))
 # Release decoded training matrices before the CLI materializes the much larger prediction.
 sources.clear();summaries.clear()
 cmd=[(__import__('shutil').which('vcc') or 'vcc'),'prep',str(out),'-g',str(d/'gene_names.csv'),'--perts',str(d/'pert_counts.csv'),'--contexts',','.join(contexts)]
 if args.synthetic:cmd+=['--expected-gene-dim',str(len(genes))]
 env=os.environ.copy();env.update(TMPDIR=str(ROOT/'artifacts/tmp'),PYTHONPATH=str(ROOT/'prep_support'),VCC_RESOURCE_LIMITS='1')
 (ROOT/'artifacts/tmp').mkdir(exist_ok=True,parents=True)
 subprocess.run(cmd+['--dry-run'],check=True,env=env)
 if args.package:subprocess.run(cmd+['-o',str(out.with_name('prediction_val.vcc' if set(contexts)=={'A','B','C'} else 'prediction_final.vcc')),'--force'],check=True,env=env)


def synthetic(directory):
 d=pathlib.Path(directory);d.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(11)
 genes=np.array([f'GENE{i:04}' for i in range(160)]);targets=genes[:6]
 pd.Series(genes).to_csv(d/'gene_names.csv',index=False,header=False);pd.DataFrame({'target_gene':targets}).to_csv(d/'pert_counts.csv',index=False)
 for c in 'ABC':
  rate=rng.lognormal(0,1,len(genes));x=rng.poisson(rate,size=(240,len(genes)))
  a=ad.AnnData(sp.csr_matrix(x,dtype=np.int32),obs=pd.DataFrame({'context':c,'target_gene':NT},index=[f'{c}_{i}' for i in range(len(x))]),var=pd.DataFrame(index=genes));a.write_h5ad(d/f'context_{c}.h5ad')

if __name__=='__main__':
 p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
 b=sub.add_parser('build');b.add_argument('--data',default='../data');b.add_argument('--contexts',default='A,B,C');b.add_argument('--output',default='out/prediction.h5ad');b.add_argument('--signatures',default='artifacts/signatures');b.add_argument('--variant',default='mean',choices=BASE_VARIANTS+[v+'+empirical' for v in BASE_VARIANTS]);b.add_argument('--means');b.add_argument('--seed',type=int,default=42);b.add_argument('--synthetic',action='store_true');b.add_argument('--package',action='store_true')
 s=sub.add_parser('synthetic');s.add_argument('--data',default='artifacts/synthetic')
 a=p.parse_args()
 if a.command=='build':build(a)
 else:synthetic(a.data)
