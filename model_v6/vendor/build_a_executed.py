"""Stream selected models, independently certify, retain one regenerable H5AD at a time."""
import runtime,argparse,time,gc,os
os.environ['VCC_ASH_WORKERS']='3'
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from model import template,heads,emit,atlas
from verify_counts import moment_certificate,audit_h5
rt=runtime.ROOT
class Writer:
 def __init__(self,path,genes):
  self.n=0;self.nnz=0;self.obs=[];self.f=h5py.File(path,'w');self.genes=genes;self.g=self.f.create_group('X');self.f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'});self.g.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0'})
  for key,dt in [('data','int32'),('indices','int32'),('indptr','int64')]:self.g.create_dataset(key,shape=(1 if key=='indptr' else 0,),maxshape=(None,),dtype=dt,chunks=(4096 if key=='indptr' else 262144,),compression='gzip',compression_opts=1,shuffle=True)
  self.g['indptr'][0]=0
 def append(self,x,c,t):
  runtime.guard();x=sparse.csr_matrix(x);assert self.nnz+x.nnz<=4750000000
  for key,value in [('data',x.data.astype(np.int32)),('indices',x.indices.astype(np.int32)),('indptr',x.indptr[1:].astype(np.int64)+self.nnz)]:
   ds=self.g[key];start=len(ds);ds.resize((start+len(value),));ds[start:]=value
  self.nnz+=x.nnz;self.n+=x.shape[0];self.obs.extend([(c,t)]*x.shape[0])
 def close(self):
  self.g.attrs['shape']=(self.n,len(self.genes));obs=pd.DataFrame(self.obs,columns=['context','target_gene'],index=[f'pred_{i}' for i in range(self.n)])
  for k in obs:obs[k]=pd.Categorical(obs[k])
  ad.io.write_elem(self.f,'obs',obs);ad.io.write_elem(self.f,'var',pd.DataFrame(index=pd.Index(self.genes)))
  for k in ['uns','obsm','varm','obsp','varp','layers']:ad.io.write_elem(self.f,k,{})
  self.f.close()

_WORKER_STATE=None
def init_worker(q,d,variant):
 global _WORKER_STATE
 _WORKER_STATE=(q,d,variant)
def emit_worker(item):
 i,t,p,b,seed=item;q,d,variant=_WORKER_STATE;x=emit(q,p,b,d,seed,variant);cert=moment_certificate(x,d,p,b);return x,cert

def main(args):
 genes=pd.read_csv(rt.parent/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(rt.parent/'data/pert_counts.csv').target_gene.to_numpy(str);assert (len(genes),len(targets))==(18533,300)
 name='prediction_v4'+args.arm;path=rt/'out'/f'{name}.h5ad';assert not path.exists();writer=Writer(path,genes);mom=[];controls=[];stats={};start=time.time()
 try:
  for ci,c in enumerate('ABC'):
   control_path=rt.parent/'data'/f'context_{c}.h5ad';a=ad.read_h5ad(control_path);assert a.n_obs==18400;controls.append(dict(context=c,file=control_path.name,sha256=runtime.sha(control_path)))
   q,d,mean,bulk=template(a,genes,seed=20260910+ci);del a;gc.collect()
   cache=rt/'artifacts'/f'{name}_{c}_heads.npz'
   if cache.exists():
    data=np.load(cache);p,b=data['p'],data['b']
   else:
    (p,b),coverage=heads(targets,genes,mean,bulk,args.variant,stats=stats);np.savez_compressed(cache,p=p,b=b);del coverage
   import multiprocessing
   from concurrent.futures import ProcessPoolExecutor
   # Three single-core numeric workers plus the single-core HDF5 writer stay within four cores.
   with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn'),initializer=init_worker,initargs=(q,d,args.variant)) as executor:
    pending={};next_submit=0
    for i,t in enumerate(targets):
     while next_submit<min(len(targets),i+5):
      k=next_submit;pending[k]=executor.submit(emit_worker,(k,targets[k],p[k],b[k],atlas.seed_for(c+':'+targets[k],20260910)));next_submit+=1
     x,cert=pending.pop(i).result();assert cert['mean_L1_upper']<.002 if not args.variant.endswith('_mean') else True
     writer.append(x,c,t);mom.append(dict(context=c,target=t,**cert));del x
     if i%10==0:
      pd.DataFrame(mom).to_csv(rt/'artifacts'/f'{name}_moments.csv',index=False);print(name,c,i,'nnz',writer.nnz,'seconds',round(time.time()-start),'free GiB',round(runtime.shutil.disk_usage('/').free/1024**3,2),flush=True)
   del p,b,q;gc.collect()
 finally:writer.close()
 pd.DataFrame(mom).to_csv(rt/'artifacts'/f'{name}_moments.csv',index=False)
 runtime.save(rt/'out'/f'{name}.manifest.json',dict(variant=args.variant,shape=[writer.n,len(genes)],nnz=writer.nnz,controls=controls,seed=20260910,adaptive_stats=stats,model_sha256=runtime.sha(rt/'model.py'),builder_sha256=runtime.sha(__file__),status='Biological quality UNVERIFIED; A selected using public H1/HEK attribution; B requested full-integration ablation'))
 verified=audit_h5(path,genes,targets);runtime.save(rt/'artifacts'/f'{name}_structure.json',verified);print('VERIFIED',name,verified['nnz'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('arm',choices=['a','b']);p.add_argument('--variant',required=True);main(p.parse_args())
