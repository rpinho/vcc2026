"""Label-free v7 dual-moment builds. No local biological fitness or submissions.
Credits: AtlasShift MIT, v5/v6 control/sources and budgets, Chen Zhang PIE padding.
"""
import runtime
import argparse,gc,time,json,os
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
import model,v5
from feasible_bulk import project_bulk
from verify_counts import audit_h5,moment_certificate
rt=runtime.ROOT

class Writer:
 def __init__(self,path,genes):
  self.n=0;self.nnz=0;self.obs=[];self.genes=genes;self.f=h5py.File(path,'w')
  self.f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
  self.g=self.f.create_group('X');self.g.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0'})
  for k,dt in [('data','int32'),('indices','int32'),('indptr','int64')]:self.g.create_dataset(k,shape=(1 if k=='indptr' else 0,),maxshape=(None,),dtype=dt,chunks=(4096 if k=='indptr' else 262144,),compression='gzip',compression_opts=1,shuffle=True)
  self.g['indptr'][0]=0
 def append(self,x,c,t):
  runtime.guard();x=sparse.csr_matrix(x)
  for k,a in [('data',x.data.astype(np.int32)),('indices',x.indices.astype(np.int32)),('indptr',x.indptr[1:].astype(np.int64)+self.nnz)]:
   z=self.g[k];j=len(z);z.resize((j+len(a),));z[j:]=a
  self.n+=x.shape[0];self.nnz+=x.nnz;self.obs.extend([(c,t)]*x.shape[0])
 def close(self):
  self.g.attrs['shape']=(self.n,len(self.genes));o=pd.DataFrame(self.obs,columns=['context','target_gene'],index=[f'pred_{i}' for i in range(self.n)])
  for k in o:o[k]=pd.Categorical(o[k])
  ad.io.write_elem(self.f,'obs',o);ad.io.write_elem(self.f,'var',pd.DataFrame(index=pd.Index(self.genes)))
  for k in ['uns','obsm','varm','obsp','varp','layers']:ad.io.write_elem(self.f,k,{})
  self.f.close()

def pie_padding(mask,delta,targets,genes,path,k_add=300):
 # Same target-centering and eligibility rule as model_v6/pipeline.py; no score/label access.
 with np.load(path) as z:
  ri=pd.Index(z['targets'].astype(str)).get_indexer(targets);ci=pd.Index(z['genes'].astype(str)).get_indexer(genes)
  rr=np.flatnonzero(ri>=0);cc=np.flatnonzero(ci>=0);P=np.full(mask.shape,np.nan,np.float32)
  P[np.ix_(rr,cc)]=z['p_de'].astype(np.float32)[np.ix_(ri[rr],ci[cc])]
 mu=np.divide(np.nansum(P,axis=0),np.isfinite(P).sum(0),out=np.zeros(len(genes),np.float32),where=np.isfinite(P).sum(0)>0);P-=mu
 out=mask.copy()
 for i in rr:
  ok=np.flatnonzero(~mask[i]&np.isfinite(P[i])&(np.abs(delta[i])>1e-9));top=ok[np.argsort(-P[i,ok],kind='stable')[:k_add]];out[i,top]=True
 return out,dict(pie_file=str(Path(path).relative_to(rt.parent)),pie_sha256=runtime.sha(path),added_mean=float((out&~mask).sum(1).mean()),pie_targets=len(rr),pie_genes=len(cc))

class Context:
 def __init__(self,a,targets,genes,label,ci,arm,pie_path=None):
  self.targets=targets;self.genes=genes;self.label=label;self.arm=arm;self.ci=ci
  print(label,arm,"template",flush=True)
  self.q4,self.d,self.mean,self.bulk=model.template(a,genes,pool_k=4,seed=20260910+ci)
  print(label,arm,"heads",flush=True)
  variant='anchor_dual' if arm=='a' else 'anchor_single_mean'
  (self.p,self.b),_=model.heads(targets,genes,self.mean,self.bulk,variant,de_amplitude=.9)
  print(label,arm,"control axes",flush=True)
  V,s=v5.control_axes(a.X,k=50,seed=ci)
  self.Lc=v5.logbulk(self.bulk)[None,:];self.L=self.Lc+v5.deflate(v5.logbulk(self.b)-self.Lc,V)
  self.q1=None;self.mask=None
  self.info=dict(context=label,arm=arm,control_cells=a.n_obs,targets=len(targets),de_amplitude=.9,bulk_head=('native_bulk_delta' if arm=='a' else 'log2fc'),bulk_deflation_k=50,pool_called=4,pool_other=4 if arm=='a' else 1)
  if arm=='b':
   self.q1,_,_,_=model.template(a,genes,pool_k=1,seed=20260910+ci)
   names=model.load_banks();src=[model.aligned_bank(n,targets,genes,'log2fc') for n in names];cd=model.cd4_family(targets,genes,'log2fc')
   if cd is not None:src.append(cd)
   delta=self.p-self.mean[None,:]
   self.mask,net=v5.de_call_mask(src,delta,k=len(genes))
   self.info['support_mask_mean']=float(self.mask.sum(1).mean());del src,cd;gc.collect()
   if pie_path is not None:
    self.mask,pinfo=pie_padding(self.mask,delta,targets,genes,pie_path);self.info.update(pinfo)
   self.info['low_variance_genes_mean']=float(self.mask.sum(1).mean())
  self.w=v5._row_weights(len(targets),len(genes),v5.excl_cols(targets,genes))
  # Controls-only cap proxy. The server cap uses hidden treated JK, which remains unknown.
  rng=np.random.default_rng(73591+ci);X=sparse.csr_matrix(a.X)
  print(label,arm,"control JK",flush=True)
  self.ctrl_jk=np.array([self.jk(X[rng.choice(a.n_obs,400,replace=False)]) for _ in range(8)])
  self.cap_jk=float(self.ctrl_jk.mean());self.info['control400_jk_samples']=self.ctrl_jk.tolist();self.info['control400_claim']=float(self.cap_jk*self.w.sum())
  self.probe=np.sort(np.random.default_rng(ci).choice(len(targets),min(8,len(targets)),replace=False));self.alpha=1.
 def jk(self,x):return float(v5.jackknife_correction(sparse.csr_matrix(x),np.zeros(x.shape[0],np.intp),1,v5.TS,chunk=32)[0])
 def emit(self,i,b):
  runtime.guard();q=self.q4
  if self.q1 is not None:
   q=self.q1.copy();q[:,self.mask[i]]=self.q4[:,self.mask[i]];q/=q.sum(1,keepdims=True)
  x=model.atlas.dual_moment_counts(q,self.p[i],b,depths=self.d,seed=model.atlas.seed_for(self.label+':'+self.targets[i],20260910),iterations=400)
  cert=moment_certificate(x,self.d,self.p[i],b)
  cert.update(float_fit_steps=int(model.atlas.dual_moment_counts.last_fit[0]),float_mean_error=float(model.atlas.dual_moment_counts.last_fit[1]),float_bulk_error=float(model.atlas.dual_moment_counts.last_fit[2]))
  if cert['mean_L1_upper']>.015 or cert['bulk_L1_upper']>.0015:
   raise RuntimeError(f"dual moment violation {self.label}/{self.targets[i]} {cert}")
  return x,cert
 def projected(self,alpha):
  raw=v5.from_logbulk(v5.rescale(self.L,alpha,self.Lc));z=self.d/self.d.mean()
  return project_bulk(raw,self.p,float(z.min()),float(z.max()))
 def solve_alpha(self,target):
  lo=0.;hi=8.
  def budget(alpha):return v5.budget(self.projected(alpha)[0],self.targets,self.genes)
  if budget(hi)<target:return hi,budget(hi)
  assert budget(lo)<=target,'feasible DE moments already exceed requested bulk budget'
  for _ in range(38):
   mid=.5*(lo+hi)
   if budget(mid)<target:lo=mid
   else:hi=mid
  alpha=.5*(lo+hi);return alpha,budget(alpha)
 def calibrate(self):
  hist=[]
  for iteration in range(5):
   bf,_=self.projected(self.alpha);jks=[]
   for i in self.probe:
    t0=time.time();jks.append(self.jk(self.emit(i,bf[i])[0]));print(self.label,self.arm,'probe',int(i),'seconds',round(time.time()-t0,2),flush=True)
   jks=np.array(jks)
   claim=float(np.mean(np.minimum(jks,self.cap_jk))*self.w.sum());alpha,B=self.solve_alpha(.95*claim)
   hist.append(dict(iteration=iteration,alpha_measured=self.alpha,claim_capped_proxy=claim,pred_claim_probe=float(jks.mean()*self.w.sum()),alpha_next=alpha,B_next=B))
   print(self.label,self.arm,'calibration',hist[-1],flush=True)
   if abs(alpha/self.alpha-1)<.005:
    self.alpha=alpha;break
   self.alpha=alpha
  self.bf,proj=self.projected(self.alpha)
  assert proj.max()<=.03, 'source bulk exceeds public feasible-projection limit'
  self.info.update(feasible_projection_L1_max=float(proj.max()),feasible_projection_L1_mean=float(proj.mean()))
  self.info.update(alpha=self.alpha,calibration=hist,B_requested=v5.budget(self.bf,self.targets,self.genes),probe_indices=self.probe.tolist(),rho=.95)

def parser():
 p=argparse.ArgumentParser();p.add_argument('arm',choices=['a','b']);p.add_argument('--fixture',action='store_true');p.add_argument('--targets',type=Path,default=rt.parent/'data/pert_counts.csv');p.add_argument('--genes',type=Path,default=rt.parent/'data/gene_names.csv');p.add_argument('--contexts',default='A='+str(rt.parent/'data/context_A.h5ad')+',B='+str(rt.parent/'data/context_B.h5ad')+',C='+str(rt.parent/'data/context_C.h5ad'));p.add_argument('--banks',type=Path,default=model.BANKS);p.add_argument('--pie-dir',type=Path,default=rt.parent/'model_v6/pie/out');p.add_argument('--pie-prefix',default='ens4');return p

def main(args):
 model.BANKS=args.banks;start=time.time();genes=pd.read_csv(args.genes).gene_name.to_numpy(str);targets=pd.read_csv(args.targets).target_gene.to_numpy(str);contexts=dict(s.split('=',1) for s in args.contexts.split(','));name='v7'+args.arm
 path=rt/'out'/f'{name}.h5ad';assert args.fixture or not path.exists(), str(path)+' already exists'
 writer=None if args.fixture else Writer(path,genes);allmom=[];recs=[]
 if args.fixture:contexts=dict(list(contexts.items())[:1])
 try:
  for ci,(c,fp) in enumerate(contexts.items()):
   runtime.guard();a=ad.read_h5ad(fp);assert np.array_equal(a.var_names,genes)
   pie=args.pie_dir/f'{args.pie_prefix}_{c}.npz' if args.arm=='b' else None
   ctx=Context(a,targets,genes,c,ci,args.arm,pie);del a;gc.collect();ctx.calibrate()
   P=[];jks=[];ids=ctx.probe[:4] if args.fixture else range(len(targets))
   for i in ids:
    x,cert=ctx.emit(i,ctx.bf[i]);jks.append(ctx.jk(x));P.append(x.sum(0,dtype=np.int64));allmom.append(dict(context=c,target=targets[i],**cert))
    if args.fixture:np.savez_compressed(rt/'artifacts'/f'fixture_{args.arm}_{i}.npz',x=x,depths=ctx.d,p=ctx.p[i],b=ctx.bf[i],mask=ctx.mask[i] if ctx.mask is not None else np.ones(len(genes),bool))
    else:writer.append(sparse.csr_matrix(x),c,targets[i])
    if i%20==0:print(name,c,i,'elapsed_s',round(time.time()-start),'nnz',writer.nnz if writer else None,flush=True)
   if not args.fixture:
    P=np.array(P);Pn=P/P.sum(1,keepdims=True);ctx.info.update(B_emitted=v5.budget(Pn,targets,genes),claim_pred_emitted=float(ctx.w@np.array(jks)),claim_control_capped_emitted=float(ctx.w@np.minimum(jks,ctx.cap_jk)),jk_range=[float(np.min(jks)),float(np.max(jks))]);ctx.info['B_over_pred_claim']=ctx.info['B_emitted']/ctx.info['claim_pred_emitted'];ctx.info['B_over_control_capped_claim']=ctx.info['B_emitted']/ctx.info['claim_control_capped_emitted']
   recs.append(ctx.info);runtime.save(rt/'artifacts'/f'{name}_progress.json',dict(contexts=recs,seconds=time.time()-start));print('CONTEXT DONE',name,json.dumps(ctx.info),flush=True)
   del ctx;gc.collect()
 finally:
  if writer is not None:writer.close()
 suffix='fixture' if args.fixture else 'build';runtime.save(rt/'artifacts'/f'{name}_{suffix}.json',dict(status='software checks only; official score UNVERIFIED',contexts=recs,moments=allmom,seconds=time.time()-start,shape=[writer.n,len(genes)] if writer else None,nnz=writer.nnz if writer else None,source_hashes={n:runtime.sha(rt/n) for n in ['build.py','model.py','v5.py','fast_fit.py','fast_round.py','feasible_bulk.py','vendor/atlas_model.py']}))
 if not args.fixture:
  cert=audit_h5(path,genes,targets,contexts=tuple(contexts));runtime.save(rt/'artifacts'/f'{name}_structure.json',cert);print('BUILD VERIFIED',name,cert['shape'],cert['nnz'],flush=True)
if __name__=='__main__':main(parser().parse_args())
