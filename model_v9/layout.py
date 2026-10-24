"""Smallest scored-column edits on the exact v5b cells.

Freeze all non-DE columns, all original row depths and each gene's count in
each of forty depth blocks. Search marginal skew/support only where its
original rank statistic conflicts with the chosen call set. No new mean head.
"""
import runtime, sys, math
from pathlib import Path
import numpy as np, numba as nb
from scipy.special import ndtri
sys.path.append(str(Path(__file__).resolve().parents[1]/'model_v8'))
from layout_refine import balance

@nb.njit(cache=True)
def finish_variable_depth_balance(x,desired,protected):
 # v8's whole-column swaps are useful here, but its last-resort test assumes
 # equal row capacities. With original v5b depths an overfull row can have
 # fewer counts in EVERY gene than an underfull row. A positive-count move
 # still repairs that feasible transport problem and preserves column sums.
 load=x.sum(1);cursor=0;steps=0;n,g=x.shape
 while True:
  delta=load-desired;a=np.argmax(delta);b=np.argmin(delta)
  if delta[a]<=0:break
  limit=min(delta[a],-delta[b]);found=False
  for protect in [False,True]:
   for off in range(g):
    j=(cursor+off)%g
    if protected[j]!=protect or x[a,j]<=0:continue
    move=min(limit,x[a,j]);x[a,j]-=move;x[b,j]+=move;load[a]-=move;load[b]+=move;cursor=(j+1)%g;found=True;break
   if found:break
  if not found:break
  steps+=1
 return load-desired,steps

@nb.njit(cache=True)
def properties(v,d,T,S,cdf,scale,n,tie0,nzero):
 u=0.;h1=0.;h2=0.;zeros=0;mean=0.;bulk=math.log1p(50000*T/S)
 for i in range(len(v)):
  count=v[i]
  if count==0:u+=nzero*.5;zeros+=1
  else:
   pos=math.log1p(count/(d[i]/1e6))/scale*(len(cdf)-1)
   if pos>=len(cdf)-1:u+=n
   else:
    k=int(pos);f=pos-k;u+=cdf[k]*(1-f)+cdf[k+1]*f
  h=math.log1p(50000*(T-count)/(S-d[i]))-bulk;h1+=h;h2+=h*h
  mean+=count/d[i]*1e6/len(v)
 extra=(nzero+zeros)**3-(nzero+zeros)-(nzero**3-nzero)
 variance=len(v)*n/12*((len(v)+n+1)-(tie0+extra)/((len(v)+n)*(len(v)+n-1)))
 z=(u-.5*len(v)*n)/math.sqrt(max(variance,1e-20))
 jk=(len(v)-1)/len(v)*max(0.,h2-h1*h1/len(v))
 return z,jk,mean

@nb.njit(cache=True)
def candidate(original,order,blocktot,frac,power,block):
 out=np.zeros(len(original),np.int32)
 for b in range(len(blocktot)):
  a=b*block;end=min(a+block,len(original));n=end-a;t=blocktot[b]
  if not t:continue
  # Stagger the one-count support decision across blocks. This gives roughly
  # one-cell resolution globally instead of a twenty-five-cell staircase.
  support=min(t,n,max(1,int(frac*n+(b+.5)/len(blocktot))))
  weights=np.empty(support);ws=0.
  for r in range(support):
   row=order[end-support+r];w=(original[row]+.25)**power;weights[r]=w;ws+=w
  remainder=t-support;used=0;residual=np.empty(support)
  for r in range(support):
   value=weights[r]/ws*remainder;iv=int(math.floor(value));out[order[end-support+r]]=iv+1;used+=iv+1;residual[r]=value-iv
  left=t-used;cum=0.;next_at=.5
  for r in range(support):
   cum+=residual[r]
   if left>0 and cum>=next_at:
    out[order[end-support+r]]+=1;left-=1;next_at+=1.
  if left:out[order[end-1]]+=left
 return out

@nb.njit(cache=True)
def scored_candidate(orig,order,bt,frac,power,block,depth,T,S,cdf,scale,n,tie,nzero,goal,target,jk0,mean0,weight_mode=0):
 if weight_mode:v=candidate(depth,order,bt,frac,power,block)
 else:v=candidate(orig,order,bt,frac,power,block)
 z,jk,mean=properties(v,depth,T,S,cdf,scale,n,tie,nzero)
 err=max(0.,abs(abs(z)-abs(target))-.10) if goal else max(0.,abs(z)-1.75)
 change=np.abs(v.astype(np.float64)-orig).sum()/max(T,1)
 cost=abs(jk-jk0)/max(jk0,1e-5)+.05*change+10*abs(mean-mean0)/max(mean0,1.)
 return v,z,err,cost

@nb.njit(cache=True)
def columns(original,depth,keep,goals,cdf,scales,n,ties,nzero,block):
 nrow,g=original.shape;out=original.copy();modified=np.zeros(g,np.bool_);fallback=0
 S=depth.sum()
 powers=np.array([-102.,-103.,-100.,-101.,-2.,-1.,-.5,0.,.25,.5,1.,2.,4.]) if block==400 else np.array([0.,.25,.5,1.,2.,4.])
 for j in range(g):
  if not keep[j]:continue
  orig=original[:,j];T=orig.sum()
  if not T:continue
  z0,jk0,mean0=properties(orig,depth,T,S,cdf[j],scales[j],n,ties[j],nzero[j])
  goal=goals[j]
  # Keep every already-neutral non-call exactly, including its expression noise.
  if goal==0 and abs(z0)<=1.75:continue
  if goal!=0 and abs(abs(z0)-abs(goal))<=.10:continue
  target=goal if goal else math.copysign(1.65,z0)
  order=np.empty(nrow,np.int64);bt=np.empty((nrow+block-1)//block,np.int64)
  for b in range(len(bt)):
   a=b*block;end=min(a+block,nrow)
   order[a:end]=a+np.argsort(orig[a:end]/depth[a:end]);bt[b]=orig[a:end].sum()
  best=1e300;besterr=1e300;bv=orig.copy();bz=z0
  original_order=order.copy()
  for code in powers:
   power=code;weight_mode=0
   if code==-102.:
    order[:]=np.argsort(-depth);power=1.;weight_mode=1
   elif code==-103.:
    order[:]=np.argsort(depth);power=1.;weight_mode=1
   elif code==-100.:
    order[:]=np.argsort(-depth);power=0.
   elif code==-101.:
    order[:]=np.argsort(depth);power=0.
   else:order[:]=original_order
   previous_q=0.;previous_z=0.
   # U(q) can peak below full support; inspect both crossing branches.
   for qi,frac in enumerate(np.array([0.,.2,.4,.6,.8,.9,1.])):
    v,z,err,cost=scored_candidate(orig,order,bt,frac,power,block,depth,T,S,cdf[j],scales[j],n,ties[j],nzero[j],goal,target,jk0,mean0,weight_mode)
    if err<besterr-1e-10 or (abs(err-besterr)<=1e-10 and cost<best):besterr=err;best=cost;bv=v.copy();bz=z
    for side in range(2):
     root=target if side==0 else -target
     # The scored direction is arithmetic LFC, not the signed MWU score.
     if qi>0 and (previous_z-root)*(z-root)<=0.:
      low=previous_q;high=frac;zlow=previous_z
      for _ in range(7):
       mid=(low+high)*.5
       mv,mz,me,mc=scored_candidate(orig,order,bt,mid,power,block,depth,T,S,cdf[j],scales[j],n,ties[j],nzero[j],goal,target,jk0,mean0,weight_mode)
       if me<besterr-1e-10 or (abs(me-besterr)<=1e-10 and mc<best):besterr=me;best=mc;bv=mv.copy();bz=mz
       if (zlow-root)*(mz-root)<=0.:high=mid
       else:low=mid;zlow=mz
    previous_q=frac;previous_z=z
  out[:,j]=bv;modified[j]=True
  if (goal and abs(abs(bz)-abs(target))>.35) or (not goal and abs(bz)>1.95):fallback+=1
 return out,modified,fallback

def selection(original,ref,confidence,nconf,own,strict=None):
 depth=original.sum(1,dtype=np.int64)
 keep=ref.keep.copy();keep[own]=False
 arith=(original/(depth[:,None]/1e6)).mean(0)
 T=original.sum(0,dtype=np.int64);direction=np.sign(arith-ref.mean)
 bulkdir=np.sign(T/T.sum()*1e6-ref.bulk)
 eligible=keep&(confidence>=1e6)&(direction==bulkdir)
 if strict is not None:eligible&=strict
 order=np.flatnonzero(eligible);order=order[np.argsort(-confidence[order],kind='stable')]
 called=order[:nconf];goals=np.zeros(len(T))
 tail=max(4.2,-ndtri(.5*.005*max(1,len(called))/keep.sum()))
 goals[called]=direction[called]*np.linspace(max(6.5,tail+1.5),tail,len(called))
 return keep,called,goals

def repair(x,original,modified,goals,block):
 depth=original.sum(1,dtype=np.int64)
 for a in range(0,400,block):
  end=min(400,a+block)
  desired=depth[a:end]-(original[a:end][:,~modified].sum(1,dtype=np.int64))
  sub=x[a:end][:,modified].copy()
  if sub.shape[1]:
   residual,steps=balance(sub,np.ones(sub.shape[1],bool),desired,goals[modified]!=0)
   if residual.any():residual,extra=finish_variable_depth_balance(sub,desired,goals[modified]!=0)
   assert not residual.any(),('row repair',residual.min(),residual.max())
   x[a:end,modified]=sub
 return x

def emit(original,ref,confidence,nconf,own,strict=None,block=10):
 original=np.asarray(original,np.int32);depth=original.sum(1,dtype=np.int64);T=original.sum(0,dtype=np.int64)
 keep,called,goals=selection(original,ref,confidence,nconf,own,strict)
 x,modified,fallback=columns(np.asfortranarray(original),depth,keep,goals,ref.cdf,ref.cdf_scale,ref.n,ref.ties,ref.n-np.diff(ref.ptr),block)
 # Only changed DE columns can move. Each block's column totals stay exact.
 x=repair(x,original,modified,goals,block)
 assert np.array_equal(x.sum(0,dtype=np.int64),T)
 assert np.array_equal(x.sum(1,dtype=np.int64),depth)
 assert np.array_equal(x[:,~keep],original[:,~keep])
 return x,called,dict(modified_columns=int(modified.sum()),infeasible_before_balance=fallback,block=block)

def restrict(original,a,ref,confidence,acalled,strict,own,block=10):
 """Only suppress calls removed by the frozen stricter source-evidence rule."""
 allowed=np.asarray([j for j in acalled if strict[j]],dtype=np.int64)
 removed=np.setdiff1d(acalled,allowed);keep=np.zeros(original.shape[1],bool);keep[removed]=True
 goals=np.zeros(original.shape[1]);depth=original.sum(1,dtype=np.int64)
 values,changed,fallback=columns(np.asfortranarray(original),depth,keep,goals,ref.cdf,ref.cdf_scale,ref.n,ref.ties,ref.n-np.diff(ref.ptr),block)
 x=a.copy();x[:,keep]=values[:,keep]
 # Freeze retained calls exactly at their v9a counts. Compensate row loads
 # using removed calls and already edited, uncalled DE columns only.
 modified=np.any(a!=original,axis=0)|keep;modified[allowed]=False
 modified[~ref.keep]=False;modified[own]=False
 x=repair(x,a,modified,goals,block)
 # BH is checked on the actual integer cells. A few columns cannot become
 # neutral inside ten-row blocks; relax only those columns' block constraint.
 rounds=[]
 for wider in [40,400,400,400,400,400]:
  _,_,q=ref.test(x);bad=(q<.05);bad[allowed]=False;bad[own]=False
  if not bad.any():break
  values,changed,extra=columns(np.asfortranarray(original),depth,bad,goals,ref.cdf,ref.cdf_scale,ref.n,ref.ties,ref.n-np.diff(ref.ptr),wider)
  x[:,bad]=values[:,bad]
  x=repair(x,a,modified|bad,goals,wider)
  rounds.append(dict(block=wider,columns=int(bad.sum()),infeasible=extra))
 _,_,finalq=ref.test(x);outside=finalq<.05;outside[allowed]=False;outside[own]=False
 assert np.array_equal(x.sum(0,dtype=np.int64),original.sum(0,dtype=np.int64))
 assert np.array_equal(x.sum(1,dtype=np.int64),depth)
 assert np.array_equal(x[:,~ref.keep],original[:,~ref.keep])
 assert np.array_equal(x[:,allowed],a[:,allowed])
 return x,allowed,dict(modified_columns=int(modified.sum()),removed_selected=len(removed),infeasible_before_balance=fallback,kept_call_counts_identical_to_v9a=True,block=block,suppression_rounds=rounds,residual_unsupported_calls=int(outside.sum()))
