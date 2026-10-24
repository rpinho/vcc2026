"""Two depth strata: fixed bulk totals, explicit mean head, calibrated rank tests.
Retrieval supplies positive-cell ranks and positive abundance patterns. Histogram
selection controls U and leave-one-out variability; row swaps preserve histograms.
"""
import runtime,math,time
import numpy as np,numba as nb
from scipy import sparse
from scipy.special import ndtri
from cell_eval2.moments import jackknife_per_gene

@nb.njit(cache=True)
def histogram(total,k,q):
 out=np.zeros(len(q),np.int64)
 if total==0:return out
 if k==0:k=1
 w=q[-k:].copy()+1e-9;w=w/w.sum()*(total-k)
 v=np.floor(w).astype(np.int64)+1;left=total-v.sum()
 frac=w-np.floor(w);cum=0.;next_at=.5;remain=left
 for a in range(k):
  cum+=frac[a]
  if remain>0 and cum>=next_at:v[a]+=1;remain-=1;next_at+=1.
 if remain:v[-1]+=remain
 out[-k:]=v
 return out

@nb.njit(cache=True)
def table(total,depth,T,S,cdf,scale,n,nz):
 ans=np.empty((total+1,2));last=len(cdf)-1;L=math.log1p(50000*T/S)
 for count in range(total+1):
  if count==0:U=nz*.5
  else:
   pos=math.log1p(count/(depth/1e6))/scale*last
   if pos>=last:U=n
   else:
    k=int(pos);f=pos-k;U=cdf[k]*(1-f)+cdf[k+1]*f
  h=math.log1p(50000*(T-count)/(S-depth))-L;ans[count,0]=U;ans[count,1]=h
 return ans

@nb.njit(cache=True)
def properties(v,tab):
 U=0.;h1=0.;h2=0.
 for count in v:
  U+=tab[count,0];h=tab[count,1];h1+=h;h2+=h*h
 return U,h1,h2

@nb.njit(cache=True)
def make_columns(T,TL,dl,dh,ql,oh,ol,cdf,scales,n,nzero,ties,keep,goals,jktarget,step_low,step_high):
 g=len(T);nl=ql.shape[0];nh=400-nl;out=np.zeros((400,g),np.int32);fallback=0;S=T.sum()
 for j in range(g):
  if T[j]==0:continue
  if not keep[j]:continue
  a=TL[j];b=T[j]-a
  maxl=min(a,nl);maxh=min(b,nh)
  # Dense k grid at small support, stride 4 otherwise; include full support.
  kl0=np.unique(np.concatenate((np.arange(1,maxl+1,step_low),np.array([maxl])))) if maxl else np.array([0])
  kl=np.empty(4*len(kl0),np.int64)
  for ii in range(len(kl)):kl[ii]=kl0[ii%len(kl0)]
  kh=np.unique(np.concatenate((np.arange(1,maxh+1,step_high),np.array([maxh])))) if maxh else np.array([0])
  vl=np.zeros((len(kl),nl),np.int64);pl=np.zeros((len(kl),3));vh=np.zeros((len(kh),nh),np.int64);ph=np.zeros((len(kh),3))
  tabl=table(a,dl,T[j],S,cdf[j],scales[j],n,nzero[j]);tabh=table(b,dh,T[j],S,cdf[j],scales[j],n,nzero[j])
  allweights=np.empty((4,nl))
  for power_index,power in enumerate((0.,.25,1.,3.)):allweights[power_index]=(ql[:,j]+1e-6)**power
  for ii,k in enumerate(kl):
   weights=allweights[ii//len(kl0)]
   v=histogram(a,k,weights);vl[ii]=v;pl[ii]=properties(v,tabl)
  for ii,k in enumerate(kh):
   v=np.zeros(nh,np.int64)
   if k:
    v[-k:]=b//k
    for t in range(b%k):v[nh-1-t]+=1
   vh[ii]=v;ph[ii]=properties(v,tabh)
  vmax=np.empty(len(kh),np.int64)
  for jj in range(len(kh)):vmax[jj]=vh[jj].max()
  sigma=np.empty(401)
  for positive in range(401):
   zeros=400-positive;nz=nzero[j];extra=(nz+zeros)**3-(nz+zeros)-(nz**3-nz)
   var=400*n/12*((400+n+1)-(ties[j]+extra)/((400+n)*(400+n-1)));sigma[positive]=math.sqrt(max(var,1e-20))
  best=1e300;bi=0;bj=0;feasible=False
  for ii,k in enumerate(kl):
   if vl[ii].max()>dl:continue
   for jj,h in enumerate(kh):
    if vmax[jj]>dh:continue
    z=(pl[ii,0]+ph[jj,0]-200*n)/sigma[k+h]
    # Neutral p > .05 even before BH; reserve .11 z for interpolation/ties.
    err=max(abs(z)-1.85,0.) if goals[j]==0 else max(abs(z-goals[j])-.12,0.)
    h1=pl[ii,1]+ph[jj,1];h2=pl[ii,2]+ph[jj,2];jk=399/400*(h2-h1*h1/400)
    cost=err*1e5+abs(jk-jktarget[j])
    if cost<best:best=cost;bi=ii;bj=jj;feasible=err==0
  if not feasible:fallback+=1
  for row in range(nl):out[ol[row,j],j]=vl[bi,row]
  for row in range(nh):out[nl+oh[row,j],j]=vh[bj,row]
 return out,fallback

@nb.njit(cache=True)
def balance(x,keep,depth,protected):
 # Swapping counts within one gene preserves its marginal distribution exactly.
 n,g=x.shape;load=x.sum(1);cursor=0;steps=0
 while steps<3000000:
  delta=load-depth;a=np.argmax(delta);b=np.argmin(delta)
  if delta[a]<=0 or delta[b]>=0:break
  limit=min(delta[a],-delta[b]);found=False
  for offset in range(g):
   j=(cursor+offset)%g
   if not keep[j]:continue
   diff=x[a,j]-x[b,j]
   if diff>0 and diff<=limit:
    temp=x[a,j];x[a,j]=x[b,j];x[b,j]=temp;load[a]-=diff;load[b]+=diff;cursor=(j+1)%g;found=True;break
  if not found:
   # An integer residue may be smaller than every available whole-histogram
   # swap. Repair only an uncalled column; exact DE is recomputed later.
   for offset in range(g):
    j=(cursor+offset)%g
    if keep[j] and not protected[j] and x[a,j]>x[b,j]:
     move=min(limit,x[a,j]-x[b,j]);x[a,j]-=move;x[b,j]+=move;load[a]-=move;load[b]+=move;cursor=(j+1)%g;found=True;break
  if not found:
   # Last resort: a called-gene outlier can exceed every single deficit.
   # Split it across deficits. This keeps both integer moment constraints;
   # final independent DE gates include the resulting histogram change.
   for offset in range(g):
    j=(cursor+offset)%g
    if keep[j] and x[a,j]>x[b,j]:
     move=min(limit,x[a,j]-x[b,j]);x[a,j]-=move;x[b,j]+=move;load[a]-=move;load[b]+=move;cursor=(j+1)%g;found=True;break
  if not found:break
  steps+=1
 return load-depth,steps

def fill(x,cols,tot,desired):
 assert np.all(desired>=0) and int(desired.sum())==int(tot.sum())
 n=len(x);q=tot//n;r=tot%n;x[:,cols]=q[None,:]
 remaining=desired-q.sum()
 for idx in np.argsort(-r,kind='stable'):
  degree=int(r[idx])
  if degree:
   rows=np.argsort(-remaining,kind='stable')[:degree];x[rows,cols[idx]]+=1;remaining[rows]-=1
 if remaining.any():
  # Unscored genes provide a capacity-constrained integer transportation
  # fallback. Counts and both stratum totals remain exact.
  x[:,cols]=0;remaining=desired.copy()
  for idx,t in enumerate(tot):
   if not t:continue
   ideal=remaining.astype(float)*(int(t)/int(remaining.sum()));v=np.floor(ideal).astype(np.int64);left=int(t-v.sum())
   rows=np.argsort(-(ideal-v),kind='stable')[:left];v[rows]+=1
   assert np.all(v<=remaining);x[:,cols[idx]]=v;remaining-=v
 assert not remaining.any(),('filler row degrees',remaining.min(),remaining.max())

def emit(totals,reference,confidence,nconf,target_index=None,cap=1.,**kwargs):
 start_time=time.monotonic()
 ref=reference;T=totals.astype(np.int64);S=int(T.sum());nl=336;nh=64
 dl=max(1,int(.40*S/400));dh=(S-nl*dl)//nh
 depths=np.r_[np.full(nl,dl),np.full(nh,dh)];depths[nl:nl+(S-int(depths.sum()))]+=1
 if not hasattr(ref,'two_depth_cache'):
  ol=np.argsort(ref.base[:nl],axis=0,kind='stable');oh=np.argsort(ref.base[nl:],axis=0,kind='stable');ql=np.take_along_axis(ref.base[:nl],ol,axis=0)
  raw=np.rint(ref.base*ref.depth[:,None]/1e6).astype(np.int32)
  jkg,_=jackknife_per_gene(sparse.csr_matrix(raw),np.zeros(400,np.intp),1,50000.,chunk=64)
  ref.two_depth_cache=(np.asfortranarray(ql),ol,oh,np.maximum(jkg[0],0))
 ql,ol,oh,jkg=ref.two_depth_cache
 keep=ref.keep.copy()
 if target_index is not None:keep[target_index]=False
 bulk=T/S*1e6;sgn=np.sign(bulk-ref.bulk);sgn[sgn==0]=1
 lo=T/(400*dh)*1e6;hi=T/(400*dl)*1e6;margin=.001+1e-4*ref.mean
 feasible=np.where(sgn>0,hi>=ref.mean+margin,lo<=np.maximum(0,ref.mean-margin))
 order=np.flatnonzero(keep&(confidence>=1e6)&(T>0)&feasible);order=order[np.argsort(-confidence[order],kind='stable')];called=order[:int(nconf)]
 goals=np.zeros(len(T));tail=-ndtri(.5*.0001*max(1,len(called))/keep.sum());head=max(8.,tail+1.)
 goals[called]=sgn[called]*np.linspace(head,tail,len(called))
 # A separate, declared DE mean head. Bulk changes are retained in integer T.
 mu=ref.mean.copy();mu[called]*=np.divide(bulk[called],ref.bulk[called],out=np.ones(len(called)),where=ref.bulk[called]>0)
 up=called[sgn[called]>0];dn=called[sgn[called]<0]
 lo[up]=np.maximum(lo[up],ref.mean[up]+.001+1e-4*ref.mean[up]);hi[dn]=np.minimum(hi[dn],np.maximum(0,ref.mean[dn]-.001-1e-4*ref.mean[dn]))
 assert np.all(lo<=hi),('infeasible declared mean',np.sum(lo>hi))
 low=0.;high=10.
 while np.clip(high*mu,lo,hi).sum()<1e6:high*=2
 for _ in range(55):
  mid=(low+high)/2
  if np.clip(mid*mu,lo,hi).sum()<1e6:low=mid
  else:high=mid
 requested=np.clip(mu*(low+high)/2,lo,hi);v=(400*requested/1e6-T/dh)/(1/dl-1/dh);v=np.clip(v,0,T)
 TL=np.floor(v).astype(np.int64);TL[up]=np.ceil(v[up]).astype(np.int64)
 need=int(nl*dl-TL.sum());free=np.ones(len(T),bool);free[called]=False
 candidates=np.flatnonzero(free&((TL<T) if need>0 else (TL>0)))
 priority=v[candidates]-TL[candidates];candidates=candidates[np.argsort(-priority if need>0 else priority,kind='stable')]
 for j in candidates:
  if need==0:break
  if need>0:TL[j]+=1;need-=1
  else:TL[j]-=1;need+=1
 assert need==0,('low stratum rounding',need)
 emit.last_means=requested.copy()
 nzero=ref.n-np.diff(ref.ptr)
 dense=cap<0 or cap>=2;rho=cap-2 if cap>=2 else max(cap,0.)
 before=time.monotonic();x,fallback=make_columns(T,TL,dl,dh,ql,oh,ol,ref.cdf,ref.cdf_scale,ref.n,nzero,ref.ties,keep,goals,jkg*rho,8 if dense else 16,2 if dense else 4);after=time.monotonic()
 fillers=np.flatnonzero(~keep)
 for start,n,d,stot in [(0,nl,dl,TL),(nl,nh,dh,T-TL)]:
  sub=x[start:start+n];f=stot[fillers];ftotal=int(f.sum());frow=np.full(n,ftotal//n);frow[:ftotal%n]+=1
  desired=depths[start:start+n]-frow
  residual,steps=balance(sub,keep,desired,goals!=0)
  assert not residual.any(),('row balancing',residual.min(),residual.max(),steps)
  fill(sub,fillers,f,depths[start:start+n]-sub.sum(1))
 assert np.array_equal(x.sum(0,dtype=np.int64),T) and np.array_equal(x.sum(1),depths)
 emit.last_fallback=fallback
 emit.last_timings=dict(setup=before-start_time,histogram=after-before,balance=time.monotonic()-after)
 return x,called,goals
