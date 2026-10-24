"""Compiled AtlasShift integer row rounding and exact margin repair (MIT algorithm).
Stable mergesort matches upstream stable ordering; RNG uniforms drawn by NumPy Generator.
"""
import numpy as np
from numba import njit
@njit(cache=True)
def round_rows(expected,depths,uniforms):
 n,g=expected.shape;x=np.empty((n,g),np.int32)
 for i in range(n):
  s=0;fracsum=0.
  for j in range(g):
   z=int(np.floor(expected[i,j]));x[i,j]=z;s+=z;fracsum+=expected[i,j]-z
  need=int(depths[i])-s
  if need:
   scale=need/fracsum;cur=0.;placed=0
   for j in range(g):
    cur+=expected[i,j]-np.floor(expected[i,j])
    while placed<need and cur*scale>placed+uniforms[i]:x[i,j]+=1;placed+=1
   assert placed==need
 return x
@njit(cache=True)
def repair(x,expected):
 n,g=x.shape;columns=np.zeros(g);target=np.empty(g,np.int64);actual=np.zeros(g,np.int64);S=0
 for i in range(n):
  for j in range(g):columns[j]+=expected[i,j];actual[j]+=x[i,j];S+=x[i,j]
 for j in range(g):target[j]=int(np.floor(columns[j]))
 rem=S-np.sum(target);assert rem>=0 and rem<=g
 if rem:
  chosen=np.argsort(columns-target,kind='mergesort')[-rem:]
  for j in chosen:target[j]+=1
 diff=actual-target;capacity=np.zeros(n,np.int64)
 for gene in np.flatnonzero(diff>0):
  excess=diff[gene];ru=x[:,gene]-np.floor(expected[:,gene]).astype(np.int64);rows=np.flatnonzero(ru>0);assert np.sum(ru[rows])>=excess
  score=-(x[rows,gene]-expected[rows,gene]);rows=rows[np.argsort(score,kind='mergesort')]
  if excess<=len(rows):
   for ri in range(excess):row=rows[ri];x[row,gene]-=1;capacity[row]+=1
  else:
   for row in rows:
    take=min(excess,ru[row]);x[row,gene]-=take;capacity[row]+=take;excess-=take
    if excess==0:break
 genes=np.flatnonzero(diff<0);genes=genes[np.argsort(diff[genes],kind='mergesort')]
 for gene in genes:
  deficit=-diff[gene]
  while deficit:
   rows=np.flatnonzero(capacity>0);assert len(rows)>0;take=min(deficit,len(rows));score=expected[rows,gene]-x[rows,gene];rows=rows[np.argsort(-score,kind='mergesort')]
   for k in range(take):row=rows[k];x[row,gene]+=1;capacity[row]-=1
   deficit-=take
 assert np.sum(capacity)==0
 return x

def quantize(expected,depths,seed):
 rng=np.random.default_rng(seed);need=depths-np.floor(expected).sum(1).astype(np.int64);u=np.zeros(len(depths));u[need!=0]=rng.random(int((need!=0).sum()))
 return repair(round_rows(expected,np.asarray(depths,np.int64),u),expected)
