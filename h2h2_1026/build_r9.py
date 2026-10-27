"""AMENDMENT_R9 arms (consult_1026_amp round 9): K250, RH20, DET885 on B0's response. Label-free: uses only B0's
response, B0's frozen confidence and the control reference mean CPM (inputs/references/ref_H1/mean.npy)."""
import json,hashlib,shutil,numpy as np
from pathlib import Path
from scipy.stats import rankdata
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
with np.load(R/'scratch/confidences.npz') as z:c=z['B0'].astype(np.float64)
with np.load(R/'results/predictions_B0.npz') as z:r=z['response'].astype(np.float64)
m=np.load(R/'inputs/references/ref_H1/mean.npy').astype(np.float64)
elig=c>=1e6
# Expected reference detectability (Poisson z of the mean shift): sqrt(control CPM) * |expm1(response)|.
zdet=np.sqrt(np.maximum(m,0))[None,:]*np.abs(np.expm1(r))
def pct(x,mask):
  out=np.zeros_like(x)
  for i in range(x.shape[0]):
    e=np.flatnonzero(mask[i]);out[i,e]=(rankdata(x[i,e])-0.5)/len(e) if len(e) else 0
  return out
pz,pc=pct(zdet,elig),pct(c,elig)
# DET885: Skadi-style mixs75 at the frozen budget; 0.75 detectability + 0.25 B0 confidence rank, eligibility unchanged.
det=np.where(elig,1e6+0.75*pz+0.25*pc,c)
# RH20: B0's own top-885 (by confidence among eligible) keeps its membership; the 20 with the highest
# pct(detectability) x pct(confidence) move to positions 1-20 (ordered by that score); the rest keep B0's order.
rh=c.copy();head_overlap=[]
for i in range(c.shape[0]):
  e=np.flatnonzero(elig[i]);top=e[np.argsort(-c[i,e],kind='stable')][:885]
  H=pz[i,top]*pc[i,top];h=top[np.argsort(-H,kind='stable')][:20]
  rh[i,h]=1e12+H[np.argsort(-H,kind='stable')][:20]*1e6
  head_overlap.append(len(np.intersect1d(h,top[:20])))
assert np.array_equal(det>=1e6,elig) and np.array_equal(rh>=1e6,elig)
np.savez_compressed(R/'scratch/confidences_R9.npz',K250=c,RH20=rh,DET885=det,B0=c)
for a in ['K250','RH20','DET885']:shutil.copyfile(R/'results/predictions_B0.npz',R/f'results/predictions_{a}.npz')
ov=[]
for i in range(c.shape[0]):
  e=np.flatnonzero(elig[i]);a=e[np.argsort(-c[i,e],kind='stable')][:885];b=e[np.argsort(-det[i,e],kind='stable')][:885]
  ov.append(len(np.intersect1d(a,b))/max(1,len(a)))
out=dict(budgets=dict(K250=250,RH20=885,DET885=885),det_top885_overlap_vs_B0=dict(min=float(min(ov)),median=float(np.median(ov)),max=float(max(ov))),
  rh20_head_overlap_with_B0_top20=dict(median=float(np.median(head_overlap)),min=int(min(head_overlap)),max=int(max(head_overlap))),
  sha256={a:sha(R/f'results/predictions_{a}.npz') for a in ['K250','RH20','DET885']},confidence_sha256=sha(R/'scratch/confidences_R9.npz'))
(R/'results/r9_arms.json').write_text(json.dumps(out,indent=1));print(json.dumps(out,indent=1))
