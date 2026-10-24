"""Bill-disclosure reconstruction, not Bill's unreleased implementation."""
import common as rt
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

class ResponseMLP(nn.Module):
    def __init__(self, genes=18533, hidden=256):
        super().__init__(); self.genes=genes
        self.net=nn.Sequential(nn.Linear(5*genes,hidden),nn.ReLU(),nn.Dropout(.1),
            nn.Linear(hidden,hidden),nn.ReLU(),nn.Dropout(.1),nn.Linear(hidden,genes))
        nn.init.zeros_(self.net[-1].weight); nn.init.zeros_(self.net[-1].bias)
    def forward(self,x): return x[:,:self.genes]+self.net(x)

def objective(pred, real, mask, target_ids, weight=.05, temperature=.1):
    # Each row has its own measured-output and own-target exclusion mask.
    p=pred*mask; y=real*mask
    mse=(((p-y)**2).sum(1)/mask.sum(1).clamp_min(1)).mean()
    logits=F.normalize(p,dim=1)@F.normalize(y,dim=1).T/temperature
    same=target_ids[:,None]==target_ids[None,:]
    eye=torch.eye(len(p),dtype=torch.bool,device=p.device)
    logits=logits.masked_fill(same&~eye,-1e4)
    labels=torch.arange(len(p),device=p.device)
    contrast=.5*(F.cross_entropy(logits,labels)+F.cross_entropy(logits.T,labels))
    return mse+weight*contrast, mse, contrast

def blend(pooled, learned, coefficient=.3):
    # A missing/degenerate learned response falls back to the donor direction.
    norm=np.linalg.norm(pooled,axis=-1,keepdims=True)
    ln=np.linalg.norm(learned,axis=-1,keepdims=True)
    unit=pooled/np.maximum(norm,1e-12)
    lu=np.where(ln>1e-12,learned/np.maximum(ln,1e-12),unit)
    mixed=(1-coefficient)*unit+coefficient*lu
    return mixed/np.maximum(np.linalg.norm(mixed,axis=-1,keepdims=True),1e-12)*norm
