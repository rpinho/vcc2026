"""Three optimizer steps only: resumed and uninterrupted CPU outputs must match."""
import common as rt
import torch,types
from train import run
folder=rt.ROOT/'artifacts/smoke/production'
resumed=torch.load(folder/'seed1019_weights.pt',map_location='cpu',weights_only=True)
(folder/'seed1019.pt').unlink()
args=types.SimpleNamespace(smoke=True,device='cpu',steps=3)
run(args,'production',1019)
fresh=torch.load(folder/'seed1019_weights.pt',map_location='cpu',weights_only=True)
assert set(fresh)==set(resumed)
assert all(torch.equal(fresh[k],resumed[k]) for k in fresh)
rt.save(rt.ROOT/'artifacts/resume_verification.json',dict(status='PASS exact CPU weights: 2+1 resumed steps equal 3 uninterrupted steps',steps=3,seed=1019,scope='optimizer checkpoint, numpy batch RNG and torch dropout RNG; CUDA reproducibility unmeasured'))
print('PASS exact checkpoint resume',flush=True)
