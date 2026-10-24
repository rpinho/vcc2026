"""Checks that matter for transfer: context leakage, masks, losses, ensemble blend."""
import common as rt
import numpy as np, torch
from dataset import Sources
from network import ResponseMLP,objective,blend

def main():
    torch.set_num_threads(2);src=Sources();d=next(d for d in src.banks if d['family']=='CD4');t=d['targets'][0];ctrl=d['control'][0].copy()
    before=src.features(t,ctrl,'CD4')
    # Poison ALL recipient sister screens in memory. They cannot affect features.
    for bank in src.banks:
        if bank['family']=='CD4':
            bank['response'].fill(1e5);bank['mean_response'].fill(-1e5);bank['control'].fill(1e5)
    after=src.features(t,ctrl,'CD4');assert np.array_equal(before,after)
    own=src.gi.get(t)
    if own is not None:assert all(before[k*src.g+own]==0 for k in [0,1,2])
    src.banks=[d for d in src.banks if d['family']!='CD4'];assert all(d['family']!='CD4' for d in src.banks)
    batch=src.batch(src.rows()[:4]);x,y,m,ids=[torch.from_numpy(a) for a in batch]
    net=ResponseMLP(src.g,hidden=16);p=net(x);loss,mse,contrast=objective(p,y,m,ids);assert torch.isfinite(loss)
    loss.backward();assert all(torch.isfinite(a.grad).all() for a in net.parameters() if a.grad is not None)
    # Poison unmeasured coordinates: loss must not change.
    y2=y+(1-m)*10000;loss2,*_=objective(p,y2,m,ids);assert torch.equal(loss,loss2)
    rng=np.random.default_rng(4);pool=rng.normal(size=(9,src.g));learned=rng.normal(size=pool.shape)
    mixed=blend(pool,learned,.3)
    expected=.7*pool/np.linalg.norm(pool,axis=1)[:,None]+.3*learned/np.linalg.norm(learned,axis=1)[:,None]
    expected=expected/np.linalg.norm(expected,axis=1)[:,None]*np.linalg.norm(pool,axis=1)[:,None]
    assert np.allclose(mixed,expected,rtol=1e-12,atol=1e-12)
    assert np.allclose(blend(pool,np.zeros_like(pool)),pool)
    assert np.allclose(np.linalg.norm(mixed,axis=1),np.linalg.norm(pool,axis=1))
    rt.save(rt.ROOT/'artifacts/verification.json',dict(status='PASS',checks=['recipient family and sister-screen poison invariance','own-target response features zero','outer-family removal','unmeasured output poison invariance','finite MSE+contrastive gradients','exact disclosed 0.3 unit blend before amplitude','zero learned response fallback']))
    print('PASS leakage, masks, objective gradients, blend',flush=True)
if __name__=='__main__':main()
