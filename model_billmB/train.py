"""Resumable three-seed fitting; explicit preflight enforces local 6-hour stop."""
import common as rt
import argparse,time,os
import numpy as np,torch
from cache import FeatureCache
from network import ResponseMLP,objective
SEEDS=[1019,2027,4093]

def run(args,fold,seed):
    excluded=[] if fold=='production' else [fold]
    src=FeatureCache(fold);rng=np.random.default_rng(seed)
    torch.manual_seed(seed);torch.set_num_threads(2)
    device=args.device;net=ResponseMLP(src.g).to(device);opt=torch.optim.AdamW(net.parameters(),lr=.001,weight_decay=.01)
    folder=rt.ROOT/('artifacts/smoke' if args.smoke else 'checkpoints')/fold;folder.mkdir(parents=True,exist_ok=True)
    cp=folder/f'seed{seed}.pt';start_step=0
    signature={p.stem:rt.sha(p) for p in (rt.ROOT/'data').glob('*.npz')}
    code_signature={n:rt.sha(rt.ROOT/n) for n in ['train.py','dataset.py','network.py','cache.py','common.py']}
    if cp.exists():
        checkpoint=torch.load(cp,map_location='cpu',weights_only=False)
        assert checkpoint['data_sha256']==signature and checkpoint['fold']==fold and checkpoint['seed']==seed
        assert checkpoint['code_sha256']==code_signature,'Resume code changed'
        net.load_state_dict(checkpoint['model']);opt.load_state_dict(checkpoint['optimizer'])
        for state in opt.state.values():
            for k,v in state.items():
                if isinstance(v,torch.Tensor):state[k]=v.to(device)
        torch.set_rng_state(checkpoint['torch_rng']);rng.bit_generator.state=checkpoint['numpy_rng'];start_step=checkpoint['step']
        if device=='cuda' and checkpoint.get('cuda_rng') is not None:torch.cuda.set_rng_state_all(checkpoint['cuda_rng'])
        if device=='mps' and checkpoint.get('mps_rng') is not None:torch.mps.set_rng_state(checkpoint['mps_rng'])
    started=time.monotonic();log=[]
    for step in range(start_step,args.steps):
        selection=src.sample_rows(rng,64);batch=src.batch(selection)
        x,y,m,ids=[torch.from_numpy(a).to(device) for a in batch]
        opt.zero_grad(set_to_none=True);pred=net(x);loss,mse,contrast=objective(pred,y,m,ids)
        if not torch.isfinite(loss):raise RuntimeError('nonfinite loss')
        loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),5);opt.step()
        if (step+1)%100==0 or step+1==args.steps:
            rec=dict(step=step+1,total=float(loss.detach().cpu()),mse=float(mse.detach().cpu()),contrast=float(contrast.detach().cpu()),seconds=time.monotonic()-started)
            log.append(rec);print(fold,seed,rec,flush=True)
        if (step+1)%2000==0 or step+1==args.steps:
            tmp=cp.with_suffix('.tmp')
            torch.save(dict(model=net.state_dict(),optimizer=opt.state_dict(),step=step+1,fold=fold,seed=seed,genes=np.load(rt.ROOT/'data/Jurkat.npz')['genes'],
                data_sha256=signature,code_sha256=code_signature,numpy_rng=rng.bit_generator.state,torch_rng=torch.get_rng_state(),
                cuda_rng=torch.cuda.get_rng_state_all() if device=='cuda' else None,
                mps_rng=torch.mps.get_rng_state() if device=='mps' else None),tmp);tmp.replace(cp)
            rt.save(folder/f'seed{seed}.json',dict(fold=fold,seed=seed,steps=step+1,log=log,excluded_families=excluded,
                donor_rule='recipient family excluded from every training example; outer family absent from all banks',data_sha256=signature,code_sha256=code_signature,
                sampling='uniform family then uniform row, distinct rows within batch; CD4 conditions share one family'))
    # Export inference-only CPU weights, avoiding optimizer payload on return.
    out=folder/f'seed{seed}_weights.pt';torch.save({k:v.detach().cpu() for k,v in net.state_dict().items()},out)
    return dict(fold=fold,seed=seed,steps=args.steps,weights=str(out.relative_to(rt.ROOT)),seconds=time.monotonic()-started)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--device',choices=['mps','cpu','cuda'],required=True)
    ap.add_argument('--steps',type=int,default=80000);ap.add_argument('--folds',default='H1,CD4,production');ap.add_argument('--seeds',default='1019,2027,4093');ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    if args.smoke:assert args.steps<=3,'Smoke is at most three steps'
    elif args.device!='cuda':
        path=rt.ROOT/f'artifacts/benchmark_{args.device}.json';assert path.exists(),'Run exact local benchmark first'
        measurement=rt.json.loads(path.read_text())
        hours=measurement['mean_seconds']*args.steps*len(args.seeds.split(','))*len(args.folds.split(','))/3600
        if hours>6:raise SystemExit(f'STOP: measured full requested training {hours:.2f}h >6h; use GPU_NEEDED.md')
    # Validate upload hashes before a paid run.
    manifest=rt.json.loads((rt.ROOT/'data/manifest.json').read_text())
    for name,rec in manifest.items():assert rt.sha(rt.ROOT/f'data/{name}.npz')==rec['sha256']
    result=[]
    for fold in args.folds.split(','):
        assert fold in ['H1','CD4','production']
        for seed in map(int,args.seeds.split(',')):
            result.append(run(args,fold,seed));rt.save(rt.ROOT/('artifacts/smoke/training.json' if args.smoke else 'artifacts/training.json'),result)
if __name__=='__main__':main()
