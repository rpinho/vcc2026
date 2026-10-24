"""Measure actual feature construction + optimizer; never launch full training."""
import common as rt
import argparse,time,platform,gc
import numpy as np,torch
from cache import FeatureCache
from network import ResponseMLP,objective

def sync(device):
    if device=='mps':torch.mps.synchronize()
    elif device=='cuda':torch.cuda.synchronize()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--device',choices=['mps','cpu','cuda'],required=True);ap.add_argument('--steps',type=int,default=24);args=ap.parse_args()
    torch.set_num_threads(2);torch.manual_seed(1019)
    src=FeatureCache();rng=np.random.default_rng(1019)
    device=args.device;net=ResponseMLP(src.g).to(device);opt=torch.optim.AdamW(net.parameters(),lr=.001,weight_decay=.01)
    times=[];feature_times=[]; losses=[]
    for i in range(args.steps+5):
        sync(device);start=time.perf_counter()
        select=rng.choice(src.rows,64,replace=False);batch=src.batch(select);ft=time.perf_counter()-start
        x,y,m,ids=[torch.from_numpy(a).to(device) for a in batch]
        opt.zero_grad(set_to_none=True);p=net(x);loss,mse,contrast=objective(p,y,m,ids);loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),5);opt.step();sync(device)
        elapsed=time.perf_counter()-start
        if i>=5:times.append(elapsed);feature_times.append(ft);losses.append([float(loss.detach().cpu()),float(mse.detach().cpu()),float(contrast.detach().cpu())])
        if i%5==0:print(device,i,'end-to-end seconds',elapsed,'features',ft,flush=True)
    median=float(np.median(times));mean=float(np.mean(times));steps=80000
    rec=dict(device=device,platform=platform.platform(),torch=torch.__version__,batch_size=64,genes=src.g,hidden=256,
        parameters=sum(p.numel() for p in net.parameters()),source_rows=src.rows,warmup_steps=5,measured_steps=len(times),
        step_seconds=times,feature_seconds=feature_times,median_seconds=median,mean_seconds=mean,
        p90_seconds=float(np.quantile(times,.9)),losses=losses,
        schedule_origin='80,000-step cap, batch64 and lr1e-3 from GenBio reference launcher; not disclosed by Bill M',
        full_production_steps=3*steps,production_hours_mean=mean*3*steps/3600,
        production_hours_median=median*3*steps/3600,
        two_loco_folds_plus_production_hours_mean=mean*9*steps/3600,
        includes='cached feature gathering, host-device transfer, forward, MSE+contrastive, backward, AdamW, device synchronize',
        excludes='validation, checkpointing, count generation and official scoring',
        mps_peak_allocated_bytes=torch.mps.driver_allocated_memory() if device=='mps' else None)
    rt.save(rt.ROOT/f'artifacts/benchmark_{device}.json',rec)
    print('RESULT',device,rec['production_hours_mean'],rec['two_loco_folds_plus_production_hours_mean'],flush=True)
    rt.status('Benchmark '+device+' complete; production estimate '+str(rec['production_hours_mean'])+' hours. Full training has NOT started. See artifacts/benchmark_'+device+'.json.')
if __name__=='__main__':main()
