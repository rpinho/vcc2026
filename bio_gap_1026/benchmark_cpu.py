"""Bounded exact-network optimizer timing; no model/checkpoint retained."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,torch,time,argparse,os,subprocess
from torch import nn
from torch.nn import functional as F
torch.set_num_threads(1);torch.set_num_interop_threads(1)
def child(label):
    ns=dict(np=np,torch=torch,nn=nn,F=F)
    rt.definitions(rt.REPO/'model_billm/network.py',['ResponseMLP','objective'],ns)
    with np.load(rt.ROOT/'scratch/benchmark_batch.npz') as z:x=z['x'];y=z['y'].astype(np.float32);mask=z['mask']
    net=ns['ResponseMLP'](18533);opt=torch.optim.AdamW(net.parameters(),lr=.001,weight_decay=.01);rng=np.random.default_rng(1019)
    target_ids=torch.arange(64);times=[];ft=[]
    for step in range(110):
        rt.guard();start=time.monotonic();ids=rng.permutation(64)
        bx=torch.from_numpy(x[ids]);by=torch.from_numpy(y[ids]);bm=torch.from_numpy(mask[ids]);ft.append(time.monotonic()-start)
        opt.zero_grad(set_to_none=True);p=net(bx);loss,*_=ns['objective'](p,by,bm,target_ids)
        loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),5);opt.step()
        if step>=10:times.append(time.monotonic()-start)
    rt.save(f'results/benchmark_{label}.json',dict(pid=os.getpid(),nice=os.getpriority(os.PRIO_PROCESS,0),torch_threads=torch.get_num_threads(),blas_threads=__import__('threadpoolctl').threadpool_info(),batch=64,parameters=sum(p.numel() for p in net.parameters()),warmup_steps=10,measured_steps=len(times),mean_seconds=float(np.mean(times)),median_seconds=float(np.median(times)),p90_seconds=float(np.quantile(times,.9)),feature_gather_seconds=float(np.mean(ft[10:])),seconds=times,scope='Real cached feature gather and optimizer, exact 28.55M architecture/loss. No full source-cache creation, validation/checkpoint/emission. Short-run extrapolation only.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--child');args=ap.parse_args()
    if args.child:child(args.child)
    else:
        rt.status('Bounded CPU benchmark: one worker, then two concurrent one-thread workers, all nice >=10.')
        sessions=[]
        for labels in [['one'],['two_a','two_b']]:
            jobs=[]
            for label in labels:
                log=(rt.ROOT/f'benchmark_{label}.log').open('w');p=subprocess.Popen([os.sys.executable,'-B',__file__,'--child',label],stdout=log,stderr=subprocess.STDOUT)
                jobs.append((p,log));sessions.append(p.pid)
            rt.save('results/benchmark_pids.json',dict(children=sessions,controller=os.getpid()))
            for p,log in jobs:
                code=p.wait();log.close();assert code==0
        one=rt.json.loads((rt.ROOT/'results/benchmark_one.json').read_text())['mean_seconds']
        a=rt.json.loads((rt.ROOT/'results/benchmark_two_a.json').read_text())['mean_seconds'];b=rt.json.loads((rt.ROOT/'results/benchmark_two_b.json').read_text())['mean_seconds'];two=max(a,b)
        rt.save('results/benchmark_projection.json',dict(one_worker_seconds=one,two_worker_seconds=[a,b],one_fit_hours=80000*one/3600,three_seeds_two_workers_hours=80000*(two+one)/3600,nine_fits_two_workers_hours=80000*(4*two+one)/3600,nine_fits_one_worker_hours=720000*one/3600,old_report_hours=16.09,excludes='cache creation, validation, checkpointing and scoring; no convergence evidence; not a sustainable-throughput guarantee'))
        (rt.ROOT/'scratch/benchmark_batch.npz').unlink()
        rt.status('CPU benchmark complete; children reaped; temporary batch removed. No training run launched.')
