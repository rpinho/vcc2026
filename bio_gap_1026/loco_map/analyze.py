"""Strict nested leave-one-recipient-context-out diagnostic, CPU only."""
import runtime as rt
rt.enforce(); rt.guard(True)
import argparse, gc, os, signal, time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_info
from models import Data, RECIPIENTS, HYPERS, train_map, predict_candidates, direct_h1, mask_own
from selection import Validation, choose, spread
from verify import evaluate

START=time.monotonic()

def materialize(candidate,data,r):
    L,W=candidate
    return mask_own(L @ W,data.recips[r]['targets'],data.genes).astype(np.float32)

def bounded(base,component,w):
    scale=np.sqrt(spread(base)/spread(component))
    return ((1-w)*base.astype(float)+w*scale*component.astype(float)).astype(np.float32)

def run_outer(data,outer):
    rt.guard(True)
    training=[r for r in RECIPIENTS if r!=outer]
    rt.status(f'OUTER {outer}: no {outer} outcome row enters any fit or source transform. Inner contexts={training}.')
    choicepath=rt.ROOT/f'results/{outer}_choices.json'
    choices=rt.json.loads(choicepath.read_text()) if choicepath.exists() else {}
    selectionpath=rt.ROOT/f'results/{outer}_selection.csv'
    allrecords=pd.read_csv(selectionpath).to_dict('records') if selectionpath.exists() else []
    fit_records=[]
    for conditional,arm in ((False,'pooled'),(True,'control')):
        if arm in choices:
            assert all(outer not in c['training_recipients'] for c in choices[arm].values())
            rt.status(f'{outer}/{arm}: reusing completed training-only selection checkpoint; no inner fits repeated.')
            fit_records.append(dict(outer=outer,arm=arm,resumed_inner_selection=True,
                selection_file=str(selectionpath),choice_file=str(choicepath)))
            continue
        records=[]
        for valid in training:
            trainers=[r for r in training if r!=valid]
            rt.status(f'{outer}/{arm}: fitting on {trainers}, validating {valid}; sources globally exclude {outer},{valid}.')
            model=train_map(data,trainers,{outer,valid},conditional)
            fit_records.append(dict(outer=outer,arm=arm,validation=valid,training=model['training'],
                retained_input_energy=model['retained_energy'],verification=model['fit_verification']))
            p,_=data.plain(valid,{outer,valid})
            candidates=predict_candidates(model,p,data.recips[valid]['weight'],data.recips[valid]['targets'],data.genes)
            val=Validation(data,valid,p)
            rr=val.grid(candidates)
            for row in rr: row.update(outer=outer,arm=arm)
            records.extend(rr)
            rt.status(f'{outer}/{arm}: inner {valid} complete; PCA retained={model["retained_energy"]:.4f}; elapsed={time.monotonic()-START:.1f}s.')
            del model,candidates,val,rr,p; gc.collect()
        choices[arm]=choose(records); allrecords.extend(records)
        rt.status(f'{outer}/{arm} TRAINING-ONLY choices: {choices[arm]}')
        pd.DataFrame(allrecords).to_csv(rt.ROOT/f'results/{outer}_selection.csv',index=False)
        rt.save(f'results/{outer}_choices.json',choices)
    # Final fits use only other recipients; no access to outer Y until evaluate calls below.
    p,coverage=data.plain(outer,{outer})
    metricpath=rt.ROOT/f'results/{outer}_metrics.json'
    exactpath=rt.ROOT/f'results/{outer}_exact.json'
    statpath=rt.ROOT/f'results/{outer}_stats.npz'
    summary=rt.json.loads(metricpath.read_text()) if metricpath.exists() else {}
    certificates=rt.json.loads(exactpath.read_text()) if exactpath.exists() else {}
    stats={}
    if statpath.exists():
        with np.load(statpath) as saved: stats={k:saved[k] for k in saved.files}
    basecache=None
    if 'plain' in summary:
        basecache={k:stats['plain__'+k] for k in ('dot','energy','pds','Bboot')}
        basecache.update(B=summary['plain']['raw_B']*float(stats['plain__truth_energy'].sum()),
                         certificate=certificates['plain'])
    target=data.recips[outer]['targets']; truth=data.recips[outer]['truth']
    predictions={'plain':p}
    def score(name,q,baseline=p,genes=data.genes,y=truth,targets=target,cache=True):
        nonlocal basecache
        if name in summary:
            rt.status(f'{outer}/{name}: reusing completed independently certified metrics; no scoring repeated.')
            return summary[name]
        rt.status(f'{outer}/{name}: independent metric/PDS/bootstrap and exact arithmetic verification starting.')
        item,certificate,small,b= evaluate(q,baseline,y,targets,genes,baseline_cache=basecache if cache else None)
        if cache: basecache=b
        summary[name]=item; certificates[name]=certificate
        for k,v in small.items(): stats[name+'__'+k]=v
        rt.save(f'results/{outer}_metrics.json',summary)
        rt.save(f'results/{outer}_exact.json',certificates)
        np.savez_compressed(rt.ROOT/f'results/{outer}_stats.npz',**stats)
        rt.status(f'{outer}/{name}: C={item["cross"]:.6f}, A={item["aligned_per_spread"]:.6f}, PDS={item["PDS"]:.6f}, deltaC CI={item["delta_cross_CI95"]}, deltaPDS={item["delta_PDS"]:.6f}.')
        return item
    score('plain',p)
    if outer=='H1':
        score('B0_common',data.B0)
        # Original axis is retained for the historical 0.1154 baseline comparison.
        score('B0_original',data.fullB0,data.fullB0,data.fullgenes,data.fullY,data.h1targets,False)
        assert abs(summary['B0_original']['cross']-.1154)<.0001
        assert abs(summary['B0_original']['PDS']-.8558953521528372)<1e-6
    for conditional,arm in ((False,'pooled'),(True,'control')):
        model=train_map(data,training,{outer},conditional)
        fit_records.append(dict(outer=outer,arm=arm,validation=None,training=model['training'],
            retained_input_energy=model['retained_energy'],verification=model['fit_verification']))
        candidates=predict_candidates(model,p,data.recips[outer]['weight'],target,data.genes)
        mapchoice=choices[arm]['map']; blendchoice=choices[arm]['bounded']
        qm=materialize(candidates[mapchoice['h']],data,outer)
        score(arm+'_map',qm)
        del qm
        comp=materialize(candidates[blendchoice['h']],data,outer)
        blend=bounded(p,comp,blendchoice['weight'])
        score(arm+'_bounded',blend)
        if outer=='H1':
            bcommon=bounded(data.B0,comp,blendchoice['weight'])
            score('B0_'+arm+'_bounded_common',bcommon)
            # Embed only the modeled gene update, preserving B0 elsewhere.
            full=data.fullB0.copy(); full[:,data.fullix]=bcommon
            score('B0_'+arm+'_bounded_original',full,data.fullB0,data.fullgenes,data.fullY,data.h1targets,False)
            del bcommon,full
        del comp,blend,candidates,model; gc.collect()
    if outer!='H1':
        q,info=direct_h1(data,outer)
        score('H1_direct',q); fit_records.append(dict(outer=outer,arm='H1_direct',**info))
        del q
    rt.save(f'results/{outer}_fits.json',fit_records)
    rt.save(f'results/{outer}_completion.json',dict(PASS=True,outer=outer,coverage=int((coverage>0).sum()),
        total=len(p),seconds=time.monotonic()-START,free_gib=rt.guard(),pid=os.getpid(),
        blas=threadpool_info(),nice=os.getpriority(os.PRIO_PROCESS,0),workers=1,bootstrap=2000))
    rt.status(f'OUTER {outer} COMPLETE. All final numerical diagnostics exactly certified; bootstrap panel/fits fixed.')
    del predictions,p,basecache,stats; gc.collect()

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--outer',choices=RECIPIENTS); args=parser.parse_args()
    signal.alarm(4*3600-15*60)
    rt.save('results/run.json',dict(pid=os.getpid(),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        free_gib=rt.guard(True),workers=1,requested_nice=10,effective_nice=os.getpriority(os.PRIO_PROCESS,0)))
    data=Data()
    order=[args.outer] if args.outer else ['H1','RPE1','Jurkat','HepG2','K562']
    for outer in order:
        if (rt.ROOT/f'results/{outer}_completion.json').exists():
            rt.status(f'OUTER {outer}: complete checkpoint reused.'); continue
        if time.monotonic()-START>3.5*3600:
            rt.status('Time-cap checkpoint before next recipient; report completed recipients only.'); break
        run_outer(data,outer)
    signal.alarm(0)
    rt.status(f'Analysis finished after {time.monotonic()-START:.1f}s. No worker pool or children started.')

if __name__=='__main__': main()
