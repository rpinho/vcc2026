import runtime_local as j
j.guard(launch=True)
import numpy as np, csv
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
from recover_cache import tidy

M='expr_mse_unbiased_capped_norm'
def read_json(p):return j.json.loads(p.read_text())
def load_means(path,targets):
    with np.load(path) as z:
        ix={t:i for i,t in enumerate(z['targets'])}
        return np.stack([z['means'][ix[t]] for t in targets])
def excluded_ss(x,own):
    ss=np.einsum('ij,ij->i',x,x)
    for i,g in own.items():ss[i]-=x[i,g]**2
    return ss

def main():
    z=np.load(j.ROOT/'reference.npz');ts=z['targets'][:-1];genes=z['genes'];real=z['means'][:-1];rctrl=z['means'][-1]
    c=z['jk'][:-1];cctrl=z['jk'][-1];c0=z['null_jk'];inp=z['input_mean']
    own={i:int(np.flatnonzero(genes==t)[0]) for i,t in enumerate(ts) if np.any(genes==t)}
    w=_row_weights(len(ts),len(genes),own)
    den=(excluded_ss(real-rctrl,own)-c-cctrl)*w;D=float(den.sum());R=float(w @ c)
    anchors=read_json(j.VCC/'h2h2_1026/results/official_B0.json')['anchors'];span=anchors['baseline'][M]-anchors['replicate'][M]
    v5=load_means(j.ROOT/'cached_v5b.npz',ts)
    amp=load_means(j.ROOT/'bulk_v8b_like.npz',ts);bill=load_means(j.ROOT/'bulk_levers_knob_none.npz',ts)
    packages={};stats={};checks={}
    arms=dict(v5b=v5,v9a_885=v5,v9a_694=v5,v8b_like=amp,levers_knob_none=bill,incumbent=bill)
    for arm,means in arms.items():
        j.guard();rows=tidy(arm);record=read_json(j.VCC/'final_push/artifacts/de_check/H1_rich_search'/arm/'score.json')
        u=np.array([rows['expr_mse_unbiased'][t] for t in ts]);cap=np.array([rows['expr_mse_unbiased_capped'][t] for t in ts])
        e=excluded_ss(means-real,own)
        # Identity recovery of whole-gene prediction jackknife from exact sums
        # and original uncapped per-target metric. No parameter estimation.
        chat=e-c-u/w
        if arm=='v5b':
            cz=np.load(j.ROOT/'cached_v5b.npz');ci={t:i for i,t in enumerate(cz['targets'])};direct=np.array([cz['jk'][ci[t]] for t in ts])
            assert np.max(abs(chat-direct))<1e-10
        assert chat.min()>0
        B=_across_pert_budget(means,own,w);J=float(w @ chat);Q=float(w @ np.minimum(chat,c));Q0=float(w @ np.minimum(chat,c0))
        rho=min(1.,B/Q);rho0=min(1.,B/Q0)
        predcap=(e-c-rho*np.minimum(chat,c))*w
        err=float(abs(predcap-cap).max())
        assert err<2e-12,(arm,err)
        b_err=abs(B-record['B_emitted']);j_err=abs(J-record['claim_emitted'])
        assert b_err<1e-12 and j_err<1e-12,(arm,B,record['B_emitted'],J,record['claim_emitted'])
        assert max(abs(den-np.array([rows['expr_distance_unbiased'][t] for t in ts])))<1e-12
        current=float(cap.sum()/D);assert abs(current-record['raw'][M])<1e-12
        E=float(w @ e)
        vals={
            'capped_clamp_lifted':current,
            'unbiased_uncapped':float(u.sum()/D),
            'trace_cap_lifted_budget_retained':(E-R-min(B,J))/D,
            'budget_labels':(E-R-rho*Q)/D,
            'budget_label_free':(E-R-rho0*Q0)/D,
            'label_free_rho_labeled_claim':(E-R-rho0*Q)/D,
            'plain_pseudobulk_delta':float(w @ excluded_ss((means-inp)-(real-rctrl),own)/D),
            'plain_expression':E/D,
        }
        scaled={k:(anchors['baseline'][M]-v)/span for k,v in vals.items()}
        scaled['current_clamped']=float(np.clip(scaled['capped_clamp_lifted'],0,1))
        vals['current_clamped']=current
        packages[arm]=dict(B=B,self_claim=J,claim_labels=Q,rho_labels=rho,claim_label_free=Q0,rho_label_free=rho0,
            E=E,real_correction=R,D=D,cap_fires=int((chat>c).sum()),budget_binds=B<Q,
            pred_jk_min=float(chat.min()),pred_jk_max=float(chat.max()),raw=vals,scaled=scaled,
            original_prediction_sha256=record['prediction_sha256'])
        stats[arm]=(means,chat)
        checks[arm]=dict(capped_row_max_error=err,B_error=b_err,self_claim_error=j_err,
            normalized_raw_error=abs(current-record['raw'][M]))
    # Scores are read after all variants have been frozen in VARIANTS.md.
    pairs={'C1':('v5b','v9a_885',-.0007),'C3':('v5b','v9a_694',.0663),
        'C4':('v8b_like','v9a_694',.0791),'C5':('levers_knob_none','v9a_885',-.0219),
        'C6':('v9a_694','v9a_885',-.0670)}
    contrasts={}
    variants=list(packages['v5b']['scaled'])
    for v in variants:
        contrasts[v]={}
        for name,(a,b,off) in pairs.items():
            delta=packages[a]['scaled'][v]-packages[b]['scaled'][v]
            if abs(delta)<1e-12:delta=0.
            ok=bool(np.sign(delta)==np.sign(off))
            contrasts[v][name]=dict(delta_scaled=delta,official_scaled=off,ratio=delta/off,sign_ok=ok,
                pair=[a,b],genuine_local_analog=name!='C6',
                loo_sign_ok=ok,loo_method='Fixed formula; no fitted parameters, identical when this contrast is held out.')
    gate={v:dict(surrogate_signs_pass=all(contrasts[v][k]['sign_ok'] for k in ['C4','C3','C6']),
        c4_ratio_pass=.3<=contrasts[v]['C4']['ratio']<=3,
        qualified=False,reason='Genuine C6 absent; also inspect direct sign/ratio failures.') for v in variants}
    out=dict(packages=packages,contrasts=contrasts,gate=gate,checks=checks,span_local=span,
        D=D,n_targets=len(ts),genes=len(genes),own_exclusions=len(own),
        no_parameters_fitted=True,official_79_used_only_for_diagnostic_contrasts=True,
        C6='Budget 694-vs-885 surrogate only; genuine v8-v9 local pair missing.',
        verdict='instrument not found',conditional_amplitude_reading='Not issued: instrument gate failed.')
    j.save('results.json',out)
    # Compact statistics retained for a separate scorer/aggregation verifier.
    np.savez_compressed(j.ROOT/'prediction_jk.npz',targets=ts,**{a:v[1] for a,v in stats.items()})
    with open(j.ROOT/'budgets.csv','w') as f:
        keys=['B','self_claim','claim_labels','rho_labels','claim_label_free','rho_label_free','cap_fires','budget_binds']
        wri=csv.DictWriter(f,fieldnames=['package']+keys);wri.writeheader()
        for arm,p in packages.items():wri.writerow({'package':arm,**{k:p[k] for k in keys}})
    with open(j.ROOT/'contrast_variants.csv','w') as f:
        wri=csv.DictWriter(f,fieldnames=['variant','contrast','delta_scaled','official_scaled','ratio','sign_ok','genuine_local_analog','loo_sign_ok']);wri.writeheader()
        for v,cs in contrasts.items():
            for k,cx in cs.items():wri.writerow(dict(variant=v,contrast=k,**{key:cx[key] for key in wri.fieldnames[2:]}))
    for v,cs in contrasts.items():print(v,' '.join(f'{k}={c["delta_scaled"]:+.6g} r={c["ratio"]:+.4g}' for k,c in cs.items()),flush=True)
    j.status('All prespecified variants computed; retained row, B, claim and denominator checks passed. Instrument not found.')
if __name__=='__main__':main()
