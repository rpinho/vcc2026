"""Cross-fitted bins from verified, compacted per-target half metrics."""
import runtime as rt
rt.enforce();resource=rt.guard()
import numpy as np,pandas as pd,json
state=json.loads((rt.ROOT/'stream_receipt.json').read_text());assert state['complete']
records=[json.loads(p.read_text()) for p in sorted((rt.ROOT/'finalized').glob('*.json'))]
assert len(records)==168
per=[row for record in records for row in record['rows']]
frame=pd.DataFrame(per)
assert len(frame)==1680 and frame.groupby(['seed','select']).size().eq(168).all()
directional=[]
for (seed,select),group in frame.groupby(['seed','select']):
    order=np.lexsort((group.target.to_numpy(),group.norm_select.to_numpy()))
    terc=np.empty(168,int);terc[order]=np.arange(168)//56
    frame.loc[group.index,'norm_bin']=terc
    bins=[('All',np.ones(168,bool))]+[(f'Norm {name}',terc==k) for k,name in enumerate(['low','middle','high'])]+[(f'Proxy {name}',group.proxy_bin.to_numpy()==k) for k,name in enumerate(['<10','10–999','≥1000'])]
    for name,mask in bins:
        part=group[mask];row=dict(seed=int(seed),select=select,bin=name,n=len(part))
        if len(part):
            row.update(ceiling_cosine=float(part.ceiling_cosine.mean()),B0_cosine=float(part.B0_cosine.mean()),ceiling_E=float(part.ceiling_sse.sum()/part.truth_energy.sum()),B0_E=float(part.B0_sse.sum()/part.truth_energy.sum()),SB_repeat_full=float(part.SB_repeat_full.mean()),latent_full_cosine=float(part.latent_full_cosine.mean()) if part.latent_full_cosine.notna().all() else None,full_B0_cosine=float(part.full_B0_cosine.mean()),invalid_SB=int(part.latent_full_cosine.isna().sum()))
        directional.append(row)
per=frame.to_dict('records')
truth_max_diff=max(r['full_truth_max_abs_error'] for r in records)
per=pd.DataFrame(per);per.to_csv(rt.ROOT/'per_target_seed_role.csv',index=False)
d=pd.DataFrame(directional);d.to_csv(rt.ROOT/'directional_bins.csv',index=False)
# Equal weight each selector role; sparse bins with an empty role are explicitly NA.
metrics=['ceiling_cosine','B0_cosine','ceiling_E','B0_E','SB_repeat_full','latent_full_cosine','full_B0_cosine']
seedrows=[];rows=[]
for name in d.bin.unique():
    grp=d[d.bin==name];sr=[]
    for seed in state['seeds']:
        pair=grp[grp.seed==seed];row=dict(bin=name,seed=seed,n=float(pair.n.mean()))
        for key in metrics:row[key]=float(pair[key].mean()) if len(pair)==2 and pair[key].notna().all() else None
        if row['ceiling_cosine'] is not None:row['cosine_ratio']=row['B0_cosine']/row['ceiling_cosine'];row['E_ratio']=row['B0_E']/row['ceiling_E']
        sr.append(row);seedrows.append(row)
    row=dict(bin=name,n=float(grp.n.mean()),n_min=int(grp.n.min()),n_max=int(grp.n.max()),valid_seed_count=sum(r['ceiling_cosine'] is not None for r in sr))
    for key in metrics:
        values=[r[key] for r in sr if r[key] is not None]
        row[key]=float(np.mean(values)) if len(values)==5 else None
        row[key+'_sd']=float(np.std(values,ddof=1)) if len(values)==5 else None
        row[key+'_range']=[min(values),max(values)] if len(values)==5 else None
    if row['ceiling_cosine'] is not None:
        row['cosine_ratio']=row['B0_cosine']/row['ceiling_cosine'];row['E_ratio']=row['B0_E']/row['ceiling_E']
        for ratio in ['cosine_ratio','E_ratio']:
            values=[r[ratio] for r in sr];row[ratio+'_seed_sd']=float(np.std(values,ddof=1));row[ratio+'_seed_range']=[min(values),max(values)]
        if row['latent_full_cosine'] is not None:row['full_B0_latent_ratio']=row['full_B0_cosine']/row['latent_full_cosine']
    rows.append(row)
pd.DataFrame(seedrows).to_csv(rt.ROOT/'seed_bins.csv',index=False)
pd.DataFrame(rows).to_csv(rt.ROOT/'summary.csv',index=False)
rt.save('results.json',dict(rows=rows,seeds=state['seeds'],resource=resource,streamed_full_truth_max_abs_error=truth_max_diff,split_cell_count_range=[int(per.n_select.min()),int(per.n_select.max())],negative_split_cosines=int((per.ceiling_cosine<0).sum()/2),invalid_SB_count=int(per.latent_full_cosine.isna().sum()/2),heldout_expression_rows_read=0,aggregation='cosine: target mean; E: ratio of summed SSE/evaluation energy, then equal A/B-role average, then seed mean; SD: five seed averages, ddof=1',SB='Per target r: repeat-full R=2r/(1+r), ideal noiseless predictor cosine sqrt(R); uncentered-cosine parallel-halves approximation; no E correction'))
rt.status('Five-seed cross-fitted analysis complete; summary.csv and results.json written, awaiting independent verification.')
print(pd.DataFrame(rows)[['bin','n','ceiling_cosine','ceiling_cosine_sd','B0_cosine','cosine_ratio','ceiling_E','B0_E','E_ratio','SB_repeat_full','latent_full_cosine']].to_string(index=False))
