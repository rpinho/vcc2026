"""Verify retained scalar receipts and independently replay selection/aggregation.
Long-double expression-vector checks run in finalize_release.py before compaction.
"""
import runtime as rt
rt.enforce();rt.guard()
import numpy as np,pandas as pd,json
records=[json.loads(p.read_text()) for p in sorted((rt.ROOT/'finalized').glob('*.json'))];assert len(records)==168
mpchecks=[r for r in records if r['mpmath_50digit_max_relative_error'] is not None]
assert len(mpchecks)==3 and max(r['mpmath_50digit_max_relative_error'] for r in mpchecks)<1e-11
state=json.loads((rt.ROOT/'stream_receipt.json').read_text());assert state['complete']
assert len(state['receipts'])==3
assert {r['url']:r['size'] for r in state['receipts']}==rt.URLS
assert sum(r['allowed_cells'] for r in state['receipts'])==234659
assert state['full_sum_checks']['count_sum']['max_abs']==0
assert sum(r['n_total'] for r in records)==196483
with np.load(rt.ROOT.parent/'results/baseline.npz') as z:targets=z['targets'].astype(str)
assert set(r['target'] for r in records)==set(targets)
df=pd.read_csv(rt.ROOT/'per_target_seed_role.csv').set_index(['target','seed','select'])
assert len(df)==1680
for record in records:
    assert record['full_sum_checks']['count_sum']['max_abs']==0
    assert record['independent_longdouble_max_relative_error']<1e-9
    assert record['proxy_BH_longdouble_checks']==10
    assert record['full_truth_max_abs_error']<1e-6
    for row in record['rows']:
        saved=df.loc[(row['target'],row['seed'],row['select'])]
        for key in ['n_select','n_evaluate','norm_select','proxy_nconf','proxy_bin','ceiling_cosine','B0_cosine','ceiling_sse','B0_sse','truth_energy','SB_repeat_full']:
            assert abs(row[key]-saved[key])<1e-10*max(1,abs(row[key]))
        assert saved.n_select+saved.n_evaluate==record['n_total']
        assert saved.control_n_select+saved.control_n_evaluate==38176
for seed in state['seeds']:
    for h in 'AB':
        rows=df.xs((seed,h),level=('seed','select')).loc[targets]
        terc=np.empty(168,int);terc[np.lexsort((targets,rows.norm_select.to_numpy()))]=np.arange(168)//56
        np.testing.assert_array_equal(rows.norm_bin,terc)
        np.testing.assert_array_equal(rows.proxy_bin,np.digitize(rows.proxy_nconf,[10,1000]))
direction=pd.read_csv(rt.ROOT/'directional_bins.csv');worst=0.
for _,row in direction.iterrows():
    frame=df.xs((int(row.seed),row['select']),level=('seed','select'))
    name=row['bin']
    if name.startswith('Norm '):frame=frame[frame.norm_bin==['low','middle','high'].index(name[5:])]
    if name.startswith('Proxy '):frame=frame[frame.proxy_bin==['<10','10–999','≥1000'].index(name[6:])]
    assert len(frame)==int(row.n)
    if len(frame):
        # Long-double aggregate replay uses independently verified scalar SSE/energies.
        e=frame.truth_energy.to_numpy(np.longdouble)
        values=dict(ceiling_cosine=frame.ceiling_cosine.to_numpy(np.longdouble).mean(),B0_cosine=frame.B0_cosine.to_numpy(np.longdouble).mean(),ceiling_E=frame.ceiling_sse.to_numpy(np.longdouble).sum()/e.sum(),B0_E=frame.B0_sse.to_numpy(np.longdouble).sum()/e.sum())
        for key,value in values.items():
            worst=max(worst,float(abs(value-row[key])));assert abs(value-row[key])<1e-10
res=json.loads((rt.ROOT/'results.json').read_text())
for row in res['rows']:
    panel=direction[direction.bin==row['bin']]
    for key in ['ceiling_cosine','B0_cosine','ceiling_E','B0_E','SB_repeat_full','latent_full_cosine','full_B0_cosine']:
        if row[key] is None:continue
        vals=panel.groupby('seed')[key].mean().to_numpy()
        assert abs(np.mean(vals)-row[key])<1e-10
        assert abs(np.std(vals,ddof=1)-row[key+'_sd'])<1e-10
    if row.get('cosine_ratio') is not None:
        assert abs(row['B0_cosine']/row['ceiling_cosine']-row['cosine_ratio'])<1e-12
        assert abs(row['B0_E']/row['ceiling_E']-row['E_ratio'])<1e-12
rt.save('verification.json',dict(status='PASS',numpy_longdouble_precision_bits=53,mpmath_50digit_targets=3,mpmath_50digit_max_relative_error=max(r['mpmath_50digit_max_relative_error'] for r in mpchecks),target_seed_role_rows=len(df),independent_longdouble_max_relative_error=max(r['independent_longdouble_max_relative_error'] for r in records),proxy_BH_longdouble_checks=sum(r['proxy_BH_longdouble_checks'] for r in records),directional_aggregate_max_abs_error=worst,selection_memberships_verified=True,seed_mean_and_sample_sd_verified=True,full_sum_checks_per_target_verified=168,all_control_release_hashes_match=len({r['control_ids_sha256'] for r in state['receipts']})==1,heldout_expression_rows_read=sum(r['heldout_expression_rows_read'] for r in state['receipts']),resource=rt.guard()))
rt.status(f'Verification PASS: 1680 target/seed/role rows independently checked before compaction (NumPy longdouble=53 bits here), plus three 50-digit target checks; selecting-half bins, role aggregation and five-seed mean/SD replayed.')
