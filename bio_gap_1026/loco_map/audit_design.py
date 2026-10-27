"""Independent checkpoint, source-exclusion and selector audit; no model imports."""
import runtime as rt
rt.enforce(); rt.guard(True)
import csv

recipients=('H1','K562','RPE1','HepG2','Jurkat')
weights=(0.,.005,.01,.02,.03,.05,.1,.2,.5,1.)
receipt=rt.json.loads((rt.ROOT/'results/input_receipt.json').read_text())
for path,digest in receipt['pins'].items():
    rt.guard()
    assert rt.sha(path)==digest,path
assert len(receipt['genes'])==6006
assert receipt['heldout79_identities_excluded_every_bank']
checks=[]
for outer in recipients:
    completion=rt.json.loads((rt.ROOT/f'results/{outer}_completion.json').read_text())
    assert completion['PASS'] and completion['workers']==1 and completion['nice']>=10
    assert all(pool['num_threads']==1 for pool in completion['blas'])
    choices=rt.json.loads((rt.ROOT/f'results/{outer}_choices.json').read_text())
    with (rt.ROOT/f'results/{outer}_selection.csv').open() as f: rows=list(csv.DictReader(f))
    for arm in ('pooled','control'):
        rr=[row for row in rows if row['arm']==arm]
        candidates=[]
        for h in range(9):
            for w in weights:
                group=[row for row in rr if int(row['h'])==h and float(row['weight'])==w]
                assert sorted(row['recipient'] for row in group)==sorted(set(recipients)-{outer})
                score=sum(float(row['cross']) for row in group)/4
                low=min(float(row['PDS_lower_normal']) for row in group)
                candidates.append(dict(h=h,weight=w,score=score,low=low))
        for kind in ('map','bounded'):
            eligible=[c for c in candidates if (c['weight']==1 if kind=='map' else c['low']>=-.005-1e-14)]
            expected=max(eligible,key=lambda c:round(c['score'],14))
            chosen=choices[arm][kind]
            assert chosen['h']==expected['h'] and chosen['weight']==expected['weight']
            assert abs(chosen['training_macro_cross']-expected['score'])<1e-12
            assert sorted(chosen['training_recipients'])==sorted(set(recipients)-{outer})
            checks.append(dict(outer=outer,arm=arm,kind=kind,h=expected['h'],weight=expected['weight']))
    fits=rt.json.loads((rt.ROOT/f'results/{outer}_fits.json').read_text())
    for fit in fits:
        if fit.get('resumed_inner_selection'):
            assert outer=='H1'
            continue
        if fit['arm']=='H1_direct':
            assert fit['training_recipient']=='H1' and fit['training_targets']==168
            assert set(fit['source_excluded'])=={'H1',outer}
            continue
        valid=fit['validation']
        expected=set(recipients)-{outer}-(set() if valid is None else {valid})
        assert {entry['recipient'] for entry in fit['training']}==expected
        for entry in fit['training']:
            required={outer,entry['recipient']}|(set() if valid is None else {valid})
            assert set(entry['source_exclusions'])==required
        v=fit['verification']
        assert v['normal_equation_relative']<1e-9
        assert v['independent_cholesky_SVD_RRR_relative']<1e-7
rt.save('results/design_verification.json',dict(PASS=True,pinned_inputs=len(receipt['pins']),
    selectors=checks,source_exclusions=True,completed_recipients=list(recipients),
    limitation='H1 inner-fit numerical checks inherited from the completed initial run; selected choices independently replayed from its saved full grid.'))
rt.status(f'PASS independent design audit: {len(checks)} selectors, input hashes, final/source-exclusion receipts and resources.')
