import runtime_local as j
r=j.json.loads((j.ROOT/'results.json').read_text())
v=j.json.loads((j.ROOT/'verification.json').read_text())
assert v['status']=='PASS'
assert v['results_sha256']==j.sha(j.ROOT/'results.json')

names=[('current_clamped','Current clamp'),('capped_clamp_lifted','Capped, clamp lifted'),
 ('unbiased_uncapped','Uncapped unbiased'),('trace_cap_lifted_budget_retained','Trace cap lifted'),
 ('budget_labels','Budget, labels'),('budget_label_free','Budget, label-free'),
 ('plain_pseudobulk_delta','Plain bulk-delta')]
table=['| Variant | C4 +.0791 | C3 +.0663 | C6* −.0670 | C5 −.0219 | C1 −.0007 |',
       '|---|---:|---:|---:|---:|---:|']
for key,name in names:
    cells=[]
    for c in ['C4','C3','C6','C5','C1']:
        q=r['contrasts'][key][c]
        cells.append('0' if q['delta_scaled']==0 else f'{q["ratio"]:+.4f}')
    table.append('| '+name+' | '+' | '.join(cells)+' |')
report='''# vcc_mse_instr_1026

**Verdict: instrument not found.** Every formulation fails C4’s sign and 0.3–3× ratio gate. No parameters were fitted. The 79 targets were used only to read contrasts.

Entries below are **local/official scaled-delta ratios**: positive means correct sign; zero means blind. The fixed ruler is the existing SEARCH anchor span, 0.95798309, without calibration to official outcomes.

'''+ '\n'.join(table)+'''

*C6 is only the 694-versus-885 budget surrogate. Calib3 contains no genuine v8/v9 layout pair; its “v8a” uses v9 code. Thus even a favorable surrogate sign cannot validate C6.

**Leave-one-contrast-out:** no fitted mapping, so held-out signs equal direct signs. In C1/C3/C4/C5 order, uncapped unbiased gives ✓/✗/✗/✗; every other listed variant gives ✗/✗/✗/✗. Genuine C6 is unavailable; only uncapped unbiased matches its surrogate. Shared package legs make these checks dependent.

**Saturation:** the real-side jackknife cap fires on 0/79 targets in every analog. B binds native v5b and both v9 analogs, making layout invisible; 1.5× amplitude clears B but increases plug-in error more than earned credit. Raw MSE rises 0.985006→0.986167; lifting the score clamp exposes C4=−0.001212 scaled. The [0,1] clamp otherwise maps both to zero because the local baseline is 0.961342.

Label-free control resamples at n_p=400 leave claims unchanged. Per-package B/claim/ρ, exact scorer locations, formulas, provenance and reproduction commands are in [DETAILS.md](DETAILS.md) and [budgets.csv](budgets.csv).

**Verification PASS:** independent installed-kernel replay matches retained rows within 2e−12; separate budget and rational aggregation checks pass. All 79 native v5b bulk reconstructions match cached moments within 8.89e−16. Missing amplitude/billm totals were reconstructed without integer cells or layout rerenders.

**A125/A150/1.38×:** no qualified instrument, so no trustworthy official MSE prediction is issued. Next prerequisite: a faithful v8 layout replay on a representative panel. No submission, download, package, promotion, commit, or external write occurred.
'''
assert len(report.split())<=400,len(report.split())
(j.ROOT/'REPORT.md').write_text(report)

details='''# MSE instrument diagnostic details

The main report is intentionally <=400 words. Full precision results: `results.json`; per-contrast ratios and fixed-formula leave-one-contrast-out records: `contrast_variants.csv`; prediction budgets: `budgets.csv`; independent verification: `verification.json`.

## Exact installed scorer trace

Installed `cell-eval2` 0.18.0 in `.venv/lib/python3.12/site-packages/cell_eval2/` (source hashes in `verification.json`):

- `moments.py:160`, `jackknife_correction`: group-sum CP50K log1p delete-one-cell correction. With group gene counts P, total S, cell count vector y_i and depth l_i, v_ig=log1p(50000*(P_g-y_ig)/(S-l_i)); C=(n-1)/n sum_ig(v_ig-mean_i v_ig)^2. This is a nonlinear bulk jackknife, not the variance of per-cell log expression. Chunk=32 changes working memory only.
- `metrics/delta.py:43`: PRED_TRACE_CAP_K=1. At `:699`, prediction credit becomes min(Chat,C_real). This cap never fires in these analogs. The metric value itself is signed and uncapped.
- `metrics/delta.py:276`, `_drop_on_target`: remove own-target squared-error coordinate; retain WHOLE-gene scalar sampling corrections. Dropping a target gene from C as well would silently change the definition. Row divisor is G-1 when resolved and G otherwise (all 79 labels resolve on this 18,077-gene panel).
- `metrics/delta.py:157`, `_across_pert_budget`: centered weighted sum of squares across predicted perturbations, each own-target coordinate omitted, with a separate weighted mean over retained rows for each gene. This is not an unweighted variance or a second moment around the control.
- `metrics/delta.py:832`: Q=sum_p w_p min(Chat_p,C_p). Lines 835–859 compute B/Q, and :860 scales every prediction credit by that ratio when below 1. Thus total credit=min(B,Q). The real correction is always subtracted separately. The warning threshold 0.7 is log-only. No hidden tolerance multiplier is present.
- `metrics/delta.py:1122`, `distance_unbiased`: denominator rows are (||real_p-real_control||²-C_p-C_control)/G_p, target coordinate omitted as above. Denominator never reads predictions.
- `catalog.py:345`: `expr_mse_unbiased_capped_norm` is derived from capped numerator and distance denominator, with ratio-of-sums aggregation. `run.py:_derived_value` at :1822 checks aligned target sets and divides summed rows; `aggregate_metrics` at :1952 invokes it. Negative individual denominator rows are not filtered or individually divided.
- `catalog.py:350` declares score clamps [0,1]. `score.py:795` delegates to `scoring.score_one` (:461); `scoring.py:561` and :575 apply the clamp to scaled scores, never to the raw normalized ratio. Calib3's `scale_local` also clips MSE below at zero. Here raw values exceed baseline 0.9613418220, so baseline/replicate scaling is negative and the floor erases all differences.

## Definitions and data provenance

`VARIANTS.md` froze the formulations before computing their contrast outcomes. No fitted parameters, context weights, official-contrast regressions, or 79-target tuning. Scaling uses the same fixed original anchor span for every D-normalized variant. Alternative variant-specific anchors or a fitted transfer factor would be different instruments and were not selected after inspecting the official outcomes.

Reference: `rsi/artifacts/tasks/H1_rich_search.h5ad`, held-out 79-target real portion, 400 cells per target; input controls are separate (3200 cells). Target identity and task-file hash are recorded in `reference.json`; gene axis and numerical moments in `reference.npz`. D and every retained denominator row were checked. Label-free null C0 is mean jackknife over 16 bootstrap samples WITH replacement at n_p=400 from input controls, seed 20261026. The score still needs labels for its residual/real subtraction; only the null cap, claim and rho are label-free. The averaged C0 exceeds every analog's maximum Chat, so this null cap is inactive. No bootstrap count or seed was selected using contrast outcomes.

Original package analogs are the existing `final_push/artifacts/de_check/H1_rich_search/{v5b,v9a_885,v9a_694,v8b_like,levers_knob_none,incumbent}/` records used by calib3. No .vcc files or DE caches were read or written. Their original prediction hashes are carried in `results.json`.

All 9551 retained `rsi/artifacts/stat_cache/*.npz` files were inspected. All 79 v5b mean/jackknife records match their stored uncapped score identities at <1e-13, and their original baseline digest equals final_push's v5b digest. These means are also the means of v9a885 and v9a694 because v9 preserves group totals. No amplitude/billm records matched. For those two only, reconstructed the frozen source profiles with the stored alpha 0.4442702372827524 (no amplitude probe/refit), and the exact prefix of `mean_counts` plus `repair_bulk_totals`' largest-remainder sum rule. Skipped integer cell allocation and all v9 layouts. This exact sum reconstruction first passed on all 79 cached v5b groups at <8.89e-16. B and self-claim for the reconstructed amplitude/billm arms then match independently recorded values to <1e-12.

The original uncapped numerator and recovered exact means imply Chat_p=||pred-real||²-C_p-G_p*N_uncapped,p. These are recovered sufficient statistics, not fitted parameters. Their independently recorded total self-claim and per-target capped numerators are verification checks. Native billm and billm-v9 share means but retain separate Chat vectors.

## Package budgets

All values are in the scorer's gene-averaged, across-target sum units, not per-target means. Q is the capped claim, J is the self-claim. In this panel Q=Q0=J for every analog.

| Package | B | J=Q=Q0 | rho=rho0 | Trace cap hits |
|---|---:|---:|---:|---:|
'''
for arm,p in r['packages'].items():
    details+=f'| {arm} | {p["B"]:.12g} | {p["self_claim"]:.12g} | {p["rho_labels"]:.12g} | {p["cap_fires"]}/79 |\n'
details+='\nD='+str(r['D'])+'; local scaling span='+str(r['span_local'])+'.\n'
details+='''
## Contrast deltas, ratios and LOO

Columns below show signed scaled delta followed by local/official ratio. C6* is a budget surrogate only. Positive ratios agree with official signs. No parameter is learned from any contrast, so leaving each out gives precisely these same predicted signs; this is not an independent generalization test or a claim of five independent contrasts. To avoid hiding missing evidence, `genuine_local_analog=false` is recorded on every C6 row. Every variant fails C4 even without considering C6.

| Variant | C1 | C3 | C4 | C5 | C6* | LOO genuine C1/C3/C4/C5 | C6* sign |
|---|---:|---:|---:|---:|---:|---|---|
'''
for key,cs in r['contrasts'].items():
    cells=[f'{cs[c]["delta_scaled"]:+.8f} / {cs[c]["ratio"]:+.6g}' for c in ['C1','C3','C4','C5','C6']]
    signs='/'.join('PASS' if cs[c]['loo_sign_ok'] else 'FAIL' for c in ['C1','C3','C4','C5'])
    details+='| '+key+' | '+' | '.join(cells)+f' | {signs} | '+('PASS' if cs['C6']['sign_ok'] else 'FAIL')+' |\n'
d=v['C4_decomposition']
details+=f'''
C4 mechanism: plug-in error rises by {d['plugin_error_increase']:.15g}; prediction credit rises by {d['prediction_credit_increase']:.15g}; their difference is positive. Normalizing gives +{d['raw_normalized_increase']:.15g} raw MSE (worse). This directly explains why lifting the score clamp still mis-signs amplitude. The more dramatic uncapped improvement assigned to v9 layouts is credit unconstrained by B and reverses C3. Label-free versus labeled caps are inactive locally, so swapping their rho cannot repair the ranking.

## Verification and reproduction

`verify.py` invokes unmodified installed `mse_unbiased`, `mse_unbiased_capped`, `distance_unbiased` and `aggregate_metrics` on recovered moments. It checks every row against original tidy scores, recomputes B using an independently weighted per-gene mean, and sums/divides native float rows as exact rational numbers. All checks pass; this validates the diagnostic computations and reproduction, not hidden-panel predictive performance. Model/source hashes are in `bulk_reconstruction.json`; scorer hashes and all numerical residuals are in `verification.json`.

From this job directory, sequentially:

```sh
nice -n 10 python -B reference.py
nice -n 10 python -B recover_cache.py
nice -n 10 python -B reconstruct_bulk.py
nice -n 10 python -B compute_variants.py
nice -n 10 python -B verify.py
nice -n 10 python -B make_report.py
```

Each numerical script launches only after `df -g /` and >=31.5 GiB free, checks the 30 GiB hard floor during work, runs CPU nice 10 with BLAS threads one, and uses one worker. A Python audit hook blocks writes outside this directory and access to both prohibited paths. No raw per-cell predictions, packages, large copies, or downloads are retained; total job footprint stays under 50 MB. No worker is left running. There is no authorized/predictive reading for A125, A150 or 1.38x from a failed instrument, so the conditional projection step was not invoked.
'''
(j.ROOT/'DETAILS.md').write_text(details)
print('REPORT words',len(report.split()))
j.status('Final REPORT.md and detailed trace written. Verdict: instrument not found; conditional amplitude forecasts not issued.')
