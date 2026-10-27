# MSE instrument diagnostic details

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
| v5b | 0.00981690755298 | 0.00982278485981 | 0.999401665932 | 0/79 |
| v9a_885 | 0.00981690755298 | 0.0200042421295 | 0.490741288244 | 0/79 |
| v9a_694 | 0.00981690755298 | 0.0194499834255 | 0.504725754167 | 0/79 |
| v8b_like | 0.0220735742559 | 0.0199237280683 | 1 | 0/79 |
| levers_knob_none | 0.00899557574224 | 0.0200358861306 | 0.448973191582 | 0/79 |
| incumbent | 0.00899557574224 | 0.00988457747962 | 0.910061736153 | 0/79 |

D=0.1590316336099593; local scaling span=0.9579830875970241.

## Contrast deltas, ratios and LOO

Columns below show signed scaled delta followed by local/official ratio. C6* is a budget surrogate only. Positive ratios agree with official signs. No parameter is learned from any contrast, so leaving each out gives precisely these same predicted signs; this is not an independent generalization test or a claim of five independent contrasts. To avoid hiding missing evidence, `genuine_local_analog=false` is recorded on every C6 row. Every variant fails C4 even without considering C6.

| Variant | C1 | C3 | C4 | C5 | C6* | LOO genuine C1/C3/C4/C5 | C6* sign |
|---|---:|---:|---:|---:|---:|---|---|
| capped_clamp_lifted | +0.00000000 / -0 | +0.00000000 / +0 | -0.00121238 / -0.0153272 | +0.00153337 / -0.070017 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| unbiased_uncapped | -0.06682956 / +95.4708 | -0.06319149 / -0.953114 | -0.06444245 / -0.814696 | +0.00713218 / -0.32567 | -0.00363807 / +0.0542996 | PASS/FAIL/FAIL/FAIL | PASS |
| trace_cap_lifted_budget_retained | +0.00000000 / -0 | +0.00000000 / +0 | -0.00121238 / -0.0153272 | +0.00153337 / -0.070017 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| budget_labels | +0.00000000 / -0 | +0.00000000 / +0 | -0.00121238 / -0.0153272 | +0.00153337 / -0.070017 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| budget_label_free | +0.00000000 / -0 | +0.00000000 / +0 | -0.00121238 / -0.0153272 | +0.00153337 / -0.070017 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| label_free_rho_labeled_claim | +0.00000000 / -0 | +0.00000000 / +0 | -0.00121238 / -0.0153272 | +0.00153337 / -0.070017 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| plain_pseudobulk_delta | +0.00000000 / -0 | +0.00000000 / +0 | -0.06780155 / -0.857162 | +0.00684806 / -0.312697 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| plain_expression | +0.00000000 / -0 | +0.00000000 / +0 | -0.06755204 / -0.854008 | +0.00692447 / -0.316186 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |
| current_clamped | +0.00000000 / -0 | +0.00000000 / +0 | +0.00000000 / +0 | +0.00000000 / -0 | +0.00000000 / -0 | FAIL/FAIL/FAIL/FAIL | FAIL |

C4 mechanism: plug-in error rises by 0.0102915265782831; prediction credit rises by 0.0101068205152891; their difference is positive. Normalizing gives +0.00116144227913106 raw MSE (worse). This directly explains why lifting the score clamp still mis-signs amplitude. The more dramatic uncapped improvement assigned to v9 layouts is credit unconstrained by B and reverses C3. Label-free versus labeled caps are inactive locally, so swapping their rho cannot repair the ranking.

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
