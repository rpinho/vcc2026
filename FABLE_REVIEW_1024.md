# Fable cloud review of rpinho/vcc2026 (commit fa2546d), 2026-10-24 (read-only; claims inferred from code + scorer algebra, NOT yet tested)
Scorer cross-checked against ArcInstitute/cell-eval2 src/cell_eval2/metrics/delta.py. Scaled MSE ≈ (0.989 − u)/0.95; billmA 0.061 → u≈0.93 (barely beats control); leaders 0.44 → u≈0.57.

## Why MSE is low
- Amplitude is set by scorer mechanics: v5.solve_alpha (model_v5/v5.py:58) bisects to B = claim exactly; callers pass rho ≤ 1 (model_v6/pipeline.py:212, final_push/de_check.py:50, final_push/emitter/build.py:106). Official r=min(1,B/claim) (delta.py:988) makes B=claim the minimum for full credit, not an optimum. Only official probe (v8b ×1.5): MSE 0.016→0.095 (undershoot) but NMAE collapsed (−0.16).
- Claim is ~1/4 of what the reference allows: model_v5/v5.py:4-6 says pool_k=1 fixed this, yet shipped builders pool 4 real cells per emitted cell (model_v6/pipeline.py:149, model_v5/build_v5.py:74, model_billm_pack/renderer_preflight.py:28-36; rsi/freeze_v5a_loco.py:81 admits "stale pool_k1 prose"). v8a swing (0.082→0.016, identical totals) consistent with the v8 layout lowering the jackknife (shaped_emit/_drain_top flattens columns).
- Direction ceiling: bulk profile zeroed on genes without cross-source sign agreement (model_v6/pipeline.py:176-185, v5.de_call_mask min_net=1) and 50 control PCs projected out (pipeline.py:169-170); v9 layout now picks calls from `confidence`, so sparsification only costs MSE/PDS. `blend` (model_billm/network.py:28-35) keeps only the MLP's direction and restores the pooled norm although the MLP is MSE-trained in the scorer's space (model_billm/prepare.py:24).
- Centering is NOT the main loss (≲0.012 scaled); don't chase it.

## Ranked tests (each via final_push/de_check.py-style arms on H1_rich_search, six official metrics, record B/claim)
1. rho sweep 1.0/1.4/2.0 (existing --rho flag, model_v6/pipeline.py:282); replicate anchor's B/claim is 1.38–1.58. Pass = overall up, NMAE not below v5b.
2. pool_k=1 (or 2) template, alpha re-solved (model_v6/pipeline.py:149, renderer_preflight.py:28; model_v5/evaluate_v5.py:65-68 has _ipf_p variants). Watch moment_certificate L1 gates and NMAE.
3. Re-solve alpha AFTER the blend: model_billmB/panel_pipeline.py:48-56 (same in final_push/main/panel_pipeline.py:84-92) solves B pre-blend, rotates 0.3 toward learned, never rechecks; if B(new) < claim, r<1 and refund leaks. One-line check.
4. Dense bulk for MSE/PDS, sparse only for calls: skip line 185, try deflate(..., eta=0); pre-screen run_selector.py recipes base vs gate.
5. Use the learned norm in blend (network.py:35) / per-target amplitude from `learned` before solve_alpha; check pred_over_true_norm, magnitude_spearman, rho_bar_oracle_amp (final_push/metrics.py) on H1/CD4 LOCO folds, then official arm.

## Flags
pool_k docstring/code contradiction; unchecked post-blend B; final_push/metrics.py:57 est_scaled_mse_fullcredit hard-codes a≈0.705 and ignores a²−c; model_v6/model.py:329-333 hard-codes `/ .6`; dangling doc refs in the public repo (docs/vcc2026-metrics-brief.md, consult_1021/FABLE.md, FINAL_PUSH.md, PUSH_1021.md cited in final_push/metric_payoff.py:4, metrics.py:5, main/panel_pipeline.py:36, calib_compare.py:10).
Uncertainty: artifacts excluded from the public repo, so gain sizes are inferred from metric algebra + the single v8b datapoint; NMAE is the binding risk for tests 1, 2, 5.
