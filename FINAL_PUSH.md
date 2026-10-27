# FINAL_PUSH.md — running log (builder: Claude Fable, from 2026-10-21 15:10 UTC)

Job directory: `final_push/`. Historical research and build checkpoints; machine supervision and transient storage observations are omitted.

## Checkpoints

- 2026-10-21 15:10 UTC — read consults (FABLE/ASTRA), PUSH_1021, OFFICIAL_SCORES, CONSTRAINTS, _codex/RULES.md, v5/v7/v8/v9/billm notes,
  v6 run_final.sh, metric brief. State: billmB source acquisition was running; `vcc datasets list` 15:07 UTC: controls-final still "unavailable — only available during the final week".
  Plan: (1) selector = weighted Frobenius cosine rho-bar in log1p(CP50K) bulk space on leave-one-family-out public folds
  (banks in billm format: model_billm/data + model_billmB/data shards), calibrated on v5-direction vs billmA-direction (H1/CD4 LOCO
  checkpoints exist); (2) emitter; (3) biology knobs on the selector; (4) final pipeline on mock D/E/F, safety package when
  controls-final lands.

- 2026-10-21 15:40 UTC — **SELECTOR BUILT AND CALIBRATED** (`final_push/{banks,metrics,recipes,run_selector,vdir,magnitude}.py`).
  Folds: leave-one-source-family-out over the 11 public screens we hold in billm format (HCT116, HEK293T, CD4 x3, H1,
  K562 GWPS, HepG2, Jurkat, RPE1, K562), 300-target panels, truth = that screen's own measured log1p(CP50K) response.
  Metrics: rho_bar (weighted Frobenius cosine, own target excluded, official 1/|G_p| weights), **rho_bar_centred**
  (both matrices centred across targets first), official pds_cosine, oracle-amplitude rho_bar, magnitude Spearman.
  * **rho_bar_centred is the selector**; plain rho_bar is not. The across-target MEAN response is a large shared
    component (truth_mean_share .02-.25 by context) and the official baseline IS that mean response, so energy spent
    on it buys rho_bar locally but nothing officially. Measured: dropping source centring ("nocentre") lifts rho_bar
    .1200 -> .1463 while rho_bar_centred barely moves (.1245 -> .1261) and pds falls .781 -> .751. Keep centring.
  * **Calibration vs official history** (`artifacts/vdir1.csv`, real model_v6 code path, matched 4 folds H1+CD4x3):
    rho_bar_centred v4a .0395 < v5b-gate .0455 < billmA-blend .0510, which is the official order
    (v4a .1438 < v5b .1528 < billmA .16438) — **3/3 pairs correct**. Plain rho_bar gets 2/3 (it puts billmA under v5b);
    local pds gets v4a/v5b right but billmA WRONG (local -.005 vs official +.026). So: trust rho_bar_centred for
    direction/magnitude; do not trust local pds at the .02 level; n = 3 directions, this is weak evidence.
  * **Negative result that kills a planned lever:** full-screen source-wide centring vs the 437-target panel cut is a
    wash (rho_bar_centred .1245 vs .1255). The "source-wide centred banks" idea is not the +.02-.05 that was hoped.
  * Small positives, pds-neutral: expression-compatibility masking +.0023, confidence-weighted amplitude +.0027.
    Hard silence gate is NEGATIVE (-.001 rho, -.014 pds): a zeroed row scores pds exactly .5. Retesting with a floor.
  * Magnitude (`artifacts/magnitude_features.csv`): predicted norm correlates with true norm .20-.68 by context;
    the destination's own control CPM of the target gene correlates with the ORACLE amplitude on every fold (.12-.36),
    so the silence idea is real but must shrink, not zero. A LOCO ridge on (pred norm, source norm, n_src, agreement,
    control CPM) does NOT beat the raw predicted norm on 7 of 11 folds — no magnitude model is adopted.

- 2026-10-21 16:05 UTC — **BIOLOGY KNOBS TESTED ON THE SELECTOR** (`artifacts/selector_round2.csv`). Best combination
  found: v4 source weights + expression-compatibility mask (control CPM >= 1) + confidence-weighted per-target
  amplitude, optionally a FLOORED silence gate. On the seven validation-like folds rho_bar_centred .0667 -> .0719
  (+7.8% relative), pds .719 -> .727; on all eleven folds .1245 -> .1291 (+3.7%). Floored silence recovers the pds
  the hard gate lost (.767 -> .788) and is worth +.000-.002. Honest size: ~+5% relative on the direction quality,
  which by the break-even-credit model is ~+5% of the MSE and PDS components, i.e. **order +0.001-0.002 overall**.
  The biology levers we can reach are small; the emitter is the larger one.

- 2026-10-21 16:10 UTC — **EMITTER BUILT AND KNOB-TESTED** (`final_push/emitter/{emitter,build_emitter}.py`, v7
  numerics copied unchanged plus `vendor/`). Independent dual-moment heads (bulk pseudobulk vs per-cell CPM, fitted
  jointly inside the depth-feasible bounds by AtlasShift `dual_moment_counts`, feasible projection applied BEFORE the
  budget solve). Switchable knobs: `scatter` (pool4 | pool1 | selective), `claim_rho`, `claim_source`, `amp_alloc`,
  `silence`, `compat_cpm`, `cpm_amplitude`, `call_budget`, `pvalue_order`; `Knobs.off()` reproduces v7a.
  **HARD GATE** `emitter.ratio_gate`: the realised B/claim must be in [0.9, 1.05] or the build RAISES (the slim_v5b
  lesson — that package shipped at 6.5 because its gate had been demoted to a diagnostic).
  Fixture on the real context A (24-target sub-panel, `artifacts/fix_*_fixture.json`, 16-18 s each):
  | scatter | own probe claim | control-400 cap | emitted B | B/claim | share of cap |
  |---|---:|---:|---:|---:|---:|
  | pool4 (v5/v7a/billmA incumbent) | 0.00881 | 0.01847 | 0.00855 | 0.970 PASS | **48%** |
  | selective (v7b rule) | 0.01865 | 0.01847 | 0.01773 | 0.970 PASS | **99%** |
  So the under-dispersed 4-cell pooled template, not the budget solve, is what caps our credit at ~half, exactly as
  model_v5/NOTES.md line 110 suspected; selective scatter doubles the credit with the ratio gate still passing and
  the CPM mean head untouched. This is lever 1's mechanism, now measured. It is NOT an official gain: credit pays
  only through the cross term, and v8b showed an amplitude move can wreck NMAE. **Tonight's v7a/v7b pair is exactly
  this contrast (claim .11-.13 vs .23-.25) and is the only official read we will ever get on it.**

- 2026-10-21 16:15 UTC — **FINAL CODE PATH VERIFIED ON D/E/F TODAY** (`final_push/artifacts/smoke/`). 12-target panel,
  mock contexts D/E/F, v6 arm: build 69 s -> installed `vcc prep` 0.2.2 exit 0 -> exact decoded CSR/axes equality
  PASS, 14,400 x 18,533, nnz 60,305,585, package sha256 f5ef9b13...cfa. Outputs deleted afterwards (regenerable).
  **Operational finding that must become a gate:** the mock panel has zero overlap with the 437-target bank cut, so
  every target was coverage-filled, alpha saturated at the solver's hi=8 bound and B/claim came out at **0.029**. A
  final run without the panel bank fetch would therefore produce a valid-but-worthless package. Pre-submit gates:
  (1) targets_with_source_row >= 0.9 x panel, (2) alpha strictly < 8, (3) B/claim in [0.9, 1.05], (4) `vcc prep`
  exit 0 + independent decode, (5) renderer/full-group receipts for the billm path.

- 2026-10-21 16:25 UTC — **STATUS / WHAT IS NOT DONE.** No package was built or submitted; `controls-final` was still
  unavailable at 15:46 UTC. billmB untouched (9 shards fetched, still running); the RSI loop untouched.
  Implemented and tested: selector (+ calibration), biology knobs, emitter scatter/claim/allocation/silence/compat
  knobs with the hard B/claim gate, final code path on D/E/F, `final_push/GATES.md`.
  **Declared but NOT implemented** in the emitter: `call_budget` and `pvalue_order` are hooks only. Calls and
  p-value ordering belong to the v9 renderer (`model_v9/layout.py`, `call_budget.py`), which billmB's
  `panel_pipeline.py` already makes panel-agnostic; the right integration is emitter -> profile/scatter/claim,
  v9 layout -> calls/ordering. Replacing v8's fixed per-context NCONF with a per-target budget ~ predicted n_real
  (the fidelity hole: fidelity is k/max(n_pred, n_real)) is a small change in `call_budget.py` that is NOT yet wired
  into any built package.
  Next, in order: (1) read tonight's v7a/v7b official result — it decides whether the main package uses selective
  scatter; (2) the moment `controls-final` lands, fetch panel banks and build the SAFETY package (v6 arm) against
  GATES.md; (3) wire the per-target call budget into the v9 layout and build the main package with the chosen
  emitter + best direction (v4 weights + compat mask + confidence allocation + floored silence); (4) billmB vs
  billmA on the selector when billmB's weights land.

- 2026-10-21 17:05 UTC — **TWO-TIER RECIPE TESTED AND REJECTED (as a bulk rule); the DE-side jackpot is not reachable.**
  Built: `final_push/strength.py` (target-strength predictor + leave-one-FAMILY-out validation),
  `source_signal.py` (per-(source,target) noise-corrected signal cache, 10,781 rows),
  `de_count.py` (official-rule reference DE counts on real cells), `metric_payoff.py` + `artifacts/tier_payoff.json`,
  and recipes `tier/tier_best/tier_hard/tier_soft/tier_sim/respgene/tier_resp` in `recipes.py` + `run_selector.py`.

  1. **Strength predictor (LOFO).** Macro Spearman with held-out true signal **0.380**, which is NOT better than its
     single best feature alone (max source signal strength, **0.389**); pooled-prediction norm gets 0.302. Silent-vs-not
     AUC 0.643 macro. Crucially it is weakest exactly where we need it: CD4 rho .14 / AUC .53, HEK293T .16 / .56,
     HCT116 .23 / .57, versus HepG2 .64 / .78 and Jurkat .51 / .72 — i.e. strength rank transfers between similar
     lines and barely transfers to dissimilar ones, which reproduces the coordinator's finding and means a NEW cell
     line (D/E/F) is the bad case. Caveat on the "~30% silent" numbers: the b/n floor estimator sets b at the 30th
     percentile of e.n, so a ~30% silent fraction is partly DEFINITIONAL, not measured.
  2. **Two-tier on the selector** (`artifacts/selector_tier1.csv`): it does lift direction quality — rho_bar_centred
     +4.6% on all folds, **+10.1% on the seven validation-like folds** (the best single direction gain I have found) —
     but it **costs 0.067 raw pds**, because the official pds_cosine scores a zero-effect row at EXACTLY 0.5 while our
     per-target pds is ~0.72. Converting both with the published anchors (pds span 0.4555; MSE via scaled ~ 1.41 rho_bar
     applied as a relative change to billmA's official 0.061):

     | recipe | d rho_centred | d pds raw | d MSE scaled | d pds scaled | **d overall** |
     |---|---:|---:|---:|---:|---:|
     | tier_best (all folds) | +4.6% | -0.067 | +0.003 | -0.147 | **-0.024** |
     | tier_best (validation-like) | +10.1% | -0.045 | +0.006 | -0.098 | **-0.015** |
     | tier_soft (gentler 20/80, weak 0.8) | +3.1% | -0.040 | +0.002 | -0.087 | **-0.014** |

     **Reject.** The pds loss is mechanical and certain; the MSE gain is small and model-dependent. Silencing the bulk
     response is strictly worse under the 6-metric average, by roughly an order of magnitude.
  3. **Metric mechanics: the empty-reference jackpot does not exist at the official power.** Reference DE counts
     recomputed on real cells with the official rules (400 cells, thinned to 20,000 UMI median, Wilcoxon two-sided
     normal approximation with tie correction, BH within perturbation, 5-CPM control filter), `artifacts/de_count_*.csv`:
     * H1_rich_search (400 cells, official count): **n_real == 0 for 0 of 79 targets**, n_real < 10 for 34%, median 107, p90 2,609.
     * CD4_depth192_search (192 cells, lower power): n_real == 0 for 2 of 79 (2.5%), median 30.
     * The official brief corroborates independently: the fidelity cohort is **295-300 of 300** and the reach cohort
       290-300, so at most ~1.7% of official targets have an empty reference set; the Jaccard cohort is always 300.
     Arithmetic: a correctly silenced empty-reference target is worth +2.70 scaled Jaccard (empty union = 1), but a
     wrongly silenced responding target costs -1.77 scaled fidelity (n_pred = 0 => F_p = 0), so the silenced set needs
     **39.5% precision** to break even against a base rate of 0-2.5%. With AUC 0.643 that is unreachable. Even at
     PERFECT identification the ceiling is 0.017 x 2.70 / 6 = **+0.008 overall**, and we cannot get near it.
     What the cohorts DO say: the NMAE cohort is 209-261 of 300, so 13-30% of targets have reference gates under 10
     genes and are dropped from NMAE entirely - amplitude error on weak targets is already free on NMAE, which is
     another reason the tier rule buys little.
  4. **Responsive-gene mask: no effect.** rho_bar_centred -0.5%, pds +0.001 (`respgene`, genes above the within-source
     median in >= 2 screens). The expression-compatibility mask already in the best recipe captures what is there.

  **Net recommendation: do not adopt the two-tier rule or the responsive-gene mask.** The current best direction
  (v4 weights + compatibility mask + confidence-weighted amplitude + floored silence) stands. The one piece worth
  keeping from this line of work is the confirmation that per-target amplitude allocation helps while per-target
  SILENCING hurts: shrink weak targets, never zero them — which is exactly what the floored silence gate already does.

- 2026-10-22 00:50 UTC — **FINAL PHASE OPEN. controls-final downloaded.** CLI upgraded uv tool vcc-cli 0.2.2 -> **0.3.0**
  0.3.0 `vcc prep` has no `--force`; removed it from `model_v6/prep.py` (backup `prep.py.bak-0.2.2`); billmB's
  `panel_prep.py` never used it. Bundle at `data/final/vcc_2026_controls_final.zip` (613.6 MiB, crc32c verified)
  unpacked to `data/final/bundle/`: **panel `vcc2026-test-1`, partition test, contexts D/E/F, 300 targets (one panel
  for all three contexts), 18,400 control cells per context (46 NTC ids), 400 cells per perturbation, ground-truth
  138,400 cells per context, gene axis 18,533 = byte-identical to the validation gene_names.csv.** Controls are
  float32-encoded counts, 100-104M nnz per context.
  **Panel coverage:** 129/300 targets overlap the validation panel, 134/300 are in the 437-target bank cut, and
  **161/300 have NO donor row in billmA's current banks** (100 have 4 families, 20 have 3). So G1 (>= 270 covered)
  requires the panel-agnostic donor fetch; started 00:48 UTC via `model_billmB/fetch_panel.py` (unchanged
  `model_v6/banks.py`, K562 GWPS + HCT116 + HEK293T + CD4 + H1 + promoters) into `model_billmB/final/banks`.
  **billmB verdict (Codex job, complete 23:13 UTC):** weights returned and verified, but frozen-SEARCH acceptance
  FAILED on CD4 (macro -0.0346 vs billmA -0.0276) -> no billmB package; **billmA weights remain the incumbent**.
  Its final-panel mock with billmA weights PASSED all 900 groups + prep + decode (probes 255 s, full 4,217 s, prep
  376 s), so the incumbent safety path is `model_billmB/RUNBOOK_FINAL.md` steps 2-3 with billmA weights; fallback is
  the v6 arm. The queued v7a/v7b validation pair was refused (phase flipped): there will be no official MSE read of
  the selective-scatter emitter, so the SAFETY package uses the incumbent emitter unchanged.

- 2026-10-22 01:05 UTC — **vcc 0.3.0 prep VERIFIED on the real D/E/F controls** (12 final targets, v6 arm; `final_push/artifacts/smoke/smoke030_*.json`):
  `vcc prep` exit 0, "targets verified against the official list", contexts D/E/F = 4,800 cells each, counts-preserved,
  exact decoded CSR/axes equality; meta now carries `cli_version 0.3.0`, `contexts`, `perts_fingerprint`. Outputs deleted.
  **SAFETY build launched** 00:47 UTC: `final_push/run_safety.sh final_billmA` (waits for the donor fetch, then
  runbook steps 2-3 with **billmA** weights -> `model_billmB/out/final_billmA.vcc`; timeline in
  `final_push/logs/safety_timeline.txt`; gates via `final_push/check_gates.py final_billmA --emitted`).
  ETA from the mock timings (fetch 1-3 h incl. pauses; responses+coverage min; probes 255 s; full 4,217 s; prep+decode
  376 s): **package ~04:00-06:00 UTC 10-22**, later if execution pauses or the HEK293T/CD4 ranges run slow.

- 2026-10-22 02:40 UTC — **Donor fetch complete (00:48 -> 02:30 UTC, 1 h 42 min, no execution pauses):** panel banks in
  `model_billmB/final/banks` cover HCT116 298, HEK293T 298, K562 GWPS 288, CD4 290 (x3 conditions + DE), H1 134 of 300
  (receipt `banks_receipt_K562_GWPS_raw_HCT116_HEK293T_CD4_H1_promoters.json`). **Donor-family coverage of the final
  panel went from 161 targets with NO donor to 2** (225 targets have 4 families, 25 have 5, 40 have 3); G1 passes with
  298/300 >= 270. Responses + coverage 02:31, **probes PASS 02:35**, full 900-group build running (mock timing 4,217 s).

- 2026-10-22 04:40 UTC — **SAFETY FINAL PACKAGE READY: `model_billmB/out/final_billmA.vcc`** (2,595,584,000 bytes,
  sha256 ce78c81c...9dc4). Incumbent billmA recipe on the real D/E/F bundle through RUNBOOK_FINAL steps 2-3 (billmA
  weights; billmB rejected by its own acceptance). Timeline: fetch 1 h 42 min, probes 4 min, full build 1 h 49 min
  (900 groups at 7.3 s), prep+decode 9 min; **3 h 47 min release-to-package**, no execution pauses. All gates in
  GATES.md pass; emitted B/claim 0.51-0.55 is the incumbent's known v9-layout regime (documented asymmetry added to G3).
  Evidence + suggested name in `rsi/NEXT_SUBMISSIONS.md`. **STOPPED before submission; the parent submits.**
  Next (main package): decide between (a) incumbent + selector-best direction knobs (compat mask + confidence
  allocation + floored silence; +5-8% direction on the selector, ~+0.001-0.002 overall, low risk) and (b) the
  selective-scatter emitter (claims 99% vs 48% of the control-400 cap; largest mechanical lever, NO official read since
  the v7 pair was refused; NMAE risk per v8b). Either needs the parent's call; both reuse the fetched panel banks.

- 2026-10-22 05:05 UTC — LOCAL DE CHECK DONE (`final_push/de_check.py`, official cell-eval2 via the frozen RSI harness on
  H1_rich_search, 79 targets x 400 cells; same billmA LOCO biology; results in `final_push/artifacts/de_check/H1_rich_search/summary.json`).
  incumbent: overall +0.048, fid -0.115, jac -0.067, nmae -0.037, raw MSE 0.9835 (reproduces billmA's stored local score bitwise).
  sel_fixed (scatter only): overall -0.056, fid -0.635, jac -0.214, raw MSE unchanged.
  sel_recal (scatter + amplitude x1.62 to B = .97 x capped claim): overall -0.016, pds +0.13 scaled, raw MSE 0.9715, fid -0.613, jac -0.195.
  Reading: scatter destroys the call count (fidelity = k/max(n_pred, n_real)); extra claim is worthless while B < claim.
  Decision under the parent's rule 3: (b) rejected; main package = (a). CD4_depth192 not run (H1 at 400 cells is decisive).

## 2026-10-22 05:23 UTC — SAFETY ENTRY SUBMITTED (parent)
Entry BKgFJX4gmteiZB88HSLe "final billmA v9 layout", final status published, partition test, panel vcc2026-test-1.
Upload 33 min @1.2 MiB/s (plan ~40 min upload time for the main package). Builder resumed on: local DE check (A/B/C)
for selective scatter vs incumbent, then (a) knobs package, then (b) only if the DE check holds. Main entry by Nov 3.
- 2026-10-22 05:25 UTC — **(a) main candidate building** as `main_a2` from `final_push/main/` (job-local copy of the billmB
  panel pipeline; writes only under final_push/main; path fixes for model_v8/gates_code imports). Knobs: compat mask
  (control CPM >= 1, ~10.8-11.2k genes kept), confidence-weighted per-target amplitude (support count, mean 1; factor
  5th-95th pct 0.84-1.23, only ANG/DPY30/C5orf22 under 0.2), floored silence (0.25; 14/23/33 targets below 1).
  Profiles: alpha 1.016/1.020/1.060, B = claim, coverage 298/300. Probes PASS 05:21; full build running (7.3 s/group).
- 2026-10-22 07:35 UTC — **MAIN CANDIDATE (a) READY: `final_push/main/out/main_a2.vcc`** (sha256 82e39463...5987, all
  gates pass; emitted B/claim 0.51-0.55 = the documented v9-layout regime). Evidence + recommendation in
  `rsi/NEXT_SUBMISSIONS.md`. (b) rejected on the local DE check. The regenerable build output was removed.
  STOPPED before submission.
- 2026-10-22 08:10 UTC — **HARNESS CALIBRATION (task 1):** `final_push/de_calib.py` + `calib_compare.py` on H1_rich_search
  (`artifacts/de_check/H1_rich_search/calibration_summary.json`). v5b arm reproduces the frozen v5b local score exactly (0.03795).
  v9a-like (v5b + v9 layout 885) vs v5b: fid down, reach up, jac down, nmae flat = official v9a vs v5b on all four (8/8 with
  the 694-budget arm as v8a). v8b-like (x1.5 amplitude) vs v9a-like: pds up OK, fid/jac flat OK, but **nmae local +0.007 vs
  official -0.161 (MISS)** and reach flat vs +0.016 (MISS). 10/12 signs; 2 of 3 contrasts ordered. Verdict: usable for
  call-set levers (budget, ordering); BLIND to amplitude overshoot on NMAE, so any DE-amplitude-head verdict from it is
  not evidence on NMAE. Local fidelity magnitudes are exaggerated (H1 median n_real 107 vs budgets 694-885).

## 2026-10-22 08:13 UTC — MAIN ENTRY main_a2 SUBMITTED (parent)
Entry LIREtFnPJLkMXZprp7Q9, published; now the last entry. (b) rejected on the H1_rich_search DE check. Builder: calibrate the harness on v5b/v8a/v8b DE metrics, then single-factor DE levers.
- 2026-10-22 09:00 UTC — **DE-side levers on main_a2 biology (task 2, H1 harness):** call budget is monotone locally
  (fid/reach/jac: 450 < 650 < 885; base_885 overall -0.027, 650 -0.045, 450 -0.067) -> budget stays 885, no lever.
  Independent CPM-head arms (0.6/0.9/1.2) could not be emitted: the public dual-moment kernel fails its tolerance at the
  4-cell pooled template without v7's feasible projection; and the harness is blind to NMAE overshoot (calibration),
  so no trustworthy local verdict exists for that lever -> not built. PIE ordering arm running (leaky on H1: upper bound).
  Open question for the LIVE entry: the main_a2-like base scores pds 0.467 vs billmA's 0.530 on this task; knob
  ablation (none/compat/conf/silence) queued to attribute it.
- 2026-10-22 09:15 UTC — PIE-ordered calls (leaky upper bound on H1): overall -0.031 vs base -0.027, reach 0.041 -> 0.021,
  jac -0.180 -> -0.186, fid flat -> dead lever even with the leak. Knob factors on H1_rich: conf factor < 0.2 for 0
  targets, < 0.5 for 2; silence < 1 for 3; compat keeps 12,708 / 18,077 genes. So the 0.467-vs-0.530 pds gap is not
  row scaling; the compat mask (30% of genes zeroed) is the remaining suspect. Ablation arms running.
- 2026-10-22 09:50 UTC — **KNOB ABLATION ON THE OFFICIAL SCORER (H1_rich_search, v9 layout 885, same biology otherwise):**
  none (= safety recipe) -0.0192 | compat -0.0193 (neutral) | conf -0.0269 (pds -0.052, nmae -0.009, fid +0.020, reach -0.006)
  | silence -0.0224 (pds -0.008, reach -0.012) | all three (= main_a2) -0.0268. **Locally main_a2 is 0.008 overall BELOW
  the safety recipe**, driven by the confidence allocation. Mechanism the selector could not see: the scorer measures
  each row's delta against ITS reference controls, so every row carries a shared input-vs-reference control offset;
  rows scaled below 1 let the offset dominate the cosine and lose pds rank (the same mechanism behind v8b's official
  pds +0.028 and sel_recal's local +0.13 when amplitude goes UP). Officially the offset is smaller (same control pool),
  so the loss should be smaller than local, but the sign is clear and there is no evidence of a gain.
  **Recommendation: make the safety package `final_billmA.vcc` the last entry again** (still on disk, sha ce78c81c...9dc4),
  unless the parent weighs the local fid +0.02 over the pds/nmae/reach losses. Lesson recorded: amplitude allocation must
  never shrink rows below 1 (boost-only), and the selector's log-bulk metrics miss emission/offset effects.

## 2026-10-22 09:50 UTC — DECISION: revert to the safety recipe (parent)
Ablation on the official scorer: main_a2 −0.008 vs safety locally (conf allocation shrinks rows below 1 → shared control offset dominates cosine → PDS loss). My 08:13 submit of main_a2 was premature (selector blind to control offset). Safety package queued for resubmission 2026-10-23 00:01:30 UTC. Rule: allocation knobs boost-only.

## 2026-10-23 ~01:37 UTC — FINAL-ENTRY HARNESS AUDIT: STOP, SAFETY STAYS LAST

Read the full log and all requested score/gate/runbook/prior-art documents.
Recomputed retained official-scorer SEARCH ranks with exact rational rank
arithmetic, independently checked with SciPy Spearman/Kendall. Original five:
rho .90, pairs 9/10 (v4a/v4b reversed), reproduces stored CLOSED gate. Expanded
seven with v5a/v5b: rho .785714 (11/14), pairs 18/21; reversals v4a/v4b,
v4b/v5a, v4b/v5b. H1 four layout/amplitude surrogates: rho 0, pairs 3/6.
Exact comparable v9a/v9b/slim_v5b/billmA replay is not established: v9a is an
885-budget surrogate, v9b/slim local replays absent, billmA local sanity score
uses native rather than final v9 emission. Anchor/reference hashes and score
arithmetic checked; no simple numerical bug found. No honest local-to-official
gain interval exists. Per user's step 1, no candidate experiments/builds run.

Read-only public refresh completed: board1,635, unchanged since22:50UTC;
em2.bio .308708 recipe still Coming soon (sitemap + three JS bundles checked).
Sidechain SER-29 .222582 now reports four->nine sources plus per-context bulk
strength, MSE .212->.299; settings omitted and transfer unverified. Public
GitHub update searches refreshed; no usable new public X recipe located.
Full report and receipts: final_push/harness_audit_20261023/REPORT.md.

Current last entry 5NeeoWn96lfm23hgExcm remains safety ce78c81c...9dc4;
no submission/upload/GPU/paid compute/HOLDOUT/model training. Appended clear
NO NEW SUBMISSION decision to rsi/NEXT_SUBMISSIONS.md; RSI gate/supervisor
unchanged. The rejected regenerable main_a2 package was removed; its sources/receipts and the safety package were retained.
