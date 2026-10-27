**The biggest immediate correction is the calendar.** The current [official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/#_5-submit) says final-phase submissions reject validation contexts; final scores are hidden. Therefore two daily *validation* experiments through November 5 are not an available resource. The queued October 22 v7 pair works only if the server still accepts validation at 00:01. Check phase before uploading; do not bypass rejection using old packaging. [Arc’s announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026) gives October 22 release and November 5 **23:59 UTC** deadline. Keep November 4 as the internal deadline. This materially changes the proposed research strategy.

**What MSE actually rewards.** Below, `M` denotes `competitors/code/cell-eval2/docs/vcc2026_metrics/vcc2026-metrics.md`; `D` denotes `competitors/code/cell-eval2/src/cell_eval2/metrics/delta.py`. These are the checked-in official implementation/docs, not inferred leaderboard mechanics.

For perturbation p, sum counts P, then form b=log(1+50,000P/ΣP). Let wₚ=1/Gₚ, with the own-target coordinate excluded from distances. Define:

\[
E=\sum_p w_p\|\hat b_p-b_p\|^2,\quad
J=\sum_p w_p C_p,\quad
Q=\sum_p w_p\min(\hat C_p,C_p),
\]
\[
B=\sum_g\sum_{p\in R_g}w_p(\hat b_{pg}-\overline{\hat b}^{\,w}_{g})^2,
\quad T=\sum_p w_p(\|b_p-b_0\|^2-C_p-C_0).
\]

The raw score is **u=(E−J−min(B,Q))/T**. C is the delete-one-cell jackknife of the log-normalized pooled profile, including the changing library-size denominator—not ordinary per-gene cell variance. The scalar corrections retain the own-target contribution even though distances exclude it. Sources: `M:190–264`; `D:698–701,835–882`; `moments.py:160–175` in the same source directory.

Each context scales as s=clip((b_anchor−u)/(b_anchor−r_anchor),0,1), then contexts/components are averaged. Using illustrative documented anchors .99/.04, .061 scaled means u≈.932; .444 means u≈.568. Closing that gap requires reducing the corrected numerator by approximately **.364T**, worth **.064 overall**. These anchors are historical examples, not recovered r4 constants (`M:77–97,293–315`). The component median is not a real competitor (`CONSTRAINTS.md:20–31`).

Matching B/Q≈1 is useful but insufficient. For fixed claim Q and a centered predicted direction scaled by a, the numerator has form K+a²V−2aH−min(a²V,Q). Below the boundary the quadratic cancels; improvement depends on alignment H with truth. Above it, excess amplitude costs quadratically. Thus .95 is a heuristic near a kink, not a universal optimum. Random spread buys approximately its own error back; it supplies no missing biological alignment. Shared offsets cost error without increasing B.

Our evidence supports **both an emission bottleneck and response error**, not an identified allocation of the .38 gap. Pool4 suppresses usable jackknife credit; coupled CPM/bulk heads force compromises. v7b roughly doubles estimated claim versus v7a (.228–.249 versus .114–.130), but also changes geometry and amplitude (`model_v7/NOTES.md:7–18,77–86`). PDS .706 measures relative identification, not calibrated distances. It cannot establish good MSE.

Use three controlled outputs: accurate pooled profiles, independently calibrated arithmetic-CPM means, and gene-selective marginal scatter. For cell composition x and depth L, pooled composition minus arithmetic composition equals Cov(L,x)/mean(L). Equal depths eliminate that freedom; differing depths permit it only within feasible bounds. This is why independent heads require a constrained dual-moment emitter, not merely two requested vectors (`model_v7/NOTES.md:62`).

Pin realized totals; fit both moments within those bounds; preserve power on confident DE genes while allowing uncertainty elsewhere. Measure emitted B and every group’s jackknife, not eight probes alone. Control-derived caps remain estimates of hidden treated caps. At fixed profiles, increasing every predicted C cannot worsen MSE; v9 already implements that certificate (`model_v9/NOTES.md:18–20`). Before clipping, the gain is Δmin(B,Q)/[T(b_anchor−r_anchor)]; there is no benefit beyond saturation. v5b→v8a’s .066 scaled loss corresponds to about .063T under the illustrative anchors—not an 80% loss of jackknife credit.

Protect NMAE by preserving arithmetic-CPM fold changes; it scores reference-significant genes whether we call them or not. Protect fidelity by preserving directional yield: **k/max(n_pred,n_real)**. A 300-call ceiling with n_real=800 caps raw fidelity at .375 even with perfect signs—below the documented ~.51 baseline. Reach instead rewards the deepest ≥90%-correct prefix of reference-significant genes, including uncalled ones. Jaccard rewards set agreement. These are different objectives (`M:334–416,458–489,523–591`).

**What stronger teams plausibly do.** Disclosures identify independent DE/bulk heads, moment-preserving count transport, selective dispersion, confidence encoded into p-value order, and broader source preparation. Bill’s disclosed .350 MSE with .697 PDS is especially informative: discrimination near ours can coexist with much better error. Our reconstruction explicitly lacks his independent DE head (`competitors/LEVERS-2026-10-18.md:9–18`; `model_billm/PACK_REPORT.md:33`). Broader family-balanced screens and full-pool centering likely improve transferable magnitudes and signs. These are hypotheses supported by descriptions, not recovered winning pipelines: the top-500 audit reproduced **zero complete entrants** (`../vcc_top500_1019/REPORT.md:3–5,23–36,48–54`).

**Priority by expected official gain per workday.** These are decision ranges, not forecasts; gains overlap and can be negative.

1. Dual moments/selective scatter on frozen billmA biology: 1–2 days, plausible +.005–.03 overall; largest quick opportunity, substantial downside.
2. Finish billmB and final-panel donor coverage: 1–2 days, plausible validation +.003–.01; essential final insurance, not optional experimentation.
3. Independent DE amplitude and source-informed ranking: ~1 day each, plausible 0–.01 each; evaluate yield, not precision alone.
4. Full-source centering/reliability filtering: 1–2 days, plausible 0–.01; reuse expanded banks.
5. New architectures or massive atlases: poor expected return before deadline.

**Day-by-day tests and final operations.** The following two columns specify single-factor contrasts against that day’s frozen parent; all other predictions, seeds and settings stay fixed. They are an official-slot contingency **only if validation remains available contrary to the guide**. Otherwise run bounded external-context diagnostics, prioritize rows 23–26, and do not manufacture official evidence from them. Avoid building thirty full packages. October 22’s already-committed pair is explicitly a multifactor exception, not causal isolation.

| UTC date | Contrast 1: hypothesis → metric | Contrast 2: hypothesis → metric |
|---|---|---|
| Oct 22 | Queued v7a: native geometry → PDS/MSE | Queued v7b: selective scatter → MSE; confounded |
| Oct 23 | billmB weights versus A, same emitter → PDS/NMAE | Independent CPM head, same bulk → NMAE |
| Oct 24 | Pool2 versus pool4 on uncertain genes, fixed moments → MSE | Pool1 versus pool4 on uncertain genes → MSE |
| Oct 25 | Bulk amplitude ×.85, fixed CPM → MSE | Bulk amplitude ×1.15, fixed CPM → MSE |
| Oct 26 | DE amplitude ×.85, fixed bulk → NMAE | DE amplitude ×1.15, fixed bulk → NMAE |
| Oct 27 | Target-centered PIE ordering only → reach | Family-consensus ordering only → fidelity/reach |
| Oct 28 | Call budget 450 only → fidelity | Call budget 800 only → yield/reach |
| Oct 29 | Full-pool source centering only → MSE | Reliability weights only → NMAE |
| Oct 30 | Learned-direction mix .15 versus .30 → MSE | Mix .45 versus .30 → PDS/MSE |
| Oct 31 | Mild source-trained responder attenuation → NMAE | Fixed 50:50 A/B profile ensemble → MSE |
| Nov 1 | Freeze; repeat seed only → stability | Leave-one-source-family-out → robustness |
| Nov 2 | Final-context control split 1 → moment stability | Split 2 → stability; no truth-side score |
| Nov 3 | Full finalist render/audit | Independent decode/axes/depth verification |
| Nov 4 | Submit verified final candidate | Reserve for failed-upload recovery |
| Nov 5 | Verify latest final entry/status | Reserve for correctness repair only |

Under the documented blind phase, leave submission slots unused except a safety entry after controls arrive and the selected replacement by November 4. Fetch the actual new panel’s donors immediately; dry-run generic D/E/F handling. Keep billmA+v9 as the incumbent; adopt B only with credible source-held-out robustness and no gross metric failures. Do not let unreliable proxy scalar wins select a complicated emitter. Freeze selection November 1. Stream one numerical worker, respect execution pauses, and use the already-budgeted ≤$20 GPU fit. Final receipt must bind package hash, panel, contexts and latest successful entry.

**Errors to stop repeating.** “billmB does not touch MSE” is false: changing direction changes E and B. The v8 loss does not refute min(B,Q): a modest credit loss can erase most of a small baseline-relative score. v9 diagnosed the lost jackknife (`PUSH_1021.md:119`; `model_v9/NOTES.md:18–20`). Drop gain4 emission, fixed 250–300-call dogma, and the claim that two confounded v7 failures would disprove emission. Stop the no-progress RSI cycle through its owner; retain its harness, not its all-pairs ranking veto (`PUSH_1021.md:78–108`). Drop further broad code hunts, scratch foundation models, and incompatible-summary conversions. The failed GenBio/GRAFT/SLIM adaptations do not isolate architecture: baseline and candidate emitters differed (`../vcc_adapt_1020/REPORT.md:23`). Package correctness is necessary; it is not biological validation.
