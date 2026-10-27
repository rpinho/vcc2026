> Historical local H1 campaign report. “Official” means the official scorer
> run locally with SEARCH anchors; these are not leaderboard scores.
> See [the afternoon/evening update](../RESEARCH_UPDATE_20261026b.md) and
> [export limits](../EXPORT_20261026b.md) for later corrections and missing inputs.

# REJECT — permanently close manifold/STRING

**Independent verification: PASS.** One out-of-sample H1 SEARCH experiment completed: train 168 H1 targets, predict and score 79 disjoint H1 targets. No submission, upload, package build/promotion, leaderboard interaction, network access/download, git mutation, other job, or external write.

B0 shipped billmA; B1 F1 STRING-only; B2 F3 STRING + projected donor direction; B3 inherited family-balanced pooled donor branch; B4 training mean; B5 50:50 shipped/training-mean response blend. All are same-run official cell_eval2 vcc2026 CPU/pdex scores with identical 79 targets, 18,077 genes, legal input-control template, renderer rule, budget 885, and local SEARCH anchors. MSE/NMAE lower is better. These are local SEARCH scores, not leaderboard or hidden-panel scores.

| Raw component | B0 | B1 | B2 | B3 | B4 | B5 |
|---|---:|---:|---:|---:|---:|---:|
| Expression MSE | 1.410715392 | 0.951898566 | 0.966530194 | 1.409413878 | 1.022118317 | 0.951298015 |
| PDS cosine | 0.850048685 | 0.554690036 | 0.655955858 | 0.858649789 | 0.500000000 | 0.676079195 |
| LFC NMAE | 0.928575728 | 0.899174641 | 0.869510082 | 0.921898388 | 0.933297856 | 0.919199267 |
| Direction fidelity | 0.448724814 | 0.472286050 | 0.467647072 | 0.449450631 | 0.467290106 | 0.455338246 |
| Direction reach | 0.126573077 | 0.137138102 | 0.150616095 | 0.118586087 | 0.082284940 | 0.139809910 |
| Jaccard | 0.042958861 | 0.042311580 | 0.044833859 | 0.043501914 | 0.041685063 | 0.042109166 |

| Scaled component | B0 | B1 | B2 | B3 | B4 | B5 |
|---|---:|---:|---:|---:|---:|---:|
| Expression MSE | 0.000000000 | 0.009857435 | 0.000000000 | 0.000000000 | 0.000000000 | 0.010484326 |
| PDS cosine | 0.731186441 | 0.114237288 | 0.325762712 | 0.749152542 | 0.000000000 | 0.367796610 |
| LFC NMAE | 0.012424440 | 0.060831695 | 0.109672742 | 0.023418308 | 0.004649719 | 0.027862262 |
| Direction fidelity | -0.479534514 | -0.394112658 | -0.410931391 | -0.476903045 | -0.412225582 | -0.455557351 |
| Direction reach | 0.045579803 | 0.057264404 | 0.072170659 | 0.036746433 | -0.003401547 | 0.060219343 |
| Jaccard | -0.180327000 | -0.182240870 | -0.174783041 | -0.178721312 | -0.184093345 | -0.182839364 |
| score_avg | 0.021554862 | -0.055693784 | -0.013018053 | 0.025615488 | -0.099178459 | -0.028672362 |

B3 comparisons and B4/B5 versus B0 are descriptive; none can rescue a failed manifold gate.

Paired target 95% percentile CIs: 20,000 draws, seed 20261026. Each cell is delta [lower, upper]. Positive lower bounds are required on all three endpoints against B0 and both artifact controls.

| Comparison | Raw fidelity | Raw PDS | Scaled score_avg | All three pass |
|---|---:|---:|---:|---|
| B1 minus B0 | +0.023561 [+0.014167, +0.032503] | -0.295359 [-0.376339, -0.213892] | -0.077249 [-0.125428, -0.026041] | False |
| B2 minus B0 | +0.018922 [+0.010487, +0.027588] | -0.194093 [-0.276534, -0.113113] | -0.034573 [-0.082797, +0.019196] | False |
| B1 minus B4 | +0.004996 [+0.001408, +0.008948] | +0.054690 [+0.032132, +0.079033] | +0.043485 [+0.023807, +0.066852] | True |
| B1 minus B5 | +0.016948 [+0.010563, +0.023453] | -0.121389 [-0.171860, -0.072055] | -0.027021 [-0.060065, +0.007519] | False |
| B2 minus B4 | +0.000357 [-0.006890, +0.007796] | +0.155956 [+0.103862, +0.209834] | +0.086160 [+0.056102, +0.121377] | False |
| B2 minus B5 | +0.012309 [+0.004716, +0.020195] | -0.020123 [-0.085524, +0.045440] | +0.015654 [-0.022785, +0.058202] | False |
| B1 minus B3 | +0.022835 [+0.013192, +0.031937] | -0.303960 [-0.383970, -0.223462] | -0.081309 [-0.129011, -0.030831] | False |
| B2 minus B3 | +0.018196 [+0.009493, +0.027212] | -0.202694 [-0.283836, -0.122687] | -0.038634 [-0.086733, +0.015364] | False |
| B4 minus B0 | +0.018565 [+0.009864, +0.026430] | -0.350049 [-0.425187, -0.274586] | -0.120733 [-0.158823, -0.082933] | False |
| B5 minus B0 | +0.006613 [-0.000044, +0.011929] | -0.173969 [-0.233690, -0.117332] | -0.050227 [-0.076281, -0.023691] | False |

Mechanical adoption gates: B1=False, B2=False. B2 does not clear either B4 or B5 on all three tests. The registered failure closes the manifold/STRING line permanently. No follow-up variant, rescue sweep, package, or promotion is authorized. A failure does not identify shrinkage as the sole biological mechanism; endpoint/control results above determine exactly what failed.

Leakage certificate: **PASS**. H1 train/test intersection 0; 168 response-basis/normalizer fitting targets; 0 official targets in fitting. Three training CV scopes are disjoint; all centering, bases, feature scales, ridge intercepts and k/alpha selection use training-only rows. No held-out truth was loaded by fit.py. All H1 donor banks excluded; other-family same-target donor measurements are legal inference features, not H1 held-out responses. Certificate is printed in logs/fit.log and retained in results/leakage_certificate.json.

Frozen hyperparameters: F1 k=64, alpha=1000, training CV SSE=11182.287737640, F3 k=32, alpha=100, training CV SSE=10702.720875423. Three-fold training-only rank/cell-count stratified CV selected summed expression SSE; final fits used all 168. B5 weight 0.5 was fixed in the preregistration before scoring. No arm or hyperparameter was selected using the 79.

Coverage: donor targets 79/79; STRING 77/79. KDR and SOX2 have no STRING embedding; both remain in every arm, with the inherited explicitly declared zero-feature fallback. Training missing STRING: KRT18/SALL4. Held-out truth has 62/146/372 cells min/median/max. Target- and matched-control CPM row-sum ranges [999999.9983398467, 1000000.0014823601] and [999999.9988807399, 1000000.0009364461]; each is explicitly normalized to one before log1p(50000*p)-log1p(50000*ctrl). Saved diagnostic truth reconstructs independently. Official six components score real held-out cells.

Gene alignment: bank 18080, manifold 18533 (18077 measured), scoring 18077. Every scoring gene maps by exact symbol to both inputs and is measured; no gene imputation. Excluded bank genes: ['TBCE-1', 'HSPA14-1', 'TMSB15B-1']. Gene map SHA256 `cc912b655d9f987e7d028e9deaaeeba0e756492273b9202d24e368d5c9c7b1b4`.

B0 reproduced the frozen shipped count digest and all six raw scores within 1e-12. Its raw direction fidelity 0.448724814 is **above its coverage-matched fair-coin expectation of approximately 0.40**. Fidelity is sign precision multiplied by a coverage penalty, so 0.5 is not its matched null. The 0.4929 historical validation value belongs to another panel and is not a chance threshold. B0's out-of-H1-training status rests on provenance, independently of its fidelity. Other-family same-target data are allowed in B0/B2/B3.

B0's scaled expression-MSE component is exactly **0.000000**: its raw MSE 1.410715392 is worse than the local anchor baseline 0.961341822, and this metric clamps at zero. B1's raw MSE 0.951898566 is better than that baseline: the manifold genuinely predicts expression better on this metric and still loses the composite.

B1 official component deltas versus B0 (composite points = scaled delta / 6):

| Component | Scaled B0 | Scaled B1 | Composite-point delta |
|---|---:|---:|---:|
| Expression MSE | +0.000000000 | +0.009857435 | +0.001642906 |
| PDS cosine | +0.731186441 | +0.114237288 | -0.102824859 |
| LFC NMAE | +0.012424440 | +0.060831695 | +0.008067876 |
| Direction fidelity | -0.479534514 | -0.394112658 | +0.014236976 |
| Direction reach | +0.045579803 | +0.057264404 | +0.001947433 |
| Jaccard | -0.180327000 | -0.182240870 | -0.000318978 |
| Total score_avg | +0.021554862 | -0.055693784 | -0.077248646 |

B1 gained +0.00164 expression-MSE points and lost -0.10282 PDS points, a 62.6x (approximately 63x) worse trade. Better expression MSE did not offset the loss of perturbation discrimination. These observations do not change any preregistered gate.

Saved response-space diagnostics are from these 79 targets only, not the earlier manifold cohort:

| Arm | Response cosine | Summed SSE |
|---|---:|---:|
| B0 | 0.045252873 | 6797.992567163 |
| B1 | 0.125553455 | 4921.177837516 |
| B2 | 0.138572124 | 5323.429965456 |
| B3 | 0.049414150 | 6770.929277297 |
| B4 | 0.092282147 | 5037.328316880 |
| B5 | 0.073545267 | 5288.838346350 |

Missing-STRING target official endpoints (these are included, not excluded subgroup tests):

| Arm / target | Fidelity | PDS |
|---|---:|---:|
| B0 / KDR | 0.5188323917137476 | 0.9615384615384616 |
| B0 / SOX2 | 0.09616074632221026 | 0.14102564102564108 |
| B1 / KDR | 0.5518672199170125 | 0.2948717948717948 |
| B1 / SOX2 | 0.09939002511661285 | 0.5 |
| B2 / KDR | 0.5349301397205589 | 0.39743589743589747 |
| B2 / SOX2 | 0.1072838177251525 | 0.9743589743589743 |
| B3 / KDR | 0.5268022181146026 | 0.9743589743589743 |
| B3 / SOX2 | 0.0966989594546107 | 0.1282051282051282 |
| B4 / KDR | 0.5493827160493827 | 0.2692307692307693 |
| B4 / SOX2 | 0.10154287764621457 | 0.5641025641025641 |
| B5 / KDR | 0.5186680121089808 | 0.6153846153846154 |
| B5 / SOX2 | 0.087908144958737 | 0.14102564102564108 |

Caveats: one H1 family and a previously used campaign panel; this tests target generalization, not a new cell line or independent batch. Bootstrap conditions on fixed fits, anchors and full-panel PDS retrieval contributions; it does not refit models or rerank duplicated retrieval competitors. Official finite metric cohorts and paired-finite MSE ratio-of-sums are preserved, with missing values saved; no zero-call targets are silently discarded. Percentile intervals have no multiplicity correction for the two prespecified candidates. Narrow local anchor spans and nonlinear clipping can magnify scaled changes; they do not forecast leaderboard points. B5 tests one predeclared 50% shrinkage weight, not an optimized shrinkage envelope. Rendering includes the frozen sign-dependent donor-confidence rule for every arm; B1 response fitting itself uses no incumbent features. KDR/SOX2 zero features are a coverage limitation. Independent verification replays saved model predictions and controls, axis/cohort identity, all per-target aggregation and bootstrap, but does not rerun all per-cell MWU or independently refit CV grids. Numeric agreement with exact binary64 gate decisions is not a formal real-arithmetic or biological certificate.

Operational interruption: worker 17892 exited 75 during B2 after 63/79 reductions at 33.240699 GB decimal free. Its inherited 31 GiB reserve was stricter than the requested 30 GB decimal. The first correction changed that operational threshold to 31 GB decimal, retaining the 30 GB hard guard. B0/B1 results and all six prediction files were hash-preserved; no model refit or scientific change occurred. B2 was deterministically replayed using the same predictions and exact DE caches, then the remaining arms ran sequentially. Original and revised policy/code receipts are retained and independently checked. results/input_pins.json stays byte-identical; its operational-resume receipt contains the six frozen prediction hashes. The verifier checks the exact recorded operational replacements against original pinned source snapshots.

Resource-floor breach: the second worker 32225 exited 75 after observing 28,284,035,072 bytes free, below the required 30 GB decimal. This transient breach is a constraint failure and must not be described as clean resource compliance. The initial correction to 31 GB decimal reserve was insufficient for the observed sudden multi-GB drop. The job initially escalated to a 35 GB launch / 34 GB stop hold. ORCHESTRATOR_AUTHORIZATION.md then explicitly superseded that self-imposed escalation with a 31 GB decimal launch / 30.5 GB decimal proactive stop. The resumed run retained the warm cache, scored one arm at a time, and checkpointed completed arms. No clean compliance claim is made for the earlier breach; resumed minimum sampled disk and launch receipts are retained separately. Predictions, scientific code, hyperparameters and the locked preregistration were preserved. Both interrupted B2 passes are operational replay, not additional fitted variants. See resource_hold.json and both resource_stop receipts.

Provenance: inherited manifold basis/ridge/features/split/pooled/embed and shipped renderer/scorer are credited in PRIOR_ART.md. No live literature/leaderboard checks were made because this task explicitly prohibits network/leaderboard use. No novelty claim. The predecessor’s zero intersection was real; the successful fit-on-168/predict-79 construction demonstrates why it was not an impossibility result. Original directories and incident artifacts were not modified or used as fitted models.

Resume/reproduce scoring from the frozen predictions only, in this directory, with >=31 GB decimal launch headroom and the retained warm DE cache. Each command runs after the preceding arm checkpoint is verified. Completed official arms are retained. Do not refit, reselect, or regenerate pins.

```sh
nice -n 10 python -B -u score_run.py B2
nice -n 10 python -B -u score_run.py B3
nice -n 10 python -B -u score_run.py B4
nice -n 10 python -B -u score_run.py B5
nice -n 10 python -B analyze.py
nice -n 10 python -B verify.py
nice -n 10 python -B make_report.py
```

For a read-only numerical recheck of the completed run, run verify.py alone; no model refit is needed. Do not reopen this panel as a tuning loop.

Preregistration SHA256 `488b3abc1209c2b3a8983fe0ad621dc306f2247be8d365caa3f6a934985d4f7d`; 112 pinned input files, 47 official package source hashes verified. Six predictions, fitted basis/weights, truth, target metric vectors, render count sums, all primary bootstrap draws and target identities are retained as NPZ. Per-target official DE tables are retained as Parquet. Input references are read-only links, rather than duplicate large data banks. Mandatory complete arrays take precedence over the older <50 MB general retention guideline; results size is logged in closeout. See results/verification.json and results/closeout.json for exact audit/resource/process receipts.
