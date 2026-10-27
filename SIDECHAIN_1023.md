# Sidechain bulk-strength and source-breadth test — 2026-10-23

Updated 2026-10-25T00:56:39.791934+00:00. Status: bounded diagnostic tests complete.

**Verdict: no candidate qualifies for a final package.** Historical ordering fails and no local-to-official noise bound is established. These results cannot show that a candidate beats the incumbent with the required margin. No combination of source and amplitude factors, refit, GPU request, upload or submission is justified.

No new package exists. The recorded current last entry remains **5NeeoWn96lfm23hgExcm**, package `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`. This job has made no changes to that entry.

The required `../RULES.md` and root `GATES.md` are absent; their applicable copies are `../_codex/RULES.md` and `final_push/GATES.md`. Prior art was recorded before this work in PRIOR_ART.md (source reference: `sidechain_1023/PRIOR_ART.md`; not included). This tests disclosed ideas; it does not reproduce SER-29.

## Reference prerequisite

The [Refara guide](https://github.com/hz424/Refara/blob/main/docs/VCC_AUDIT.md) was used only for reference roles, calibration bindings and complete-panel scoring. Model-input controls are disjoint from scoring cells; observed and predicted effects use the same measured scoring controls (`control_source=real`). The frozen scorer and measured from-replicate anchors were retained. Each 79-target panel enters global PDS/MSE together. Distinct cell IDs do not establish independent acquisitions.

| Retained history | Spearman | Correct pairs | Gate |
|---|---:|---:|---|
| Five recipes | 0.9 | 9/10 | FAIL |
| Seven recipes | 11/14 | 18/21 | FAIL |
| Four H1 layout/amplitude surrogates | 0 | 3/6 | FAIL |

Exact recent quartet coverage remains incomplete. The retained H1 incumbent replay is bitwise identical to the earlier same-adaptation result; numerical repeatability is not transfer uncertainty. See audit.json (source reference: `sidechain_1023/artifacts/audit.json`; not included), per-task reference receipts and verification.json (source reference: `sidechain_1023/artifacts/verification.json`; not included).

## Fixed experiment contract

Frozen billmA LOCO weights, original four-cell template, v9 layout and generic 885-call budget. Bulk arms change only the multiplier on log bulk effects about input controls. Source arms retain incumbent amplitude and weights. These are public SEARCH adaptations, not exact validation-package replays or hidden-final scores. Every scored arm emits 79 × 400 integer cells and records independent spread/jackknife calculations and moment certificates. Both measured-task and input-control caps are checked against the hard emitted upper bound 1.05. No over-spread arm can become eligible.

The 1.05× follow-up was added after the wide sweep showed a narrow band passing the requested ratio check. This is an adaptive diagnostic sweep; no final-context setting is inferred from it.

Moment certificates in this SEARCH adaptation compare the v9 output to the corresponding pre-layout integer emission. The original emitter separately asserts a converged bulk fit. These receipts do not replace final-panel requested-head certificates or full renderer verification.

The requested-ratio lower bound and emitted-ratio upper bound have different roles: v9 can increase claim after preserving integer totals. The retained CD4 incumbent requested ratio is 0.88598, below 0.90; that exception is reported as a failed check, not silently waived. Final-panel coverage, packaging, decoding and full renderer gates were not run because promotion is blocked.

## Global bulk multiplier

Equal-context macro deltas from the emitted incumbent. All values are scaled; higher is better.

| Multiplier | Δoverall | ΔPDS | ΔMSE | ΔNMAE | Δfidelity | Δreach | ΔJaccard | Requested / capped checks |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1.0 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | FAIL / PASS |
| 0.75 | +0.019390 | -0.046112 | +0.000000 | -0.002454 | +0.201231 | -0.036835 | +0.000509 | FAIL / PASS |
| 1.025 | -0.000401 | +0.005424 | +0.000000 | -0.002576 | +0.005050 | -0.010130 | -0.000174 | PASS / PASS |
| 1.05 | -0.005916 | +0.007625 | +0.000000 | +0.000349 | -0.030322 | -0.013261 | +0.000114 | PASS / PASS |
| 1.125 | -0.010553 | +0.017459 | +0.000000 | -0.008698 | -0.041105 | -0.031265 | +0.000290 | FAIL / PASS |
| 1.25 | -0.010728 | +0.030847 | +0.000000 | -0.009907 | -0.083074 | -0.002925 | +0.000688 | FAIL / PASS |
| 1.5 | -0.010293 | +0.055935 | +0.000000 | -0.008677 | -0.096077 | -0.013345 | +0.000409 | FAIL / FAIL |

Other completed SEARCH count checks (G2, G6, G7 and conditional JK dominance): all PASS.

## Per-context multiplier grid

Each cell is the equal-context Δoverall from two already-scored count emissions. Rows are H1 multipliers; columns are CD4 multipliers. A dagger (†) marks a failed requested-ratio or emitted-count check. All cells fail the historical-order/noise gate. This is a descriptive reuse of the sweep, not independent validation or a rule fitted for D/E/F.

| H1 / CD4 | 0.75 | 1.0 | 1.025 | 1.05 | 1.125 | 1.25 | 1.5 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.75 | +0.019390† | -0.007852† | -0.009136† | -0.014588† | -0.021269† | -0.024471† | -0.027208† |
| 1.0 | +0.027242† | +0.000000† | -0.001284 | -0.006737 | -0.013417† | -0.016619† | -0.019356† |
| 1.025 | +0.028125† | +0.000883† | -0.000401 | -0.005853 | -0.012534† | -0.015736† | -0.018473† |
| 1.05 | +0.028062† | +0.000821† | -0.000463 | -0.005916 | -0.012597† | -0.015798† | -0.018535† |
| 1.125 | +0.030106† | +0.002864† | +0.001580† | -0.003872† | -0.010553† | -0.013755† | -0.016492† |
| 1.25 | +0.033132† | +0.005890† | +0.004606† | -0.000846† | -0.007527† | -0.010728† | -0.013465† |
| 1.5 | +0.036305† | +0.009063† | +0.007779† | +0.002327† | -0.004354† | -0.007555† | -0.010293† |

Largest completed grid delta passing these numerical checks: H1 1.025×, CD4 1.025×, Δoverall -0.000401. Promotion remains blocked.

## Per-context and standalone source results

Deltas use the same task’s emitted incumbent. MSE scaling can clip at zero; raw metrics and raw deltas are retained in summary.json (source reference: `sidechain_1023/artifacts/summary.json`; not included). `Req` is the requested B/claim check; `Cap` is the hard emitted upper bound with both caps. Every row is diagnostic and fails the overall promotion gate.

### H1_rich_search

| Arm | Overall | Δoverall | ΔPDS | ΔMSE | ΔNMAE | Δfidelity | Δreach | ΔJaccard | B/capped | Req / Cap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| base | -0.019164 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | 0.448973 | PASS / PASS |
| bulk_0.75 | -0.034868 | -0.015704 | -0.087797 | +0.000000 | -0.008016 | +0.008928 | -0.007523 | +0.000185 | 0.251738 | FAIL / PASS |
| bulk_1.025 | -0.017398 | +0.001766 | +0.010847 | +0.000000 | -0.000403 | -0.001674 | +0.001843 | -0.000016 | 0.471806 | PASS / PASS |
| bulk_1.05 | -0.017523 | +0.001641 | +0.016271 | +0.000000 | +0.000341 | -0.001943 | -0.004761 | -0.000060 | 0.495122 | PASS / PASS |
| bulk_1.125 | -0.013436 | +0.005728 | +0.034237 | +0.000000 | +0.002573 | -0.002745 | +0.000432 | -0.000127 | 0.568695 | FAIL / PASS |
| bulk_1.25 | -0.007383 | +0.011781 | +0.061695 | +0.000000 | +0.004074 | -0.007365 | +0.012510 | -0.000230 | 0.702516 | FAIL / PASS |
| bulk_1.5 | -0.001037 | +0.018127 | +0.110847 | +0.000000 | +0.009162 | -0.012648 | +0.001466 | -0.000068 | 1.009841 | FAIL / PASS |
| without_cd4 | -0.019210 | -0.000046 | -0.001356 | +0.000000 | -0.001093 | -0.000557 | +0.003130 | -0.000397 | 0.453549 | PASS / PASS |
| source_feng | -0.018563 | +0.000601 | +0.001356 | +0.000000 | +0.000250 | +0.000247 | +0.001978 | -0.000227 | 0.448707 | PASS / PASS |

Descriptive within-task maximum: `bulk_1.5`, Δoverall +0.018127. This uses the sweep data and is not a validated per-context rule. No D/E/F setting was fitted.

Raw error metrics (lower is better), retained because scaled MSE can clip at zero:

| Arm | Raw normalized capped MSE | Δraw MSE | Raw NMAE | Δraw NMAE |
|---|---:|---:|---:|---:|
| base | 0.983537089 | +0.000000000 | 0.958685448 | +0.000000000 |
| bulk_0.75 | 0.989800813 | +0.006263724 | 0.963554192 | +0.004868744 |
| bulk_1.025 | 0.982961444 | -0.000575645 | 0.958930492 | +0.000245044 |
| bulk_1.05 | 0.982392067 | -0.001145022 | 0.958478476 | -0.000206972 |
| bulk_1.125 | 0.980717363 | -0.002819726 | 0.957122552 | -0.001562896 |
| bulk_1.25 | 0.978050443 | -0.005486646 | 0.956210832 | -0.002474615 |
| bulk_1.5 | 0.974549786 | -0.008987302 | 0.953120603 | -0.005564845 |
| without_cd4 | 0.983007497 | -0.000529592 | 0.959349432 | +0.000663984 |
| source_feng | 0.983408685 | -0.000128404 | 0.958533894 | -0.000151553 |

### CD4_depth192_search

| Arm | Overall | Δoverall | ΔPDS | ΔMSE | ΔNMAE | Δfidelity | Δreach | ΔJaccard | B/capped | Req / Cap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| base | -0.862876 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000000 | 0.639857 | FAIL / PASS |
| bulk_0.75 | -0.808393 | +0.054483 | -0.004427 | +0.000000 | +0.003108 | +0.393533 | -0.066147 | +0.000832 | 0.357094 | FAIL / PASS |
| bulk_1.025 | -0.865444 | -0.002568 | +0.000000 | +0.000000 | -0.004748 | +0.011774 | -0.022104 | -0.000333 | 0.672461 | PASS / PASS |
| bulk_1.05 | -0.876349 | -0.013473 | -0.001022 | +0.000000 | +0.000358 | -0.058701 | -0.021761 | +0.000288 | 0.705893 | PASS / PASS |
| bulk_1.125 | -0.889711 | -0.026835 | +0.000681 | +0.000000 | -0.019969 | -0.079465 | -0.062962 | +0.000706 | 0.810736 | FAIL / PASS |
| bulk_1.25 | -0.896114 | -0.033238 | +0.000000 | +0.000000 | -0.023889 | -0.158783 | -0.018360 | +0.001607 | 1.000431 | FAIL / PASS |
| bulk_1.5 | -0.901588 | -0.038712 | +0.001022 | +0.000000 | -0.026516 | -0.179505 | -0.028156 | +0.000885 | 1.438206 | FAIL / FAIL |
| source_h1 | -0.863309 | -0.000433 | +0.000000 | +0.000000 | +0.001393 | -0.003990 | +0.000000 | -0.000001 | 0.639674 | FAIL / PASS |
| source_pan | -0.863970 | -0.001094 | +0.000681 | +0.000000 | +0.000662 | +0.011668 | -0.019616 | +0.000043 | 0.639080 | FAIL / PASS |
| source_feng | -0.863469 | -0.000593 | +0.000681 | +0.000000 | -0.001361 | -0.002879 | +0.000000 | -0.000002 | 0.639696 | FAIL / PASS |

Descriptive within-task maximum: `bulk_0.75`, Δoverall +0.054483. This uses the sweep data and is not a validated per-context rule. No D/E/F setting was fitted.

Raw error metrics (lower is better), retained because scaled MSE can clip at zero:

| Arm | Raw normalized capped MSE | Δraw MSE | Raw NMAE | Δraw NMAE |
|---|---:|---:|---:|---:|
| base | 1.119290374 | +0.000000000 | 0.935981026 | +0.000000000 |
| bulk_0.75 | 1.112913643 | -0.006376731 | 0.935110383 | -0.000870644 |
| bulk_1.025 | 1.120039392 | +0.000749018 | 0.937310872 | +0.001329846 |
| bulk_1.05 | 1.120840450 | +0.001550076 | 0.935880721 | -0.000100306 |
| bulk_1.125 | 1.123442346 | +0.004151971 | 0.941574482 | +0.005593455 |
| bulk_1.25 | 1.128692107 | +0.009401732 | 0.942672491 | +0.006691465 |
| bulk_1.5 | 1.273960050 | +0.154669676 | 0.943408295 | +0.007427268 |
| source_h1 | 1.119353442 | +0.000063068 | 0.935590724 | -0.000390302 |
| source_pan | 1.119338687 | +0.000048313 | 0.935795471 | -0.000185555 |
| source_feng | 1.119305381 | +0.000015006 | 0.936362223 | +0.000381197 |

## Source-family scope and acceptance

| Factor | Test scope | Gate result |
|---|---|---|
| CD4 | Already present. `without_cd4` removes only learned-feature CD4 evidence on H1; the frozen v5b head remains. CD4 stays excluded on CD4. | Limited inclusion contrast; not a new-family test; promotion blocked |
| Existing H1 row breadth | `source_h1` appends one disjoint row relevant to CD4 SEARCH; old masks, means and weights remain fixed. | Separate from second-screen result; promotion blocked |
| Independent second H1 screen | `source_pan`: streamed Pan2026 lanes, assigned within-lane controls, native library depths, ≥30 source cells; shares the existing H1 family vote. H1 destination excluded. | Inference-only comparison; no training refit; promotion blocked |
| Pooled stem-cell evidence | `source_feng`: eligible published adjusted-coefficient directions only, fixed 1/8 vote with incumbent row norm restored. Magnitudes discarded; no invented raw source moments. Cached raw KOLF banks have zero SEARCH overlap. | Limited direction adaptation; not exact Sidechain source reproduction; promotion blocked |
| Calu-3 | No compatible verified local source bank; retained metadata has ambiguous controls and a different condition. | UNTESTED; not a negative result |
| Glioblastoma | Native control assignment and raw-count/readout preservation remain unverified. | UNTESTED; not a negative result |
| Combination | None: no factor passes the required evidence gate alone. | NOT BUILT |

Source eligibility and provenance: SOURCE_PLAN.md (source reference: `sidechain_1023/SOURCE_PLAN.md`; not included) and source_inventory.json (source reference: `sidechain_1023/artifacts/source_inventory.json`; not included). Missing-family experiments remain incomplete; they are not evidence that those source families hurt.

Observed standalone source deltas in scaled overall (positive means a local improvement; metric tradeoffs remain in the six-metric tables):

| Source contrast | Context | Δoverall | Local direction | Promotion |
|---|---|---:|---|---|
| CD4 learned-feature inclusion (incumbent minus removal) | H1_rich_search | +0.000046 | helps locally | BLOCKED |
| source_feng | H1_rich_search | +0.000601 | helps locally | BLOCKED |
| source_h1 | CD4_depth192_search | -0.000433 | hurts locally | BLOCKED |
| source_pan | CD4_depth192_search | -0.001094 | hurts locally | BLOCKED |
| source_feng | CD4_depth192_search | -0.000593 | hurts locally | BLOCKED |

## Reproduction and retention

From the repository root, completed scores are reused:

```sh
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/run.py --task H1_rich_search
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/run.py --task CD4_depth192_search
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/run.py --task H1_rich_search --arms without_cd4 --scratch-tag cd4_removal
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/run.py --task CD4_depth192_search --arms source_h1 --scratch-tag h1_rows
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/acquire_pan.py
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/sources.py --task H1_rich_search --source feng
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/sources.py --task CD4_depth192_search --source feng
nice -n 10 ../.venv/bin/python -Bu sidechain_1023/sources.py --task CD4_depth192_search --source pan
nice -n 10 ../.venv/bin/python -B sidechain_1023/verify.py
nice -n 10 ../.venv/bin/python -B sidechain_1023/report.py --complete
```

Exact raw/scaled scores, per-target count hashes, moment/jackknife certificates and source receipts are under `sidechain_1023/artifacts/`; logs are under `sidechain_1023/logs/`. The separate verifier checks exact stored-number scaling with rational arithmetic, group identities, emitted checks and the closed promotion decision. It does not independently recompute the official raw metrics.


The next useful evidence is representative native SEARCH data and faithful recent historical replays. No upload command is supplied because no package meets the required gate.

