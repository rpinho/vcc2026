# DECISION REQUEST — Virtual Cell Challenge 2026, final phase

You are one of three independent advisors. Ricardo (the competitor) has explicitly asked that the
next move be decided by the models, not by him or by the orchestrating agent. The orchestrator has
written this brief to be neutral and states no preference. **Your job is to make the call.**

Answer in English. Be decisive and concrete. <=900 words. Structure your answer as:
  (A) VERDICT: one of [PURSUE A NEW PACKAGE | STOP AND SHIP CURRENT | HYBRID - specify]
  (B) WHAT TO DO, ranked, with concrete first action and expected hours
  (C) WHAT WOULD HAVE TO BE TRUE for your verdict to be wrong, and the cheapest test of it
  (D) RISK of your recommendation, stated plainly

## Hard constraints (not up for debate)

- Deadline **Nov 05 23:59 UTC**. Today is **Oct 25**. 11 days.
- **Only the last final entry counts.** No score is shown for final-phase entries.
- Current last entry `5NeeoWn96lfm23hgExcm`, package SHA `ce78c81c…9dc4`, status published/FINAL.
  Account reports `can_submit: true`, phase `final`, panel `vcc2026-test-1`, contexts D/E/F.
- A submission can be replaced, but replacing is a one-way bet: if the new package scores worse, the
  better one is gone (no revert to a prior final entry — the last one counts).
- Local CPU only. No GPU/pod budget. One part-time operator, no human in the loop overnight.
- The agent is forbidden from submitting; Ricardo executes any submission himself.

## Scoreboard facts

Our official **validation** scores (panel vcc2026-val-1, scores visible; NOT the final blind panel):

| entry | overall | PDS | MSE | NMAE | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|---:|
| **billmA** (best) | **0.16438** | 0.706 | 0.061 | 0.116 | -0.066 | 0.186 | -0.016 |
| v9a | 0.16135 | 0.680 | 0.083 | 0.108 | -0.068 | 0.181 | -0.016 |
| v5b | 0.15279 | 0.680 | 0.082 | 0.109 | -0.031 | 0.077 | -0.000 |
| v5a | 0.14990 | 0.678 | 0.078 | 0.094 | -0.020 | 0.071 | -0.001 |
| v4a | 0.14380 | 0.678 | 0 | 0.129 | -0.019 | 0.076 | -0.001 |
| v8a | 0.14700 | 0.680 | 0.016 | 0.102 | -0.064 | 0.162 | -0.015 |
| v8b | 0.14120 | 0.708 | 0.095 | -0.059 | -0.060 | 0.178 | -0.015 |
| slim_v5b | 0.08833 | 0.667 | 0 | -0.276 | -0.048 | 0.200 | -0.013 |
| v9b | 0.12049 | - | - | - | - | - | - |

Public reference points: leader Illumina/iCell **0.4716**; ranks 10-50 median **~0.31** (PDS 0.812,
MSE **0.444**, NMAE 0.244, fid 0.046, reach 0.273, jac 0.012); ranks 150-300 median ~0.22 (MSE 0.226).
AtlasShift 0.1546; Sidechain 0.1657 (MSE 0.151). Roughly 1,620 entries.

**The MSE gap dominates.** billmA 0.061 vs the 150-300 band's 0.226 and the 10-50 band's 0.444.

The current final entry is billmA's *recipe* re-run panel-agnostically on the real D/E/F bundle
(donor banks refetched for the 300 final targets, same source fusion / PC-50 deflation / cross-source
gate / B = 1.0x probe claim / 0.3 learned blend / v9 layout at 885 calls). All package gates
(G1-G3e, G6-G8) PASS; independent decode 900/900; closeout verifier audited PASS. Its blind score is
unknown. It is NOT identical to the scored 0.16438 entry — it is that recipe applied to new contexts.

## What has been tried and failed (all measured, all negative)

Since 10-21 (billmA) every attempted improvement has lost:
- v9b 0.12049 (worse than v9a 0.16135).
- slim_v5b 0.08833 (4x amplitude, MSE clamped to 0, NMAE -0.276).
- main_a2 (10-22, replacing the safety on D/E/F): scored **-0.008 vs safety** on the official-scorer
  ablation, reverted on 10-23.
- Local frozen-scorer ablations on a 79-target public H1 SEARCH panel: rho 1.4/2.0 (local +0.007/+0.016
  but violate the amplitude gate G3); pool_k=1 (+0.0017); post-blend alpha re-solve (+0.0011);
  dense bulk (drop sign gate + PC deflation) **-0.0111 overall, PDS -0.041, raw error worse**;
  eta0 (drop only deflation) running now; learned-norm substitution running now.
- Sidechain source/amplitude arms: all blocked, `transfer_noise_bound: null`.
- Family-quality H1 closure: NOT PROMOTED (budget gate failed, delta 0.051 > 0.020 tolerance).
- Refara control-reference audit, sigRecon comparator, dense-bulk, learned-norm: no candidate emerged.

## Three measured facts about the local harness (this is the crux)

1. **It is structurally blind to scaled MSE.** Its anchor pair is baseline u = 0.961342,
   replicate u = 0.003359, linear between and clamped at 0. Every one of the seven historical recipes
   sits *above* 0.961342, so all seven scale to MSE exactly 0.000000 — including v5a (0.985618) and
   v5b (0.985006), which officially earn 0.078 and 0.082. So the local harness has had **no signal**
   on the component that carries ~0.064 of the overall gap. Raw u (lower better) does order the
   officially non-tied pairs 11/11, but 8 of those 11 are "earns credit vs tied at zero".
2. **It misorders.** Excluding the v3b outlier, on the H1 panel it gets **6/15** recipe pairings
   backwards. The errors are directional: every inversion ranks an officially-better later recipe
   below an earlier worse one (e.g. v1 is locally 0.055 above v4b while officially 0.037 below).
   Local PDS for the v4b/v5a/v5b family is anti-correlated with official PDS.
3. **Our effect sizes are below its resolution.** Measured candidate gains ~0.002; smallest observed
   inversion 0.014. So local deltas of the size we can produce carry no ordering information.

Related: our primary local panel (`H1_rich_search`, 79 targets) overlaps the official 300-target
panel by only **8 targets (10.1%)** — the least representative of any fold we hold. Same-line source
evidence is worth ~4.8x the response correlation of cross-line evidence, and the shipped fusion sits
at cross-line level.

## Also true, and relevant

- A matched local/official inventory found **0 exact-matched pairs** out of 67 records, so no
  transfer-aware promotion rule can be derived from existing data. (Bounded scope, not proof of
  absence.)
- The anchors are *measurements*, not knobs. Re-deriving them to make our arms score positive would
  be fitting the yardstick to the answer. A legitimate re-anchor would need its own held-out
  construction and a check that it preserves the historical pairings.
- Rebuilding a better local MSE harness is real work (est. 6-12h) with no guarantee it produces a
  packagable gain in the remaining time, and its output would still need a leap of faith on transfer.

## The question

Given 11 days, no GPU, no official feedback on the final panel, a one-way replacement bet, a local
harness that is blind to the dominant error component and inverts small differences, and a history of
every recent attempt measuring worse than the incumbent:

**Should we (a) pursue a new package for the dominant MSE gap, (b) stop and ship the current entry,
or (c) something else?** If (a) or (c), say exactly what, in what order, and what you would abandon
first when it does not work. Consider explicitly whether any *high-variance* strategy is rational
here — e.g. a deliberately aggressive amplitude or layout bet that is likely worse but might be much
better, given that the current entry is already submitted and any replacement destroys it.
