# DECISION 2026-10-25 03:07 ET — removing the expression mask HURTS: the B=claim solve neutralises it

## First arm of the 2x2 (H1_rich_search, paired, same run, official anchors)
| arm | raw PDS | raw MSE norm | raw NMAE | alpha | B/claim | G2/G3e/G6/G7/G8 |
|---|---:|---:|---:|---:|---:|---|
| `incumbent` | 0.753651 | 0.983537 | 0.958685 | 0.4443 | 0.449 | all pass |
| `no_mask` | 0.722650 | 0.991230 | 0.962857 | 0.2611 | 0.428 | all pass |

Officially-anchored deltas for `no_mask`:
`PDS −0.0688, MSE −0.0081, NMAE −0.0069, fidelity +0.0047, reach −0.0258, Jaccard −0.0012`,
sum **−0.1060 scaled → −0.0177 official overall**. Not promoted; it fails the priced NMAE rule outright
(the NMAE move is small, but there is no MSE/PDS gain to pay for it).

## Why this is a mechanism result and not noise
The prediction was that unmasking restores suppressed response direction. **The opposite happened, and the
reason is visible in the frozen calibration rule.** `no_mask` uses a larger response surface, but the
label-free B=claim bisection solved **alpha 0.2611 instead of 0.4443** — it scales the whole profile down to
hit the same budget. So unmasking does not add usable amplitude; it adds structure that then gets divided
away, and what survives is noisier:
- PDS falls 0.031 raw (0.7537 → 0.7227): the masked profile was **better discriminated**, not suppressed.
- Raw MSE rises 0.0077 (worse) even though MSE is the component unmasking was supposed to help.
- Reach falls 0.0258 scaled: fewer of the added structure survives into confident DE calls.

## What this closes and what it does not
- It closes the claim that the **expression mask** is what suppresses direction in this pipeline. At fixed
  budget, the mask is protective: it removes structure that does not survive normalisation.
- It does **not** yet speak to the 50-PC deflation (`no_pc`) or the interaction (`unmasked`). Those are the
  remaining two arms; the pre-registered rule requires all four before the lever is closed.
- It strengthens the standing lesson from `decisions/20261025_0105`: at fixed B=claim, **the amplitude solve
  is a normaliser that silently absorbs profile changes**, so any "more structure" arm must be read at fixed
  alpha as well as fixed B, or the comparison measures the solver, not the biology.

## Next
`no_pc` is emitting; `unmasked` follows. Apply the frozen rule on receipt:
- If both remaining arms are also null/negative → close the expression-mask/deflation line and launch the
  queued budget-matched `family_quality` arm (`_codex/vcc_family_quality_matched_1025.md`).
- If `no_pc` (mask kept, PC removed) helps → the deflation, not the mask, is the suppressor, and that is a
  single, well-identified change worth a second recipient.
No package change, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.
