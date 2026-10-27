# 2026-10-25 14:25Z — pair-ordering error rate of the local harness, per scoreboard

## Why

Both consultants (`decisions/consult_1025_dense/`) independently said the binding constraint is
not which lever to try next, but whether a local delta can justify replacing the package.
Astra: "the paired reversal disproves reliable ordering for this contrast; it does not
calibrate a correction." So measure how often the ordering is already wrong, and at what
delta size.

## Method

`final_push/harness_audit_20261023/audit.json` retains a local frozen-scorer overall and an
official validation overall for 10 recipes. Its four scoreboards are **not independent** and
must be reported separately:

- `by_task.H1_rich_search` — 7 recipes scored on the public H1 rich panel. **This is the panel
  every candidate we build is scored on tonight.**
- `by_task.CD4_depth192_search` — the same 7 recipes on a CD4 panel.
- `historical_seven_macro` — the average of the two above; a re-expression of them, not new data.
- `approximate_layout_amplitude_four` — 4 layout/amplitude surrogates, from a different run.

I enumerated every unordered pair within each board and asked whether the sign of the local
delta matches the sign of the official delta. Read-only; no scoring, no new data.

## Result

| scoreboard | pairs | sign-correct | excluding the v3b outlier |
|---|---:|---:|---:|
| **H1_rich (our panel)** | 21 | 15 (71%) | **9/15 (60%)** |
| CD4_depth192 | 21 | 18 (86%) | 12/15 (80%) |
| macro (H1,CD4) | 21 | 18 (86%) | 12/15 (80%) |
| layout/amplitude (4) | 6 | 3 (50%) | 3/6 (50%) |

The H1 row is the one that governs tonight's work, and it is the weakest: **6 of 15 pairs wrong
once the catastrophic v3b arm is set aside** — barely better than a coin flip.

Every wrong H1 pair, with the local and official deltas (a negative local delta with a positive
official delta means the harness ranked a better recipe lower):

| pair | local delta | official delta |
|---|---:|---:|
| v1 → v4b | **-0.0550** | **+0.0366** |
| v4a → v5a | **-0.0519** | **+0.0061** |
| v4a → v5b | **-0.0426** | **+0.0090** |
| v2b → v4b | **-0.0304** | **+0.0626** |
| v1 → v5a | **-0.0236** | **+0.0622** |
| v1 → v5b | **-0.0144** | **+0.0651** |

The pattern is systematic, not noise: every one of these six inversions ranks a *later,
officially better* recipe below an earlier, worse one. The harness does not merely add
symmetric noise around the truth; on H1 it carries a **pessimistic drift against the recipe
direction we have been improving in**. v1 is locally 0.014-0.055 *above* v5a/v5b, while
officially it is 0.062-0.065 *below* them.

Wrong pairs on the other boards are also all near-ties plus the two v4b inversions
(v4a→v4b local +0.147/official -0.020 and v4b→v5a local -0.122/official +0.026), which appear
only on CD4 and the macro where the CD4 panel's much larger local scale amplifies them.

## What this means for tonight

Every candidate gain we have measured or can plausibly measure tonight is **smaller than the
smallest H1 inversion (0.014)**:

- pool1 +0.0017, postblend +0.0011 (previous job);
- anything `final_push/dense_bulk_1025` produces will be in that range at best.

So a local +0.0017 cannot distinguish candidate from incumbent. Worse, because the H1 drift is
*directional* — biased against exactly the kind of later-recipe change we are testing — a small
negative local result here is not evidence against a candidate either. The panel's resolution
is coarser than the effects, and its bias points the wrong way.

This is a sharper statement than the existing `STOP_NO_TRUSTWORTHY_LOCAL_ORDER` verdict:
the harness is not uninformative, it is **misleading in a known direction at the scale we
operate on**.

## Corrections to earlier framing

- Earlier tonight I repeated "9/10 pairwise, 18/21 macro, 3/6 surrogates" as accuracy figures.
  Those mix overlapping boards and include v3b. The per-board table above is the honest
  summary and should be quoted from now on.
- My first draft of this record pooled all three sets into 27 pairs and reported 15/21. That
  still mixed a derived board (the macro) with its own inputs. The per-board table replaces it;
  the corrected H1 number is worse, not better.
- Astra's incumbent-vs-v5b pair (local +0.057 for v5b, official -0.0116) is arithmetically
  correct and is an additional inversion at a delta larger than anything we can produce.

## Consequence

1. No local-only result from `dense_bulk_1025` can justify a replacement, whatever its sign.
   A pass is a mechanism diagnostic, exactly as both consultants said; a fail is not a
   refutation either.
2. The only path to a *better package* before Nov 2 runs through an official-scored
   observation, not a local one.
3. The right use of the local harness from here is **mechanism diagnosis** — did removing the
   gate change the emitted DE sign structure, the amplitude, the call set — not ranking.
4. Do not spend further overnight hours on local arm sweeps hoping for a large delta. The
   measured effects are ~0.002 and the harness inverts at 0.014.

Consultation records: `decisions/consult_1025_dense/{brief,kimi,astra,astra_prompt}.md`,
consulted before the arm numbers were known.
