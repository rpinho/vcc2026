# FREEZE MEMO — VCC 2026 final entry

**Issued:** 2026-10-25 17:20Z
**Authority:** unanimous decision of Astra (`openai/gpt-6-astra`), Kimi (`moonshotai/kimi-k3.1-pro`) and
Fable (`anthropic/claude-fable-5.1`), 3-0 at 15:25Z, reaffirmed 3-0 at 15:40Z on new evidence.
Records: `decisions/20261025_1525_THREE_MODEL_DECISION.md`,
`decisions/20261025_1540_second_consultation.md`, `decisions/consult_1025_decide/`,
`decisions/consult_1025_harness/`.

## Rule

**Do not replace the final entry.** No speculative performance improvement justifies a swap. Ricardo
alone executes any submission; the agent is forbidden to submit.

## The frozen artifact

| field | value |
|---|---|
| entry | `5NeeoWn96lfm23hgExcm` |
| name | final billmA v9 layout r2 |
| package | `model_billmB/out/final_billmA.vcc` |
| size | 2,595,584,000 bytes |
| SHA-256 | `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4` |
| status | published, FINAL, current last entry |
| panel | `vcc2026-test-1`, contexts D, E, F, 300 targets |
| deadline | 2026-11-05 23:59 UTC |

Verified unchanged 2026-10-25 16:33Z, after the reconciliation job.

## Why (the asymmetry)

Only the last final entry counts, and official validation scoring is closed (last public entry
2026-10-22; zero entries 10-23 through 10-25). So:

- No candidate built now can be scored before the deadline.
- Replacing is one-way; a worse entry destroys the better one with no revert.
- Fable: "the closed server dominates every one of these facts."

Our base rate for "local harness says better → officially better" is **0 for 7** since 2026-10-21.

## When this rule may be broken

Only for a **concrete submission-invalidating defect**, e.g. the account shows a different last entry,
the artifact is truncated or hash-mismatched, or the panel/context mapping is proven wrong. In that case
the remedy is to restore **the same audited package**, not a variant. A disappointing score is not a
defect.

## Known-bad reasons to reopen (all already rejected)

| tempting | why not |
|---|---|
| "Amplitude/layout sweep might hit" | Measured twice (v8b, slim_v5b); negative-sum trade on NMAE |
| "Dense/unmasked response raises rho" | Untested, and unscoreable before the deadline |
| "We have 11 days, do something" | Occupancy is not evidence; Astra and Kimi rejected it explicitly |
| "Local harness now agrees" | It does not rank official-panel recipes; verdict (b) |
| "Just try one package" | That is precisely the one-way bet the advisors forbade |

## Standing corrections (do not regress)

- `learnednorm` is **VOID**, not negative — aborted at 50/79 on its own 30 GB disk guard.
- `u_min = 1 - rho^2` is **not** an identity for this scorer. See
  `decisions/20261025_1645_harness_reconcile_result.md`.
- Harness ordering is per-fold, not uniform: H1_rich 15/21, CD4_depth192 18/21, pooled 9/10.

## Remaining authorized work

1. Read-only status checks on Nov 4-5. Confirm the entry above is still last. Intervene only for a
   concrete defect.
2. Learning-only, firewalled from submission: representative folds, controlled recipe pairs, and the
   four unimplemented rho-raising mechanisms in `decisions/20261025_1540_second_consultation.md`.
3. Nothing that touches the package.

## Correction (2026-10-26 23:50Z)
The "0 for 7" base rate above overstates its source. Fable's line (decisions/20261025_1525_THREE_MODEL_DECISION.md)
was "the last seven attempts measured worse or null", i.e. seven candidates that came out worse or flat locally. It
was NOT seven cases where the local harness said "better" and the official board disagreed; no official reading
exists after 10-22. The verified local-vs-official record is calib3_1026/REPORT.md §1: on 4 checkable contrasts the
local composite sign was right 2/4 (e.g. v9 layout: local +0.0655, official −0.0086). The freeze decision itself
does not depend on the wording.
