# 2026-10-25 14:55Z — my own experiment spec was wrong; caught by a reproduction gate

## What happened

I launched `vcc_dense_bulk_1025` at 10:12 EDT with a spec that said:

> "`confidence` must be computed consistently with the profile each arm actually emits (the v9
> call order uses the sign of the emitted response). Recompute it per arm…"

That instruction is wrong, and it made the job's `incumbent` arm **not the incumbent**.

Production computes the v9 confidence array from `ctx.L` — the sign-gated, deflated fused
response — *before* the 0.3 learned blend and *before* alpha rescaling
(`model_billmB/panel_pipeline.py:58-66`; `final_push/main/panel_pipeline.py:94`). The emitted
profile used for expression is the blended one, but the DE call ordering is keyed on the gated
array. My spec told the job to key it on the blended profile instead.

## Why it mattered

`sign(gated L - Lc)` and `sign(blended profile - Lc)` disagree on **82.25% of gene entries**
(1,174,636 of 1,428,083), because the gate zeroes genes that the blend then reintroduces at
small magnitude. Measured consequences in attempt 1:

| | attempt 1 | production-faithful |
|---|---|---|
| prediction digest | `a9f90b1a…` | `2cbf35e2…` |
| overall | **-0.018159** | **-0.019164** |
| fidelity | 0.455575 | 0.453333 |
| reach | 0.113989 | 0.117822 |
| Jaccard | 0.043544 | 0.042820 |
| PDS, raw u | identical (≤1e-16) | identical |

The overall shift is **+1.005e-3** — the same size as pool1's entire measured gain (+1.7e-3).
Had I let it run, every arm would have been compared against a comparator that was itself off by
more than the effect being hunted.

## How it was caught, and why that matters

Not by review — by the *frozen reference*. The previous job's incumbent receipt
(`final_push/fable_tests_1024/artifacts/H1_rich_search/incumbent/score.json`) carries a
prediction digest. Attempt 1's digest did not match it, and the DE components moved while PDS
and raw u did not. That asymmetry is what pointed at the call-ordering path rather than the
expression path.

I did not have that check in the spec. I have now added it as a **mandatory reproduction gate**:
`incumbent` must reproduce digest `2cbf35e2ac60d3d6…` and overall `-0.01916405775558278` within
1e-12 before any other arm is scored. The lesson generalizes: any ablation that reuses a frozen
baseline must gate on that baseline's *exact* reproduction, not on "it looks about right".

## What I did

- Stopped the run at ~16/79 of the second arm; no arm's result was used.
- Archived attempt 1's receipts to `final_push/dense_bulk_1025_attempt1_notes/` (72 KB), deleted
  484 MB of regenerable arrays, and recorded that its conclusions must not be reused.
- Carried forward only the frozen biology/source arrays (`profiles.npz`), which are independent
  of the confidence bug and hash-checked against `profiles.json` on load.
- Relaunched in a clean directory with the corrected spec: production-faithful confidence for
  `incumbent`/`eta0`/`learnednorm`; for `dense` — which has no gate to key on — confidence from
  the ungated response, stated explicitly as a necessary deviation rather than a silent one.
- Safety package hash verified before and after; unchanged (`ce78c81c…9dc4`).

## Also amended in the spec

The pre-registered "scaled MSE > 0" clause is **void by construction** and must not be evaluated
(see `decisions/20261025_1445_raw_u_recovers_mse.md`): the H1 anchor makes it unreachable for
every arm in this family. Raw u replaces it as a mechanism readout, not as a pass condition.
This was amended before I saw any arm's numbers from either attempt.
