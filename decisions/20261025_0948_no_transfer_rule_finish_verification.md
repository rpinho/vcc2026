# 2026-10-25 09:48Z — no transfer-aware rule; finish verification, hand off

## Question

Can the existing matched local/official candidate records support a
transfer-aware acceptance rule for a future candidate? If not, what is the
minimum missing observation?

## Independent answers (no shared context)

- **Kimi** (`orq best`, `consult_1025_0920/kimi.md`): **No.** Seven recipes are
  not fourteen independent labels; the 21 pair contrasts share candidates. The
  four recent surrogates have Spearman 0 and 3/6 direction agreement. A linear
  fit's slope moves from 0.105 to 0.270 when v3b is excluded, and leave-one-out
  median absolute error is 0.0313, which is diagnostic, not a transfer bound.
  Anchors do not repair v9a-surrogate minus v5b (-0.05958 locally vs +0.00855
  officially). Finish verification first; reopen calibration only with a
  faithful matched native contrast.
- **Astra** (`codex exec -m gpt-6-astra`, `consult_1025_0920/astra.md`): **No.**
  The records support exploratory calibration and falsification, not a positive
  acceptance threshold. Seven recipe-level comparisons, observed reversals, and
  no independent transfer-error estimate. One faithful contrast could falsify
  or begin validating a mapping; a reusable 95% bound needs several or a
  validated error model.
- **Opus 5.5** (`orq anthropic/claude-opus-5.5`, `consult_1025_0920/opus.md`):
  **No.** Same estimand and same missing observation. Recovered anchors give
  15/21 agreement with MSE floored or 20/21 unclipped, but still predict
  v4a→v4b +0.04493 locally against -0.01950 officially. Choosing the
  better-looking convention after the fact would fit the check.

All three independently ranked **finish verification and hand off** first and a
bounded reconciliation of already-retained matches second. None endorsed a new
parameter, blend, source, or expression sweep. None would accept a target
bootstrap as a transfer-error bound.

## Decision

1. Do not launch a new model fit, prediction sweep, or unregistered arm.
2. Finish the audits bound to the exact current revisions and record their
   actual verdicts. AUDIT13 has now passed rebuild revision 13.
3. Keep the safety package unchanged. Do not submit.
4. Reopen calibration only if a faithful, recipe-matched local contrast on a
   frozen representative native panel becomes available under existing
   constraints. A candidate rule would still require independent transfer
   errors, a held-out forecast test, and a pre-frozen one-sided 95% predictive
   lower bound above zero, plus every existing gate.

## Supporting evidence

Kimi's compact check is retained in `transfer_acceptance_1025/DECISION.md` and
`transfer_acceptance_1025/evidence.json`; all 18 named historical inputs still
match their recorded hashes. Opus and Astra cite
`final_push/harness_audit_20261023/audit.json` and
`final_push/artifacts/de_check/H1_rich_search/calibration_summary.json`.

## Safety

Nothing submitted. Safety entry `5NeeoWn96lfm23hgExcm`; package
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`
(2,595,584,000 bytes) unchanged.
