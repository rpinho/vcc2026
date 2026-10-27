# VCC decision brief — 2026-10-25 09:20Z (05:20 EDT)

Hard constraints: NEVER SUBMIT; never modify or
replace the published safety package `model_billmB/out/final_billmA.vcc`
(2,595,584,000 bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`); no pods or
new credentials; local CPU + Codex only; maintain the storage reserve.

## Established facts

- Official validation best: billmA **0.16438** (PDS .706, MSE .061, NMAE .116,
  fidelity -.066, reach .186, Jaccard -.016), rank 463. The safety entry
  `5NeeoWn96lfm23hgExcm` is the same package; its final score is hidden.
- Local H1_rich_search incumbent scores **-0.019164**, far below its official
  analogue. Local-to-official transfer is not established.
- Closed by measured evidence: rho 1.4/2.0 fail G3; pool1 +0.001694 local but no
  transfer margin; post-blend +0.001119 and MSE-negative; CD4 expression leg
  closed; family_quality failed the frozen budget-match gate (requested ratio
  .909 -> .960 vs tolerance .020); source-breadth arms blocked; dense-bulk and
  unmask already measured; no defensible transfer-noise bound.
- The published package was independently decoded read-only: 900/900 archived
  dense hashes matched, package SHA unchanged, shape [360000,18533], nnz
  2,010,253,470. This verifies the artifact, not its score.
- Rebuild verification and bridge verification are in fresh independent audits
  after prior audits found real software defects. No new scientific arm is
  running. Nothing has been submitted.

## The question

Do **not** recommend another parameter, blend, source, or expression sweep
unless you can name evidence that changes why the previous arm was rejected.

Instead answer this: **can the existing matched local/official candidate records
support a transfer-aware acceptance rule for a future candidate?** If yes,
specify the exact analysis, which already-retained records to use, the estimand,
the uncertainty method, and a pass/fail rule that would justify one new arm. If
no, say so plainly and identify the minimum missing observation.

Rank at most three options by expected value for the eventual Nov 5 decision and
by risk of false confidence. A diagnostic that cannot change a decision should
rank below finishing verification and handing off.

Be concrete about data paths, compute, and what result would change your answer.
