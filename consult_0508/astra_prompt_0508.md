You are an independent scientific reviewer. Read the brief below and answer in under 900 words. Do not run experiments or modify files.

# VCC decision brief — 2026-10-25 05:08 EDT

We are Ricardo's VCC research agent. The decision window ends at 08:00 ET.
Hard constraints: NEVER SUBMIT; never modify the published safety package
model_billmB/out/final_billmA.vcc (2,595,584,000 bytes, SHA-256
ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4); no pods/new
credentials; local CPU + Codex only; maintain the storage reserve.

Current verified facts:
- Official validation best: billmA 0.16438 (PDS .706, MSE .061, NMAE .116, fidelity
  -.066, reach .186, Jaccard -.016), rank 463. Safety final entry
  5NeeoWn96lfm23hgExcm is the same package; final score hidden.
- The local H1_rich_search ablation incumbent scores -0.019164, far below its
  official analogue; local-to-official transfer is not established.
- Closed by measured evidence: rho 1.4/2.0 fail G3 budget; pool1 tiny local gain
  but no transfer margin; post-blend re-solve tiny and MSE-negative; CD4
  expression leg closed; family_quality/equal arm failed the frozen budget-match
  gate (requested ratio moved .909->.960, tolerance .020); source-breadth arms
  all blocked; unmask and dense-bulk arms already measured; no defensible
  transfer-noise margin.
- Two independent read-only software audits are running now (rebuild
  reproducibility revision 11 and bridge receipt revision 17). The actual
  production evidence check PASSes: 900/900 groups, identical package SHA, etc.
- The rebuild was resealed at 08:59Z with a strict verifier PASS. Package bytes
  were unchanged.

Question: given those closed levers, what is the single most valuable NEW
scientific experiment we can run locally before 08:00 ET that has a real chance
to improve the hidden final score (not merely produce another diagnostic), and
what exact pass/fail gate and minimal run should we use? If the honest answer is
"none; finish verification and hand off", say so plainly and explain what
evidence would justify reopening a lever. Rank at most 3 alternatives with
expected value and risk. Be concrete about data and compute.
