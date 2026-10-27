# DECISION 2026-10-25 04:05 ET — the safety package is byte-identical on disk; provenance is intact

## Why this mattered
Both overnight consults (Kimi `consult_1025/kimi_0402.md`, Fable `consult_1025/fable_0402.md`) said to stop
hunting response-transform arms and spend the remaining CPU verifying the retained final package. That advice
is only actionable if the published artifact can be *identified*. The obvious-looking file is the wrong one:
`model_billm/billmA.vcc` is 2,700,052,480 bytes, SHA256 `24e6dfdff30421311035d66b31b402b8fee68cee81844b12707ac16fb81ade39`.

Comparing *that* file to the published digest would have produced a false "package missing/mismatched"
finding. It is a different build, not the submitted one.

## The verified chain
The packaged, published artifact is **`model_billmB/out/final_billmA.vcc`**:

| property | measured now | recorded |
|---|---|---|
| path | `model_billmB/out/final_billmA.vcc` | same, per `artifacts/panels/final_billmA/decode_receipt.json:command` |
| size | **2,595,584,000 bytes** | 2,595,584,000 (`decode_receipt.json:bytes`) |
| SHA256 | **`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`** | identical |
| shape | per receipt 360,000 × 18,533 | 360000 × 18533 |
| nnz | per receipt 2,010,253,470 | same |
| contexts | D,E,F; pert list `data/final/bundle/pert_counts.csv` | same |

The digest on disk **matches the recorded published digest exactly**, so the safety entry
`5NeeoWn96lfm23hgExcm` is backed by a byte-identical artifact. Nothing has been rebuilt, moved, or replaced.

## What this closes and what it does not
- **Closed:** the provenance question "can we point at the exact bytes that were published, and are they
  intact?" Yes. This is the single most important operational fact for the rest of the night: the fallback is
  real and untouched.
- **Not closed:** a *rebuild* reproducibility test (does re-running the pipeline regenerate these bytes?). That
  is a separate, much more expensive claim and is **not** implied by an on-disk hash match. Do not conflate them.
- **Not closed:** any claim about the hidden final score. The package is published with a hidden score and we
  have no evidence about it.

## Consequence for tonight
Any future candidate must be compared against **this** artifact's identity, not `model_billm/billmA.vcc`. Any
report that cites `billmA.vcc` as "the submitted package" is wrong by 104 MB and a different digest.

No package was modified, rebuilt, uploaded or replaced. No submission. Safety entry unchanged.
