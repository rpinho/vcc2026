# DECISION 2026-10-24 22:15 ET — the "acquire a same-line screen" lever is devalued ~5x; keep it closed

Source: `context_identity_1025/REPORT.md:123,125`. Derived this session; recorded so it is not
re-litigated or over-promised.

## The measurement
On the **94 official-overlap targets**, holding fold and scorer fixed, **cell-disjoint, same-line**
transfer is:

| comparison | rho_opt |
|---|---:|
| same-line, cell-disjoint | **0.143471** |
| pooled cross-line | **0.120482** |
| **difference** | **+0.022989 (1.191x)** |

Full 300-target difference: +0.022898. Descriptive split diagnostic — no claimed significance, no
independent-screen replication, no final-panel forecast.

## What this kills
The identity result (D = MOLM-13-like, E = U-87MG ATCC-like) plus 139/300 direct cross-line
exact-target coverage *looked* like it reopened the MSE lever via a public MOLM-13 / U-87MG screen.
**It does not.** Going from "different cell line" to "same cell line, different screen" buys
**1.19x**, not the 4.8x that motivated the search.

## Why the 4.8x was never ours to claim
The 4.8x (0.73 vs 0.15) came from **H1 = same experiment** — same cell line, *same screen*, same
guides, disjoint cells. `REPORT.md:125` explicitly states the H1 multiplier is **"not a justified
multiplier for the final panel."** Same-experiment is not same-line. Nobody can buy same-experiment
for D/E/F; the final panel is a cross-line problem by construction.

## The contamination caveat (limits even the 1.19x)
`REPORT.md:125`: the **unamplified input-control mean alone** gives rho_opt **0.135225** (B=0.275108).
These oracle correlations contain a **shared input/reference-control component**, so 0.135 of the
0.143 "same-line signal" may be a control artifact rather than target-specific transfer. Do not read
each point of rho_opt as target-specific signal; high untouched-bank scores are contaminated.

## Binding consequences
1. **Do not resurrect the data-ceiling claim** as if perturbation screens lack information — it is a
   statement about *cross-line transfer of target-specific direction*, not about screens in general.
2. **Do not over-promise a same-line bank either.** A MOLM-13 or U-87MG public screen is worth
   ≈1.19x on the transfer term, and possibly less after the control component is removed. That is
   not a reason to spend the remaining runway on data acquisition.
3. The open levers are the ones that act on *layout/call-density* and *scoring geometry*, which is
   why `layout_cd4_1025` is the live test. Cross-recipient generality is the question that decides
   whether layout work is package-level or H1-specific.
