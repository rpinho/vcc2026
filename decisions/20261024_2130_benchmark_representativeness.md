# DECISION 2026-10-24 21:30 ET — the local benchmark is not representative of the final panel

## Measured (all from retained artifacts, reproducible)

Target overlap between the official 300-target panel (`data/final/bundle/pert_counts.csv`) and
each public signature bank (`model_billm/data/*.npz`):

| bank | targets | overlap with official 300 |
|---|---:|---:|
| HCT116 | 408 | **124** |
| HEK293T | 421 | **131** |
| K562_GWPS_raw | 675 | **126** |
| CD4_Stim8hr | 404 | **128** |
| H1 | 132 | **11** |

Overlap of the local benchmark *tasks* with the official panel:

| task | n | overlap | frac |
|---|---:|---:|---:|
| H1_rich_search (our primary SEARCH) | 79 | 8 | 10.1% |
| HCT116_search | 300 | 89 | 29.7% |
| HEK293T_holdout | 300 | 94 | 31.3% |

So our primary calibration fold is the *least* representative of the final panel of any fold we hold.
This is a separate issue from any model choice and was not recorded before tonight.

## Measured — same-line vs cross-line response correlation (H1_rich_search, amplitude-matched)

Recomputed in `ceiling_audit_1024/` (see that directory's scripts). Every candidate is
amplitude-matched to the shipped budget B=0.142779 via the emitter's own scaling, then scored in the
scorer's `bulk_lognorm` space with own-target exclusion and 1/|G_p| weights.

| candidate | rho_opt | u | bank-vs-truth delta corr (mean) |
|---|---:|---:|---:|
| H1 bank (same screen/line) | **0.610** | 0.628 | 0.402 |
| HCT116 bank (different line) | 0.138 | 0.981 | 0.042 |
| HEK293T bank (different line) | 0.129 | 0.983 | 0.025 |
| K562_GWPS_raw | 0.121 | 0.985 | 0.055 |
| CD4_Stim8hr | 0.115 | 0.987 | 0.032 |
| **shipped incumbent** | **0.128** | 0.984 | — |

Same-line evidence is worth ~4.8x the response correlation. The shipped fusion sits at the
cross-line level, which is what the earlier "data ceiling" finding was actually measuring.

## What this corrects
`mse_direction_1024/FINDING_per_source_ceiling.md` and `FINDING_ceiling.md` concluded "the ceiling is
DATA". That is true only for *cross-line transfer*. The scripts excluded the H1 bank as
"leave-one-context-out", which is correct for estimating cross-line transfer but wrong as a statement
about the information available on a panel whose contexts can be identified.

Two further corrections to those scripts, both applied in `ceiling_audit_1024/audit_units.py`:
1. Units: bank `response` is `log1p(50000*p) - log1p(50000*control)` (`model_billm/prepare.py:24`),
   a log1p-CPM fold change. The ceiling scripts added it directly to raw counts (`cp + fac*D`), where
   `cp ~ 1e-2` and `D ~ O(1)`. Re-scoring under the correct multiplicative CPM mapping changes
   individual rho_opt values by +0.004..+0.010 but does **not** change the ordering or the conclusion.
   The units bug is real; its effect on the headline was small. Recorded rather than buried.
2. Leakage check: the H1 bank is not the same cells as the benchmark (delta correlation 0.40, not
   >0.95), so 0.610 is same-line transfer from a disjoint screen, not in-sample leakage.

## Consequence
The single highest-value open question is now: **what are official contexts D/E/F?**
`data/final/bundle/context_{D,E,F}.h5ad` contain only 18,400 control cells each. Contexts A/B/C were
identified from exactly this kind of data at CCLE Spearman 0.88-0.90 (`data/CONTEXT_IDENTITY.md`).
If D/E/F match a line we hold, the final package can use same-line signatures and the MSE lever reopens
with a measured 4.8x. If they do not, the final panel is a cross-line problem and the ceiling finding
stands for it.

Job launched: `_codex/vcc_context_identity_1025.md` -> `vcc2026/context_identity_1025/`.

## Also recorded: what this does NOT license
- No claim that the incumbent is worse than believed; its 0.061 official MSE is measured.
- No tuning on D/E/F (no truth exists locally).
- No submission.
