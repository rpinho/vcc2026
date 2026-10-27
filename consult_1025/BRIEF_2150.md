# VCC final-phase consultation, 2026-10-24 ~22:00 ET

You are advising on the Virtual Cell Challenge 2026 final phase. Read-only; answer in <=1200 words.
Do not run jobs, do not download, do not submit.

## Hard facts (measured; do not re-derive, do not contradict without a citation to a file)
- Official validation scores: best is **billmA 0.16438 (rank 463/1580)**: pds .706, mse .061, nmae .116,
  fid -.066, reach .186, jac -.016. Our earlier v5b 0.1528, v9a 0.16135, v8a 0.1470.
- FINAL PHASE IS BLIND: new contexts D/E/F, new 300-target panel, only the LAST entry counts, due 2026-11-05.
  Our safety entry 5NeeoWn96lfm23hgExcm (package ce78c81c...) is already published with a hidden score.
- Local proxy scores have INVERTED official ordering twice. The local overall composite is NOT a decision quantity.

## What is closed (measured, do not propose re-testing)
Amplitude/rho sweep, pool_k, global bulk multipliers, B/claim ratio, PC projection, learned-blend coefficient,
source breadth/family fusion variants, per-target call-budget shrinkage (donor counts predict true n only weakly:
Spearman 0.64, median est 65 vs true 107), call-density sweep (H1 model is monotone toward the dense endpoint,
yet the official dense endpoint scored WORSE: 0.1528 vs 0.16135).

## The central measured fact
Response-direction transfer, amplitude-matched, in the scorer's bulk-lognorm space at bulk_target_sum=50000:
| regime | rho_opt |
|---|---:|
| same experiment, disjoint cells (H1 bank -> H1 fold) | **0.732** |
| same cell line, cell-disjoint (HEK holdout) | **0.1435** |
| cross-line pooled (HCT116/HEK293T/K562/CD4 donors) | 0.1205 |
| shipped incumbent | 0.128 |

Contexts identified from control expression (CCLE Spearman 0.86-0.91, HPA cross-check agrees):
- **D = MOLM-13-like** myeloid leukemia (MOLM-14 close alternative; HVG gap only 0.0136)
- **E = U-87MG ATCC-like** glioblastoma (large margin)
- **F = unresolved fibroblast-like/mesenchymal** (no confident clone; HA1E kidney alternative at gap 0.0006)
We hold 11 perturbation banks: H1, HCT116, HEK293T, K562 (raw + GWPS), CD4 (Rest/Stim8hr/Stim48hr), HepG2, Jurkat, RPE1.
**Same-line coverage for D/E/F is 0/300. Direct cross-line exact-target coverage is 139/300.**

## Officially measured layout trade (H1, one recipient, paired, same head/totals/B)
v9 confidence-ordered layout (885 calls): fid 0.4533, jac 0.04282, reach 0.11782
native dense layout (~4475 calls):          fid 0.5492, jac 0.08133, reach 0.06638
A job is running RIGHT NOW to replicate this on a second recipient (CD4). Verdict frozen in advance:
REPLICATES only if fid delta > +0.020, jac > +0.010, reach < -0.005, and the paired target-bootstrap 95%
fidelity interval excludes zero.

## Constraints
- Local Mac Mini CPU only; no pods, no new credentials, no uploads, nothing leaves the machine except reading
  public sources. Public data downloads ARE permitted (reading public archives).
- ~9 days to the deadline. Two validation submissions/day are available but the final panel is blind and the last
  entry is the only one that counts.
- Every claim needs an independently written exact verifier; negative results are successes.

## The question
Given D/E/F are cross-line targets at rho_opt ~0.12-0.14, and given ~9 days of local CPU plus the 2/day validation
slots during the remaining visible period:

1. **Rank the levers by expected official-score gain per day of work**, and say for each whether it is measurable
   locally or only on the official scorer. Be concrete about the mechanism that would move which metric component.
2. Is **acquiring a public perturbation screen in MOLM-13, U-87MG, or a fibroblast line** the best use of the next
   36 hours? If yes, name specific candidate public datasets/resources worth checking (dataset names, accessions,
   repositories) and the minimum bar that would justify incorporating one. If no, say what beats it and why.
3. Which single metric component is most under-exploited given the table above, and what is the concrete
   experiment that would move it without a local-proxy trap?
4. Name the traps specific to this situation: places where our own history shows a local measurement misled us.
