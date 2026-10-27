You are deciding the next move for a team in the Virtual Cell Challenge 2026.
Be brutal. We want YOUR decision. Do not give a menu.

## Position
Official best billmA **0.16438** (rank ~463/1635). Leader Illumina/iCell **0.47203** (proprietary
Billion Cell Atlas; no public code). Top-10 median 0.321, top-100 median 0.287.
Frozen safety entry stays last; we CANNOT submit without Ricardo’s authorization.

## Two lines closed tonight (both real, both negative)
1. **Recipient-control conditioning** (dirmap2): H1 official fold, four arms in one run —
   incumbent 0.021555, A-shared 0.008394, B-real 0.011010, B-shuffled 0.009601. Incumbent wins all.
2. **Effect-size gated calling** (gatecall): emit only |log2FC|>0.8 genes. H1 raw fidelity
   A0 0.4487 -> A1 0.1511, A2 0.1509. Both CI-exclude-zero NEGATIVE. Also fails on CD4.
   Conclusion: our high-|LFC| genes are not the same genes as the reference DE set.

## The one positive signal, and its exact weakness
**Manifold experiment** (manifold_1026): per family, 5-fold over targets, SVD on training-target
responses, ridge from 64-d STRING target embeddings -> PC scores, reconstruct held-out response.
Verified numerically against a separately written verifier.

| family | cand cosine | pooled-branch cosine | training-mean cosine | MSE/pooled |
|---|---:|---:|---:|---:|
| H1 | 0.226 | 0.083 | 0.181 | 0.466 |
| K562 | 0.344 | 0.181 | 0.198 | 0.534 |
| Jurkat | 0.326 | 0.164 | 0.209 | 0.665 |
| HepG2 | 0.375 | 0.203 | 0.286 | 0.722 |
| RPE1 | 0.474 | 0.202 | 0.425 | 0.712 |
| HCT116 | 0.184 | 0.032 | 0.183 | 0.744 |
| HEK293T | 0.153 | 0.031 | 0.138 | 0.651 |

STRING-only arm (F1, no incumbent stacking) gives 7/7 over pooled (0.216/0.284/0.260/0.313/0.440/
0.183/0.139) and 6/7 over training-mean.

**Its preregistered gate FAILED, for a statistical reason the report identifies honestly:** after
response centering with an unpenalized intercept, a heavily shrunk shuffled ridge predicts the
training response mean, so "shuffled must not beat do-nothing" is really testing the generic response
mean, not leakage. The verifier found no scope overlap.

**The weakness that matters:** the comparator was the POOLED BRANCH, not our shipped MLP ensemble.
Our shipped model's official PDS raw is 0.8187, implying a true cosine far above 0.083. So
"beats incumbent 7/7" may only mean "beats a weak baseline". A head-to-head against the ACTUAL shipped
prediction is running now (arms: shipped B0, manifold F1, manifold F3, pooled B3; paired bootstrap).

## Our remaining levers, with status
- Expression MSE is 45% of the leaderboard gap; official raw 0.973 (=2.7% skill). Everything tried on
  it has failed except this manifold signal.
- DE metrics: fidelity raw 0.4929 (chance), reach 0.2438, jaccard 0.0245, nmae 0.9301.
- Untried from the audit: response-reliability source weighting; combined STRING/sequence embedding;
  priors for the 161 targets with no source row; saturation routing; responder gating; 250-300/target
  DE budget; Sidechain staged bulk-strength source progression (SER-21 .1757 -> SER-27 .1875 -> SER-29 .2226).
- We have 16 GB CPU Mac Mini, Modal credentials, a documented prior RunPod path (~$3 for 3 h on a
  4090), and ~10 days to the Nov 5 deadline.

## Answer decisively
1. If the head-to-head shows the manifold does NOT beat the shipped model — what is the single next
   experiment? Name one, with its protocol and a kill test.
2. If it DOES beat the shipped model — what is the single next step, and what could still invalidate it
   before a package is worth building?
3. Given the evidence, is the realistic objective (a) move up from rank 463, (b) reach ~0.25, or
   (c) contest a specific component? Pick one and say what it implies we should stop doing.
4. Should we rent a GPU for any remaining work, or is CPU sufficient? Say plainly yes or no and why.

FORMAT: five labelled lines first, nothing else on them:
VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
OBJECTIVE:
Then justification, max 240 words. Cite only numbers above; invent nothing.
