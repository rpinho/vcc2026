You are deciding the next move for a team in the Virtual Cell Challenge 2026.
Be brutal. We want YOUR decision. Do not give a menu.

## An important correction to what you were told earlier
You (Fable, Astra) were asked ~40 minutes ago to rule on the manifold experiment. Astra said "if the
manifold wins, run a frozen head-to-head through the official evaluation pipeline." Fable said "the
manifold is unproven until F3 beats shipped B0 with a paired-bootstrap CI excluding zero; otherwise
kill STRING and move to Sidechain staged source progression."

Then a head-to-head job (`head2head_1026`) tried to run that test and stopped, declaring it
IMPOSSIBLE: "manifold H1 has 168 targets, the official H1 SEARCH panel has 79, intersection is ZERO,
so no same-target comparison exists."

**That stop is now known to be wrong.** The zero intersection is real but it is the anti-leakage
design, not a blocker. I verified directly:
- The 168 are H1 *training* targets with cached responses (manifold_1026/inputs/data/H1_full.npz).
- The 79 official targets are a DIFFERENT set, and all 79 have real local ground truth in
  gatecall_1026/inputs/banks/H1.npz (62 / 146 / 372 cells at min / median / max).
- Donor families (all non-H1 banks) cover 79/79 of the official targets.
- STRING embeddings cover 77/79.
- The shipped billmA prediction for those 79 also exists (inference_H1.npz).
- Truth reconstructs exactly as log1p(50000*p) - log1p(50000*ctrl).

So the honest test is: FIT ON THE 168, PREDICT THE 79, score against held-out truth. Nothing about
the 79 touches fitting. That is a genuine out-of-sample comparison and it is running now as
`h2h2_1026`, with arms B0 (shipped), B1 (manifold F1 STRING-only), B2 (manifold F3), B3 (pooled
branch), B4 (training-mean), B5 (shipped shrunk toward training-mean).

## Where we actually stand
Official best billmA **0.16438** (rank ~463/1635). Leader 0.47203. Top-100 median 0.287.
Frozen safety entry stays last; we CANNOT submit without Ricardo’s authorization.
Three lines are CLOSED with real negative results: recipient-control conditioning (dirmap2),
effect-size gated DE calling (gatecall), and earlier amplitude/shape/transfer work.
Gap decomposition measured from the live board: expression MSE is 40-52% of the gap to every target
(leader, em2.bio, top-10 median, top-100 median). We already MATCH the top-100 median on DE
fidelity/reach. So DE work cannot move us much; expression MSE is where the remaining points are.
Manifold is the only live expression signal: H1 cosine 0.226 (F3) / 0.216 (F1) vs pooled 0.083,
training-mean 0.181; beats pooled 7/7 families.

## The thing that worries us most
Shipped B0's raw direction fidelity on these 79 official targets is 0.4487 — BELOW chance 0.4929.
Interpretation: our shipped model is essentially not predicting DE direction on this panel. It also
means B0 is itself out-of-sample here, so the comparison is fair. But it raises a question: is our
whole score carried by something other than DE skill, and if so what?

## Answer decisively — pick ONE
1. If h2h2 shows the manifold does NOT beat shipped B0: what is the single next experiment
   (protocol + kill test)? Note Fable previously said Sidechain staged source progression and Astra
   said response-reliability source weighting. They disagree. Say which one we run Monday morning,
   and why the other waits.
2. If h2h2 shows it DOES beat B0: what must be true before a package is worth building, and what is
   the first thing that would invalidate it?
3. B0's below-chance direction fidelity: is this a bug worth investigating before anything else, or
   is it expected for a model that never fit these targets? Be specific.
4. Is the realistic objective (a) move up from rank 463, (b) reach ~0.25, or (c) contest a specific
   component? Pick one and say what it implies we should STOP doing.
5. Is a GPU rental warranted for anything in the next 48 hours? Prior evidence says no.
