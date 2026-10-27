You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 2026-11-05 23:59 UTC). The
orchestrator executes what the advisors rule; it does not make the science call. Give ONE decision, not
a menu. Be specific enough that it can be executed mechanically.

## What you already ruled (round consult_1026_pds, 08:15-09:05 EDT today)
Unanimous: PDS is the protected asset (it is +0.1219 of our +0.0216 composite); expression MSE is a
dead objective (B0's scaled MSE is already clamped at 0.000, raw 1.4107 vs anchor 0.9613); the next
modelling experiment is amplitude calibration of B0; no GPU; packaging requires official score_avg to
strictly beat 0.021555, PDS >= 0.845, positive paired-bootstrap lower bound, confirmation fold.

But you gave THREE INCOMPATIBLE FORMS of the amplitude experiment:
- Astra: pure scaling of the control-relative delta, delta' = alpha_f * delta_B0, alpha per family, LOO on
  the 168. "Do not substitute mu_train + alpha(B0 - mu_train) - that rotates the effect."
- Kimi: exactly y = mu_train + alpha(y_B0 - mu_train), per family, LOO on the 168.
- Fable: GLOBAL scalar ladder k in {1.0, 0.5, 0.25, 0.1} on the delta, plus renderer dispersion
  {0.5x, 1x, 2x} at k=1.0. Explicitly no per-family scalars.

Only one family (H1) is in the official measurement loop, so "per family" vs "global" is moot for this run.

## Measurements made since, that none of you had (proxy = response-space cosine retrieval rank of own
## truth among the 79, self-gene excluded; rank->PDS calibrated on official B0/B1, slope -0.014812/rank;
## validated: the constant arm B4 returns rank 40.0 = chance, PDS 0.505)

1. Pure scaling alpha*B0 is EXACTLY PDS-invariant in the proxy (rank identical for alpha 0.5..3.0). It
   cannot buy PDS. Any effect would come only through the renderer's DE calls (unmeasured).
2. A scalar has ALREADY been run officially (mse_1025, CD4 fold, integer v9): source-fitted alpha=0.343
   gave MSE -13.7% but raw DE fidelity 0.1767 -> 0.1251 (-29%). Verdict then: do not promote.
   So Fable's shrink ladder (0.5/0.25/0.1) is moving toward a region already measured to cost fidelity.
3. Kimi's form mu + alpha(B0 - mu) with alpha < 1 is catastrophic for PDS: with the legitimate 168-target
   training mean, alpha=0.5 -> proxy PDS 0.647 (-0.07 composite). B5 (the 0.5 shrink toward the training
   mean, officially scored below) is exactly this form at alpha=0.5.
4. NEW, and it reverses the predecessor's "oracle artifact" claim: Kimi's form with alpha > 1 is
   EXTRAPOLATION AWAY from the shared training mean (it amplifies the target-specific part of B0 and
   suppresses the pan-perturbation shared component). With the legitimate 168 mean (no held-out data):

   | alpha | mu-interp rank | est PDS | scale rank | est PDS |
   |---|---:|---:|---:|---:|
   | 0.50 | 30.46 | 0.647 | 16.57 | 0.853 |
   | 0.75 | 23.48 | 0.750 | 16.57 | 0.853 |
   | 1.00 | 16.57 | 0.853 | 16.57 | 0.853 |
   | 1.10 | 15.42 | 0.870 | 16.57 | 0.853 |
   | 1.25 | **14.73** | **0.880** | 16.57 | 0.853 |
   | 1.50 | 15.05 | 0.875 | 16.57 | 0.853 |
   | 2.00 | 17.38 | 0.841 | 16.57 | 0.853 |
   | B0 - mu (full centring) | 23.23 | 0.754 | | |

   (My proxy reproduces the predecessor's B0 rank as 16.57 vs 16.75 - a small implementation difference,
   tie/mask handling; relative comparisons within one proxy are what matter.)
   Est. gain at alpha=1.25: +0.027 PDS ~= +0.0095 composite points from PDS alone. MSE gets worse but is
   already clamped at 0, so it is free. DE effects are unmodelled.

5. TRAINING-ONLY replication (no held-out data loaded, so alpha can legally be chosen here): on the 168
   fitting targets, using the inherited pooled() donor response as B0's surrogate (B0 itself exists only
   for the 79; the pooled arm B3 matches B0 on the 79 proxy, rank 16.80 vs 16.75), leave-one-out mean,
   retrieval rank of own truth among 168:

   | alpha | mu-interp (LOO mu) | pure scale |
   |---|---:|---:|
   | 0.50 | 52.29 | 25.10 |
   | 1.00 | 25.10 | 25.10 |
   | 1.10 | 22.35 | 25.10 |
   | 1.25 | **21.35** | 25.10 |
   | 1.50 | 22.17 | 25.10 |
   | 2.00 | 24.37 | 25.10 |
   | 3.00 | 26.92 | 25.10 |

   The training-only optimum is alpha = 1.25, the same as on the 79. A fixed grid, not tuned on the 79.

So: the conflict between Astra and Kimi is real, and the measurement says each is half right. Astra is
right that pure scaling preserves the effect geometry - and the proxy says it is therefore inert for PDS.
Kimi's form rotates the effect - and the proxy says the rotation is HARMFUL below 1 and HELPFUL above 1,
because extrapolating away from the shared mean is a de-meaning / discrimination-sharpening operation.
Nobody proposed alpha > 1 in that form; the LOO-on-168 fit selects it.

Caveats: proxy only. The official renderer exponentiates, clips, renormalises and emits integers; DE
fidelity/jaccard/reach/nmae are not modelled. B0 surrogate on the 168 is pooled(), not the shipped model.

## External evidence found today (competitors/DAILY-2026-10-26.md)
- Team "Res" (T3V200, validation 0.2612, rank 114; repo published 10-25, no licence - ideas only) gets
  its lead over our 0.16438 mostly from MSE: scaled MSE 0.391 vs our 0.061, PDS 0.801 vs 0.706 (MSE is
  ~0.055 of the 0.097 gap). Their lever "P1a" is an amplitude calibration: a magnitude lookup table
  indexed by (predicted direction rank x control-mean expression), weighted medians, fitted on PUBLIC truth
  DE (H1, HCT116, HEK293T, K562, KOLF), never on target labels; plus "P1N" (2x5 table for non-called genes,
  +NMAE/reach at flat PDS/MSE). So "MSE is dead" holds only for direction-changing arms like B1: a
  PDS-preserving, gene-wise (not scalar) amplitude calibration earned a top-120 team real MSE points.
  Note our H1 B0 raw MSE 1.41 vs anchor 0.96 - a scalar cannot fix a per-gene magnitude error profile.
- bioRxiv 10.64898/2026.10.20.760514 (10-25): centering an independent source on known training targets
  raised PDS on every fold; transferring signs + magnitude RANKS with isotonic rendering beat the 2025
  champion's DE score. This independently corroborates measurement 4 (de-meaning sharpens PDS).

## Official B2-B5 results (the locked h2h2 head-to-head; filled in when scoring finished)
| arm | score_avg | raw PDS | raw fidelity | raw MSE | d mse | d pds_cosine | d lfc_nmae | d direction_fide | d direction_reac | d sig_jaccard |
|---|---|---|---|---|---|---|---|---|---|---|
| B0 | +0.021555 | 0.8500 | 0.4487 | 1.4107 | +0.00000 | +0.00000 | +0.00000 | +0.00000 | +0.00000 | +0.00000 |
| B1 | -0.055694 | 0.5547 | 0.4723 | 0.9519 | +0.00164 | -0.10282 | +0.00807 | +0.01424 | +0.00195 | -0.00032 |
| B2 | -0.013018 | 0.6560 | 0.4676 | 0.9665 | +0.00000 | -0.06757 | +0.01621 | +0.01143 | +0.00443 | +0.00092 |
| B3 | +0.025615 | 0.8586 | 0.4495 | 1.4094 | +0.00000 | +0.00299 | +0.00183 | +0.00044 | -0.00147 | +0.00027 |
| B4 | -0.099178 | 0.5000 | 0.4673 | 1.0221 | +0.00000 | -0.12186 | -0.00130 | +0.01122 | -0.00816 | -0.00063 |

(d columns = composite points vs B0, scaled delta / 6.) Proxy predictions made BEFORE these landed: B2 PDS 0.668,
B3 0.849, B4 0.505 - compare to the raw PDS column. B5 (= mu + 0.5(B0 - mu), i.e. Kimi's form at alpha=0.5)
was stopped at 55/79 by a disk-space guard (machine-wide swap growth, not this job); it relaunches
automatically when disk recovers. Proxy prediction for B5: PDS 0.644. Locked verdict: B1 and B2 both fail the
primary PDS test vs B0 on the point estimate by a wide margin, so the preregistered manifold ADOPT condition
cannot be met regardless of B5; the formal REJECT will be written once B5 + verify.py complete.

**SURPRISE IN THE LEDGER: B3 (the inherited pooled() donor response, a continuity comparator, NOT a
manifold arm) scored +0.025615 officially vs B0 +0.021555 (+0.00406), with higher PDS (0.8586 vs 0.8500)
and gains on nmae/fidelity/jaccard, a small reach loss.** Under the preregistration B3 is descriptive and
cannot change the manifold verdict. But it is the first arm in this campaign to beat the shipped direction
on the official pipeline. It has now been SEEN on the 79, so promoting it would be post-hoc selection on the
acceptance set; the paired bootstrap CI for B3-B0 is not yet computed (analysis waits for B5). The gain is
below Fable's earlier +0.010 margin. The mu-extrapolation alpha=1.25 can be applied to B3 as well as B0.

## Your ruling (answer all four, briefly)
1. FORM: which arm(s) get an official run on the 79? Choose from: (a) mu-extrapolation at the
   training-fitted alpha = 1.25 (y = mu_168 + 1.25 (B0 - mu_168)); (b) pure scaling alpha*B0 at some
   alpha you name and why; (c) Fable's dispersion ladder; (d) none - freeze B0. Max 3 arms (~12 min each).
2. Is choosing alpha = 1.25 from the 168-target training-only LOO proxy admissible under our rule
   "never fit/select/tune on the 79"? If not, what fitting protocol instead?
3. Gate: is the converged packaging gate (score_avg > 0.021555 strictly, PDS >= 0.845, composite paired
   lower bound > 0, then an untouched confirmation fold) still right for this arm? What is the
   confirmation fold, concretely, given every H1 target is now either in the 168 or the 79?
4. B3 beat B0 by +0.0041. What do we do with that? Is it a candidate, and under what confirmation
   protocol, given it was observed on the 79? Should the alpha arm be built on B0, on B3, or both?
5. Should a gene-wise magnitude calibration in the P1a style (fitted on public/training truth only, holding
   direction and therefore PDS) replace or follow the scalar arm? It is a bigger build (~hours, CPU).
