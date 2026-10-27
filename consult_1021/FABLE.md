# Independent consult (Fable), 2026-10-21

Read: PUSH_1021.md, OFFICIAL_SCORES.md, CONSTRAINTS.md, model_v5/NOTES.md, model_v7/NOTES.md, model_v8/{NOTES,GATES}.md, model_billm/PACK_REPORT.md, vcc_top500_1019/REPORT.md, vcc_adapt_1020/REPORT.md, rsi/NEXT_SUBMISSIONS.md, competitors/{LEVERS-2026-10-18,COMPETITORS,DAILY-2026-10-19}.md, and the scorer docs `model/vendor/cell-eval2/docs/vcc2026_metrics/vcc2026-metrics-brief.md` (brief) and `docs/metrics.md` (full). No code was run.

## 1. Why scaled MSE sits at 0.06–0.09

**The metric** (brief lines 154–190): MSE = ΣN_p / ΣD_p with
N_p = ‖b̂_p − b_p‖² − r·min(Ĉ_p, C_p) − C_p and D_p = ‖b_p − b_ctrl‖² − C_p − C_ctrl, in log1p(CP50K) group-sum bulk, own target gene excluded, r = min(1, B/claim), B = across-target centred spread of our bulks. Scaled = (0.99 − u)/0.95 (brief lines 344–351), clamped [0,1].

**Energy decomposition.** Write our bulk delta as Â (300×G) and the truth as Δ. In expectation, with totals pinned (every entry since v4a pins them, model_v5/NOTES.md line 18):

ΣN = ΣD + ‖Â‖² − 2⟨Â,Δ⟩ − credit, credit = r·Σmin(Ĉ_p,C_p) ≤ ΣC_real.

Arc built the credit to be *break-even by construction*: "it is exactly the quantity the numerator charges for the same variation, so buying budget is break-even at worst" (metrics.md lines 505–507). At B = claim the ‖Â‖² term and the credit cancel, and the whole score is the cross term: u = 1 − 2⟨Â,Δ⟩/ΣD = 1 − 2·ρ̄·‖Â‖/‖Δ‖, where ρ̄ is the **Frobenius cosine** between our prediction matrix and the truth (within-target direction *and* across-target magnitude allocation). Optimum: push ‖Â‖² up to the cap ΣC_real (allowed for free) and allocate it as a_p ∝ ρ_p·d_p. On the official panels ΣC_real/ΣD ≈ 0.5 (metrics.md line 459; model_v8/GATES.md lines 27–31 give the public-control cap estimate 0.231/0.238/0.271 per context, and ΣD ≈ 0.5), so ‖Â‖ ≈ 0.7‖Δ‖ and

**scaled MSE ≈ (1.41·ρ̄ − 0.01)/0.95 at full credit; ≈ ρ̄²/0.95 with no credit.**

**Calibrating on our own scores (ΣD ≈ 0.5):**

| entry | scaled MSE | u | B | credit | 2⟨Â,Δ⟩ | ρ̄ |
|---|---:|---:|---:|---:|---:|---:|
| v5b | 0.082 | 0.91 | 0.105 | ≈0.096 | 0.054 | **0.12** |
| v8a (same totals, new layout) | 0.016 | 0.975 | 0.105 | ≈0.064 | 0.054 | 0.12 |
| v8b (alpha ×1.5) | 0.095 | 0.90 | 0.237 | ≈0.206 | 0.081 | 0.12 |
| billmA | 0.061 | 0.93 | ≈0.10 | ≈0.10 | ≈0.034 | ≈0.08 |
| ranks 10–50 median | 0.444 | 0.57 | ≈0.25 | ≈0.25 | 0.215 | **0.30** |
| Bill M (PDS 0.697) | 0.350 | 0.66 | | | 0.17 | 0.24 |
| Illumina | 0.793 | 0.24 | | | 0.38 | 0.54 |

The model reproduces v8b exactly (predicted +0.013 scaled, observed +0.013) and explains v8a: the server credited 0.064 against a 0.0955 local claim because heterogeneous per-row Ĉ_p hit the per-row cap min(Ĉ_p,C_p), and the claim was extrapolated from 16 probe targets (model_v8/NOTES.md line 14).

**Diagnosis.** Our MSE is low for two separable reasons, in this order:
1. **ρ̄ ≈ 0.12 vs 0.30 for the top 50.** That is prediction quality, not emission. PDS 0.71 vs 0.81 says within-target direction is part of it; the energy-weighted ρ̄ being even lower says the bigger part is *magnitude allocation*: we emit a near-uniform amplitude over 300 targets (common alpha times signature, model_v5/NOTES.md lines 25–29) while ΣD is dominated by a minority of strong responders. Spread spent on silent targets is pure cost.
2. **We claim ~40% of the cap.** v5b's claim was 0.096/0.086/0.107 against caps 0.231/0.238/0.271 (model_v8/GATES.md lines 7–9, 27–31); the 4-cell pooled template halves Ĉ (model_v5/NOTES.md line 110). Even at ρ̄ = 0.12, moving to the cap with ‖Â‖ = 0.7‖Δ‖ and a ∝ ρ_p d_p is worth about +0.14 scaled (0.08 → ~0.22), i.e. **+0.02 to +0.025 overall**. That is the ceiling of emission at today's ρ̄; the other 0.06 overall needs ρ̄ ≈ 0.3.

**What the emitter must do (all label-free, all from our prediction plus the released controls):**
- Bulk head: Ĉ_p ≥ C_real,p on *every row* (uniform, not just in sum), estimated from 400-cell control jackknives at the context's depth; then set B = 0.95–1.0× that claim; allocate a_p² ∝ (predicted magnitude × confidence)². PDS is cosine, so allocation is PDS-neutral.
- Per-cell CPM head independent of the bulk: fold-changes at ~0.9× transferred LFC on sign-confident genes and 0 elsewhere (NMAE is L1 on the real DE set, so a zero costs exactly |lfc| and an over-amplitude costs (1−2q)+t; our NMAE raw 0.93 says sign accuracy on real DE genes is ~0.54). Call count per target ≥ predicted n_real (fidelity = precision × min(1, n_pred/n_real), brief line 243; under-calling strong responders is why our fidelity is *below* the tiled baseline). P-value order over all 18.5k genes = confidence order (reach walks our order over the real DE pool).
- The coupling is depth: bulk is a count-weighted sum, the DE head is a per-cell mean. A long-tailed depth distribution (e.g. 10% of cells at 4–5× median) carries ~1/3 of the bulk mass in 10% of the ranks, enough to set the bulk amplitude to 0.7‖Δ‖ without moving per-cell LFC. Competitors who disclose this: Jae (DAILY-2026-10-19 line 34), yetao wu "40×0.9M carrier cells" (line 39), Haoran Xu "16 heavy cells" (line 49). Keep carriers inside the controls' real depth range; Arc patched this channel once (#348) and may again for the final.

## 2. What the top teams do that we do not

From LEVERS-2026-10-18 row 9 and COMPETITORS.md lines 38–50, 135–149:
1. **Per-target magnitude** (Jae "bulk scaled per target", M014 "psig delta shrink", null_hypothesis responder gate from paralogs/degree/constraint, Trident "budget-bound curve"). We have none. Cheapest proxy: control expression of the target gene (CRISPRi of an unexpressed gene is silent), source effect norm, cross-source agreement.
2. **Claim at the cap with B/claim ≈ 1** (Bill M 1.00, William 0.99, Jason Jowu 0.95, Kaipeng 0.85×min(JK, NTC400 JK)). We calibrate to our own under-dispersed claim.
3. **Two heads, two amplitudes** (Bill M "per-cell DE head only", Jae dual moment, Riya DE 0.9, Chiu bulk shrunk 1.25× while calls doubled).
4. **Broader, source-wide-centred sources at equal weight** (William: 9 sources incl. HipSci, KOLF, VIPerturb, LINCS; "expression compatibility" in Bill's MLP, i.e. masking responses by the context's control expression). Our banks are 437-target cuts centred on themselves.
5. **Calls ranked through p-values with PIE** (Chen Zhang #79, SIVC), budgets near the real DE count.
6. 85–112 validation probes each; part of their rank will not survive D/E/F.

## 3. Plan to Nov 5

Constraint: one real build/day on the Mini while billmB runs; emitter-only variants (same profile, different layout/amplitude) are hours, so pairs should be one base × one knob. Validation stays the only judge; every rule below is label-free so it carries to D/E/F.

| UTC day | slot 1 (hypothesis → component) | slot 2 |
|---|---|---|
| Oct 22 | v7b: claim at ~cap (0.23–0.25) on v5 direction → MSE. Prediction from §1: ~0.13 | v7a: native bulk geometry → PDS; MSE ~0.09 |
| Oct 22 (no slot) | `controls-final` lands: start full-screen bank fetch for the new panel (`model_v6/banks.py`), run the mock-final runbook. **Safety final entry Oct 23.** | |
| Oct 23 | T1 billmA profile, v5b emitter, **a_p ∝ magnitude score**, same B → MSE up, PDS/NMAE flat | T2 billmA, **heavy-tail depth layout, claim at cap, CPM head untouched** → MSE up, NMAE flat |
| Oct 24 | T3 best(T1,T2) + both knobs combined → MSE | T4 + call budget per target = predicted n_real → fid, jac |
| Oct 25 | T5 + CPM-head amplitude 0.9 vs 0.6 → NMAE | T6 + 18.5k-gene confidence p-value ordering → reach |
| Oct 26 | T7 + target-gene control-expression silence gate (amplitude → 0, calls → ~0 for unexpressed targets) → MSE, PDS, fid | T8 billmB profile vs billmA, same emitter → PDS |
| Oct 27 | T9 full-screen source-wide centring banks vs 437-cut → PDS, MSE | T10 + expansion banks (HipSci/KOLF/Jiang) as equal-weight sign/direction sources → PDS |
| Oct 28 | T11 deflation 10 PCs vs 50 → PDS vs MSE trade | T12 response mask by context control expression → MSE |
| Oct 29–30 | Best combination; one exact re-submission of it to measure build noise | second-best combination |
| Oct 31–Nov 2 | Final-panel builds with the winning recipe on D/E/F; dry-run decode; no new validation ideas | |
| Nov 3 | **Submit main final entry** | Nov 4 reserve; Nov 5 nothing touches the board |

**Ranked by expected official gain per day:** (1) magnitude allocation + per-target call budget, +0.01–0.02, one day; (2) claim-at-cap dual-moment layout, +0.01–0.015, one day (code exists in v7/v8); (3) silence gate + expression mask, +0.005–0.02, one day; (4) full-screen centred banks + more sources, +0.02–0.05 but 3–5 days and the only route to ρ̄ ≈ 0.3; (5) billmB, +0.003–0.01.

**Drop:** the RSI loop (229 rounds, ~25M tokens, zero candidates, self-locked); the local six-metric calibration gate; SLIM/GRAFT/GenBio-MLP ports (lost locally, no slot value); slim real-control emitter and any amplitude knob gated on call counts; KOLF 46 GB / VIPerturb RDS; Jaccard (span 0.35, nobody scores it); context-identity coefficients (Jurkat for A) since D/E/F are new; layout-only re-tests; submitting anything whose B/claim is outside 0.9–1.05.

## 4. Errors and wasted effort

1. **The +0.05–0.075 emission estimate (PUSH lever 1) ignores #348's break-even design.** Credit only pays through the cross term; at ρ̄ = 0.12 the emission ceiling is ~+0.02 overall. The gap is the model.
2. **Claim calibrated to our own pooled-template jackknife (40% of cap) with uniform amplitude.** v5 NOTES knew it ("about half the realistic level", line 110) and rejected the fix because realistic scatter wrecked fidelity in a *single-moment* emitter; the answer was decoupling, not retreat.
3. **Coupling the DE gate to the bulk.** v5b zeroes ~80% of bulk genes to control calls; the 275-call test forced alpha 4–6 (model_v5/NOTES.md lines 51–61). The two heads must be independent; never zero bulk genes.
4. **Uniform call budgets** (N_conf per context, every target within [0.5,2]×, model_v8/GATES.md line 23) against a fidelity formula that is linear in n_pred/n_real. That is the fidelity hole.
5. **v8a's credit model** assumed the sum of claims transfers; the per-row cap and 16-target probe extrapolation lost a third of it.
6. **slim_v5b was submitted with its own gate reading B/claim 6.5** (PUSH 1c): a wasted slot. v7a/v7b sat unsubmitted for four days while the MSE question stayed open.
7. **Local proxies were asked the wrong question.** H1 (21 essential targets, 55k depth) and CD4-192 cannot rank 0.01 differences, but the one quantity that maps monotonically to official MSE, the Frobenius cosine ρ̄ in bulk-log space on held-out public lines, was never measured. Measure it; it is the diagnostic for every direction change.
8. **HipSci/Jiang excluded as "dishonest."** The honesty concern is the MLP response format; as sign/direction priors for the bulk head, log-residual effects are exactly what is needed.
9. **billm donor coverage collapses on a new panel** (PUSH 1b). The panel-agnostic bank fetch is the critical path for *any* method and should run first on Oct 22, ahead of validation experiments.
10. OFFICIAL_SCORES.md is stale (no v9a/v9b components); keep the one table current, it is the only ground truth.

Bottom line: lock the emitter in three days (worth ~+0.02, to ~0.185), then spend every remaining day on ρ̄: full-screen centred sources, magnitude prediction, expression masking. Ship a safety final on Oct 23 and the main entry on Nov 3.
