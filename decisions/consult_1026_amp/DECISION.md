> Historical AI decision record. “Official runs” below are local scorer runs,
> not new validation submissions. Premises are corrected in later rounds;
> use the final rulings and linked research update for current interpretation.

# DECISION — consult_1026_amp (amplitude-form conflict). Round 1: 2026-10-26 14:16 UTC. Round 2: 14:22 UTC. Round 3: pending.

Advisors (policy 2026-10-26): Fable 5.1 and Opus 5.5 via `claude -p` (read-only tools), Astra via codex
(read-only sandbox). Kimi was not consulted in this series. Briefs: brief_final.md (R1), brief2.md (R2).

## Round 1 — unanimous
- The amplitude experiment is **mean-extrapolation at alpha > 1** (`mu_168 + alpha(B - mu_168)`), alpha chosen
  training-only on the 168 (pooled() surrogate, LOO mean). Admissible as a fitted parameter; the 79 are now a
  development panel, not independent confirmation (measurement 4 was viewed first — disclosed).
- **Dropped:** pure scaling (PDS-inert), Fable's shrink ladder and dispersion ladder (Fable withdrew both), any alpha < 1.
- 79-gate retained: score_avg > 0.021555 strictly, raw PDS >= 0.845, paired composite bootstrap lower bound > 0.
- **P1a-style gene-wise magnitude calibration follows and is built in parallel now**, fitted on 168 + public non-H1
  only, never the 79; PDS must be re-measured (gene-wise rescaling rotates the vector). -> launched as p1a_1026/.
- B3-B0 paired bootstrap to be computed (descriptive).

## Round 1 — split
| question | Fable | Opus | Astra |
|---|---|---|---|
| arms on 79 | A=mu+a*(B3-mu), B=mu+a*(B0-mu) | N1 = norm-matched mu+1.25(B0-mu) (only candidate) + raw A1 descriptive | E0=1.25B0-.25mu, E3=1.25B3-.25mu |
| alpha | re-fit on fine grid WITH leaderboard MSE priced; "argmax 1.00 => alpha line dead" | 1.25 frozen | 1.25 frozen |
| B3 | candidate (fallback / base of A) | not a candidate (+0.004 = noise floor) | candidate |
| confirmation | CD4 fold, then a leaderboard submission (Ricardo only) | cross-fitted official run on the 168 | none exists in the 247 -> no promotion, B0 ships |

## Measurement that contradicts the majority (put back in Round 2)
Fable's MSE-priced objective on the 168 (results/amp_fable_objective168.json): **argmax = 1.00** for raw
extrapolation — a 5% response-MSE rise zeroes our leaderboard scaled MSE (0.061, ~0.0102 pts), more than the
best PDS gain (+0.0087). Opus's norm-matched N1 rises only 1.1% in MSE at 1.25 (results/amp_n1_mse168.json):
Obj(N1, 1.25) = **+0.0070**, positive and the argmax of {1.10, 1.25, 1.50}. Caveat: 0.061 is a validation-panel
figure; the final-panel MSE term is hidden.

## Executed / blocked
- Leak check (Fable): pooled() excludes all H1 banks -> PASS.
- Local scorer runs temporarily blocked by the disk reserve (29.7 GB < 30 GB guard).
  Unrelated machine-workload details are omitted from this publication.
- No official arm launched pending Round 2 ruling on the form.

## Round 2 — 14:22 UTC (fable_r2.md, opus_r2.md, astra_r2.md)
**Unanimous:**
- **N1 only**, alpha frozen 1.25: `A = mu168 + 1.25(B0 - mu168)`, `N1_i = A_i * ||B0_i||/||A_i||`. Raw E0/E3 and Fable's
  A/B are dead (Fable holds to its pre-declared argmax-1.00 rule). No B3-based N1, no grid widening, no re-optimising.
- Kill rule (Opus, agreed by Astra/Fable): response-MSE ratio N1/B0 < 1.055 on the 79. **Measured 1.00805: ALIVE.**
- N1 79-gate: score_avg > 0.021555; PDS >= 0.845; paired composite LB > 0; **official raw-MSE ratio N1/B0 <= 1.03** (Fable/Opus).
- **P1a is the primary development lever**; N1 then P1a in fixed order; fit P1a on the leading base; MSE-priced
  objective incl. identity; raw MSE ratio <= 0.97 vs base (Fable) / r < 1 and Obj > 0 (Opus); rank/168 must not worsen
  by > 0.5 (Opus); fit set 168 + public non-H1 only; Nov 2 cutoff. (Constraints forwarded to the p1a_1026 build.)
- B3 is a candidate only with a paired composite LB > 0 (Fable/Opus; Astra: full gate incl. LB).

**Executed:** N1 built and hash-locked (h2h2_1026/AMENDMENT_N1.md, results/n1_build.json, score_n1.py — runner is
score_run.py with exactly two textual substitutions). B3-B0 bootstrap: composite +0.00406, CI [-0.00197, +0.01048]
-> **B3 DROPPED**. PDS CI [+0.00097, +0.01688] (descriptive). N1 official run queued (local queue) behind
the disk guard; then N1-B0 bootstrap; then B5 and the locked chain.

**Split -> Round 3:** confirmation. Fable: CD4 fold (b). Opus, Astra: Ricardo leaderboard submission (d).
Correction put back: (d) is infeasible — final-phase scores are hidden and the validation board is frozen.

## Round 3 — confirmation (fable_r3.md, opus_r3.md, astra_r3.md)
**Unanimous:** route (b) CD4 fold check with the locked N1 transform, no retuning; fallback (a) only if no legal CD4
training mean: 79-gate + composite >= B0 + 0.008 + raw-MSE ratio <= 1.03. No improvised mean from other cell types.
Final-phase rule verified: **last entry counts** (CONSTRAINTS.md) -> N1 would replace B0 blind (Fable's branch).
Pass rule (majority-composed) registered in h2h2_1026/AMENDMENT_N1_CD4.md before any CD4 build: composite point > 0,
raw-MSE ratio <= 1.03, PDS >= B0 - 0.005; Astra's LB > 0 reported. CD4 runs only after N1 passes on the 79 (Opus).
Legal mu_CD4 looks available: CD4 banks hold ~1,690 targets per condition, 15 overlap the CD4 search fold (excluded).

## RESULT — N1 official (2026-10-26 15:27 UTC): **FAIL** (mechanical, AMENDMENT_N1.md; results/n1_gate.json)
| | B0 | N1 |
|---|---:|---:|
| score_avg | 0.021555 | **0.012746** |
| raw PDS | 0.85005 | **0.83966** (< 0.845 floor) |
| raw fidelity | 0.44872 | 0.44518 |
| raw MSE | 1.41072 | 1.42998 (ratio 1.0137, <= 1.03 ok) |
Paired N1-B0: composite -0.00881 [-0.01790, +0.00067]; PDS -0.01039 [-0.02840, +0.00779]; fidelity -0.00354
[-0.00529, -0.00182]. Fails score, PDS floor and LB. CD4 not run (AMENDMENT_N1_CD4: only after a 79 pass).
(The first gate run crashed on a SyntaxError in n1_gate.py; fixed and re-run; verdict unchanged.)

**Why — the proxy was right about the prediction and wrong about the score: the renderer.**
Official-form PDS (cosine, all 79 target genes masked, official tie-rank; results/diag_n1_pds.json):
| arm | pre-render (response vs truth) | official (rendered) | renderer adds |
|---|---:|---:|---:|
| B0 | 0.7960 | 0.8500 | +0.054 |
| B3 | 0.7965 | 0.8586 | +0.062 |
| N1 | **0.8229** | 0.8397 | +0.017 |
| B5 | 0.6204 | 0.6761 | +0.056 |
| P1aB0 | 0.8103 | (pending) | |
| P1aN1 | 0.8322 | (pending) | |
N1 genuinely sharpens the predicted direction (+0.027 pre-render) but loses most of the renderer's PDS bonus.
Per-target N1-B0 deltas correlate only 0.55 between pre-render and official. Response-space proxies are NOT
reliable for small direction changes; only official runs count. The renderer (anchor emission + confidence-ordered
885-call layout, which uses donor z-scores agreeing with the predicted SIGN; N1 flips ~4% of signs) is itself a
~0.05-0.06 PDS lever that nobody has mapped. -> for the advisors, with the P1a official results.

Note (Codex lane, 2026-10-26): read B3 as PDS evidence only — paired PDS LB +0.000974 > 0, composite LB -0.001967 < 0;
the +0.004 composite is noise, not signal. Consistent with "B3 DROPPED" above.

## RESULT — P1a official (2026-10-26 ~16:50 UTC): **FAIL, both arms** (AMENDMENT_P1A.md; results/p1a_gate.json)
| arm | score_avg | PDS | fidelity | raw MSE | composite vs base [95% CI] |
|---|---:|---:|---:|---:|---|
| B0 | 0.021555 | 0.8500 | 0.4487 | 1.4107 | — |
| P1aB0 | 0.005789 | 0.7970 | 0.4648 | 0.9628 (ratio 0.682) | -0.0158 [-0.0311, +0.0015] |
| P1aN1 | -0.002919 | 0.7842 | 0.4632 | 0.9681 (ratio 0.677) | -0.0157 [-0.0289, -0.0018] vs N1 |
MSE ratio criterion passes; score, PDS floor and LB fail. Scaled MSE still 0: raw 0.9628 missed the 0.9613 anchor.

**Causal finding (overturns the "pure scaling is PDS-inert" premise of R1-R2):** P1aB0 has B0's signs exactly (100%),
so the renderer's sign-dependent confidence/layout is identical to B0's; only amplitude changed (0.28x norm).
Pre-render official-form PDS ROSE (0.810 vs 0.796) but rendered PDS FELL 0.053. => **Rendered PDS depends on
predicted amplitude** (renderer signal-to-sampling-noise); the response-space invariance does not hold in the
official metric. Shrinking also RAISED fidelity (+0.016), opposite to mse_1025's CD4 scalar. Put back in round 4.

## Round 4 — 16:45 UTC (fable_r4.md, opus_r4.md, astra_r4.md)
**Unanimous:** freeze B0 as the submission incumbent/default now; research continues, bounded; pure-scaling ladder on
B0 runs; renderer becomes the research target. Astra explicitly withdrew the official scale-invariance claim.
**Majority / median composition (executed):**
- Ladder is an exploration MAP, no promotion (Fable, Astra). Alphas = per-slot median {0.75, 1.25, 1.5}.
- First renderer experiment = render-seed replicate of B0 (Fable, Opus): 2 new seeds -> noise floor for every
  verdict so far (bootstrap CIs resample targets, not renders). Astra's uniform-confidence ablation and Fable's 2x
  budget are the next renderer experiments, after the noise floor.
- Cutoff 11-02 23:59 UTC (Opus). Fable's "validation-leaderboard check" is infeasible (board frozen; R3 correction).
Registered: h2h2_1026/AMENDMENT_RENDER.md. Queued locally; machine scheduling details omitted.
results/render_summary.md). Preflight 16:46Z: 35.7 GB free, swap 2.5/4.0 GB.

## RESULT — AMENDMENT_RENDER official (2026-10-26 ~18:00 UTC; h2h2_1026/results/render_summary.md)
| arm | score_avg | PDS | fidelity | raw MSE | reach | composite vs B0 [95% CI] |
|---|---:|---:|---:|---:|---:|---|
| B0 | +0.021555 | 0.8500 | 0.4487 | 1.4107 | 0.1266 | — |
| B0s1 (seed 20261101) | +0.022318 | 0.8500 | 0.4479 | 1.4107 | 0.1340 | +0.00076 [-0.00219, +0.00538] |
| B0s2 (seed 20261102) | +0.022064 | 0.8500 | 0.4485 | 1.4107 | 0.1296 | +0.00051 [-0.00286, +0.00544] |
| A075 | +0.018126 | 0.8476 | 0.4473 | 1.1491 | 0.1207 | -0.00343 [-0.00869, +0.00255] |
| A125 | +0.021829 | 0.8523 | 0.4513 | 1.7593 | 0.1327 | +0.00027 [-0.00568, +0.00654] |
| A150 | +0.021379 | 0.8541 | 0.4547 | 2.1934 | 0.1497 | -0.00018 [-0.00948, +0.00885] |

**Noise floor:** the rendered pseudobulk (count_sums, depth) is bit-identical across render seeds, so official PDS
and MSE are deterministic given the prediction; the seed only resamples cells, moving the DE metrics (per-target
reach jumps up to 0.83, fidelity <= 0.017). Seed |dPDS| = 0 < 0.01, so no earlier verdict is flagged as render noise
(N1's -0.0104 PDS loss is real). Composite seed noise ~0.0008, all from DE.
**Amplitude map:** PDS rises monotonically with alpha (0.8476 / 0.8500 / 0.8523 / 0.8541), fidelity and reach too,
but composite is flat (all CIs straddle 0) and raw MSE rises 1.41 -> 1.76 -> 2.19. Opus's alpha>1 gate fails for
A125 and A150 (composite LB < 0; raw MSE > B0's). Scaling is not a lever. B0 stays the incumbent. -> round 5.

## Round 5 — 17:54 UTC (fable_r5.md, opus_r5.md, astra_r5.md)
**All three:** the layout preserves per-gene column sums and depth, so budget and confidence can't move PDS or raw
MSE; only reach, fidelity and jaccard can move. Renderer changes are legitimate if the spec allows them, but
legitimacy is not evidence. B0 stays frozen. Astra retracted the sampling-SNR explanation.
**Q1, majority (Fable, Astra):** run the uniform-confidence ablation at 885 calls as a diagnostic. Opus (minority)
wanted 2x budget; it runs only if the ablation shows the layout matters.
**Q3, composition:** one diagnostic, then stop unless it is actionable. Fable: stop after it. Opus: one experiment
with a stop rule. Astra: continue bounded, but stop early if nothing is actionable.
**Promotion bar, if it ever arises:** Astra's untouched confirmation targets plus Opus's Q2 list (LB > 0, two more
render seeds each beating B0 by more than 0.0008, PDS and MSE bit-identical, gain spread across targets, format
validation, Ricardo signs off).
Registered: h2h2_1026/AMENDMENT_RENDER_U.md (Fable's ±0.002 preregistered reading; Opus's budget-provenance fix).
Queued locally; machine scheduling details omitted.

## RESULT — RENDER_U official (2026-10-26 18:11 UTC; h2h2_1026/results/u_gate.md)
U885 score_avg 0.016127 vs B0 0.021555; composite d -0.00543 [-0.00982, -0.00056]. Fidelity -0.0050, jaccard -0.0037,
reach -0.0034, nmae ~0; PDS bit-identical (0.8500487); rendered pseudobulk bit-identical; counts changed on 79/79
targets; median called-set overlap with B0 0.540.
**Literal preregistered reading: INVALID**, because raw MSE is not bit-identical (1.40920 vs 1.41072, -0.11%).
**Correction (orchestrator's premise error, put back in round 6):** the official expr_mse_unbiased subtracts a
per-cell variance term, so the cell layout moves it even with identical pseudobulk. The seed replicates moved it too
(B0s1 -2.6e-5 norm). "PDS and MSE deterministic" (R5 brief) holds for PDS only. Scaled MSE is 0 for both arms
(clamped), so the composite delta is entirely DE. If validity is re-based on pseudobulk + PDS bit-identity, the
reading is d < -0.002: the confidence allocation contributes. That explains B0 but gives no entry.

## Round 6 — 18:12 UTC (fable_r6.md, opus_r6.md, astra_r6.md) — unanimous
Q1: accept U885 under a corrected check. This is a post-result correction of a mistaken invariant; the original
clause is kept in the audit trail. The branch d < -0.002 fired: confidence allocation contributes at 885 calls.
Size uncertain: the CI's upper end, -0.00056, is inside ±0.002 (Opus).
Q2: run the 2x budget once (D1770, B0 confidence) under Opus's binding stop rule (point < +0.002 or LB <= 0 ->
CLOSE, B0 final). No further arms; a pass is diagnostic only.
Q3, per-item median: pseudobulk and PDS bit-identical; scaled-MSE contribution unchanged; raw MSE reported and
no more than 0.1% worse (Opus); raw and scaled deltas reported separately (Astra).
Registered: h2h2_1026/AMENDMENT_RENDER_D.md. Queued locally; machine scheduling details omitted.
~18:16Z, 35.4 GB free, swap 2.1/3.0 GB.

## RESULT — RENDER_D official (2026-10-26 18:29 UTC; h2h2_1026/results/d_gate.md): **PASS (diagnostic only)**
D1770 score_avg 0.056045 vs B0 0.021555; composite d +0.03449 [+0.02194, +0.04770]. Validity holds: pseudobulk and PDS
bit-identical, scaled MSE 0 for both, raw MSE ratio 0.9942. Calls per target min/median/max 378/1655/1770.
Raw components: fidelity 0.4902 vs 0.4487, jaccard 0.0592 vs 0.0430, reach 0.1343 vs 0.1266, nmae flat.
Fidelity per target: 39 better and 38 worse (median 0); the top-5 targets give 0.38 of the summed gain; reach's
gain comes from a few targets.
**Conflicting prior evidence (orchestrator, from the record):**
- VCC_OVERNIGHT_LOG 21:45 (10-25): the local H1 model of 1200->3500 calls is monotone and points at the dense
  endpoint, "which is already an official measured loss (v5b native ~4,475 = 0.1528 vs v9a 885 = 0.16135)".
- The local harness reversed the sign of the official v5b-vs-v9a overall comparison.
- layout_cd4_1025: CD4 replicated the local native-layout gain, so the second recipient shares the bias.
- FINAL_PUSH 10-22: "Local fidelity magnitudes are exaggerated (H1 median n_real 107 vs budgets 694-885)".
Put to the advisors in round 7.

## Round 7 — 18:31 UTC (fable_r7.md, opus_r7.md, astra_r7.md) — **UNANIMOUS: B0 FINAL; D1770 not promoted**
Q1: D1770's local +0.0345 carries no evidential weight for the hidden test. The harness mis-signed the only official
density comparison (v5b dense lost 0.1528 vs v9a 885 at 0.16135), and D1770 sits on that same axis. The gain comes
from a few large under-called targets, which is exactly what changes if the hidden test's n_real distribution differs.
Q2: no honest confirmation exists before 11-05. The board is frozen. CD4 shares the instrument. Calibrating against
the single v5b/v9a pair would fit one point (Opus); at most it could disqualify (Fable).
Q3: keep B0 final. Do not present D1770 as a candidate; keep it as a diagnostic that passed and was not promoted.
Caveat if Ricardo asks (Opus): "The only official density comparison we have says more calls lost, and this local
gain comes from the same harness that got that comparison backwards."
**Status:** renderer and density research CLOSED. Scaling CLOSED. Every lever since 10-26 is either REJECTED or a
diagnostic that was not promoted. Frozen entry 5NeeoWn96lfm23hgExcm stays last; a local integrity check reported matching hash and size (17:24Z).

## CORRECTION to the R7 brief (orchestrator, 2026-10-26 ~19:00 UTC, after Ricardo asked for the harness calibration)
The R7 brief said the local harness "REVERSED the sign" of the official v5b-vs-v9a comparison. That is true for the
overall composite only. Per component (official PUSH_1021.md / OFFICIAL_SCORES.md vs local
final_push/artifacts/de_check/H1_rich_search/calibration_summary.json), v5b minus v9a(885):
| component | official scaled d | local scaled d | ratio local/official |
|---|---:|---:|---:|
| fidelity | +0.037 | +0.343 | ~9x |
| jaccard | +0.016 | +0.113 | ~7x |
| reach | -0.104 | -0.064 | ~0.6x |
| MSE | +0.001 | 0 (clamped locally) | blind |
| overall | -0.0085 | +0.065 | sign flipped |
All four component signs agree. The overall flip comes from magnitudes: fidelity and jaccard are exaggerated, reach
is understated. FINAL_PUSH 10-22 already said the harness is "usable for call-set levers (budget, ordering)" on
signs. D1770 moves reach UP as well (unlike the dense native layout), so the R7 inference "same axis -> same
official loss" does not follow automatically. A component-level calibration is running; it goes back to the
advisors with the result.

## Round 8 — 18:51 UTC (fable_r8.md, opus_r8.md, astra_r8.md), after the R7 correction and calib3_1026
Q1 (unanimous): the ruling is unchanged. B0 is final and D1770 is not a candidate. The R7 rationale is withdrawn:
D1770 is now labelled "indeterminate, at the seed-noise floor (central about +0.002, ~60% positive)", not a
known loss.
Q2 (majority, Fable + Opus; Astra dissents): do NOT run the CD4 D1770 test. No outcome changes the action, the
CD4 transfer is uncalibrated, and there is no CD4 seed floor. If it is ever run, the rule needs a CD4 seed
replicate, a count of targets whose selection actually changes, official slopes plus composite weights, and
reach/NMAE/MSE as veto evidence.
Q3: only an official score could justify presenting D1770. The validation board is frozen, so none is possible
before 11-05. Ricardo gets it as a one-line note, not as an option (Opus): "D1770: about 60% positive, central
+0.002, inside seed noise. Not recommended. B0 final."

## Round 9 — 18:57 UTC (fable_r9.md, opus_r9.md, astra_r9.md): progress is plausible; three call-side arms
Arms, by composition:
- K250 (budget 250): all three.
- RH20 (top-20 reach head, call set unchanged): all three.
- DET885 (Skadi-style detectability selection at 885): Fable and Opus. Astra's K100 is the minority and not run.
Instrument: per-stratum paired raw deltas, official slopes (fid 3.30, jac 2.88), reach and nmae at local scale
(conservative), mix 31.5/65/3.5 plus a mixture sweep (Astra), seed floor transformed the same way, and CD4 as the
second panel for any survivor. It is a screen, not a decision (Fable, Astra).
Bar: screen PASS = DE-only validity + central >= +0.005 and > 2x seed floor; then CD4 same sign; then "B0 + arm" to
Ricardo with an honest probability (Opus). Fable's official probe is impossible (board frozen).
Registered: h2h2_1026/AMENDMENT_R9.md. Instrument dry run: D1770 central +0.0025 (calib3: +0.0026), U885 -0.0040.
Queued locally; machine scheduling details omitted.

## RESULT — R9 official (2026-10-26 19:40 UTC; h2h2_1026/results/r9_gate.md): all three FAIL / VOID
| arm | local d [CI] | central proj | conservative | weak+mid fid better/worse | MSE ratio | screen |
|---|---|---:|---:|---|---:|---|
| RH20 | -0.00002 [-0.00017, +0.00013] | -0.00000 | -0.00002 | 0.53/0.38 | 1.00001 | FAIL (inert) |
| K250 | -0.0701 [-0.0884, -0.0524] | -0.0379 | -0.0503 | 0.45/0.55 | 1.00543 | VOID / large loss |
| DET885 | -0.0028 [-0.0062, +0.0006] | -0.0060 | -0.0068 | 0.13/0.87 | 1.00420 | VOID / loss |
Seed floor +0.0028. K250 per stratum (weak/mid/strong): fid +0.010/-0.074/-0.164, reach -0.091/-0.019/-0.002.
DET885: fid -0.0135/-0.0124/+0.0037. The sign-audit inference ("fewer, better-ranked calls") is refuted: B0's
confidence order beats detectability, and trimming costs mid-target coverage and weak-target reach.
Orchestrator check (not a lever): local raw MSE does not track official MSE. v8b_like vs v9a_694 local raw
0.98617 vs 0.98501 (worse), while official v8b vs v8a gained +0.079 scaled. Local MSE is unusable, both clamped and
mis-signed on amplitude.

## Round 10 — 19:44 UTC (fable_r10.md, opus_r10.md, astra_r10.md): UNANIMOUS yes, one bounded bet on decoupled amplitude
The only official datum above +0.01 is v8a->v8b: PDS + MSE worth about +0.018 overall before the NMAE collapse.
Composition (orchestrator; flagged for R11):
- Step 0, Fable's kill test (rho >= 0.9 -> dead): answered from the record. The emitted B/claim of the official
  billmA/v9a packages is 0.57-0.62 (rsi/NEXT_SUBMISSIONS.md), and v9a_885 is 0.49 locally. B binds, so the lever
  is alive.
- Mechanism: Astra's, chosen because it alone answers the risks Fable and Opus each named (library-size drift
  moving NMAE; changed within-group dispersion moving the Wilcoxon ranks). Keep B0's per-cell CPM vectors fixed
  and change only the cells' relative depths. Exact form: integer multipliers x'_c = k_c * x_c, which leaves the
  CPM vectors bit-identical, so NMAE and every DE input are unchanged by construction; only the pooled profile
  (PDS, MSE) moves.
- Gate (per-item median): DE inputs bit-identical (stronger than everyone's tolerance); PDS change >= 0; label-free
  B/claim in [0.875, 0.97]; projected MSE = 0.079 x (B/claim achieved - B0's)/0.5 (Opus; Fable's halved and
  Astra's zero-credit cases reported); central projected overall >= +0.010 (Opus, Astra; Fable +0.008).
- Timeline: 48 h feasibility (Astra) by 10-28 ~20:00 UTC; hard stop 11-01 23:59 UTC (Opus). A pass earns further
  validation, not promotion. Any submission is Ricardo's decision.

## Round 11 — 19:53 UTC (fable_r11.md, opus_r11.md, astra_r11.md): STOP (2-1). B0 FINAL.
Majority (Fable, Opus): do not build the depth-tilt arm. Opus's arithmetic from the R10 gate itself, with emitted
B/claim 0.57-0.62 and the 0.97 cap: max MSE credit 0.079 x (0.97 - 0.60)/0.5 ~ 0.058 scaled ~ +0.0097 overall;
PDS from A125/A150 adds about +0.001. The central ceiling is about +0.011 against a +0.010 bar, and only with zero
DE cost; Fable's halved case is +0.006, Astra's zero-credit case +0.001. Fable: the tilt is mean-preserving spread
on the called genes, which pushes the Wilcoxon toward 0.5 (the 10-22 kernel failure again, larger). Neither a pass
nor a fail could be checked before 11-05 (board frozen; local MSE mis-signed on amplitude).
Astra (minority): build one arm at s ~ 1.38 with a 48 h stop.
If Ricardo overrides: s ~ 1.38 (Opus, Astra; 1.25 cannot clear the gate); do not exclude called genes wholesale
(all three); DE within +/-0.003 on two seeds.
**Closed lines (2026-10-26):** manifold/STRING, magnitude table (P1a), N1 direction, scaling, renderer confidence,
call density (D1770, K250), reach head (RH20), detectability (DET885), decoupled amplitude (dm_1026). B0 = billmA
0.16438, entry 5NeeoWn96lfm23hgExcm stays last; local integrity checks remain in place.

## R12 (2026-10-26 22:28 UTC) — biology line: unanimous close, B0 final
Brief: brief12.md (four Codex diagnostics in bio_gap_1026/: decomposition, shared response, in-context alignment, LOCO map; plus h2h2_1026/scratch/mse_instr: no local MSE instrument).
- **Fable:** arithmetic right in form but the 0.811 conversion is unanchored (no instrument). Ceiling +0.0007 to +0.0015, under the seed floor; PDS downside asymmetric (CI −0.0052 → ~−0.0065 official). No missed lever: in a blind context, direction can come only from control structure (fails) or source priors (already used). B0 final.
- **Opus:** conclusion holds; +0.0044 scaled is an upper bound, because cross is normalised by delta energy, not the per-cell raw MSE. control_bounded is ~4× below the seed floor. GO/network priors would be another family-correction-style shrinkage. Close; B0 final. Note finding 4 in any write-up.
- **Astra:** +0.0004 is a scenario, not an estimate or a bound (the change in mean energy is not cancelled by matched spread). control_bounded is a small transferable diagnostic gain, not worth production work. Better donors and control-conditioned pathway interactions remain post-challenge hypotheses; the large gains were demonstrated with in-context truth, not proven to require it. B0 final.
**Ruling (3-0): close the biology line for this entry; no new arm before 11-03; B0 = billmA (validation 0.16438) final.**
Closed lines now: manifold/STRING, P1a, N1, scaling, renderer confidence, call density (D1770, K250), reach head, detectability, decoupled amplitude, biology (family correction, shared response, in-context map, LOCO map, coverage).
