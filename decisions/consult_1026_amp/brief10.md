Round 10 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = billmA 0.16438; final phase blind; nothing is submitted without Ricardo). <= 300 words.
Read h2h2_1026/results/r9_gate.md, competitors/IDEAS-2026-10-26.md (idea 6) and OFFICIAL_SCORES.md (v8a/v8b rows).

## R9 results: all three call-side arms failed
RH20 inert (0.0000). K250: local -0.070, projected -0.038 (trimming kills mid-target fidelity and weak-target
reach). DET885: projected -0.006, fidelity worse on 87% of weak+mid targets. B0's confidence order is better than
detectability. The DE call side looks exhausted at the +/-0.005 scale.

## Where the gap actually is
billmA official scaled: pds 0.706, mse 0.061, nmae 0.116, fid -0.066, reach 0.186, jac -0.016. Board: Skadi scaled
MSE .434 (raw normalised top-60 0.46-0.61; ours ~0.91). T3V200 "Res" (0.261) beats us mainly on MSE. An MSE gap of
~0.3 scaled is ~0.05 overall: ten times any DE lever.

## The one official datum on this axis
v8b vs v8a (1.5x amplitude, B=claim): pds +0.028, mse +0.079, nmae -0.161 (collapse), fid +0.004, reach +0.016;
overall -0.006. cell-eval2 docs (IDEAS #6): MSE and PDS read the PSEUDOBULK (group count sums); NMAE and DE read
per-cell CPM arithmetic means and Wilcoxon ranks. The MSE credit is budget-limited: rho = min(1, B / sum min(C^,C)),
with B = OUR across-perturbation pseudobulk spread; emitted B/claim ~0.5 under the v9 layout (known since 10-22).
Hypothesis: raise pseudobulk spread (toward B ~ 0.95 claim) while holding per-cell CPM means, i.e. decouple
the two channels (depth-correlated count placement, "dual-moment"). If NMAE holds, v8b's PDS+MSE gain alone is
about +0.018 overall. The public dual-moment kernel failed its tolerance on the 4-cell pooled template (10-22).
Locally: PDS and NMAE are readable; MSE is NOT (clamped, and local raw MSE moved the wrong way on v8b_like).

## Rule
1. Is the decoupled-amplitude lever worth the remaining 9 days? Yes/no, and the main risk.
2. If yes: the first concrete experiment, and the local readout (e.g. NMAE/DE bit-stable, PDS up, MSE projected
   from the official v8a->v8b slope plus the label-free B/claim computation). Name the gate.
3. If no: is there any other lever with an official datum showing more than +0.01? Otherwise say B0 final.
