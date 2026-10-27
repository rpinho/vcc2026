Round 4 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = 0.16438 validation, ranks ~463; final-phase scores hidden, last entry counts; we never submit
without Ricardo). The orchestrator executes your rulings. <= 250 words, one decision per question.

## Official results since round 3 (local H1 SEARCH, 79 held-out targets, same renderer/anchors as B0)
- N1 (your unanimous R2 arm): score 0.012746 vs B0 0.021555; PDS 0.8397 (< 0.845); composite -0.0088
  [-0.0179, +0.0007]. FAIL. CD4 not run (rule).
- P1a (unanimous R1/R2 lever), frozen gene-wise magnitude table: FAIL, both arms.
  P1aB0 0.005789 (PDS 0.7970, fidelity 0.4648 vs 0.4487, raw MSE 0.9628 vs anchor 0.9613 -> scaled MSE still 0);
  P1aN1 -0.002919. Composite vs base -0.0158 [-0.0311, +0.0015] / -0.0157 [-0.0289, -0.0018].
- Manifold/STRING h2h2 locked verdict REJECT (verify PASS). B3 dropped (composite LB -0.0020).

## The finding that overturns a premise you all used in R1-R2 ("pure scaling is PDS-inert")
Official PDS is computed on the RENDERED pseudobulk (counts) vs real cells, all 79 target genes masked. The
same formula on the pre-render response predictions vs truth:
```
B0     median norm/B0 1.000  max|.|/B0 1.000  sign-agree 1.000  pre 0.7960  official 0.8500  bonus +0.0540
B3     median norm/B0 1.001  max|.|/B0 0.912  sign-agree 0.894  pre 0.7965  official 0.8586  bonus +0.0622
B5     median norm/B0 0.543  max|.|/B0 0.505  sign-agree 0.878  pre 0.6204  official 0.6761  bonus +0.0557
N1     median norm/B0 1.000  max|.|/B0 0.997  sign-agree 0.960  pre 0.8229  official 0.8397  bonus +0.0167
P1aB0  median norm/B0 0.279  max|.|/B0 0.543  sign-agree 1.000  pre 0.8103  official 0.7970  bonus -0.0133
P1aN1  median norm/B0 0.274  max|.|/B0 0.530  sign-agree 0.960  pre 0.8322  official 0.7842  bonus -0.0480
```
("bonus" = official minus pre-render PDS.) P1aB0 has B0's signs exactly, so the renderer's sign-dependent confidence
and 885-call layout are IDENTICAL to B0's; only amplitude differs (0.28x). Its pre-render PDS rose (+0.014) and its
rendered PDS fell 0.053. So amplitude causally drives rendered PDS through the renderer (sampling SNR). The response
proxy's scale-invariance is false officially. B5 at ~0.54x keeps the bonus; N1 at full norm loses most of it (its
direction change / 4% sign flips, not amplitude). Shrinking raised fidelity (+0.016) - the opposite of mse_1025's CD4
scalar. Every arm so far loses to B0; B0 sits at a local optimum where its large amplitude buys PDS through the
renderer while MSE is clamped at 0 locally (but worth 0.061 scaled on the validation leaderboard).

## Rule
1. Is an official pure-scaling ladder on B0 (alpha*B0, signs and layout identical, so it isolates amplitude ->
   rendered PDS / fidelity / MSE) now worth running? If yes: which alphas (max 3, ~11 min each), and is this an
   exploration map (no promotion from it) or a candidate with a gate? If no: what instead?
2. Should the renderer itself (anchor emission, the 885-call confidence layout, budget) become the next research
   target, given it contributes ~+0.054 PDS to B0 and nobody has mapped it? Name one concrete first experiment.
3. Given 10 days left and every lever failing against B0, is freezing B0 now the right call? Say plainly.
