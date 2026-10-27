Round 5 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = 0.16438 validation; final-phase scores hidden, last entry counts; we never submit without
Ricardo). Your R4 rulings (unanimous): B0 frozen as incumbent; bounded research to 11-02 23:59 UTC; renderer is the
research target; seed replicate first, then a scaling MAP (no promotion). The orchestrator executed both. <= 250
words, one decision per question. You may read h2h2_1026/score_run.py and h2h2_1026/results/render_summary.md.

## Official results (local H1 SEARCH, 79 held-out targets, paired target bootstrap 20k vs B0)
| arm | score_avg | PDS | fidelity | raw MSE | reach | composite vs B0 [95% CI] |
|---|---:|---:|---:|---:|---:|---|
| B0 | +0.021555 | 0.8500 | 0.4487 | 1.4107 | 0.1266 | - |
| B0 seed 2 | +0.022318 | 0.8500 | 0.4479 | 1.4107 | 0.1340 | +0.00076 [-0.00219, +0.00538] |
| B0 seed 3 | +0.022064 | 0.8500 | 0.4485 | 1.4107 | 0.1296 | +0.00051 [-0.00286, +0.00544] |
| 0.75*B0 | +0.018126 | 0.8476 | 0.4473 | 1.1491 | 0.1207 | -0.00343 [-0.00869, +0.00255] |
| 1.25*B0 | +0.021829 | 0.8523 | 0.4513 | 1.7593 | 0.1327 | +0.00027 [-0.00568, +0.00654] |
| 1.5*B0 | +0.021379 | 0.8541 | 0.4547 | 2.1934 | 0.1497 | -0.00018 [-0.00948, +0.00885] |

## Correction to a premise in your R4 answers
The rendered pseudobulk (per-target gene count sums and depth) is BIT-IDENTICAL across render seeds. The seed only
redistributes counts among cells. So official PDS and MSE are deterministic functions of the prediction; render seed
noise is zero for them and ~0.0008 on the composite (all from cell-level DE: per-target reach moves up to 0.83,
fidelity <= 0.017). No earlier verdict is inside render noise; N1's -0.0104 PDS loss is real. The "sampling SNR"
explanation of the +0.054 renderer bonus is therefore wrong as stated: the bonus is a deterministic transform
(anchor_single_mean emission -> count layout with 885 confidence-weighted calls -> log1p(50000*counts/sum)).
Amplitude buys PDS monotonically (+0.004 at 1.5x) but composite is flat and raw MSE rises 55%, so no alpha passes
Opus's gate (composite LB > 0, raw MSE <= B0's). Scaling is closed.

## Rule
1. Given the bonus is deterministic, is Astra's uniform-confidence ablation (same 885 calls, same eligible genes,
   uniform weights; B0 response) still the right next renderer experiment? Or Fable's 2x budget (1770 calls)? Or
   neither? Pick one (one ~11 min official arm), or name a better single experiment.
2. Is it legitimate to change renderer settings (budget, confidence) in a new entry if the frozen submission spec
   allows them? Say what evidence would be needed before any renderer change could replace B0.
3. With seven days of research left and nine levers failed against B0, should the bounded research continue, or
   stop now and leave B0 as final? Say plainly.
