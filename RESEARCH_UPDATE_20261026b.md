# Research update — October 26, afternoon and evening

**billmA 0.16438 remains our best recorded official validation result.** Validation
closed on October 22; final scores are hidden. Nothing below is a new leaderboard
result. The existing final safety recipe remains the reference. The longer
[retrospective](RETROSPECTIVE.md) explains the complete submission history.

This follows the [morning update](RESEARCH_UPDATE_20261026.md). Names such as
`official_N1` in historical scripts mean the official scorer **run locally**, on
public H1 tasks with local SEARCH anchors. B0 is the locally replayed incumbent
recipe, not a reproduction of its hidden-context score.

## Direction and magnitude failed the rendered six-component test

The experiment used 168 H1 fitting targets and 79 disjoint evaluation targets.
After repeated inspection, those 79 are a development panel, not untouched
confirmation. N1 extrapolates away from the training mean at 1.25× and then
restores each B0 response norm. P1a fits gene-wise magnitude calibration.

| Local arm | Composite | Raw PDS | Raw expression MSE | Outcome |
|---|---:|---:|---:|---|
| B0 | 0.021555 | 0.8500 | 1.4107 | Reference |
| N1 | 0.012746 | 0.8397 | 1.4300 | Fails composite and PDS gates |
| P1aB0 | 0.005789 | 0.7970 | 0.9628 | Better distance, worse composite |
| P1aN1 | -0.002919 | 0.7842 | 0.9681 | Fails versus its N1 base |

N1's paired composite change was -0.00881, 95% CI [-0.01790, +0.00067]. P1aB0
reduced raw MSE by about 32%, but its scaled MSE stayed zero: it missed the local
baseline anchor, 0.9613. Its composite change was -0.0158 [-0.0311, +0.0015].
B3's +0.00406 composite change also failed confirmation: CI [-0.00197, +0.01048].
See the [decision trail, rounds 1–4](decisions/consult_1026_amp/DECISION.md) and
[campaign report](h2h2_1026/REPORT.md).

## Renderer effects were real, but did not establish a replacement

Two render seeds moved the composite by +0.00076 and +0.00051. Their PDS was
identical; the DE metrics moved. Scaling B0 by 0.75/1.25/1.50 changed the
composite by -0.00343/+0.00027/-0.00018; every paired interval contained zero.
Raw MSE increased to 1.7593 and 2.1934 for the two larger amplitudes. Response-space
cosine invariance did not imply rendered PDS invariance.

Uniform confidence at 885 calls, U885, lost -0.00543 [-0.00982, -0.00056].
Doubling the call budget, D1770, gained **+0.03449 locally** [+0.02194, +0.04770].
It passed a diagnostic gate, not a promotion gate. Its top five targets supplied
38% of the summed gain. [Renderer results](h2h2_1026/results/render_summary.md),
[uniform-confidence gate](h2h2_1026/results/u_gate.md) and
[double-budget gate](h2h2_1026/results/d_gate.md) retain the numbers.

Why reject such a large local point gain? The
[component calibration](calib3_1026/REPORT.md) found that local layout changes
inflate fidelity by **9–12×** and Jaccard by **7–9×**, while understating reach
changes at about 0.6–0.75× the official magnitude. Component signs can agree
while the composite sign reverses. A stratum model put D1770 near +0.0026 in a
hypothetical official projection, comparable to transformed seed changes of
about +0.002–0.003. This is an assumption-sensitive diagnostic, not a forecast.
Calibration outside the layout axis got only 1 of 4 held-out contrast signs
right, or 2 of 4 without MSE. D1770 is **indeterminate**, not a known official loss.

## Three follow-up call-side screens failed

| R9 arm | Local composite change | Central mapped diagnostic | Gate |
|---|---:|---:|---|
| RH20: reorder the top reach head | -0.00002 | approximately 0 | FAIL: inert |
| K250: trim the budget to 250 | -0.07007 | -0.03794 | VOID; large observed loss |
| DET885: detectability-based selection | -0.00280 | -0.00601 | VOID; observed loss |

The latter two also violated the DE-only validity check; their negative readings
do not make them valid isolated experiments. The mapped seed floor was +0.00280.
No arm reached the bar for a CD4 follow-up. See the
[registered R9 design](h2h2_1026/AMENDMENT_R9.md),
[build](h2h2_1026/build_r9.py), [scoring](h2h2_1026/score_r9.py) and
[gate results](h2h2_1026/results/r9_gate.md).

The [sign audit](signaudit_1026/REPORT.md) found 80,044 called genes, with sign
precision 0.567 overall. Renderer sign flips affected 11.6% of calls, but those
flipped calls agreed with truth **more** often (0.602 versus 0.563). The intended
head itself was weak. Reference-null genes were negative 50.006% of the time;
the proposed null down-sign skew did not appear locally. This corrects hypotheses
in the dated [public-method idea review](competitors/IDEAS-2026-10-26.md).

## No usable MSE instrument or decoupled-amplitude mechanism emerged

The [MSE instrumentation study](h2h2_1026/scratch/mse_instr/REPORT.md) tried the
current clamp, lifted clamp, uncapped unbiased distance, lifted trace cap,
label-based and label-free budgets, and plain bulk distance. Every formulation
failed the C4 sign/ratio requirement. Lifting the clamp gave -0.001212 locally
where the official amplitude contrast gained +0.0791 scaled MSE. The genuine
v8/v9 layout replay is still missing; C6 is only a surrogate.

Historical numerical checks reproduced installed-kernel rows within 2e-12 and
79 reconstructed native pseudobulks within 8.89e-16. They validate the calculation,
not its usefulness as an official instrument. [Details and formulas](h2h2_1026/scratch/mse_instr/DETAILS.md)
and the [contrast table](h2h2_1026/scratch/mse_instr/contrast_variants.csv) are included.

The [decoupled-amplitude probe](dm_1026/REPORT.md) tried whole-cell integer depth
multipliers, preserving CPM vectors. Across five targets, simple weights reached
only 1.03× amplitude with 0.20 orthogonal error; a continuous NNLS upper-bound
probe reached roughly 1.00–1.01× toward 1.2/1.38× goals. That mechanism failed.
Gene-wise depth-tilted placement was not built: it could disturb Wilcoxon ranks,
and no qualified instrument could price its official MSE benefit.

## The biology bottleneck is mostly transfer on already-covered targets

These follow-ups used training168 only, with target or recipient holdouts; they
did not use the 79 outcomes for fitting. Their metrics are vector diagnostics,
not the six-component composite.

- [Decomposition](bio_gap_1026/REPORT.md): B0 mean cosine **0.084**; normalized
  SSE **1.642**; error split **58.3% direction / 41.7% magnitude**. Covered targets
  carry **99.5%** of error. Perfectly fixing only missing rows buys at most 0.5%
  on this panel; better donors for existing rows remain a different hypothesis.
- [Shared-response test](bio_gap_1026/shared_resp/REPORT.md): a training mean plus
  shrunken B0 lowers normalized SSE to **0.943**, but PDS collapses **0.856→0.523**
  and centered spread falls to **1.46%** of B0. Family correction's roughly 10.5%
  SSE gain likewise reduces spread; unrescaled aligned cross falls
  **0.1154→0.1086**. Large distance gains were dominated by shrinkage.
- [In-context alignment](bio_gap_1026/align/REPORT.md): after matching centered
  spread, kNN using other targets' **true H1 responses** raises normalized cross
  **0.1154→0.4230**; ridge/RRR reaches **0.3984**. PDS falls to **0.5782/0.6619**.
  This demonstrates signal with recipient truth, which D/E/F do not provide.
- [Leave-one-context-out maps](bio_gap_1026/loco_map/REPORT.md): full maps fail
  PDS in all five recipients. Bounded blends pass a weaker cross-CI/mean-PDS rule
  in five of five. Against full-axis H1 B0, the best control-conditioned blend
  raises cross only **0.115447→0.118116** (about **+0.0027**); PDS changes
  **-0.002281**, CI **[-0.005204, +0.000428]**, missing stricter noninferiority.

Cross is normalized by truth energy; matching centered spread does not fix
uncentered energy. Converting +0.0027 cross into official MSE or quoting a
net +0.0004 overall gain is unsupported. Large gains were demonstrated with
in-context truth; these experiments do not prove such truth is necessary for
every future method. [Advisor round 12](decisions/consult_1026_amp/DECISION.md)
closed this production line while preserving the small transferable diagnostic.

## What this export verifies

This is a curated source-and-summary release. CSVs contain diagnostic/selection
summaries, not cells or predicted expression matrices. Fitted coefficient arrays,
datasets, model weights, bulk arrays, caches, logs and operational snapshots are
omitted. Many historical inputs are absent; scripts are research evidence, not
a complete reproduction package. See [export scope](EXPORT_20261026b.md).

Historical PASS receipts are retained as historical results, with paths normalized;
they do not certify the edited export. This publication checks syntax, structure,
secrets, privacy, links and checksums locally. It does not rerun the campaigns or
reproduce an official score. Source was prepared by AI agents under Ricardo
Pinho's direction; agreement among advisors is not independent biological review.
