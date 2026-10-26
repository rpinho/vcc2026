# Research update — 2026-10-26

**Best recorded official validation remains billmA 0.16438, historical rank 463
of about 1,580.** The final-test score is hidden. The existing final safety entry
is retained; this publication does not submit or replace an entry.

The work since the previous release concerns local diagnostics using the
installed official scorer. These scores use local SEARCH anchors and are **not
leaderboard scores**. They do not establish a better final entry.

## An expression improvement can lose the composite

One completed H1 experiment fits on 168 targets and evaluates 79 disjoint targets.
B0 is the shipped response recipe; B1 uses STRING features; B2 combines those
features with a projected donor direction. B3 is the inherited pooled donor
branch, B4 the training mean, and B5 a predeclared half blend of B0 and B4.

| Arm | Local six-component composite |
| --- | ---: |
| B0 | 0.021554862 |
| B1 | -0.055693784 |
| B2 | -0.013018053 |
| B3 | 0.025615488 |
| B4 | -0.099178459 |
| B5 | -0.028672362 |

B1 improves raw expression MSE (0.951898566 versus B0's 1.410715392) but loses
perturbation discrimination and the composite. The two fitted candidates fail
the preregistered adoption rule. The saved-result verifier reports PASS for
prediction replay, alignment, aggregation, scaling and paired bootstrap.
It does not independently refit every model or rerun every per-cell rank test.
B3's higher point estimate alone is not a validated replacement decision.

This is one previously used H1 campaign panel, not a new cell-line test.
Intervals condition on fixed fits, anchors and full-panel retrieval contributions;
they do not refit models or rerank duplicated retrieval competitors. Neither the
local score values nor their deltas forecast hidden-test points.

## Renderer diagnostics

Seed replicates keep PDS at 0.8500 while the local composite varies by up to
0.00076 in absolute difference from B0. Amplitude multipliers 0.75, 1.25 and 1.50
produce composite differences -0.00343, +0.00027 and -0.00018; their paired
95% intervals all contain zero. This experiment supplies no accepted composite
improvement. PDS is seed-invariant in these runs, not amplitude-invariant.

## Selected source release and limits

The new files are selected analysis/scoring source from the completed workbench:

- `h2h2_1026/analyze.py`: finite-cohort aggregation, paired target bootstrap and adoption gates.
- `h2h2_1026/paired_boot.py`: descriptive pair comparisons, separate from the registered verdict.
- `h2h2_1026/score_official.py` and `local_runtime.py`: scorer configuration, parity/provenance checks and confined local IO.
- `mse_1025/layout_bound.py`: a fixed-group-total lower bound on raw expression MSE, conditional on validated scorer moments and correction caps.

These are historical research sources, not a complete download-and-run campaign.
No predictions, fitted model, target metrics, cell matrices, data banks, cached
rank tests, journals or raw leaderboard responses are distributed. The scripts
require independently regenerated task files, local baseline/replicate anchors,
per-target metric arrays, prediction summaries, pinned runtime sources and
provenance receipts at the paths named in the code. The layout-bound script
requires the installed scorer source snapshot and its validated hash. Missing
inputs must fail; do not remove integrity checks to make a historical run start.

Use the existing requirements with the reference cell-eval2 0.18.0 environment.
`score_official.py --help` describes its local operations; even `--validate`
requires prepared research inputs. Running the historical fit, full per-cell
scorer or bootstrap is outside this source publication's verification scope.
Syntax, source extraction and the standalone arithmetic/aggregation routines were
checked with synthetic inputs. No production package was opened for testing.

Analysis and code were prepared by AI agents under Ricardo Pinho's direction.
Numerical replay is not independent human expert review, a formal proof, or
biological validation. Existing upstream credits and licenses remain in
[THIRD_PARTY.md](THIRD_PARTY.md).
