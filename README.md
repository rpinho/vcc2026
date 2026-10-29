# VCC 2026: public perturbation prediction pipeline

Our research code for the [Virtual Cell Challenge 2026](https://virtualcellchallenge.org/):
source-response transfer, control-conditioned integer-count generation, DE-aware
emission and local evaluation. We are publishing the implementation so others can
inspect it, reproduce components and build on it. This is an evolving research
workbench, with historical experiments retained as source.

**Best official validation: billmA 0.16438**, about rank **463 of 1,580**
when recorded. We are well behind the leading entries and would welcome help.
Start with [lessons and concrete open questions](LESSONS_AND_HELP_WANTED.md)
and the [builder and reproduction guide](REPRODUCING.md). Issues and PRs welcome.

Start with the [short retro](RETRO.md): what we did, what broke, what I would do
differently and where I think this is going. Then the [full retrospective](RETROSPECTIVE.md): the official submission
arc, failed local improvements, harness calibration and biology diagnostics.
The latest [October 26 afternoon/evening update](RESEARCH_UPDATE_20261026b.md)
follows the [morning update](RESEARCH_UPDATE_20261026.md), with results and links
for direction, magnitude, rendering and context-transfer tests. There is no new
official validation score.

`model_v4/` contains AtlasShift-anchor transfer and dual-summary experiments;
`model_v5/` adds spread calibration, cross-source DE gating and control-axis removal;
`model_v6/` generalizes the builder to arbitrary panels/contexts and adds coverage
fill. `proxy/` contains the local scoring/calibration harness; `model_de/` contains
DE experiments and exact count-margin checks. `rsi/` contains acquisition, experiment
launching, baseline reproduction and verification code. `competitors/monitor.py`
monitors the public leaderboard. Machine-specific agent/launchd supervision,
journals, credentials, datasets, weights and predictions are excluded.

The refresh also includes `model_billm/` training and inference,
`model_billm_pack/` and `model_billmB/` checked packaging, `model_v7/` and
`model_v9/` emission experiments, and `final_push/` final builders and
official-scorer ablations. See [REPRODUCING.md](REPRODUCING.md) for the source map,
required local inputs and historical reproduction limits.

The October 26 additions include `calib3_1026/` (local/official component
calibration), `signaudit_1026/`, `dm_1026/`, R9 call-side experiments in
`h2h2_1026/`, MSE instrumentation in `h2h2_1026/scratch/mse_instr/`, and
`bio_gap_1026/` with shared-response, alignment and leave-one-context-out studies.
The [decisions index](DECISIONS_INDEX.md) lists the main rulings and how the AI
advisors voted. See [export scope](EXPORT_20261026b.md)
for the selected result tables and deliberately omitted inputs. Historical
“official” scorer runs in these directories are local experiments with SEARCH
anchors, not submitted validation scores. The scripts are a partial research
release, not a complete clean-checkout reproduction of these campaigns.

## Official validation scores

Recorded on vcc2026-val-1, through the end of validation on 2026-10-22. All metrics below are scaled and
rounded; higher is better. These are submitted scores, not local proxy estimates.

| entry | overall | pds | mse | nmae | fid | reach | jac | what it was |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| v1 transfer | 0.0877 | 0.517 | 0 | -0.006 | -0.067 | 0.099 | -0.017 | public same-target signatures on the control mean (model/) |
| v2b posterior transport | 0.0617 | 0.446 | 0 | 0.042 | -0.188 | 0.092 | -0.020 | 1x signatures, transport of real control cells (model_v2/) |
| v3b mean-shift + 0.3x signatures | **-0.2482** | 0.050 | 0 | 0.014 | **-1.539** | 0.028 | -0.042 | shared average-perturbation shift on every target + shrunk signatures (model_v3/) |
| v4a AtlasShift-anchor transfer | **0.1438** | 0.678 | 0 | 0.129 | -0.019 | 0.076 | -0.001 | AtlasShift anchors + transfer (model_v4/)|
| v4b dual-summary + ash shrinkage + guide | 0.1243 | 0.658 | 0 | 0.070 | -0.025 | 0.042 | -0.000 | dual-summary anchors + ash shrinkage + guide (model_v4/)|
| v5a spread-calibrated v4a | **0.1499** | 0.678 | 0.078 | 0.094 | -0.020 | 0.071 | -0.001 | v4a anchors, per-context spread scaled to 1.0x the MSE noise credit (model_v5/)|
| v5b DE-gated + spread-calibrated + control axes removed | **0.1528** | 0.680 | 0.082 | 0.109 | -0.031 | 0.077 | -0.000 | v5a + DE calls gated by cross-source sign agreement + top-50 control axes removed (model_v5/)|
| v6 | — | — | — | — | — | — | — | arbitrary-panel builder and coverage fill; no official score recorded |
| v7a / v7b | — | — | — | — | — | — | — | native bulk / selective-scatter experiments; validation attempts refused after the phase changed |
| v8a layout-only | 0.1470 | 0.680 | 0.016 | 0.102 | -0.064 | 0.162 | -0.015 | layout changes on the exact v5b totals |
| v8b amplitude | 0.1412 | 0.708 | 0.095 | -0.059 | -0.060 | 0.178 | -0.015 | increased response amplitude |
| v9a | 0.16135 | — | — | — | — | — | — | v5b totals with confidence-ordered DE layout |
| v9b | 0.12049 | — | — | — | — | — | — | alternative v9 layout |
| slim_v5b | 0.08833 | 0.667 | 0 | -0.276 | -0.048 | 0.200 | -0.013 | real-control sampling with gain and excessive spread |
| **billmA** | **0.16438** | **0.706** | **0.061** | **0.116** | **-0.066** | **0.186** | **-0.016** | response MLP blend + v5b amplitude + v9 layout |

Only numbers present in the project's official score record are copied here;
missing components stay blank (—). Historical ranks are approximate snapshots,
not the current position of the team's latest entry on the live board.
The older v1–v3 entries provide history; their complete original directories are
not part of this export. Local numerical parity does not reproduce an official
score, and this export has not been rerun on the full official panel.

## FINAL PHASE

The [final phase](https://virtualcellchallenge.org/) uses vcc2026-test-1 and new
D/E/F contexts. Scores stay hidden until the organizers publish them. Final
submission deadline: **November 5, 2026, 23:59 UTC**.

| Submitted (UTC) | Entry | Recipe / status |
|---|---|---|
| 2026-10-22 05:23 | BKgFJX4gmteiZB88HSLe | billmA v9 layout safety recipe; submitted, score hidden |
| 2026-10-22 08:13 | LIREtFnPJLkMXZprp7Q9 | direction-knob variant; superseded, score hidden |
| **2026-10-23 01:03** | **5NeeoWn96lfm23hgExcm** | **current last entry: billmA v9 layout safety recipe**, same package as the first entry; score hidden |

The last entry was restored after the direction-knob variant lost about 0.008
locally on an official-scorer ablation. That comparison is local evidence, not a
final-test score. This repository contains local builders and packagers; publishing
it does not submit or change an entry.

## Setup and checks

Python 3.12 is the reference environment. Create an environment at the repository
root; Unix shell examples below use it explicitly.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python tools/check_export.py
.venv/bin/python model_v6/pipeline.py --help
```

The requirements preserve the core versions used by the research environment;
optional integrations for the upstream Sidechain code are in requirements-optional.txt.
The normal v4–v6 builder needs no GPU. Large builds need substantial free disk and
RAM; the runtime has disk floors and conservative thread limits.

Run self-contained synthetic checks and regenerate verification receipts locally:

```sh
(cd rsi && ../.venv/bin/python reproduce_atlas.py)
(cd rsi && ../.venv/bin/python verify_bounded.py)
(cd rsi && ../.venv/bin/python verify_bounded_score.py)
.venv/bin/python model_de/verify.py
.venv/bin/python proxy/check_full_adapter.py
```

These tests compare released numerical kernels, exact integer invariants and the
bounded dispatcher against the installed official scorer. Receipts and synthetic
scratch files stay ignored. Run the source/privacy check on a clean checkout;
it intentionally flags generated data and runtime directories.

## Inputs and running a model

Obtain challenge controls, target panel and gene axis yourself through the official
participant download flow. Keep them locally in `data/`: `gene_names.csv`,
`pert_counts.csv`, and one raw-count `context_LABEL.h5ad` per context (at least 1,600 distinct controls
for the default four-cell pooling template). Use only
inputs you have rights to use. No organizer dataset or derived cells are distributed.

Source-response banks must also be prepared locally. The `.npz` schema and moment
construction are in `model_v6/banks.py` (`Acc.bank` and `local_bank`). Source
files must retain raw counts, gene symbols and perturbation/control annotations.
For an initial **H1-only** bank, independently obtain the public 2025 H1 training
file from the [Arc atlas data instructions](https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/virtual-cell-challenge)
and place it at `public_data/H1.h5ad`:

```sh
mkdir -p model_v6/out model_v6/artifacts public_data
(cd model_v6 && ../.venv/bin/python banks.py \
  --targets ../data/pert_counts.csv --sources H1 --out artifacts/banks/example)
```

This is a limited source configuration, not a reproduction of the scored multi-source
runs. Full historical runs additionally use K562/GWPS, X-Atlas-Orion, CD4 and promoter
inputs. Their readers/preparation code is in model_v4/, model_v6/banks.py and rsi/;
some require source metadata/annotations prepared by earlier steps. They are research
scripts, not an automatic download-and-reproduce command. Historical bank rows,
selection metadata and caches are deliberately absent. Rebuilding with different
source coverage or panel centering can change predictions.

Build using your prepared banks, released controls and the official ordered gene
axis (replace the example labels and file names with your own):

```sh
(cd model_v6 && ../.venv/bin/python pipeline.py \
  --arm v5b --targets ../data/pert_counts.csv --genes ../data/gene_names.csv \
  --contexts A=../data/context_A.h5ad,B=../data/context_B.h5ad,C=../data/context_C.h5ad \
  --banks artifacts/banks/example --out example --workers 1)
```

Choose `--arm v4a`, `v4b`, `v5a`, `v5b` or `v6`; `--exclude H1` holds out that entire
source family for an H1 evaluation. Supply `--promoters path/to/promoter_pairs.csv`
when you have generated promoter annotations. Context labels and panels are inputs,
not assumptions in the v6 builder. Predictions appear in `model_v6/out/` with local
count/shape verification receipts in `model_v6/artifacts/`.

If you separately install the official participant CLI and put `vcc` on PATH,
`model_v6/prep.py` validates, packages and independently decodes a prediction locally:

```sh
(cd model_v6 && ../.venv/bin/python prep.py example \
  --perts ../data/pert_counts.csv --genes ../data/gene_names.csv --contexts A,B,C)
```

This step does not submit. Generated `.vcc` packages and every raw or derived
cell-level file must stay out of the repository. The optional PIE inference source,
checkpoints and embeddings are excluded; `--pie-dir` requires separately produced,
appropriately licensed inputs.

For DE experiments use `model_de/harness.py --help` after preparing public tasks and
source banks. The optional `--fast` rank kernel is unavailable because its upstream
license was not established; use the default installed official backend. See
[proxy/README_V2.md](proxy/README_V2.md) for local evaluation prerequisites. Missing
verification receipts deliberately prevent unverified calibration.

## Credits and license

MIT for our original code and included MIT upstream code, with original notices
retained. **AtlasShift** supplies the response-transfer and integer-count foundation;
**Sidechain** supplies shrinkage and public emission ideas/kernels; **Arc cell-eval2**
supplies the official evaluation definitions. We also credit **stage1-lab** for ideas
about DE detectability, call-set size, coexpression and replicate anchors; none of
its unlicensed source code is included. The RSI baseline harness includes the MIT
zero-shot-perturbation release. See [THIRD_PARTY.md](THIRD_PARTY.md) for URLs,
revisions, notices and exclusions. External dependencies and data keep their own terms.

[Publication check](PUBLISH_CHECK.md) records the current rules review and scans.
Open publication does not change team membership rules or permit integrity-undermining
cross-team collaboration. Forks and general improvements are welcome within those rules.
