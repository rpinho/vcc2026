# Builders, inputs and reproduction limits

This is a source release. It contains no source banks, target metadata, fitted
weights, predictions, packages, historical caches or private experiment journals.
Small numerical checks can be run without challenge inputs; full score
reproduction requires obtaining inputs and rebuilding the matching historical
receipts. A valid package is not evidence of a better biological prediction.

## Code map

| Directory | What it implements |
|---|---|
| model_v4–model_v6 | source fusion, calibrated response profiles, arbitrary-panel controls and integer-count emission |
| model_v7 | separate bulk/mean heads, feasible moment projection and scatter experiments; no official validation score |
| model_v8 | layout balancing and public real-cell call-budget estimation used by the later builders |
| model_v9 | confidence-ordered DE layout, fixed-profile MSE checks, independent count audits and local packaging |
| model_billm | billmA response preparation, source-family exclusion, three-seed training/inference and local evaluation |
| model_billm_pack | original billmA renderer and fail-closed historical packaging path, plus exact receipt recheck |
| model_billmB | expanded-source experiment and panel-general builder/packager; the safety run selects the original billmA weights |
| final_push | selector, official-scorer ablations, final safety script and rejected direction-knob/emitter experiments |
| validation | independent moment, budget and rank-test checks used by the count auditor |

The reconstruction has five full-gene input channels (pooled response,
compatibility-weighted response, disagreement, destination basal expression and
coverage), a two-hidden-layer residual MLP, MSE plus contrastive loss, and seeds
1019/2027/4093. Its learned unit direction is mixed at 0.3 into the pooled
direction before emission. CD4 sister screens form one excluded family. This is
our implementation of Bill M's public method description, not his source code.

## Local setup

Follow [README setup](README.md#setup-and-checks). The neural model also requires
PyTorch (`pip install torch` in your environment; select a build for your CPU/GPU).
Historical receipt verification additionally uses Node.js 22 or newer. The local
packagers discover the separately installed participant `vcc` CLI on PATH; the
final phase used **vcc 0.3.0**. The scorer remains pinned to the historical
**cell-eval2 0.18.0**. Treat a scorer upgrade as a separate parity experiment.

## Rebuilding billmA inputs and weights

Obtain the canonical ordered gene list through the organizers and source screens
under their own terms. `model_billm/prepare.py` shows the native normalization,
duplicate-gene aggregation, measured-gene masks and response-bank format.
Its raw readers expect authorized inputs in `public_data/`; the cached-bank
readers require locally rebuilt v4/v6 moment banks. Inspect those prerequisites
before running the script. No automatic fetch of all historical inputs is promised.

```sh
mkdir -p model_billm/data model_billm/artifacts
.venv/bin/python -B model_billm/prepare.py
.venv/bin/python -B model_billm/verify.py
.venv/bin/python -B model_billm/train.py --help
```

The training CLI supports production and held-family folds. Rebuild its feature
caches and production/H1/CD4 ensembles locally; no returned GPU archive is
distributed. Keep recipient families out of both fitting and donor-feature
construction. Do not use challenge treated observations for training.

The original billmA count construction is in
`model_billm_pack/renderer_preflight.py`: it requires the original frozen v5b
profile/amplitude receipts, source/promoter inputs and v8 per-target call budgets.
Run its probes, then `--full`, then `package_checked.py`; that path checks every
group and independently decodes the archive. Those historical receipts are not
published and must be regenerated. Publication edits change source hashes: do
not copy a private receipt or weaken a hash check to make the export pass.

The expanded billmB experiment failed its frozen local acceptance on CD4 and
was not packaged. Its `acceptance.py` and `verify_local_results.mjs` preserve
that decision boundary; the safety recipe below uses original billmA weights.

## Final safety build

`model_billmB/panel_pipeline.py` generalizes the same blend and v9 layout to an
explicit target panel and context paths. It recalibrates v5b amplitude on the
destination controls and uses a context-independent public call prior. It does
not use the expanded experiment's weights when `--weights` selects model_billm.

Before running, obtain final controls under `data/final/bundle/`, locally prepare
panel moment banks in `model_billmB/final/banks/`, and rebuild response banks from
your licensed source backbone. `model_billmB/panel_responses.py`, `donor_coverage.py`
and the bank readers implement those steps. The public scalar-prior generator
`model_v8/public_scalars.py` requires public Jurkat/HepG2/RPE1 raw cells and the
locally prepared v8 control references in `model_v8/artifacts/ref_{A,B,C}/`.
Its output `model_v8/artifacts/public_scalars.json` is an input to the final
builder; no per-cell/target records are copied into this release.

With those inputs and locally trained original billmA weights present:

```sh
# Builds and validates locally; never submits.
./final_push/run_safety.sh example_safety
```

The script runs coverage accounting, six renderer probes, all-group count checks,
official prep and independent decoded equality. It refuses existing package
outputs. Use a new name after changing inputs. Packages stay ignored.
`final_push/main/panel_pipeline.py` contains the direction-knob variant; it lost
the later local ablation and is included as a failed experiment.

## Official-scorer ablations

```sh
.venv/bin/python -B final_push/de_check.py --help
.venv/bin/python -B final_push/de_calib.py --help
.venv/bin/python -B final_push/de_levers.py --help
```

These reuse `rsi/official_scoring.py`, `bounded_score.py` and
`v5b_loco_baseline.py`. They require locally rebuilt frozen public SEARCH task
contracts, anchors, source banks and fold weights. They refuse sealed HOLDOUT
reads. `de_check.py` compares incumbent and selective scatter; `de_calib.py`
checks historical layout/amplitude contrasts; `de_levers.py` tests DE budgets,
mean-head amplitude, PIE ordering and direction knobs. Surrogates are labeled in
the code. Optional PIE inputs must be produced separately under their license.

The harness failed the later ranking audit; its outputs are diagnostics, not a
permission to claim a final-score gain. See [lessons](LESSONS_AND_HELP_WANTED.md).
