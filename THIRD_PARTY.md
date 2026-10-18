# Third-party code and ideas

The top-level MIT license applies to our original code. The original notices below
also apply to copied/adapted upstream code; they are retained verbatim.

| Work | Source | Included code and notice | Reviewed revision |
|---|---|---|---|
| AtlasShift | https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 | model_v4..v6/vendor/atlas_model.py and atlas_predict.py; rsi/vendor/atlas_model.py and atlas_release/; each ATLAS_LICENSE or LICENSE (MIT) | d24ce4fdae0cd7cba3ba29546cd8737094eae98a |
| Sidechain | https://github.com/saberhq/sidechain | adaptive_shrink.py kernels and rsi/vendor/sidechain_release/src/; SIDECHAIN_LICENSE and sidechain_release/LICENSE (MIT, Saber HQ) | pinned RSI release 2d8f3d17042030ac866c4d568584dc8e5b2bbcc5 |
| zero-shot-perturbation | https://github.com/MMatinGerami/zero-shot-perturbation | rsi/vendor/zsp_release/src/ and LICENSE (MIT, Matin Gerami) | 2a198d5fa165179633ed0c0becc0d2b050b27187 |
| Arc cell-eval2 | https://github.com/ArcInstitute/cell-eval2 | Installed dependency, not vendored; third_party/CELL_EVAL2_LICENSE (MIT, Arc Research Institute) | package 0.18.0; inspected repository a2630dcbda95fc0b819bf841e56d9561d9c1e1da |

Stage1-lab (https://github.com/danielchen26/vcc2026-stage1-lab, reviewed
3a5ea5e9f640a6dd652516fd1492711618829870) informed analysis of DE detectability,
call-set size, control coexpression and replicate anchors. No license was found;
we include no stage1-lab source code or copied implementation. Credit acknowledges
ideas, not permission to redistribute that repository.

Two local rank-walk optimizations derived from numba-mwu source were excluded:
proxy/presorted_mwu.py and model_de/mwu_kernel.py. Their upstream license was
unspecified in the inspected source/wheel. The export uses the installed official
scorer, and model_de --fast is disabled. pdex and numba-mwu remain external
requirements under their own terms; our MIT license does not relicense them.

Arc PIE source (CC BY-NC-SA 4.0), checkpoints and embeddings are excluded. The
builder's optional --pie-dir accepts separately generated summaries, whose source
terms must be respected. No State/Stack code, pretrained weights or external
training data are included or relicensed.

Upstream source distributions were narrowed to Python source needed by the
baseline harness, and example paths were made relative. Private upstream journals,
site/leaderboard dumps, source data and weights are omitted. Frozen working-tree
hash receipts are not transferred; regenerate verification for this export.
