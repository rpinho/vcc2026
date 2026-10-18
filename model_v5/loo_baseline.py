"""Sensitivity only: a leave-one-PERTURBATION-out generic baseline for the 11-target local panel.
The official baseline averages 300 constructs (self-share 1/300); on 11 targets each target's own
off-target response is 1/11 of the profile, which makes the local 0-point unrealistically strong.
Same official functions (generic_response_profile, build_baseline_prediction, emit='tile')."""
import runtime, json, sys
from dataclasses import replace
import numpy as np, pandas as pd, anndata as ad
from cell_eval2 import EvalConfig, compute_metrics, aggregate_metrics
from cell_eval2.config import DEParams
from cell_eval2.baseline import generic_response_profile, build_baseline_prediction
import evaluate_v5 as E
view = sys.argv[1]
ctrl, real, anchors = E.load_fold(view)
labels = real.obs.target_gene.astype(str).to_numpy(); targets = anchors['targets']; NT = 'non-targeting'
cfg = replace(EvalConfig.from_preset('vcc2026'), pert_col='target_gene', device='cpu', num_threads=2, gather_threads=2,
              de=DEParams(backend='pdex'), cache_real=None, allow_fractional_counts=True)
parts = []
for t in targets:
    sub = real[labels != t].copy()
    prof = generic_response_profile(sub, pert_col='target_gene', control=NT)
    full = build_baseline_prediction(prof, real, pert_col='target_gene', control=NT, emit='tile')
    parts.append(full[full.obs.target_gene.astype(str) == t].copy())
parts.append(real[labels == NT].copy())
pred = ad.concat(parts); pred.obs_names_make_unique()
raw = aggregate_metrics(compute_metrics(pred, real, config=cfg)).to_pandas().set_index('metric')['mean'].to_dict()
raw['pds_cosine'] = 0.5
runtime.save(E.OUT / 'loo_baseline.json', raw); print(json.dumps(raw, indent=1))
