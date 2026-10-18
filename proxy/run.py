"""Score recipes on a LOCO-300 proxy task with rsi's bitwise-verified official scorer (cell-eval2 0.18.0), read-only.

Scorer: proxy/scorer/{bounded_score,bounded_de,official_scoring,baseline_protocol}.py are byte copies of rsi's; this
script re-checks their SHA-256 against rsi's PASS fixtures before activating them (bounded memory scheduling only;
every metric/DE/aggregation/scaling function is official). Scaling: Arc-documented mean-response baseline over the
destination's own panel (official generic_response_profile -> 1/n self-share, n = panel size) tiled to 400 cells, and
official replicate anchors (compute_replicate_anchor, 3 splits; rsi uses 5 for H1 / 3 for CD4).
Predictions: model_v6 pipeline.Context (v4a/v4b/v5a/v5b numerics; destination family's banks excluded), 400 cells
per target, streamed to the scorer (no prediction file on disk).

python run.py --dest K562 --recipes v4a,v4b,v5a,v5b
"""
import sys, os, json, time, argparse, hashlib, gc
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['NUMBA_NUM_THREADS'] = os.environ.get('PROXY_DE_THREADS', '4')  # before numba loads (runtime would cap it at 2)
sys.path.insert(0, os.path.join(HERE, '..', 'model_v6')); sys.path.insert(1, os.path.join(HERE, 'scorer'))
import runtime
import model, pipeline, v5
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from dataclasses import replace
P = Path(HERE); RSI = runtime.VCC / 'rsi'
runtime.ROOT = P  # bounded_de's DE cache lives in proxy/artifacts/de_cache (model_v6 paths were bound at import)
import official_scoring as scoring, bounded_score, bounded_de
import cell_eval2.anchor as anchor
from cell_eval2.anchor import compute_replicate_anchor
from cell_eval2 import aggregate_metrics

FAMILY = {'K562': 'K562', 'HCT116': 'HCT116', 'HEK293T': 'HEK293T', 'CD4': 'CD4', 'H1': 'H1'}
NT = 'non-targeting'
DE_THREADS = int(os.environ.get('PROXY_DE_THREADS', '4'))  # operational only: bounded_de keys its cache with num_threads=1


def activate_scorer():
    b = json.loads((RSI / 'artifacts/bounded_fixture.json').read_text()); s = json.loads((RSI / 'artifacts/bounded_score_fixture.json').read_text())
    assert b['status'] == 'PASS bitwise equality' and runtime.sha(P / 'scorer/bounded_de.py') == b['adapter_sha256']
    assert s['status'] == 'PASS bitwise equality' and runtime.sha(P / 'scorer/bounded_score.py') == s['driver_sha256']
    assert scoring.provenance()['official_files'] == b['cell_eval2']['official_files'], 'installed cell-eval2 differs from fixture'
    for f in ['official_scoring.py', 'baseline_protocol.py']:
        assert runtime.sha(P / 'scorer' / f) == runtime.sha(RSI / f)
    bounded_de.install_unchecked(); anchor.compute_metrics = bounded_score.resident
    return dict(bounded_de=b['adapter_sha256'], bounded_score=s['driver_sha256'], cell_eval2=scoring.provenance()['version'])


def aggregate(tidy):
    f = aggregate_metrics(tidy); return dict(zip(f['metric'].to_list(), f['mean'].to_list()))


RSI_H1 = RSI / 'artifacts/tasks/H1_rich_search.h5ad'


def load_task(dest):
    if dest == 'H1':
        # rsi's H1-rich SEARCH task, read-only: 79 targets x 400 real cells, exact batch quotas. Both rsi control arms
        # (3,200 input + 3,200 reference) become the shared control set, as on the server.
        a = ad.read_h5ad(RSI_H1)
        role = a.obs.evaluation_role.astype(str).to_numpy()
        a.obs['target_gene'] = np.where(np.isin(role, ['input_control', 'reference_control']), NT, a.obs.target_gene.astype(str).to_numpy())
        panel = np.array(sorted(set(a.obs.target_gene) - {NT}))
        if not (P / 'tasks/H1_panel.csv').exists():
            pd.DataFrame({'target_gene': panel}).to_csv(P / 'tasks/H1_panel.csv', index=False)
            runtime.save(P / 'tasks/H1.json', dict(destination='H1', source=str(RSI_H1), sha256=runtime.sha(RSI_H1), n=len(panel), cells_per_target=400,
                         controls=int((a.obs.target_gene == NT).sum()), genes=int(a.n_vars), note='rsi H1_rich_search read-only; HOLDOUT not used'))
        return a, a[a.obs.target_gene.to_numpy() == NT].copy(), panel
    a = ad.read_h5ad(P / 'tasks' / f'{dest}.h5ad')
    a.obs['target_gene'] = a.obs.target_gene.astype(str)
    lab = a.obs.target_gene.to_numpy()
    ctrl = a[lab == NT].copy()
    panel = pd.read_csv(P / 'tasks' / f'{dest}_panel.csv').target_gene.astype(str).to_numpy()
    panel = np.array([t for t in panel if (lab == t).any()])
    return a, ctrl, panel


def anchors(dest, real, folder, n_splits=3):
    f = folder / 'anchors.json'; ident = scoring.matrix_id(real)
    if f.exists():
        r = json.loads(f.read_text()); assert r['real_sha256'] == ident; return r
    from baseline_protocol import mean_response
    t0 = time.time(); cfg = scoring.config(cache_real=str(folder / 'real_cache'))
    mean, desc = mean_response(real)
    tidy, _ = bounded_score.score_iter(real, bounded_score.baseline_iter(real, mean, cells=400), replace(cfg, allow_fractional_counts=True))
    base = aggregate(tidy); tb = time.time() - t0
    print(dest, 'baseline', json.dumps(base), round(tb), 's', flush=True)
    splits, frame = compute_replicate_anchor(real, config=replace(cfg, cache_real=None), n_splits=n_splits)
    rep = dict(zip(frame['metric'].to_list(), frame['replicate'].to_list()))
    scoring.scale(base, base, rep)
    r = dict(baseline=base, replicate=rep, real_sha256=ident, n_splits=n_splits, baseline_construction=desc, baseline_seconds=tb, anchor_seconds=time.time() - t0 - tb)
    runtime.save(f, r); print(dest, 'anchors', json.dumps(rep), flush=True); return r


def score(dest, recipe, real, ctrl, panel, anc, folder, fill=False, workers=3):
    out = folder / f'{recipe}{"_fill" if fill else ""}.json'
    if out.exists():
        return json.loads(out.read_text())
    t0 = time.time()
    model.configure(model.SHARED_BANKS, [FAMILY[dest]], runtime.VCC / 'model_v6/artifacts/promoter_pairs_v4.csv')
    genes = real.var_names.to_numpy(str)
    if recipe in ('v1', 'v2b', 'v3b'):
        # Frozen legacy v2 imports its own pipeline alias. Restore the modern module
        # after loading so a later Context emitter remains spawn/pickle-compatible.
        saved_pipeline = sys.modules.get('pipeline')
        sys.path.insert(0, str(P / 'legacy')); import legacy_recipes
        if saved_pipeline is not None:
            sys.modules['pipeline'] = saved_pipeline
        ctx = None; gen = legacy_recipes.predict(ctrl, list(panel), dest, recipe)
    else:
        arm, pie = recipe, None
        if recipe == 'v5b_pie':  # = v5c recipe: v5b gate + 300 gene-centred PIE p_de calls per target (4-fold ensemble)
            arm = 'v5b'; pie = pipeline.load_pie(runtime.VCC / f'model_v6/pie/out/ens4_{dest}.npz')
        ctx = pipeline.Context(arm, panel, genes, ctrl, dest, 0, fill=fill, workers=workers, pie=pie); gen = ctx.stream()
    seen = []

    def it():
        for t, x in gen:
            runtime.guard()  # yields (pauses the whole streamed scoring) while a kids' job runs or disk is low
            seen.append(t)
            if len(seen) % 50 == 1:
                print(dest, recipe, 'emitted', len(seen), '/', len(panel), round(time.time() - t0), 's', flush=True)
            yield t, x
    # restrict the reference to the panel + controls (targets with too few cells were dropped from the panel)
    tidy, pred_sha = bounded_score.score_iter(real, it(), replace(scoring.config(cache_real=str(folder / 'real_cache')), num_threads=DE_THREADS))
    tidy.write_parquet(out.with_suffix('.tidy.parquet'))
    raw = aggregate(tidy); scaled, overall = scoring.scale(raw, anc['baseline'], anc['replicate'])
    rec = dict(dest=dest, recipe=recipe, fill=fill, overall=overall, scaled=scaled, raw=raw, n_targets=len(panel), seconds=round(time.time() - t0, 1),
               prediction_sha256=pred_sha, build_info={k: v for k, v in ctx.info.items() if k not in ('targets_without_source_row',)} if ctx else 'legacy rsi recipe port',
               targets_without_source_row=ctx.info.get('targets_without_source_row') if ctx else None, status='Public-context proxy estimate; not an official score')
    runtime.save(out, rec); print('SCORE', dest, recipe, round(overall, 4), json.dumps({k: round(v, 3) for k, v in scaled.items()}), rec['seconds'], 's', flush=True)
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--dest', required=True); ap.add_argument('--recipes', default='v4a,v4b,v5a,v5b')
    ap.add_argument('--fill', action='store_true'); ap.add_argument('--splits', type=int, default=3)
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args()
    prov = activate_scorer()
    folder = P / 'artifacts/scores' / a.dest; folder.mkdir(parents=True, exist_ok=True)
    real, ctrl, panel = load_task(a.dest)
    keep = np.isin(real.obs.target_gene.to_numpy(), list(panel) + [NT]); real = real[keep].copy()
    runtime.guard(); anc = anchors(a.dest, real, folder, a.splits)
    for r in a.recipes.split(','):
        runtime.guard(2 * 1024**3)
        score(a.dest, r, real, ctrl, panel, anc, folder, fill=a.fill, workers=a.workers); gc.collect()
