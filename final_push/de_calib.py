"""Calibrate the H1_rich_search harness against official DE history before trusting it for DE-side levers.

Arms (v5b biology, no learned blend, all on the same frozen task; official cell-eval2 scoring):
  v5b        : the frozen v5b LOCO recipe (4-cell pooled template, one-moment pinned totals)      [official .1528]
  v9a_885    : v5b emission, then the v9 confidence-ordered call layout at the final-package budget 885   [~v9a .1613 / v8a .1470]
  v9a_694    : same at v8A's context-A budget 694
  v8b_like   : profile amplitude x1.5 about the control (v8b's move), v5b emission, layout at 885       [official v8b .1412]
Official validation contrasts to reproduce: v8a/v9a vs v5b = reach UP, fidelity DOWN, jaccard down, nmae ~flat;
v8b vs v8a = NMAE COLLAPSE (-0.16), pds up, fidelity ~flat. Writes only under final_push/artifacts/de_check/.
"""
import os, sys, types, gc, time, json, argparse
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; VCC = HERE.parent
sys.path.insert(0, str(HERE))
import de_check                                   # reuses the RSI harness context (env, thread caps, OUT dir)
OUT = de_check.OUT


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--task', default='H1_rich_search')
    ap.add_argument('--arms', default='v5b,v9a_885,v9a_694,v8b_like'); a = ap.parse_args(); task = a.task
    benchmark, scoring, bounded_score, baseline, rsi = de_check.rsi_context()
    helpers = baseline.helpers
    from scipy import sparse
    contract = baseline.read_spec(); model, ctrl, real = benchmark.setup(task); assert model['split'] == 'search'
    genes = ctrl.var_names.to_numpy(str); targets = np.array(model['targets']); T = len(targets)
    q4, depth, bf, amp = baseline.calibrated_profile(ctrl, task, contract)
    _, _, L, Lc, support, pinfo = baseline.profile(ctrl, task, contract)
    # v9 confidence key exactly as the final packages compute it (net cross-source sign agreement, z tie-break)
    with baseline.local_v4(model['context']) as v4:
        names = v4.load_banks(None, True)
        sources = [v4.aligned_bank(n, targets, genes, 'log2fc') for n in names]
        cd = v4.cd4_family(targets, genes, 'log2fc')
        if cd is not None:
            sources.append(cd)
    sign = np.sign(L - Lc); net = np.zeros(L.shape); zs = np.zeros(L.shape)
    for e, m, v in sources:
        z = np.where(m & np.isfinite(v), e / np.sqrt(np.maximum(v, 1e-12)), 0); sig = np.abs(z) > 1.96
        net += (sig * (np.sign(z) == sign)).astype(float) - (sig * (np.sign(z) == -sign)).astype(float)
        zs += np.minimum(np.abs(z), 100) * (np.sign(z) == sign)
    confidence = net * 1e6 + zs; del sources, cd; gc.collect()
    # v9 layout modules with their OWN runtime stub (reference files under final_push, never under rsi)
    rsi_rt = sys.modules['runtime']
    stub = types.ModuleType('runtime'); stub.ROOT = OUT; stub.json = rsi_rt.json; stub.sha = rsi_rt.sha; stub.save = rsi_rt.save
    stub.guard = lambda *a, **k: None; sys.modules['runtime'] = stub
    sys.path.insert(0, str(VCC / 'model_v8')); sys.path.insert(0, str(VCC / 'model_v9'))
    import de as v9de, layout as v9layout
    sys.modules['runtime'] = rsi_rt
    refdir = OUT / 'artifacts' / f'ref_{task}'
    if not (refdir / 'cdf_scale.npy').exists():
        t0 = time.time(); v9de.prepare(ctrl.X, refdir, np.random.default_rng(20261018)); print('reference prepared', round(time.time() - t0), 's', flush=True)
    sys.modules['runtime'] = stub; ref = v9de.Reference(task); sys.modules['runtime'] = rsi_rt
    gi = {g: i for i, g in enumerate(genes)}
    w = helpers._row_weights(T, len(genes), helpers.excl_cols(targets, genes))
    atlas = baseline.recipes.v4.atlas

    def jk(x):
        return float(helpers.jackknife_correction(sparse.csr_matrix(x), np.zeros(x.shape[0], np.intp), 1, helpers.TS, chunk=32)[0])

    def own_of(t):
        return int(gi[t]) if t in gi else np.array([], dtype=np.int64)

    def layered(i, prof, budget):
        original = baseline.emit(q4, depth, prof[i], model['context'], targets[i]).toarray().astype(np.int32)
        x, called, info = v9layout.emit(original, ref, confidence[i], budget, own_of(targets[i]))
        return sparse.csr_matrix(x), len(called)

    anchors = rsi_rt.json.loads((rsi / f'artifacts/scores/{task}/anchors.json').read_text())
    assert anchors['real_sha256'] == scoring.matrix_id(real)
    bf15 = helpers.from_logbulk(helpers.rescale(helpers.logbulk(bf), 1.5, Lc))
    results = {}
    for arm in a.arms.split(','):
        t0 = time.time(); out = OUT / task / arm; out.mkdir(parents=True, exist_ok=True)
        prof = bf15 if arm == 'v8b_like' else bf
        budget = {'v9a_885': 885, 'v9a_694': 694, 'v8b_like': 885}.get(arm)
        P = np.zeros((T, len(genes))); jks = np.zeros(T); calls = np.zeros(T, int)

        def iterator():
            for i, t in enumerate(targets):
                if budget is None:
                    x = baseline.emit(q4, depth, prof[i], model['context'], t)
                else:
                    x, calls[i] = layered(i, prof, budget)
                jks[i] = jk(x); P[i] = np.asarray(x.sum(0)).ravel()
                if i % 10 == 0:
                    print(arm, 'EMIT', i, '/', T, 'elapsed', round(time.time() - t0), flush=True)
                yield t, x
        cfg = scoring.config(cache_real=str(OUT / task / 'real_cache'))
        tidy, digest = bounded_score.score_iter(real, iterator(), cfg)
        raw = benchmark.aggregate(tidy); scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
        Pn = P / P.sum(1, keepdims=True); B = float(helpers.budget(Pn, targets, genes)); claim = float((w * jks).sum())
        rec = dict(task=task, arm=arm, raw=raw, scaled=scaled, overall=overall, prediction_sha256=digest, B_emitted=B, claim_emitted=claim,
                   ratio_B_over_claim=B / claim, budget=budget, calls_mean=float(calls.mean()) if budget else None,
                   alpha_v5b=amp['alpha'], amplitude_factor=1.5 if arm == 'v8b_like' else 1.0, seconds=round(time.time() - t0),
                   status='official scorer on a public LOCO SEARCH task; calibration of the harness against official history')
        tidy.write_csv(out / 'tidy.csv'); rsi_rt.save(out / 'score.json', rec); results[arm] = rec
        print('SCORED', task, arm, 'overall', round(overall, 5), {k: round(v, 4) for k, v in scaled.items()},
              'rawMSE', round(raw['expr_mse_unbiased_capped_norm'], 4), 'B/claim', round(B / claim, 3), flush=True)
        gc.collect()
    rsi_rt.save(OUT / task / 'calibration_summary.json', results)
    print('DONE', task, flush=True)


if __name__ == '__main__':
    main()
