"""DE-side levers on the main_a2 biology, one factor at a time, scored by the official cell-eval2 code paths on
H1_rich_search (harness calibrated in de_calib.py: trustworthy for call-set levers, blind to NMAE overshoot).

Base = main_a2 biology rebuilt here: v5b profile -> direction knobs (compat mask, confidence allocation, floored
silence) -> alpha re-solved (B = probe claim) -> 0.3 learned blend (billmA H1-fold weights) -> v5b emission ->
v9 layout at budget 885. Arms change exactly one thing:
  budget_450 / budget_650 : layout budget
  deamp_0.6 / 0.9 / 1.2   : an independent per-cell CPM head at that log2FC amplitude (bulk head unchanged), fitted by
                            AtlasShift dual_moment_counts at the pooled template, then the same layout
  pie_885                 : call ORDER by Arc PIE centred p_de among supported genes (LEAKY on H1: PIE saw Arc's 2025 H1
                            data; treat as an upper bound only)
Writes only under final_push/artifacts/de_check/<task>/levers_<arm>/.
"""
import os, sys, types, gc, time, json, argparse
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; VCC = HERE.parent
sys.path.insert(0, str(HERE))
import de_check
OUT = de_check.OUT
from infer import Ensemble
from network import blend


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--task', default='H1_rich_search')
    ap.add_argument('--arms', default='base_885,budget_450,budget_650,deamp_0.6,deamp_0.9,deamp_1.2,pie_885'); a = ap.parse_args(); task = a.task
    benchmark, scoring, bounded_score, baseline, rsi = de_check.rsi_context(); helpers = baseline.helpers
    from scipy import sparse
    contract = baseline.read_spec(); model, ctrl, real = benchmark.setup(task); fold = model['context']
    genes = ctrl.var_names.to_numpy(str); targets = np.array(model['targets']); T = len(targets); gi = {g: i for i, g in enumerate(genes)}
    ensemble = Ensemble(fold); assert fold not in {b['family'] for b in ensemble.src.banks}
    canonical = ensemble.src.genes; lookup = {g: i for i, g in enumerate(canonical)}; indices = np.array([lookup[g] for g in genes])
    bulk = np.asarray(ctrl.X.sum(0, dtype=np.float64)).ravel(); bulk /= bulk.sum()
    cc = np.zeros(len(canonical)); cc[indices] = bulk
    learned = ensemble.predict(targets, cc)['learned'][:, indices]; del ensemble; gc.collect()
    q4, depth, bf0, amp = baseline.calibrated_profile(ctrl, task, contract)
    _, _, L, Lc, support, pinfo = baseline.profile(ctrl, task, contract)
    # sources: confidence key (v9), support count (allocation knob), and the fused log2FC effect for the CPM head
    with baseline.local_v4(fold) as v4:
        names = v4.load_banks(None, True); sources = [v4.aligned_bank(n, targets, genes, 'log2fc') for n in names]
        cd = v4.cd4_family(targets, genes, 'log2fc')
        if cd is not None:
            sources.append(cd)
        _, _, mean, _ = v4.template(ctrl, genes, pool_k=4, seed=baseline.SEED)
        (p06, _), _ = v4.heads(targets, genes, mean, bulk, baseline.VARIANT, exclude=fold, stats={})
    atlas = baseline.recipes.v4.atlas
    sign = np.sign(L - Lc); net = np.zeros(L.shape); zs = np.zeros(L.shape)
    for e, m, v in sources:
        z = np.where(m & np.isfinite(v), e / np.sqrt(np.maximum(v, 1e-12)), 0); sig = np.abs(z) > 1.96
        net += (sig * (np.sign(z) == sign)).astype(float) - (sig * (np.sign(z) == -sign)).astype(float)
        zs += np.minimum(np.abs(z), 100) * (np.sign(z) == sign)
    confidence = net * 1e6 + zs; del sources, cd; gc.collect()
    # --- main_a2 direction knobs, exactly as final_push/main/panel_pipeline.py applies them
    sup = (net >= 1).sum(1).astype(float); conf = sup / max(sup.max(), 1.0); f_conf = np.maximum(conf, 1e-3); f_conf = f_conf / f_conf.mean()
    fs = np.ones(T)
    for i, t in enumerate(targets):
        j = gi.get(t)
        if j is not None:
            cpm = float(bulk[j] * 1e6); fs[i] = float(np.clip((np.log(max(cpm, 1e-9)) - np.log(1.0)) / np.log(10.0), 0, 1))
    f_sil = 0.25 + 0.75 * fs; keep = bulk * 1e6 >= 1.0
    print('knob factors: conf <0.2:', int((f_conf < 0.2).sum()), '<0.5:', int((f_conf < 0.5).sum()), '| silence <1:', int((f_sil < 1).sum()), '| compat genes kept', int(keep.sum()), '/', len(genes), flush=True)

    def profile_for(compat, conf_on, sil_on):
        delta = L - Lc
        if compat:
            delta = delta * keep[None, :]
        fac = (f_conf if conf_on else 1.0) * (f_sil if sil_on else 1.0)
        Lk = Lc + delta * (np.asarray(fac) * np.ones(T))[:, None]
        al, Bk = helpers.solve_alpha(Lk, targets, genes, amp['claim_est_probe'], center=Lc)
        bfk = helpers.from_logbulk(helpers.rescale(Lk, al, Lc))
        return helpers.from_logbulk(Lc + blend(helpers.logbulk(bfk) - Lc, learned, .3)), al
    KNOBS = {'knob_none': (False, False, False), 'knob_compat': (True, False, False), 'knob_conf': (False, True, False), 'knob_silence': (False, False, True)}
    new, alpha = profile_for(True, True, True)
    print('main_a2-like profile: alpha', round(alpha, 4), 'claim', round(amp['claim_est_probe'], 5), flush=True)
    # CPM head at amplitude a from the fused effect recovered out of the v4 .6-amplitude head
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.where((p06 > 0) & (mean[None, :] > 0), np.log2(p06 / mean[None, :]) / 0.6, 0.0)
    eff = np.clip(np.nan_to_num(eff), -5, 5)

    def cpm_head(aa):
        return atlas.desired_mean(mean, eff, space='log2fc', amplitude=aa, clip=3)
    # v9 layout with its own runtime stub
    rsi_rt = sys.modules['runtime']; stub = types.ModuleType('runtime'); stub.ROOT = OUT; stub.json = rsi_rt.json; stub.sha = rsi_rt.sha
    stub.save = rsi_rt.save; stub.guard = lambda *a, **k: None; sys.modules['runtime'] = stub
    sys.path.insert(0, str(VCC / 'model_v8')); sys.path.insert(0, str(VCC / 'model_v9'))
    import de as v9de, layout as v9layout
    refdir = OUT / 'artifacts' / f'ref_{task}'
    if not (refdir / 'cdf_scale.npy').exists():
        v9de.prepare(ctrl.X, refdir, np.random.default_rng(20261018))
    ref = v9de.Reference(task); sys.modules['runtime'] = rsi_rt
    # PIE ordering key (leaky on H1, see docstring): supported genes ordered by centred p_de
    pie_key = None
    pie_path = VCC / 'model_v6/pie/out' / f'ens4_{fold}.npz'
    if pie_path.exists():
        z = np.load(pie_path); ri = {t: i for i, t in enumerate(z['targets'].astype(str))}; ci = {g: i for i, g in enumerate(z['genes'].astype(str))}
        P = np.full((T, len(genes)), np.nan, np.float32); rows = [i for i, t in enumerate(targets) if t in ri]; cols = [j for j, g in enumerate(genes) if g in ci]
        P[np.ix_(rows, cols)] = z['p_de'].astype(np.float32)[np.ix_([ri[targets[i]] for i in rows], [ci[genes[j]] for j in cols])]
        C = P - np.nanmean(P, 0, keepdims=True); C = np.nan_to_num(C, nan=-1.0); Cn = (C - C.min()) / max(C.max() - C.min(), 1e-9)
        pie_key = np.where(net >= 1, 1e6 + 1e5 * Cn, 0.0)
    w = helpers._row_weights(T, len(genes), helpers.excl_cols(targets, genes))

    def jk(x):
        return float(helpers.jackknife_correction(sparse.csr_matrix(x), np.zeros(x.shape[0], np.intp), 1, helpers.TS, chunk=32)[0])

    def own_of(t):
        return int(gi[t]) if t in gi else np.array([], dtype=np.int64)
    anchors = rsi_rt.json.loads((rsi / f'artifacts/scores/{task}/anchors.json').read_text()); assert anchors['real_sha256'] == scoring.matrix_id(real)
    results = {}
    for arm in a.arms.split(','):
        t0 = time.time(); out = OUT / task / f'levers_{arm}'; out.mkdir(parents=True, exist_ok=True)
        budget = {'budget_450': 450, 'budget_650': 650}.get(arm, 885)
        key = pie_key if arm == 'pie_885' else confidence
        prof_arm, alpha_arm = (profile_for(*KNOBS[arm]) if arm in KNOBS else (new, alpha))
        pa = cpm_head(float(arm.split('_')[1])) if arm.startswith('deamp') else None
        P = np.zeros((T, len(genes))); jks = np.zeros(T); calls = np.zeros(T, int); projection = []

        def emit(i):
            t = targets[i]; seed = atlas.seed_for(fold + ':' + t, baseline.SEED)
            if pa is None:
                original = baseline.emit(q4, depth, prof_arm[i], fold, t).toarray().astype(np.int32)
            else:
                original = np.asarray(atlas.dual_moment_counts(q4, pa[i], prof_arm[i], depths=depth, seed=seed)).astype(np.int32)
            x, called, info = v9layout.emit(original, ref, key[i], budget, own_of(t)); calls[i] = len(called)
            return sparse.csr_matrix(x)

        def iterator():
            for i, t in enumerate(targets):
                x = emit(i); jks[i] = jk(x); P[i] = np.asarray(x.sum(0)).ravel()
                if i % 10 == 0:
                    print(arm, 'EMIT', i, '/', T, 'elapsed', round(time.time() - t0), flush=True)
                yield t, x
        try:
            cfg = scoring.config(cache_real=str(OUT / task / 'real_cache'))
            tidy, digest = bounded_score.score_iter(real, iterator(), cfg)
        except ValueError as err:
            print('ARM FAILED', arm, err, flush=True); results[arm] = dict(arm=arm, failed=str(err)); continue
        raw = benchmark.aggregate(tidy); scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
        Pn = P / P.sum(1, keepdims=True); Bm = float(helpers.budget(Pn, targets, genes)); claim = float((w * jks).sum())
        rec = dict(task=task, arm=arm, raw=raw, scaled=scaled, overall=overall, prediction_sha256=digest, B_emitted=Bm, claim_emitted=claim,
                   ratio_B_over_claim=Bm / claim, budget=budget, calls_mean=float(calls.mean()), alpha=alpha_arm, seconds=round(time.time() - t0),
                   pie_leaky=(arm == 'pie_885'), status='official scorer on a public LOCO SEARCH task; not an official challenge score')
        tidy.write_csv(out / 'tidy.csv'); rsi_rt.save(out / 'score.json', rec); results[arm] = rec
        print('SCORED', task, arm, 'overall', round(overall, 5), {k: round(v, 4) for k, v in scaled.items()}, 'rawMSE', round(raw['expr_mse_unbiased_capped_norm'], 4), 'calls', round(calls.mean()), flush=True)
        gc.collect()
    rsi_rt.save(OUT / task / 'levers_summary.json', results); print('DONE', task, flush=True)


if __name__ == '__main__':
    main()
