"""Leave-one-family-out selector runs. Usage:
  python run_selector.py --tag base --recipes base,eq,panel ... [--contexts HCT116,HEK293T,...] [--panel 300]
Each context = one public screen held out with its whole family; the panel = up to --panel of its targets (seeded draw)
that at least one other family measures; truth = its measured log-bulk response. Writes artifacts/selector_<tag>.csv/json.
"""
import argparse, json, time, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import banks as B, recipes as R, metrics as M

ART = Path(__file__).resolve().parent / 'artifacts'

RECIPES = {
    # name: knobs. 'v5dir' approximates the scored v5 direction family (v4 weights, bank-internal centring, sign gate).
    'base':        dict(center='source', weights='equal'),
    'baseline':    dict(center='source', weights='equal', shared_mean=1),   # the official baseline construction: one mean response for every target
    'v4w':         dict(center='source', weights='v4'),
    'panelc':      dict(center='panel', weights='equal'),
    'nocentre':    dict(center='none', weights='equal'),
    'v5dir':       dict(center='panel', weights='v4', gate=1),
    'gate':        dict(center='source', weights='equal', gate=1),
    'compat':      dict(center='source', weights='equal', compat=1.0),
    'silence':     dict(center='source', weights='equal', silence=(1.0, 10.0)),
    'conf_nsrc':   dict(center='source', weights='equal', conf='nsrc'),
    'conf_agree':  dict(center='source', weights='equal', conf='agree'),
    'all_soft':    dict(center='source', weights='equal', compat=1.0, silence=(1.0, 10.0), conf='agree'),
    # second round: the silence gate with a FLOOR (a zeroed row scores pds 0.5, which is why the hard gate cost pds)
    'silence_fl':  dict(center='source', weights='equal', silence=(1.0, 10.0), silence_floor=0.25),
    'compat5':     dict(center='source', weights='equal', compat=5.0),
    'best':        dict(center='source', weights='v4', compat=1.0, conf='agree'),
    'best_sil':    dict(center='source', weights='v4', compat=1.0, conf='agree', silence=(1.0, 10.0), silence_floor=0.25),
    'best_gate':   dict(center='source', weights='v4', compat=1.0, conf='agree', gate=1),
    # two-tier ("not all genes are equal"): silent tier -> zero centred response (= the baseline), strong -> full
    'tier':        dict(center='source', weights='v4', tier=(0.33, 0.67), tier_weak=0.6),
    'tier_best':   dict(center='source', weights='v4', compat=1.0, conf='agree', tier=(0.33, 0.67), tier_weak=0.6),
    'tier_hard':   dict(center='source', weights='v4', compat=1.0, conf='agree', tier=(0.33, 0.67), tier_weak=0.3),
    'tier_soft':   dict(center='source', weights='v4', compat=1.0, conf='agree', tier=(0.20, 0.80), tier_weak=0.8),
    'tier_sim':    dict(center='source', weights='v4', compat=1.0, conf='agree', tier=(0.33, 0.67), tier_weak=0.6, tier_sim=1),
    'respgene':    dict(center='source', weights='v4', compat=1.0, conf='agree', respgene=0.5),
    'tier_resp':   dict(center='source', weights='v4', compat=1.0, conf='agree', tier=(0.33, 0.67), tier_weak=0.6, respgene=0.5),
}


def predict(all_banks, targets, genes, family, knobs, dest_control, learned=None):
    entries = R.gather(all_banks, targets, family, min_cells=knobs.get('min_cells', 30))
    entries = R.centre(entries, targets, genes, knobs.get('center', 'source'))
    pred, nsrc, nfam, cov = R.pool(entries, targets, genes, weights=knobs.get('weights', 'equal'))
    info = dict(n_sources_mean=float(nsrc.mean()), n_families_mean=float(nfam.mean()), sources=[b.source for b, _, _ in entries])
    agree = R.agreement(entries, pred, targets, genes) if knobs.get('conf', 'none').startswith('agree') else np.full(len(targets), np.nan)
    if knobs.get('respgene'):
        keep, rinfo = R.responsive_gene_mask(entries, genes, quantile=float(knobs['respgene']))
        pred = pred * keep[None, :]; info['respgene'] = rinfo
    if knobs.get('tier'):
        sim = None
        if knobs.get('tier_sim'):
            dst = np.log1p(50000 * dest_control)
            sim = {b.source: float(np.exp(-np.mean(np.abs(dst[b.measured] - np.log1p(50000 * b.source_mean()['control'][b.measured])))))
                   for b, _, _ in entries}
        score = R.strength_score(entries, targets, sim)
        q_silent, q_strong = knobs['tier']
        tier = R.tiers(score, q_silent, q_strong)
        pred, f = R.apply_tiers(pred, tier, weak=knobs.get('tier_weak', 0.6))
        info['tier'] = dict(silent=int((tier == 0).sum()), weak=int((tier == 1).sum()), strong=int((tier == 2).sum()),
                            unknown=int((~np.isfinite(score)).sum()), factor_mean=float(f.mean()))
    if knobs.get('shared_mean'):
        pred = np.repeat(pred.mean(0, keepdims=True), len(targets), 0)
    if knobs.get('gate'):
        pred, mask = R.sign_gate(entries, pred, targets, genes, min_net=int(knobs['gate']), tau=knobs.get('gate_tau', 0.05))
        info['gate_calls_mean'] = float(mask.sum(1).mean())
    if learned is not None:
        pred = R.blend_learned(pred, learned, knobs.get('mix', 0.3))
    if knobs.get('compat'):
        pred = R.compat_mask(pred, dest_control, genes, tau_cpm=float(knobs['compat']))
    if knobs.get('silence'):
        lo, hi = knobs['silence']; pred, f = R.silence_gate(pred, targets, genes, dest_control, lo, hi, floor=knobs.get('silence_floor', 0.0))
        info['silence_factor_mean'] = float(f.mean()); info['silenced_targets'] = int((f == 0).sum())
    if knobs.get('conf', 'none') != 'none':
        pred, f = R.confidence_scale(pred, nsrc, agree, knobs['conf']); info['conf_factor_mean'] = float(f.mean())
    return pred, info


def context_panel(bank, all_banks, size, seed, min_cells=30):
    cand = [t for t in bank.targets if bank.n_cells[t] >= min_cells and
            any(b.family != bank.family and t in b.index and b.n_cells[t] >= min_cells for b in all_banks.values())]
    cand = np.array(sorted(cand))
    if size and len(cand) > size:
        cand = np.sort(np.random.default_rng(seed).choice(cand, size, replace=False))
    return cand


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True); ap.add_argument('--recipes', default='base')
    ap.add_argument('--contexts', default='HCT116,HEK293T,CD4_Rest,CD4_Stim8hr,CD4_Stim48hr,H1,K562_GWPS,HepG2,Jurkat,RPE1,K562')
    ap.add_argument('--panel', type=int, default=300); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--no-shards', action='store_true'); ap.add_argument('--learned', action='store_true', help='billmA LOCO blend on H1/CD4 contexts')
    a = ap.parse_args()
    all_banks = B.discover(include_shards=not a.no_shards)
    print('banks:', {k: (len(v.targets), v.family) for k, v in all_banks.items()}, flush=True)
    rows = []; t0 = time.time()
    for ctx in a.contexts.split(','):
        bank = all_banks[ctx]; genes = bank.genes
        targets = context_panel(bank, all_banks, a.panel, a.seed)
        if len(targets) < 10:
            print('skip', ctx, len(targets), flush=True); continue
        truth = np.nan_to_num(bank.rows(targets)); ctrl_rows = bank.rows(targets, 'control')
        dest_control = np.nanmean(ctrl_rows, 0).astype(np.float64); dest_control /= dest_control.sum()
        learned = None
        if a.learned and bank.family in ('H1', 'CD4'):
            from learned import billm_learned
            learned = billm_learned(bank.family, targets, dest_control, genes)
        for name in a.recipes.split(','):
            knobs = RECIPES[name]; t1 = time.time()
            pred, info = predict(all_banks, targets, genes, bank.family, knobs, dest_control, learned=learned)
            m = M.compute(pred, truth, targets, genes, bank.measured)
            rec = dict(context=ctx, family=bank.family, recipe=name, knobs=json.dumps(knobs), seconds=round(time.time() - t1, 1),
                       **{k: v for k, v in m.items() if not isinstance(v, np.ndarray)}, **{k: v for k, v in info.items() if k != 'sources'})
            rows.append(rec)
            print(f"{ctx:12s} {name:12s} rho_bar {m['rho_bar']:.4f} oracle {m['rho_bar_oracle_amp']:.4f} pds {m['pds']:.4f} "
                  f"cos {m['cos_mean']:.4f} mag {m['magnitude_spearman']:.3f} n {m['n_targets']} src {info['n_sources_mean']:.2f} ({rec['seconds']}s)", flush=True)
    df = pd.DataFrame(rows); ART.mkdir(exist_ok=True)
    df.to_csv(ART / f'selector_{a.tag}.csv', index=False)
    macro = df.groupby('recipe')[['rho_bar', 'rho_bar_centred', 'rho_bar_oracle_amp', 'pds', 'cos_mean', 'magnitude_spearman', 'truth_mean_share', 'pred_mean_share']].mean()
    macro.to_csv(ART / f'selector_{a.tag}_macro.csv')
    print(macro.round(4).to_string()); print('total seconds', round(time.time() - t0))


if __name__ == '__main__':
    main()
