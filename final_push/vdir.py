"""Calibration: score the ACTUAL scored directions (v4a fusion, v5b gate, billmA learned blend) on held-out public
contexts with the real model_v6 code path, so the selector can be checked against our official history.

v4a   = model_v6 heads 'anchor_single_mean' (AtlasShift source-centred fusion + CD4 family), held-out family excluded.
v5bg  = v4a + the v5b cross-source DE call gate (v5.de_call_mask, uncapped) applied to the log-bulk response.
billm = v4a (or v5bg) + Bill M learned-direction blend at 0.3, using the LOCO checkpoint for that family.
Deflation is NOT applied (it needs the destination's own control CELLS, which these summary banks do not carry) --
it is a within-context operation that does not change the transfer direction being compared here.

Official references: v4a 0.1438 (pds .678), v5b 0.1528 (pds .680), billmA 0.16438 (pds .706).
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
VCC = HERE.parent
sys.path.insert(0, str(HERE))
import banks as B, metrics as M

ART = HERE / 'artifacts'


def load_v6():
    sys.path.insert(0, str(VCC / 'model_v6'))
    import model, v5
    return model, v5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--contexts', default='H1,CD4_Rest,CD4_Stim8hr,CD4_Stim48hr,HCT116,HEK293T')
    ap.add_argument('--panel', type=int, default=300)
    ap.add_argument('--tag', default='vdir1')
    a = ap.parse_args()
    model, v5 = load_v6()
    genes = pd.read_csv(VCC / 'data/gene_names.csv').iloc[:, 0].astype(str).to_numpy()
    if genes[0] in ('gene_name', 'gene'):
        genes = genes[1:]
    all_banks = B.discover(include_shards=False)
    rows = []
    for ctx in a.contexts.split(','):
        bank = all_banks[ctx]
        assert np.array_equal(bank.genes, genes), f'{ctx} is not on the canonical axis'
        # panel: this screen's targets that the v4 bank universe also knows (the 437-target cut), >= 30 cells
        universe = set(np.load(VCC / 'public_data/model_v4_signature_banks/HCT116.npz', allow_pickle=False)['targets'].astype(str))
        cand = np.array(sorted(t for t in bank.targets if bank.n_cells[t] >= 30 and t in universe and t in set(genes)))
        if len(cand) > a.panel:
            cand = np.sort(np.random.default_rng(0).choice(cand, a.panel, replace=False))
        if len(cand) < 20:
            print('skip', ctx, len(cand), flush=True); continue
        truth = np.nan_to_num(bank.rows(cand))
        ctrl = np.nanmean(bank.rows(cand, 'control'), 0).astype(np.float64); ctrl = ctrl / ctrl.sum()
        model.configure(VCC / 'public_data/model_v4_signature_banks', [bank.family], VCC / 'model_v6/artifacts/promoter_pairs_v4.csv')
        t0 = time.time()
        (p, b), cov = model.heads(cand, genes, ctrl, ctrl, 'anchor_single_mean')
        Lc = v5.logbulk(ctrl)[None, :]; L = v5.logbulk(b)
        base = L - Lc
        variants = {'v4a': base}
        src = [model.aligned_bank(n, cand, genes, 'log2fc') for n in model.load_banks(None, True)]
        cd = model.cd4_family(cand, genes, 'log2fc')
        if cd is not None:
            src.append(cd)
        if src:
            mask, net = v5.de_call_mask(src, base, k=len(genes))
            variants['v5bg'] = np.where(mask, base, 0.0)
        del src, cd
        if bank.family in ('H1', 'CD4'):
            from learned import billm_learned
            learned = billm_learned(bank.family, cand, ctrl, genes)
            import importlib; R = importlib.import_module('recipes')
            variants['billm_on_v4a'] = R.blend_learned(variants['v4a'], learned, 0.3)
            if 'v5bg' in variants:
                variants['billm_on_v5bg'] = R.blend_learned(variants['v5bg'], learned, 0.3)
        for name, pred in variants.items():
            m = M.compute(np.asarray(pred, np.float32), truth, cand, genes, bank.measured)
            rec = dict(context=ctx, family=bank.family, variant=name,
                       **{k: v for k, v in m.items() if not isinstance(v, np.ndarray)})
            rows.append(rec)
            print(f"{ctx:12s} {name:14s} rho_bar {m['rho_bar']:.4f} centred {m['rho_bar_centred']:.4f} pds {m['pds']:.4f} "
                  f"cos {m['cos_mean']:.4f} n {len(cand)}", flush=True)
        print(ctx, 'seconds', round(time.time() - t0, 1), flush=True)
    df = pd.DataFrame(rows); df.to_csv(ART / f'{a.tag}.csv', index=False)
    print(df.groupby('variant')[['rho_bar', 'rho_bar_centred', 'pds', 'cos_mean']].mean().round(4).to_string())


if __name__ == '__main__':
    main()
