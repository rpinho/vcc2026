"""Evaluate final_push/GATES.md on a panel built by model_billmB/panel_pipeline.py + panel_prep.py.
Usage: python check_gates.py NAME [--emitted]   (--emitted recomputes B from the retained H5 group sums; minutes)
Prints the gate table and writes final_push/artifacts/gates_<NAME>.json. Exit 0 only if every gate passes."""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd

VCC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VCC / 'model_v6'))


def main():
    name = sys.argv[1]; emitted = '--emitted' in sys.argv
    root = VCC / 'model_billmB'
    for arg in sys.argv[2:]:
        if arg.startswith('--root='):
            root = Path(arg.split('=', 1)[1]).resolve()
    out = root / 'artifacts/panels' / name
    gates = []; ok_all = True

    def gate(g, desc, value, passed, note=''):
        nonlocal ok_all; ok_all &= bool(passed)
        gates.append(dict(gate=g, desc=desc, value=value, passed=bool(passed), note=note))

    profiles = sorted(out.glob('profile_*.json'))
    panel = pd.read_csv(VCC / 'data/final/bundle/pert_counts.csv'); T = int((panel.target_gene != 'non-targeting').sum())
    for p in profiles:
        r = json.loads(p.read_text()); info = r['profile_info']; c = r['context']
        cov = info.get('targets_with_source_row'); gate('G1', f'{c}: targets with a source row / {T}', cov, cov is not None and cov >= 0.9 * T)
        gate('G2', f'{c}: alpha < 8', r['alpha'], 0 < r['alpha'] < 8)
        ratio = r['B'] / r['claim']; gate('G3', f'{c}: requested B/claim', round(ratio, 4), 0.9 <= ratio <= 1.05, 'requested profile; emitted below if --emitted')
    fg = out / 'full_gates.json'
    if fg.exists():
        g = json.loads(fg.read_text()); recs = g['records']
        gate('G8', 'all groups all_gates_pass', f"{sum(r['all_gates_pass'] for r in recs)}/{len(recs)}", g['status'] == 'PASS' and all(r['all_gates_pass'] for r in recs))
        ml1 = max(r['mean_certificate']['mean_L1_upper'] for r in recs); bl1 = max(r['mean_certificate']['bulk_L1_upper'] for r in recs)
        gate('G6', 'max mean L1 / bulk L1 over groups', [round(ml1, 5), round(bl1, 6)], ml1 < .015 and bl1 < .0015)
        if emitted:
            import h5py
            from scipy import sparse
            from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
            genes = pd.read_csv(VCC / 'data/gene_names.csv').gene_name.to_numpy(str)
            h5 = root / 'out' / f'{name}.h5ad'
            with h5py.File(h5) as f:
                X = f['X']; ptr = X['indptr'][...]
                ctxs = sorted({r['context'] for r in recs}); per = {}
                for ci, c in enumerate(ctxs):
                    rows = [r for r in recs if r['context'] == c]; rows.sort(key=lambda r: r['index'])
                    targets = [r['target'] for r in rows]; P = np.zeros((len(rows), len(genes)))
                    for i, r in enumerate(rows):
                        lo = (ci * len(rows) + i) * 400; a, b = int(ptr[lo]), int(ptr[lo + 400])
                        m = sparse.csr_matrix((X['data'][a:b], X['indices'][a:b], ptr[lo:lo + 401] - a), shape=(400, len(genes)))
                        P[i] = np.asarray(m.sum(0)).ravel()
                    L = np.log1p(50000 * P / P.sum(1, keepdims=True))
                    gi = {g_: k for k, g_ in enumerate(genes)}; cols = {k: gi[t] for k, t in enumerate(targets) if t in gi}
                    w = _row_weights(len(rows), len(genes), cols); B = float(_across_pert_budget(L, cols, w))
                    claim = float(sum(w[i] * r['jk'] for i, r in enumerate(rows)))
                    per[c] = dict(B_emitted=B, claim_emitted=claim, ratio=B / claim)
                    # GATES.md asymmetry: the upper bound is the hard gate; below 0.9 is unclaimed credit (v9-layout regime)
                    gate('G3e', f'{c}: EMITTED B/claim (<= 1.05 hard)', round(B / claim, 4), B / claim <= 1.05,
                         'below 0.9 = unclaimed credit, expected for the v9 layout' if B / claim < 0.9 else '')
    dr = out / 'decode_receipt.json'
    if dr.exists():
        d = json.loads(dr.read_text())
        gate('G4', 'vcc prep exit code', d['exit_code'], d['exit_code'] == 0)
        gate('G5', 'independent decode', d['status'], d['status'].startswith('PASS'))
        s = d['structure']; gate('G7', 'shape', s['shape'], tuple(s['shape']) == (400 * T * 3, 18533))
        gates.append(dict(gate='pkg', desc='package', value=dict(bytes=d['bytes'], sha256=d['sha256']), passed=True))
    else:
        gate('G4', 'vcc prep', 'no decode_receipt.json yet', False)
    df = pd.DataFrame(gates); print(df.to_string(index=False))
    print('ALL GATES PASS' if ok_all else 'GATES FAILED')
    (VCC / 'final_push/artifacts' / f'gates_{name}.json').write_text(json.dumps(dict(name=name, all_pass=ok_all, gates=gates), indent=2, default=str))
    sys.exit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
