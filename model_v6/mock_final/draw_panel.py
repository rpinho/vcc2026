"""Mock-final panel: 300 FRESH targets drawn the way the validation panel looks (PIE_AND_EVIDENCE_PLAN B0):
non-essential (no Replogle-essential target), outside the 437-target v4 bank universe (so every source must be
fetched fresh), ~91% K562-GWPS-covered, ~55% GWPS responders (energy p < 0.05) among covered, in the official
gene list. Fixed SHA-256 order; no model score used."""
import sys, json, hashlib, re
sys.path.insert(0, '..')
import runtime
import numpy as np, pandas as pd, anndata as ad
V = runtime.VCC
genes = set(pd.read_csv(V / 'data/gene_names.csv').gene_name.astype(str))
universe = set(np.load(V / 'public_data/model_v4_signature_banks/HCT116.npz')['targets'].astype(str))
ess = set()
for n in ['K562', 'RPE1', 'HepG2', 'Jurkat']:
    a = ad.read_h5ad(V / f'public_data/{n}.h5ad', backed='r'); col = 'gene' if 'gene' in a.obs else 'target_gene'; ess |= set(a.obs[col].astype(str)); a.file.close()
b = ad.read_h5ad(V / 'public_data/K562_GWPS_bulk.h5ad', backed='r'); obs = b.obs.copy(); b.file.close()
obs['t'] = [re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$', s).group(1) if 'non-targeting' not in s else 'non-targeting' for s in obs.index]
obs = obs[(obs.t != 'non-targeting') & (obs.num_cells_filtered >= 50)].drop_duplicates('t')
key = lambda t: hashlib.sha256(('mockfinal:' + t).encode()).hexdigest()
ok = lambda t: t in genes and t not in ess and t not in universe
cov = obs[[ok(t) for t in obs.t]]
resp = sorted(cov.t[cov.energy_test_p_value < .05], key=key); non = sorted(cov.t[cov.energy_test_p_value >= .05], key=key)
n_cov = 272; n_resp = round(.55 * n_cov)
gw = set(obs.t)
uncovered = sorted([g for g in genes if ok(g) and g not in gw and not g.startswith(('MT-', 'RPL', 'RPS')) and re.match(r'^[A-Z][A-Z0-9]+$', g)], key=key)[:300 - n_cov]
panel = resp[:n_resp] + non[:n_cov - n_resp] + uncovered
assert len(set(panel)) == 300
pd.DataFrame({'target_gene': panel}).to_csv('panel.csv', index=False)
runtime.save('panel_receipt.json', dict(n=300, gwps_covered=n_cov, gwps_responders=n_resp, uncovered_by_gwps=len(uncovered),
             essential_overlap=0, overlap_v4_universe=0, rule=__doc__))
print(len(panel), panel[:10], uncovered[:10])
