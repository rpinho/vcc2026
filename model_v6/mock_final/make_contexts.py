"""Mock-final contexts D/E/F = control (non-targeting) cells of three public held-out lines (Replogle 2022
HepG2 / Jurkat / RPE1, raw counts), put on the official 18,533-gene axis exactly as a released context file
would be (genes the screen did not measure are zero columns: MOCK ONLY, the real D/E/F files carry all genes).
Up to 18,400 cells per context like the official files."""
import sys; sys.path.insert(0, '..')
import runtime
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
V = runtime.VCC
genes = pd.read_csv(V / 'data/gene_names.csv').gene_name.astype(str).to_numpy()
out = {}
for label, line in zip('DEF', ['HepG2', 'Jurkat', 'RPE1']):
    a = ad.read_h5ad(V / f'public_data/{line}.h5ad', backed='r')
    col = 'gene' if 'gene' in a.obs else 'target_gene'
    rows = np.flatnonzero(a.obs[col].astype(str).str.startswith('non-targeting').to_numpy())
    rows = np.sort(np.random.default_rng(7).choice(rows, min(18400, len(rows)), replace=False))
    names = a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
    X = sparse.vstack([sparse.csr_matrix(a.X[rows[i:i + 2048]]) for i in range(0, len(rows), 2048)], format='csr'); a.file.close()
    uniq, first = np.unique(names, return_index=True)
    ix = pd.Index(names[np.sort(first)]).get_indexer(genes)
    X = X[:, np.sort(first)]
    valid = np.flatnonzero(ix >= 0)
    M = sparse.csr_matrix((np.ones(len(valid)), (ix[valid], valid)), shape=(X.shape[1], len(genes)))
    Y = (X @ M).tocsr(); Y.data = np.rint(Y.data); Y = Y.astype(np.int32); Y.eliminate_zeros()
    keep = np.asarray(Y.sum(1)).ravel() > 0
    c = ad.AnnData(Y[keep], obs=pd.DataFrame({'target_gene': 'non-targeting', 'context': label}, index=[f'{line}_{r}' for r in rows[keep]]), var=pd.DataFrame(index=genes))
    c.write_h5ad(f'context_{label}.h5ad', compression='gzip')
    out[label] = dict(line=line, cells=int(c.n_obs), measured_genes=int((ix >= 0).sum()), median_umi=float(np.median(np.asarray(Y[keep].sum(1)).ravel())))
    print(label, out[label], flush=True)
runtime.save('contexts_receipt.json', out)
